#!/usr/bin/env python3
"""
refineGT.py – Refine SV genotypes from a VCF using read‑backed evidence.
Compares the original genotype (from the VCF) with a recomputed genotype
and reports/fixes mismatches as requested:
    - if 1/1 -> 0/0 : report & assign new GT
    - if 0/0 -> 1/1 : report & assign new GT
    - if 0/1 -> 0/0 : report & assign new GT

Additionally, outputs a consensus VCF where only these hard conflicts
are masked to './.' ; all other discordant calls (e.g. 1/1 -> 0/1) keep
the new genotype.
"""

from os import makedirs, path
import subprocess
import multiprocessing
import argparse
from time import time

import pandas as pd
import pysam
from tqdm import tqdm

# Import the existing genotype callers (must be in your PYTHONPATH)
from sub_srSVGT import delGT, insGT, traGT, breaks2GT


# ----------------------------------------------------------------------
# Utility functions
# ----------------------------------------------------------------------
def makedir(dirpath):
    if not path.exists(dirpath):
        makedirs(dirpath)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def split_region(start, end):
    step = (end - start) / 7
    points = [start + i * step for i in range(8)]
    del points[0]
    del points[-1]
    return points

def usepysam_count_localreads(opened_sam, chrom, start, end, maq):
    reads = opened_sam.fetch(reference=chrom, start=start, end=end)
    high_quality_reads = [r for r in reads if r.mapping_quality >= maq]
    return len(high_quality_reads)


# ----------------------------------------------------------------------
# Core function to compute new genotype for a single SV record
# ----------------------------------------------------------------------
def compute_new_gt(vcf_line, opened_sam, name, min_maq, homo_rate, ref_rate, shift, span):
    """
    Takes a VCF line (string) and returns a dict with SV info and the new genotype.
    """
    fields = vcf_line.strip().split("\t")
    if len(fields) < 10:
        return None

    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]           # e.g., <DEL>, <INS>, <INV>, <DUP>, <TRA>
    info = fields[7]
    fmt = fields[8] if len(fields) > 8 else "GT"
    sample_col = fields[9] if len(fields) > 9 else ""

    # Extract original genotype from the sample column (first field of FORMAT:GT)
    orig_gt = sample_col.split(":")[0] if sample_col else "./."

    # Determine SV type from ALT or INFO (INFO has SVTYPE=...)
    sv_type = ""
    if "<" in alt and ">" in alt:
        sv_type = alt.strip("<>").upper()
    elif "SVTYPE=" in info:
        sv_type = info.split("SVTYPE=")[1].split(";")[0]
    else:
        # Fallback: look for keywords in info
        for kw in ["DEL", "INS", "INV", "DUP", "TRA"]:
            if kw in info:
                sv_type = kw
                break

    # Parse start/end/size from VCF columns and INFO/ID
    sv_start = pos
    sv_end = None
    sv_size = 0

    # END is usually in INFO
    if "END=" in info:
        sv_end = int(info.split("END=")[1].split(";")[0])
    elif ":" in fields[2]:  # ID like "X:152928643-152928978_335"
        parts = fields[2].split(":")[1].split("_")[0]  # "152928643-152928978"
        if "-" in parts:
            sv_end = int(parts.split("-")[1])
    if sv_end is None:
        # For TRA, the mate breakpoint is in the ALT or INFO
        if sv_type == "TRA":
            # format: ALT like "chrXI:443829" or INFO gives END=chrXI:443829
            if ":" in alt and not alt.startswith("<"):
                sv_end_chr = alt.split(":")[0]
                sv_end = int(alt.split(":")[1])
            elif "END=" in info:
                end_val = info.split("END=")[1].split(";")[0]
                if ":" in end_val:
                    sv_end_chr, sv_end_str = end_val.split(":")
                    sv_end = int(sv_end_str)
        else:
            # should not happen for other types
            sv_end = sv_start

    # Size
    if "SVLEN=" in info:
        sv_size = abs(int(info.split("SVLEN=")[1].split(";")[0]))
    elif sv_end is not None:
        sv_size = abs(sv_end - sv_start)
    else:
        sv_size = 0

    # Second chromosome for TRA/translocation
    chrome2 = None
    bp2 = None
    if sv_type == "TRA":
        if "END=" in info and ":" in info.split("END=")[1]:
            end_val = info.split("END=")[1].split(";")[0]
            chrome2, bp2_str = end_val.split(":")
            bp2 = int(bp2_str)
        elif ":" in alt:
            chrome2, bp2_str = alt.split(":")
            bp2 = int(bp2_str)
        else:
            # fallback: parse from ID
            id_parts = fields[2].split(":")[1].split("_")[0]
            chrome2, bp2_str = id_parts.split(":")
            bp2 = int(bp2_str)

    # Now call the appropriate genotype function (reusing your original logic)
    genotype = None
    if sv_type == "DEL":
        if sv_end is None:
            return None
        left_sam = opened_sam.fetch(reference=chrom, start=sv_start - shift, end=sv_start + shift)
        right_sam = opened_sam.fetch(reference=chrom, start=sv_end - shift, end=sv_end + shift)
        genotype = delGT(name, left_sam, right_sam, chrom, sv_start, sv_end, sv_size,
                         min_maq, homo_rate, ref_rate, shift, span)

    elif sv_type == "INS":
        if sv_end is None:
            return None
        region_sam = opened_sam.fetch(reference=chrom, start=sv_start - shift, end=sv_end + shift)
        genotype = insGT(name, region_sam, chrom, sv_start, sv_end, sv_size,
                         min_maq, homo_rate, ref_rate, shift, span)

    elif sv_type == "TRA":
        if chrome2 is None or bp2 is None:
            return None
        bp1_left = max(sv_start - shift, 0)
        bp2_left = max(bp2 - shift, 0)
        bp1_sam = opened_sam.fetch(reference=chrom, start=bp1_left, end=sv_start + shift)
        bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2_left, end=bp2 + shift)
        genotype = traGT(name, bp1_sam, bp2_sam, chrom, chrome2, sv_start, bp2, sv_size,
                         min_maq, "TRA", shift)
        if genotype[0] != "1/1":
            bp1_sam = opened_sam.fetch(reference=chrom, start=sv_start - shift, end=sv_start + shift)
            bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2 - shift, end=bp2 + shift)
            genotype = breaks2GT(name, bp1_sam, bp2_sam, chrom, chrome2, sv_start, bp2, sv_size,
                                 min_maq, "TRA", shift)

    elif sv_type in ("DUP", "INV"):
        if sv_end is None:
            return None
        bp1_sam = opened_sam.fetch(reference=chrom, start=sv_start - shift, end=sv_start + shift)
        bp2_sam = opened_sam.fetch(reference=chrom, start=sv_end - shift, end=sv_end + shift)
        genotype = breaks2GT(name, bp1_sam, bp2_sam, chrom, chrom, sv_start, sv_end, sv_size,
                             min_maq, sv_type, shift)

    if genotype is None:
        # Return original GT if computation failed
        genotype = [orig_gt, "0", "0"]

    return {
        "chrom": chrom,
        "pos": sv_start,
        "end": sv_end,
        "size": sv_size,
        "sv_type": sv_type,
        "orig_gt": orig_gt,
        "new_gt": genotype[0],          # e.g., "0/0", "0/1", "1/1"
        "total_map_reads": genotype[1] if len(genotype) > 1 else ".",
        "sv_support": genotype[2] if len(genotype) > 2 else "."
    }


def refine_vcf(args):
    """
    Main refinement pipeline.
    """
    makedir(args.dir)
    with open(args.sv_vcf, "r") as f:
        vcf_lines = f.readlines()

    header_lines = [line for line in vcf_lines if line.startswith("#")]
    data_lines = [line for line in vcf_lines if not line.startswith("#")]

    # Group by chromosome for parallelisation
    chrom_groups = {}
    for line in data_lines:
        chrom = line.split("\t")[0]
        chrom_groups.setdefault(chrom, []).append(line)

    tasks = []
    for chrom_lines in chrom_groups.values():
        tasks.append((chrom_lines, args.mapf, args.ACC, args.maq,
                      args.homo_rate, args.ref_rate, args.shift, args.span))

    # Process in parallel
    nproc = min(args.cpu, len(tasks)) if args.cpu > 0 else 1
    results = []
    with multiprocessing.Pool(processes=nproc) as pool:
        for res_list in tqdm(pool.starmap(process_chromosome_lines, tasks), total=len(tasks)):
            results.extend(res_list)

    # Build refined VCF, consensus VCF, and report
    refined_vcf_path = path.join(args.dir, f"{args.ACC}_refined.vcf")
    consensus_vcf_path = path.join(args.dir, f"{args.ACC}_consensus.vcf")
    report_path = path.join(args.dir, f"{args.ACC}_gt_changes.txt")

    refined_lines = []
    consensus_lines = []
    report_entries = []

    # Map results by (chrom, pos, sv_type) for quick lookup
    res_dict = {}
    for r in results:
        key = (r["chrom"], r["pos"], r["sv_type"])
        res_dict[key] = r

    # Pairs of (old_gt, new_gt) that are considered *hard conflicts* and
    # will be masked to "./." in the consensus VCF.
    MASK_PAIRS = {("1/1", "0/0"), ("0/0", "1/1"), ("0/1", "0/0")}

    for line in data_lines:
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        # Determine sv_type for matching
        alt = fields[4]
        info = fields[7]
        sv_type = ""
        if "<" in alt and ">" in alt:
            sv_type = alt.strip("<>").upper()
        elif "SVTYPE=" in info:
            sv_type = info.split("SVTYPE=")[1].split(";")[0]
        key = (chrom, pos, sv_type)

        # --- Refined VCF (always uses new GT if available) ---
        if key in res_dict:
            r = res_dict[key]
            old_gt = r["orig_gt"]
            new_gt = r["new_gt"]

            # Report certain hard changes (for logging)
            report_flag = False
            if old_gt == "1/1" and new_gt == "0/0":
                report_flag = True
            elif old_gt == "0/0" and new_gt == "1/1":
                report_flag = True
            elif old_gt == "0/1" and new_gt == "0/0":
                report_flag = True

            if report_flag:
                print(f"CHANGE: {chrom}:{pos} {sv_type} {old_gt} -> {new_gt}")
                report_entries.append(f"{chrom}\t{pos}\t{sv_type}\t{old_gt}\t{new_gt}")

            # Update GT for refined VCF
            sample_parts = fields[9].split(":")
            sample_parts[0] = new_gt
            fields[9] = ":".join(sample_parts)
            refined_lines.append("\t".join(fields))

            # --- Consensus VCF: mask only hard conflicts, otherwise keep new GT ---
            consensus_fields = fields.copy()  # work on a copy
            if (old_gt, new_gt) in MASK_PAIRS:
                # Mask genotype to missing, keep the rest of the sample fields
                consensus_parts = consensus_fields[9].split(":")
                consensus_parts[0] = "./."
                consensus_fields[9] = ":".join(consensus_parts)
                consensus_lines.append("\t".join(consensus_fields))
            else:
                # For all other cases (including 1/1->0/1, 0/0->0/1, etc.)
                # keep the new genotype as consensus.
                consensus_lines.append("\t".join(consensus_fields))
        else:
            # No recomputed genotype: keep original line for both outputs
            refined_lines.append(line.strip())
            consensus_lines.append(line.strip())

    # Write refined VCF
    with open(refined_vcf_path, "w") as out:
        out.writelines(header_lines)
        for l in refined_lines:
            out.write(l + "\n")

    # Write consensus VCF
    with open(consensus_vcf_path, "w") as out:
        out.writelines(header_lines)
        for l in consensus_lines:
            out.write(l + "\n")

    # Write report
    with open(report_path, "w") as rep:
        rep.write("#CHROM\tPOS\tSVTYPE\tOLD_GT\tNEW_GT\n")
        for entry in report_entries:
            rep.write(entry + "\n")

    print(f"Refined VCF written to {refined_vcf_path}")
    print(f"Consensus VCF written to {consensus_vcf_path}")
    print(f"Genotype changes written to {report_path}")


def process_chromosome_lines(chromosome_lines, mapf_path, name, min_maq, homo_rate, ref_rate, shift, span):
    """
    Worker function called by multiprocessing.
    Returns a list of result dicts for all lines of a chromosome.
    """
    output = []
    try:
        sam = pysam.AlignmentFile(mapf_path)
        for line in chromosome_lines:
            res = compute_new_gt(line, sam, name, min_maq, homo_rate, ref_rate, shift, span)
            if res is not None:
                output.append(res)
    except Exception as e:
        print(f"Error processing chromosome batch: {e}")
    finally:
        if 'sam' in locals():
            sam.close()
    return output


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Refine SV genotypes using read-backed evidence and report changes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", dest="sv_vcf", required=True,
                        help="Input SV VCF (already genotyped, e.g., from lr GT)")
    parser.add_argument("-mapf", dest="mapf", required=True,
                        help="BAM/CRAM file with aligned reads")
    parser.add_argument("-n", dest="ACC", required=True,
                        help="Sample/accession name")
    parser.add_argument("-o", dest="dir", required=True,
                        help="Output directory")
    parser.add_argument("-maq", dest="maq", type=int, default=1,
                        help="Minimum mapping quality for reads")
    parser.add_argument("-homo_rate", dest="homo_rate", type=float, default=0.65,
                        help="Homozygous ratio threshold")
    parser.add_argument("-ref_rate", dest="ref_rate", type=float, default=0.10,
                        help="Reference (0/0) ratio threshold")
    parser.add_argument("-shift", dest="shift", type=int, default=30,
                        help="Breakpoint shift range")
    parser.add_argument("-span", dest="span", type=int, default=50,
                        help="Span parameter for heterozygous evidence")
    parser.add_argument("-cpu", dest="cpu", type=int, default=10,
                        help="Number of CPUs to use")
    args = parser.parse_args()

    start_time = time()
    refine_vcf(args)
    print(f"******************* Cost time {time() - start_time:.2f}s *********************")
