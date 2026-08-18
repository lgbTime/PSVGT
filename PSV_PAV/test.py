#!/usr/bin/env python3
"""
Genomic PAV Matrix Read-Backed Genotype Rescuer
Assemble INS in short reads maybe difficult in low coverage datasets, but those SV has breakpoint evidence in short reads
This script is aim to rescue such case SV from absent 0 of PAV matrix to present 1

Description:
  Reads a collapsed PAV matrix, locates matching sample BAM files in a target directory, 
  and calls read-backed evidence functions to rescue false-negative 0s to 1s.
  Targeted filtering active: Only evaluates INS variants with size > 100 bp.
"""

import os
import sys
import argparse
from tqdm import tqdm
import pysam

# Import your native core validation functions directly
try:
    from srGT_refiner_for_asm_V1 import compute_new_gt
except ImportError:
    print("[-] Error: Could not import 'compute_new_gt' from srGT_refiner_for_asm_V1.py.", file=sys.stderr)
    print("[*] Please ensure this script is run from the directory containing your module.", file=sys.stderr)
    sys.exit(1)

def match_bam_file(sample_header, bam_dir):
    """
    Pairs sample text headers directly with file prefixes inside your BWA output directory.
    Handles headers like 'B09.contigs.fa', 'C12-2', 'S1', etc.
    """
    prefix = sample_header.split(".")[0]

    possible_names = [
        f"{prefix}_sorted.bam",
        f"{prefix.replace('-', '_')}_sorted.bam"
    ]

    for name in possible_names:
        full_path = os.path.join(bam_dir, name)
        if os.path.isfile(full_path):
            return full_path

    return None


def load_sample_table(table_path):
    """Parse a tab-delimited sample→BAM mapping file.

    Each line:  sample_name<TAB>/path/to/file.bam
    Lines starting with '#' are skipped.
    Returns {sample_name: bam_path}.
    """
    mapping = {}
    with open(table_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping

def generate_pseudo_vcf_line(row, sample_val="0/0"):
    """
    Formulates a virtual VCF data string to feed natively into compute_new_gt().
    """
    chrom = row["#Target_name"]
    pos = row["Target_start"]
    svid = row["SVID"]
    sv_type = row["SVType"]
    end = row["Target_end"]
    size = row["SVlen"]
    
    # Construct an exact standard VCF line layout representation
    ref = "N"
    alt = f"<{sv_type}>"
    qual = "60"
    filter_status = "PASS"
    info = f"SVTYPE={sv_type};END={end};SVLEN={size}"
    fmt = "GT"
    
    vcf_line = f"{chrom}\t{pos}\t{svid}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}\t{fmt}\t{sample_val}\n"
    return vcf_line

def main():
    parser = argparse.ArgumentParser(description="Rescue false-negatives in a PAV matrix using read-backed VCF verification.")
    parser.add_argument("-m", "--matrix", required=True, help="Path to input merged PAV matrix (uniq.pav.merge.txt)")
    parser.add_argument("-b", "--bam-dir", default=None, help="Directory containing your sample *_sorted.bam files")
    parser.add_argument("-t", "--sample-table", default=None, help="Optional tab-delimited file: sample_name<TAB>bam_path (one per line)")
    parser.add_argument("-o", "--output", default="refined.pav.merge.txt", help="Output filepath for refined matrix")
    parser.add_argument("-l", "--log-changes", default=None, help="Optional output filepath to log rescued (0 -> 1) variants")
    parser.add_argument("-maq", dest="maq", type=int, default=1, help="Minimum mapping quality for reads")
    parser.add_argument("-homo_rate", dest="homo_rate", type=float, default=0.65, help="Homozygous ratio threshold")
    parser.add_argument("-ref_rate", dest="ref_rate", type=float, default=0.15, help="Reference (0/0) ratio threshold")
    parser.add_argument("-shift", dest="shift", type=int, default=30, help="Breakpoint shift range")
    parser.add_argument("-span", dest="span", type=int, default=50, help="Span parameter for heterozygous evidence")
    
    args = parser.parse_args()

    print(f"[*] Reading source matrix: '{args.matrix}'...")
    with open(args.matrix, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]

    if not lines:
        print("[-] Error: Empty matrix file.")
        return

    headers = lines[0]
    base_headers = headers[:10]
    sample_headers = headers[10:]
    
    print(f"[+] Total sample columns identified: {len(sample_headers)}")

    # --- Resolve BAM paths: sample table > auto-guess from bam-dir ---
    sample_table = {}
    if args.sample_table:
        print(f"[*] Loading sample→BAM table: '{args.sample_table}'...")
        sample_table = load_sample_table(args.sample_table)
        print(f"[+] {len(sample_table)} entries loaded from table")

    bam_handles = {}
    missing_samples = []

    for sample in sample_headers:
        # 1) explicit table entry
        if sample in sample_table:
            bam_path = sample_table[sample]
        # 2) auto-guess from bam-dir (if provided)
        elif args.bam_dir:
            bam_path = match_bam_file(sample, args.bam_dir)
        else:
            missing_samples.append(sample)
            continue

        if bam_path and os.path.isfile(bam_path):
            try:
                bam_handles[sample] = pysam.AlignmentFile(bam_path, "rb")
            except Exception as e:
                print(f"[-] Warning: Error opening '{bam_path}': {e}. Skipping sample.")
                missing_samples.append(sample)
        else:
            missing_samples.append(sample)

    if missing_samples:
        print(f"[-] Missing/Unmatched BAM files for {len(missing_samples)} sample headers.")
        print(f"    Examples missing: {', '.join(missing_samples[:5])}...")
    
    print(f"[+] Active verification handles established for {len(bam_handles)} samples.")
    print("[*] Commencing read-backed verification across calls (INS > 80bp only)...")

    rescued_count = 0
    refined_rows = []
    log_entries = []

    # Process variant rows with progress visibility
    for data in tqdm(lines[1:], desc="Processing structural profiles"):
        row_dict = dict(zip(headers, data))
        updated_genotypes = []
        
        # Extract features for filtering constraints
        sv_type = row_dict["SVType"].strip().upper()
        try:
            sv_size = abs(int(row_dict["SVlen"]))
        except ValueError:
            sv_size = 0
            
        # Hard check: Only evaluate row if it's an insertion variant larger than 100bp
        is_target_ins = (sv_type == "INS" and sv_size > 80)
        
        for sample in sample_headers:
            current_pav_val = row_dict[sample].strip()
            
            if is_target_ins and current_pav_val == "0" and sample in bam_handles:
                pseudo_line = generate_pseudo_vcf_line(row_dict, sample_val="0/0")
                
                try:
                    res = compute_new_gt(
                        vcf_line=pseudo_line, 
                        opened_sam=bam_handles[sample], 
                        name=sample.split(".")[0], 
                        min_maq=args.maq, 
                        homo_rate=args.homo_rate, 
                        ref_rate=args.ref_rate, 
                        shift=args.shift, 
                        span=args.span
                    )
                    
                    if res and res["new_gt"] in ("1/1", "0/1"):
                        updated_genotypes.append("1")
                        rescued_count += 1
                        
                        # Cache the modification details for the log file
                        log_entries.append(
                            f"{row_dict['#Target_name']}\t{row_dict['Target_start']}\t"
                            f"{row_dict['SVID']}\t{row_dict['SVType']}\t{row_dict['SVlen']}\t{sample}\t{res['new_gt']}"
                        )
                    else:
                        updated_genotypes.append("0")
                except Exception:
                    updated_genotypes.append("0")
            else:
                # Keep original value if it's a DEL, a small INS (<=100), or already marked 1
                updated_genotypes.append(current_pav_val)
                
        new_row = data[:10] + updated_genotypes
        refined_rows.append(new_row)

    for name, handle in bam_handles.items():
        handle.close()

    print(f"[+] Matrix verification completed. Total large INS variants rescued to '1': {rescued_count}")
    
    # Save modification log if parameter flag is turned on
    if args.log_changes and log_entries:
        print(f"[*] Saving rescue modification history logs to: '{args.log_changes}'...")
        with open(args.log_changes, 'w') as log_out:
            log_out.write("#CHROM\tPOS\tSVID\tSVTYPE\tSVLEN\tSAMPLE_HEADER\tRESCUED_GT\n")
            for entry in log_entries:
                log_out.write(entry + "\n")

    print(f"[*] Writing updated matrix out to: '{args.output}'...")
    with open(args.output, 'w') as out:
        out.write("\t".join(headers) + "\n")
        for r in refined_rows:
            out.write("\t".join(r) + "\n")

    print("[+] Refined PAV Matrix processing completed successfully.")

if __name__ == "__main__":
    main()
