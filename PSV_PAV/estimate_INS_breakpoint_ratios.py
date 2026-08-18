#!/usr/bin/env python3
"""
estimate_INS_breakpoint_ratios.py
Estimate breakpoint ratios distribution for INS SV calls marked as "1" in a PAV matrix.

Description:
  Reads a collapsed PAV matrix, locates matching sample BAM files via a bam.path file,
  and for each INS variant (>100 bp) where the PAV call is "1", computes the
  breakpoint ratio (clipped/breakpoint reads / total mapped reads) from the BAM.

  Collects all ratios and outputs:
    1. A raw data file with per-variant per-sample ratios
    2. A summary statistics file with distribution metrics
    3. A histogram-like bin count file

Unlike refine_INS_PAV.py which rescues 0→1 calls, this script characterizes the
read-backed evidence strength for already-called INS variants.

Usage:
  python estimate_INS_breakpoint_ratios.py \
      -m 2.PSV.uniq.SV50.pav.rename \
      -p ../bam.path \
      -o ins_breakpoint_ratios \
      --max-variants 10000
"""

import os
import sys
import argparse
import random
import re
from collections import defaultdict
from tqdm import tqdm
import pysam

# Import the core validation function
try:
    from srGT_refiner_for_asm_V1 import compute_new_gt
except ImportError:
    print("[-] Error: Could not import 'compute_new_gt' from srGT_refiner_for_asm_V1.py.", file=sys.stderr)
    print("[*] Please ensure this script is run from the directory containing your module.", file=sys.stderr)
    sys.exit(1)


def load_bam_paths(bam_path_file):
    """
    Load BAM file paths from a bam.path file.
    Format: <sample_id>\t<relative/path/to/bam>

    Paths are resolved relative to the bam.path file's directory.
    """
    bam_map = {}
    bam_dir = os.path.dirname(os.path.abspath(bam_path_file))

    with open(bam_path_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                sample_id = parts[0]
                bam_rel_path = parts[1]
                bam_abs_path = os.path.normpath(os.path.join(bam_dir, bam_rel_path))
                bam_map[sample_id] = bam_abs_path

    return bam_map


def generate_pseudo_vcf_line(row_dict, sample_val="0/0"):
    """
    Formulates a virtual VCF data string to feed into compute_new_gt().
    """
    chrom = row_dict["#Target_name"]
    pos = row_dict["Target_start"]
    svid = row_dict["SVID"]
    sv_type = row_dict["SVType"]
    end = row_dict["Target_end"]
    size = row_dict["SVlen"]

    ref = "N"
    alt = f"<{sv_type}>"
    qual = "60"
    filter_status = "PASS"
    info = f"SVTYPE={sv_type};END={end};SVLEN={size}"
    fmt = "GT"

    vcf_line = f"{chrom}\t{pos}\t{svid}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}\t{fmt}\t{sample_val}\n"
    return vcf_line


def parse_ins_rate(sv_support_str):
    """
    Parse the INS_rate value from the sv_support string returned by compute_new_gt.
    Format: "INS_rate=0.123;INS"
    Returns float ratio or None if parsing fails.
    """
    if not sv_support_str:
        return None
    match = re.search(r'INS_rate=([\d.]+)', sv_support_str)
    if match:
        return float(match.group(1))
    return None


def parse_total_reads(info_str):
    """
    Parse total_map_reads from the info string.
    Format: "total_map_reads=50;effective_spans=10"
    Returns int or None.
    """
    if not info_str:
        return None
    match = re.search(r'total_map_reads=(\d+)', info_str)
    if match:
        return int(match.group(1))
    return None


def compute_bins(ratios, num_bins=20):
    """
    Compute histogram bins for ratio distribution.
    """
    if not ratios:
        return []

    bin_width = 1.0 / num_bins
    bins = [0] * num_bins

    for r in ratios:
        idx = min(int(r / bin_width), num_bins - 1)
        bins[idx] += 1

    bin_labels = [f"{i*bin_width:.2f}-{(i+1)*bin_width:.2f}" for i in range(num_bins)]
    return list(zip(bin_labels, bins))


def main():
    parser = argparse.ArgumentParser(
        description="Estimate breakpoint ratios distribution for INS SV calls marked as '1' in a PAV matrix."
    )
    parser.add_argument("-m", "--matrix", required=True,
                        help="Path to input PAV matrix (e.g., 2.PSV.uniq.SV50.pav.rename)")
    parser.add_argument("-p", "--bam-path-file", required=True,
                        help="Path to bam.path file mapping samples to BAM files")
    parser.add_argument("-o", "--output-prefix", default="ins_breakpoint_ratios",
                        help="Output prefix for result files (default: ins_breakpoint_ratios)")

    # Filtering / sampling
    parser.add_argument("--max-variants", type=int, default=None,
                        help="Maximum number of INS variants to process (random sample; default: all)")
    parser.add_argument("--max-samples-per-variant", type=int, default=None,
                        help="Maximum samples per variant where PAV=1 (random; default: all)")
    parser.add_argument("--random-seed", type=int, default=42,
                        help="Random seed for sampling (default: 42)")
    parser.add_argument("--min-map-reads", type=int, default=5,
                        help="Minimum total mapped reads to include ratio (default: 5)")

    # GT calling parameters (passed through to compute_new_gt / insGT)
    parser.add_argument("-maq", dest="maq", type=int, default=1,
                        help="Minimum mapping quality for reads (default: 1)")
    parser.add_argument("-homo_rate", dest="homo_rate", type=float, default=0.65,
                        help="Homozygous ratio threshold (default: 0.65)")
    parser.add_argument("-ref_rate", dest="ref_rate", type=float, default=0.30,
                        help="Reference (0/0) ratio threshold (default: 0.30)")
    parser.add_argument("-shift", dest="shift", type=int, default=25,
                        help="Breakpoint shift range (default: 25)")
    parser.add_argument("-span", dest="span", type=int, default=50,
                        help="Span parameter for heterozygous evidence (default: 50)")

    args = parser.parse_args()

    random.seed(args.random_seed)

    # ------------------------------------------------------------------
    # 1. Load BAM path mappings
    # ------------------------------------------------------------------
    print(f"[*] Loading BAM path mappings from: '{args.bam_path_file}'...")
    bam_path_map = load_bam_paths(args.bam_path_file)
    print(f"[+] Loaded {len(bam_path_map)} sample-to-BAM mappings.")

    # ------------------------------------------------------------------
    # 2. Read PAV matrix
    # ------------------------------------------------------------------
    print(f"[*] Reading PAV matrix: '{args.matrix}'...")
    with open(args.matrix, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]

    if not lines:
        print("[-] Error: Empty matrix file.")
        return

    headers = lines[0]
    base_headers = headers[:10]
    sample_headers = headers[10:]

    print(f"[+] Total sample columns: {len(sample_headers)}")

    # ------------------------------------------------------------------
    # 3. Identify INS variants (>100bp) with at least one "1" call
    # Build a list of (row_idx, row_dict, samples_with_one)
    # ------------------------------------------------------------------
    print("[*] Scanning for INS variants (>100bp) with PAV=1 calls...")

    ins_variants = []

    for row_idx, data in enumerate(tqdm(lines[1:], desc="Scanning PAV matrix")):
        row_dict = dict(zip(headers, data))

        sv_type = row_dict["SVType"].strip().upper()
        try:
            sv_size = abs(int(row_dict["SVlen"]))
        except ValueError:
            sv_size = 0

        # Only INS > 100bp
        if sv_type != "INS" or sv_size <= 100:
            continue

        # Find samples where PAV = "1"
        samples_with_one = []
        for sample in sample_headers:
            if sample in bam_path_map and row_dict[sample].strip() == "1":
                samples_with_one.append(sample)

        if samples_with_one:
            ins_variants.append((row_dict, samples_with_one))

    print(f"[+] Found {len(ins_variants)} INS variants (>100bp) with >=1 PAV=1 call.")

    if not ins_variants:
        print("[-] No INS variants to process. Exiting.")
        return

    # ------------------------------------------------------------------
    # 4. Sample if requested
    # ------------------------------------------------------------------
    if args.max_variants and len(ins_variants) > args.max_variants:
        ins_variants = random.sample(ins_variants, args.max_variants)
        print(f"[*] Randomly sampled {args.max_variants} variants for processing.")

    # ------------------------------------------------------------------
    # 5. Open BAM handles (only for samples that appear in selected variants)
    # ------------------------------------------------------------------
    needed_samples = set()
    for row_dict, samples_with_one in ins_variants:
        for s in samples_with_one:
            needed_samples.add(s)

    print(f"[*] Opening BAM handles for {len(needed_samples)} needed samples...")

    bam_handles = {}
    missing_samples = []

    for sample in needed_samples:
        bam_path = bam_path_map.get(sample)
        if bam_path and os.path.isfile(bam_path):
            try:
                bam_handles[sample] = pysam.AlignmentFile(bam_path, "rb")
            except Exception as e:
                print(f"[-] Warning: Error opening '{bam_path}': {e}. Skipping sample {sample}.")
                missing_samples.append(sample)
        else:
            missing_samples.append(sample)

    if missing_samples:
        print(f"[-] Missing/Unmatched BAM files for {len(missing_samples)} samples: {', '.join(missing_samples[:10])}...")

    print(f"[+] Active BAM handles: {len(bam_handles)}")

    # ------------------------------------------------------------------
    # 6. Process: for each INS variant, for each sample where PAV=1,
    #    compute breakpoint ratio
    # ------------------------------------------------------------------
    print("[*] Computing breakpoint ratios for INS calls...")

    # Results storage
    all_records = []  # list of dicts with ratio data

    total_processed = 0
    total_skipped_low_reads = 0
    total_errors = 0

    for row_dict, samples_with_one in tqdm(ins_variants, desc="Processing INS variants"):
        # Optionally subsample the samples for this variant
        samples_to_process = samples_with_one
        if args.max_samples_per_variant and len(samples_to_process) > args.max_samples_per_variant:
            samples_to_process = random.sample(samples_to_process, args.max_samples_per_variant)

        pseudo_line = generate_pseudo_vcf_line(row_dict, sample_val="0/0")

        for sample in samples_to_process:
            if sample not in bam_handles:
                continue

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

                if res:
                    ins_ratio = parse_ins_rate(res.get("sv_support", ""))
                    total_reads = parse_total_reads(res.get("total_map_reads", ""))
                    new_gt = res.get("new_gt", "./.")

                    if ins_ratio is not None and total_reads is not None:
                        if total_reads >= args.min_map_reads:
                            all_records.append({
                                "chrom": row_dict["#Target_name"],
                                "start": row_dict["Target_start"],
                                "end": row_dict["Target_end"],
                                "svid": row_dict["SVID"],
                                "svlen": row_dict["SVlen"],
                                "sample": sample,
                                "ins_ratio": ins_ratio,
                                "total_map_reads": total_reads,
                                "computed_gt": new_gt
                            })
                            total_processed += 1
                        else:
                            total_skipped_low_reads += 1
                    else:
                        total_skipped_low_reads += 1
                else:
                    total_errors += 1
            except Exception as e:
                total_errors += 1

    # Close BAM handles
    for handle in bam_handles.values():
        handle.close()

    print(f"\n[+] Processing complete:")
    print(f"    Total breakpoint ratios collected: {total_processed}")
    print(f"    Skipped (low reads / no data):    {total_skipped_low_reads}")
    print(f"    Errors:                            {total_errors}")

    if not all_records:
        print("[-] No ratio data collected. Exiting.")
        return

    # ------------------------------------------------------------------
    # 7. Compute distribution statistics
    # ------------------------------------------------------------------
    ratios = [r["ins_ratio"] for r in all_records]
    total_reads_list = [r["total_map_reads"] for r in all_records]

    ratios_sorted = sorted(ratios)
    n = len(ratios_sorted)

    def percentile(data_sorted, p):
        """Compute percentile from sorted data."""
        if not data_sorted:
            return 0.0
        idx = int((p / 100.0) * (len(data_sorted) - 1))
        return data_sorted[idx]

    mean_ratio = sum(ratios) / n
    median_ratio = percentile(ratios_sorted, 50)

    # Compute bin distribution
    bins = compute_bins(ratios, num_bins=20)

    # Compute per-sample statistics
    sample_ratios = defaultdict(list)
    for r in all_records:
        sample_ratios[r["sample"]].append(r["ins_ratio"])

    # ------------------------------------------------------------------
    # 8. Write output files
    # ------------------------------------------------------------------

    # 8a. Raw data file
    raw_file = f"{args.output_prefix}.raw.tsv"
    print(f"[*] Writing raw ratio data to: '{raw_file}'...")
    with open(raw_file, 'w') as f:
        f.write("#CHROM\tSTART\tEND\tSVID\tSVLEN\tSAMPLE\tINS_RATIO\tTOTAL_MAP_READS\tCOMPUTED_GT\n")
        for rec in all_records:
            f.write(f"{rec['chrom']}\t{rec['start']}\t{rec['end']}\t{rec['svid']}\t"
                    f"{rec['svlen']}\t{rec['sample']}\t{rec['ins_ratio']}\t"
                    f"{rec['total_map_reads']}\t{rec['computed_gt']}\n")

    # 8b. Summary statistics file
    summary_file = f"{args.output_prefix}.summary.txt"
    print(f"[*] Writing summary statistics to: '{summary_file}'...")
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("INS Breakpoint Ratio Distribution Summary\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Source matrix:     {args.matrix}\n")
        f.write(f"Total INS variants processed: {len(ins_variants)}\n")
        f.write(f"Total ratio measurements:     {n}\n")
        f.write(f"Samples represented:          {len(sample_ratios)}\n\n")

        f.write("-" * 50 + "\n")
        f.write("Overall Distribution Statistics\n")
        f.write("-" * 50 + "\n")
        f.write(f"  Mean:                {mean_ratio:.4f}\n")
        f.write(f"  Median:              {median_ratio:.4f}\n")
        f.write(f"  Std Dev:             {sum((x-mean_ratio)**2 for x in ratios)/n:.4f}\n")
        f.write(f"  Min:                 {min(ratios):.4f}\n")
        f.write(f"  Max:                 {max(ratios):.4f}\n")
        f.write(f"  P5:                  {percentile(ratios_sorted, 5):.4f}\n")
        f.write(f"  P10:                 {percentile(ratios_sorted, 10):.4f}\n")
        f.write(f"  P25:                 {percentile(ratios_sorted, 25):.4f}\n")
        f.write(f"  P50 (Median):        {percentile(ratios_sorted, 50):.4f}\n")
        f.write(f"  P75:                 {percentile(ratios_sorted, 75):.4f}\n")
        f.write(f"  P90:                 {percentile(ratios_sorted, 90):.4f}\n")
        f.write(f"  P95:                 {percentile(ratios_sorted, 95):.4f}\n")
        f.write(f"  P99:                 {percentile(ratios_sorted, 99):.4f}\n\n")

        f.write("-" * 50 + "\n")
        f.write("Total Mapped Reads Statistics\n")
        f.write("-" * 50 + "\n")
        f.write(f"  Mean:                {sum(total_reads_list)/len(total_reads_list):.1f}\n")
        f.write(f"  Median:              {sorted(total_reads_list)[len(total_reads_list)//2]:.1f}\n")
        f.write(f"  Min:                 {min(total_reads_list)}\n")
        f.write(f"  Max:                 {max(total_reads_list)}\n\n")

        f.write("-" * 50 + "\n")
        f.write("Ratio Binned Distribution (20 bins)\n")
        f.write("-" * 50 + "\n")
        f.write(f"{'Bin Range':<16} {'Count':>8} {'Pct':>8}  Histogram\n")
        f.write("-" * 50 + "\n")
        for bin_label, count in bins:
            pct = (count / n) * 100
            bar = "#" * max(1, int(pct * 2))
            f.write(f"{bin_label:<16} {count:>8} {pct:>7.1f}%  {bar}\n")
        f.write("\n")

        # Per-sample summary (top 20 by count)
        f.write("-" * 50 + "\n")
        f.write("Per-Sample Mean Ratio (top 30 by count)\n")
        f.write("-" * 50 + "\n")
        f.write(f"{'Sample':<12} {'Count':>8} {'Mean_Ratio':>12} {'Median_Ratio':>14}\n")
        sample_stats = [(s, len(v), sum(v)/len(v), sorted(v)[len(v)//2])
                        for s, v in sample_ratios.items()]
        sample_stats.sort(key=lambda x: x[1], reverse=True)
        for s, cnt, mean_r, med_r in sample_stats[:30]:
            f.write(f"{s:<12} {cnt:>8} {mean_r:>12.4f} {med_r:>14.4f}\n")

    # 8c. GT breakdown
    gt_file = f"{args.output_prefix}.gt_breakdown.txt"
    print(f"[*] Writing GT breakdown to: '{gt_file}'...")
    gt_counts = defaultdict(int)
    for rec in all_records:
        gt_counts[rec["computed_gt"]] += 1

    with open(gt_file, 'w') as f:
        f.write("Computed Genotype Distribution for INS PAV=1 Calls\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'GT':<8} {'Count':>8} {'Pct':>8}\n")
        for gt in ["1/1", "0/1", "0/0", "./."]:
            cnt = gt_counts.get(gt, 0)
            pct = (cnt / n) * 100 if n > 0 else 0
            f.write(f"{gt:<8} {cnt:>8} {pct:>7.1f}%\n")
        # Any other GTs
        for gt, cnt in sorted(gt_counts.items()):
            if gt not in ["1/1", "0/1", "0/0", "./."]:
                pct = (cnt / n) * 100 if n > 0 else 0
                f.write(f"{gt:<8} {cnt:>8} {pct:>7.1f}%\n")

    print(f"\n[+] All output files written with prefix '{args.output_prefix}':")
    print(f"    Raw data:      {raw_file}")
    print(f"    Summary stats: {summary_file}")
    print(f"    GT breakdown:  {gt_file}")
    print("[+] Done.")


if __name__ == "__main__":
    main()
