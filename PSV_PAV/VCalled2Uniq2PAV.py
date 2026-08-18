#!/usr/bin/env python3
"""
Genomic Record Formatter and SVID-Based PAV Matrix Generator
Author: lgb & DeepSeek
Year: 2026

Description:
  Processes variant files, groups duplicate lines based STRICTLY on the SVID 
  column, and creates a wide-format Presence-Absence Variant (PAV) matrix.
"""

import os
import re
import glob
import sys
import argparse

def extract_sample_name(file_path, chr_pattern=r'_[A-Za-z0-9]+[A-Za-z0-9.]*\.record\.txt|_Clustered_Record\.txt',
                        strip_prefix='0_tmp_', ref_name=None):
    """Extract sample name from an output file path.

    Splits on the chromosome/cluster boundary to separate sample+ref from
    chr/cluster+suffix, then strips optional prefix and reference name.

    Parameters
    ----------
    file_path : str
        Path like 'stem_Lsat_v11_1.record.txt' or
        'LA1589_ont.sorted_Clustered_Record.txt'
    chr_pattern : str
        Regex marking the chromosome/cluster boundary.
        Default matches alphanumeric chr names followed by .record.txt:
          _1.record.txt            (bare number)
          _01.record.txt           (zero-padded)
          _chr1.record.txt         (chr + digit)
          _chrI.record.txt         (chr + Roman numeral)
          _SL4.0ch00.record.txt    (alphanumeric with dots)
          _Clustered_Record.txt    (pre-clustered file)
    strip_prefix : str or None
        Prefix to remove from the start (e.g. '0_tmp_').  Pass None to skip.
    ref_name : str or None
        Reference name to strip from the end of the sample+ref portion.

    Returns
    -------
    str  - the extracted sample name.
    """
    base = os.path.basename(file_path)
    if strip_prefix and base.startswith(strip_prefix):
        base = base[len(strip_prefix):]
    parts = re.split(chr_pattern, base)
    sample_ref = parts[0]
    if ref_name and sample_ref.endswith(f'_{ref_name}'):
        sample_ref = sample_ref[:-(len(ref_name) + 1)]
    return sample_ref

def _parse_record_file(file_path, sample_name, variant_registry, is_clustered=False):
    """Parse one record file and register variants by SVID.

    Handles two column layouts:
      Clustered:  #target start end svlen SVID svtype seq maq cluster_size sv_rate
      Regular:    target query start end svlen maq SVID svtype seq

    Returns (record_count, sv_type_counts) for this file.
    """
    record_count = 0
    sv_type_counts = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')

            if is_clustered:
                if len(parts) < 10:
                    continue
                sv_id    = parts[4]
                sv_type  = parts[5]
                seq      = parts[6]
                maq      = parts[7]
                sv_start = parts[1]
                sv_end   = parts[2]
                sv_len   = parts[3]
            else:
                if len(parts) < 8:
                    continue
                sv_id    = parts[6]
                sv_type  = parts[7]
                seq      = parts[8] if len(parts) > 8 else "*"
                maq      = parts[5]
                sv_start = parts[2]
                sv_end   = parts[3]
                sv_len   = parts[4]

            if sv_id not in variant_registry:
                variant_registry[sv_id] = {
                    "target_name":  parts[0],
                    "target_start": sv_start,
                    "target_end":   sv_end,
                    "sv_len":       sv_len,
                    "sv_type":      sv_type,
                    "seq":          seq,
                    "maq":          maq,
                    "count":        0,
                    "samples":      set()
                }

            variant_registry[sv_id]["count"] += 1
            variant_registry[sv_id]["samples"].add(sample_name)
            record_count += 1
            sv_type_counts[sv_type] = sv_type_counts.get(sv_type, 0) + 1

    return record_count, sv_type_counts


def parse_and_cluster_variants(input_dir):
    """Scan input directory, preferring *_Clustered_Record.txt over per-chr files.

    For each sample, if a Clustered_Record.txt exists the per-chromosome
    *record.txt and *suppAlign files for that sample are skipped entirely.

    Returns (variant_registry, sorted_samples, total_records_before, sv_type_counts_before).
    """
    all_samples = set()
    variant_registry = {}
    clustered_samples = set()
    total_records_before = 0
    sv_type_counts_before = {}

    def _merge_counts(src, dst):
        for k, v in src.items():
            dst[k] = dst.get(k, 0) + v

    # --- Pass 1: Clustered_Record.txt (preferred) ---
    clustered_files = glob.glob(os.path.join(input_dir, "*_Clustered_Record.txt"))
    if clustered_files:
        print(f"[*] Found {len(clustered_files)} Clustered_Record.txt file(s):")
    for file_path in clustered_files:
        if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
            continue
        sample_name = extract_sample_name(file_path,
                                          chr_pattern=r'_Clustered_Record\.txt')
        print(f"    + {sample_name}  ({os.path.basename(file_path)})")
        clustered_samples.add(sample_name)
        all_samples.add(sample_name)
        n, counts = _parse_record_file(file_path, sample_name, variant_registry,
                                       is_clustered=True)
        total_records_before += n
        _merge_counts(counts, sv_type_counts_before)

    # --- Pass 2: per-chromosome *record.txt + *suppAlign ---
    skipped_count = 0
    for file_path in (glob.glob(os.path.join(input_dir, "*record*txt"))
                      + glob.glob(os.path.join(input_dir, "*suppAlign"))):
        if '_Clustered_Record.txt' in file_path:
            continue  # already done in pass 1
        if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
            continue
        sample_name = extract_sample_name(file_path)
        if sample_name in clustered_samples:
            skipped_count += 1
            continue  # skip - Clustered file takes precedence
        all_samples.add(sample_name)
        n, counts = _parse_record_file(file_path, sample_name, variant_registry,
                                       is_clustered=False)
        total_records_before += n
        _merge_counts(counts, sv_type_counts_before)

    if skipped_count:
        print(f"[*] Skipped {skipped_count} per-chromosome file(s) in favor of Clustered versions")

    if not all_samples:
        print(f"[-] Warning: No matching files found in '{input_dir}'",
              file=sys.stderr)

    return variant_registry, sorted(all_samples), total_records_before, sv_type_counts_before

def main():
    parser = argparse.ArgumentParser(description="Convert genomic pipelines into an SVID-deduplicated wide PAV matrix.")
    parser.add_argument("-i", "--input_dir", default="cr_bams_out", help="Directory containing variant files")
    parser.add_argument("-o", "--output", default="pav_matrix_output.txt", help="Output matrix filepath")
    args = parser.parse_args()

    print(f"[*] Scanning for data files in '{args.input_dir}'...")
    
    # Process files using SVID clustering logic
    variant_registry, sorted_samples, total_before, sv_before = parse_and_cluster_variants(args.input_dir)

    if not sorted_samples:
        print("[-] Process ended. No valid samples or variants found.")
        return

    print(f"[+] Identified {len(sorted_samples)} unique sample columns: {', '.join(sorted_samples)}")

    # --- Summary: before/after dedup & SV-type breakdown ---
    sv_after = {}
    for vdata in variant_registry.values():
        svt = vdata["sv_type"]
        sv_after[svt] = sv_after.get(svt, 0) + 1

    all_types = sorted(set(list(sv_before.keys()) + list(sv_after.keys())))

    print(f"\n{'='*80}")
    print(f"  Total records parsed   (before dedup): {total_before:>8}")
    print(f"  Unique SVIDs           (after  dedup): {len(variant_registry):>8}")
    print(f"  Reduction: {total_before - len(variant_registry)} duplicates removed "
          f"({(1 - len(variant_registry)/max(total_before,1))*100:.1f}%)")
    print(f"\n  {'SV type':<10s} {'Before':>10s} {'After':>10s} {'Retained':>10s} {'Redundancy':>12s}")
    print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
    for svt in all_types:
        before = sv_before.get(svt, 0)
        after = sv_after.get(svt, 0)
        pct = (after / before * 100) if before else 0
        ratio = f"{before/after:.1f}x" if after else "-"
        print(f"  {svt:<10s} {before:>10,} {after:>10,} {pct:>9.1f}% {ratio:>12}")
    print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
    print(f"  {'TOTAL':<10s} {total_before:>10,} {len(variant_registry):>10,} "
          f"{len(variant_registry)/max(total_before,1)*100:>9.1f}% "
          f"{total_before/max(len(variant_registry),1):.1f}x".rjust(12))
    print(f"{'='*80}\n")

    print(f"[+] Generating PAV matrix across {len(variant_registry)} unique SVIDs...")

    with open(args.output, 'w') as out:
        # Define the custom headers requested
        base_headers = [
            "#Target_name", "Target_start", "Target_end", "SVlen", 
            "SVID", "SVType", "seq", "maq", 
            "cluster_size_prevalent", "sv_rate_prevalent"
        ]
        
        # Append dynamic sample names onto the header row
        full_header = "\t".join(base_headers + sorted_samples)
        out.write(full_header + "\n")
        
        # Output sorted by SVID coordinates or alphanumeric order
        for sv_id in sorted(variant_registry.keys()):
            vdata = variant_registry[sv_id]
            
            # Formulate individual binary genotypes (1 or 0) for each sample column
            genotypes = []
            for sample in sorted_samples:
                if sample in vdata["samples"]:
                    genotypes.append("1")
                else:
                    genotypes.append("0")
            genotype_string = "\t".join(genotypes)
            
            # Reconstruct columns exactly to match header definitions
            final_row = (
                f"{vdata['target_name']}\t{vdata['target_start']}\t{vdata['target_end']}\t{vdata['sv_len']}\t"
                f"{sv_id}\t{vdata['sv_type']}\t{vdata['seq']}\t{vdata['maq']}\t"
                f"{vdata['count']}\t{vdata['count']}\t{genotype_string}"
            )
            out.write(final_row + "\n")

    print(f"[+] PAV matrix successfully saved to: '{args.output}'")

if __name__ == "__main__":
    main()
