"""Worker functions for chromosome-level SV detection (BAM + SAM).

Uses sv_engine for all SV detection logic.  Accepts either:
  - A BAM path (worker opens and fetches the chromosome)
  - A list of pre-parsed dict records (SAM mode, where fetch() is unavailable)
"""

import pysam
from sub_Signaling_Module import (
    process_record,
    call_segment_svs_cross_read,
    compute_depth,
    record_to_dict,
    _r,
)


def _process_records(records, chromosome, chrom_size, chromosome_list,
                     min_len, max_len, min_maq, dtype, msv, out_prefix):
    """Core worker logic — process a list of records for one chromosome.

    Records may be pysam AlignedSegment objects or picklable dicts.
    _r() handles both transparently.
    """
    sv_path    = f'{out_prefix}_{chromosome}.record.txt'
    cov_path   = f'{out_prefix}_{chromosome}.record.txt.cov'
    supp_path  = f'{out_prefix}_{chromosome}.record.txt.suppAlign'
    depth_path = f'{out_prefix}_{chromosome}.record.txt.depth'
    break_path = f'{out_prefix}_{chromosome}.record.txt.breakpoint'

    with open(sv_path, 'w') as sv_out, \
         open(cov_path, 'w') as cov_out, \
         open(supp_path, 'w') as supp_out, \
         open(depth_path, 'w') as depth_out, \
         open(break_path, 'w') as break_out:

        total_map_len = 0

        for rec in records:
            cigar_svs, cov_line, supp_svs, breakpoints = process_record(
                rec, min_len, max_len, min_maq, msv, dtype, chromosome_list,
                chrom_size)

            if cigar_svs:
                sv_out.writelines(cigar_svs)
            if cov_line:
                parts = cov_line.strip().split('\t')
                if len(parts) >= 5:
                    total_map_len += (int(parts[4]) - int(parts[3]))
                cov_out.write(cov_line)
            if supp_svs:
                for sv in supp_svs:
                    supp_out.write(sv + '\n')
            if breakpoints:
                break_out.writelines(breakpoints)

        depth = compute_depth(total_map_len, chrom_size)
        depth_out.write(
            f'{chromosome}\t{chrom_size}\t{total_map_len}\t{depth}\n')

    # --- Second pass: cross-read segment SV (cr mode only) ---
    if dtype == 'cr' and msv == 'yes':
        cross_svs = call_segment_svs_cross_read(
            records, min_len, max_len, chromosome_list)
        if cross_svs:
            with open(supp_path, 'a') as supp_out:
                for sv in cross_svs:
                    supp_out.write(sv + '\n')

    # Count lines for reporting
    def _count_lines(path):
        try:
            with open(path) as f:
                return sum(1 for _ in f)
        except FileNotFoundError:
            return 0

    sv_count   = _count_lines(sv_path)
    supp_count = _count_lines(supp_path)
    return chromosome, sv_count, supp_count


# ---------------------------------------------------------------------------
# BAM mode — worker opens its own BAM handle and fetches the chromosome
# ---------------------------------------------------------------------------

def process_chromosome_bam(chromosome, chrom_size, chromosome_list,
                           bamfile_path, min_len, max_len, min_maq,
                           out_prefix, dtype, msv):
    """Worker entry for BAM mode — opens the BAM and fetches records."""
    samfile = pysam.AlignmentFile(bamfile_path, 'rb')
    records = list(samfile.fetch(chromosome))
    samfile.close()
    return _process_records(records, chromosome, chrom_size, chromosome_list,
                            min_len, max_len, min_maq, dtype, msv, out_prefix)


# ---------------------------------------------------------------------------
# SAM mode — worker receives pre-parsed dict records
# ---------------------------------------------------------------------------

def process_chromosome_dicts(chromosome, chrom_size, chromosome_list,
                             dict_records, min_len, max_len, min_maq,
                             out_prefix, dtype, msv):
    """Worker entry for SAM mode — receives a list of picklable dict records."""
    return _process_records(dict_records, chromosome, chrom_size, chromosome_list,
                            min_len, max_len, min_maq, dtype, msv, out_prefix)
