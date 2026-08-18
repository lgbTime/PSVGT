#!/usr/bin/env python
"""
sv_engine.py — Consolidated SV detection engine (pysam + SAM dict compatible).

Two detection layers:
  Layer 1 — CIGAR-based DEL/INS from primary alignment CIGAR.
  Layer 2 — Segment-pair geometric analysis (DeBreak-inspired) for
            INS, DEL, DUP, INV, TRA from SA-tag or cross-read segments.

All functions are pure (no file I/O). They accept either pysam
AlignedSegment objects or picklable dicts (from record_to_dict or
sam_line_to_dict), accessed through the universal _r() accessor.
"""

import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Record converters (pysam <-> dict)
# ---------------------------------------------------------------------------

def record_to_dict(rec):
    """Convert pysam AlignedSegment -> picklable dict (pass-through if already dict)."""
    if isinstance(rec, dict):
        return rec
    return {
        'flag':    rec.flag,
        'mapq':    rec.mapping_quality,
        'rname':   rec.reference_name,
        'rstart':  rec.reference_start,
        'rend':    rec.reference_end,
        'qname':   rec.query_name,
        'qseq':    rec.query_sequence or None,
        'cigar':   rec.cigarstring,
        'cigars':  [(op, length) for op, length in rec.cigar] if rec.cigar else [],
        'qalen':   rec.query_alignment_length,
        'sa_tag':  rec.get_tag('SA') if rec.has_tag('SA') else None,
    }


def sam_line_to_dict(line):
    """Parse a SAM text line -> picklable dict (0-based coords, matches pysam).

    Returns None for headers or unparseable lines.
    """
    if line.startswith('@'):
        return None

    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None

    qname = fields[0]
    flag  = int(fields[1])
    rname = fields[2]
    mapq  = int(fields[4])
    cigar_str = fields[5]
    seq   = fields[9]

    # Parse CIGAR -> [(code_int, length), ...]
    CIGAR_CODE = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5,
                  'P': 6, '=': 7, 'X': 8}
    cnums   = [int(x) for x in re.findall(r'\d+', cigar_str)]
    ccodes  = re.findall(r'[A-Z]', cigar_str)
    cigars  = [(CIGAR_CODE.get(c, 0), length)
               for c, length in zip(ccodes, cnums)]

    # 0-based rstart (SAM is 1-based)
    rstart = int(fields[3]) - 1

    # qalen = aligned query length (M,I,=,X only — excl S/H, matches pysam query_alignment_length)
    qalen = sum(l for op, l in cigars if op in (0, 1, 7, 8))

    # rend = 0-based exclusive reference end
    rend = rstart + sum(l for op, l in cigars if op in (0, 2, 3, 7, 8))

    # SA tag
    sa_tag = None
    for f in fields[11:]:
        if f.startswith('SA:Z:'):
            sa_tag = f[5:]
            break

    return {
        'flag':    flag,
        'mapq':    mapq,
        'rname':   rname,
        'rstart':  rstart,
        'rend':    rend,
        'qname':   qname,
        'qseq':    seq if seq != '*' else None,
        'cigar':   cigar_str,
        'cigars':  cigars,
        'qalen':   qalen,
        'sa_tag':  sa_tag,
    }


# ---------------------------------------------------------------------------
# Universal record accessor
# ---------------------------------------------------------------------------

def _r(rec, key):
    """Get a field from a pysam AlignedSegment or a dict record."""
    if isinstance(rec, dict):
        return rec[key]

    # pysam AlignedSegment
    attr_map = {
        'flag':   'flag',
        'mapq':   'mapping_quality',
        'rname':  'reference_name',
        'rstart': 'reference_start',
        'rend':   'reference_end',
        'qname':  'query_name',
        'qseq':   'query_sequence',
        'cigar':  'cigarstring',
        'cigars': 'cigar',
        'qalen':  'query_alignment_length',
        'sa_tag': None,
    }
    attr = attr_map.get(key)
    if attr is None:
        if key == 'sa_tag':
            return rec.get_tag('SA') if rec.has_tag('SA') else None
        raise KeyError(key)
    return getattr(rec, attr)


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

def parse_cigar2clipinfo(cigar_str_or_tuples):
    """Return [left_clip, mapped_length, right_clip] from a CIGAR.

    Handles both string ("100M50S") and tuple-list (((0,100),(4,50))).
    mapped_length = M + I + = + X only (excl S/H/N/P, matches pysam qalen).
    """
    if isinstance(cigar_str_or_tuples, str):
        lengths = [int(n) for n in re.findall(r'\d+', cigar_str_or_tuples)]
        ops     = re.findall(r'[A-Z]', cigar_str_or_tuples)
        tuples  = list(zip(ops, lengths))
    else:
        code_map = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 7: '=', 8: 'X'}
        tuples = [(code_map.get(op, 'M'), length) for op, length in cigar_str_or_tuples]

    left_clip  = 0
    right_clip = 0
    mapped_len = sum(l for op, l in tuples if op in ('M', '=', 'X', 'I'))

    if tuples and tuples[0][0] in ('S', 'H'):
        left_clip = tuples[0][1]
    if tuples and tuples[-1][0] in ('S', 'H'):
        right_clip = tuples[-1][1]

    return [left_clip, mapped_len, right_clip]


# ---------------------------------------------------------------------------
# Layer 1 — CIGAR-based DEL / INS
# ---------------------------------------------------------------------------

def call_cigar_del_ins(record, min_len, max_len):
    """Extract DEL and INS from CIGAR of a primary-like alignment."""
    if _r(record, 'flag') == 4 or _r(record, 'mapq') < 1:
        return []

    results = []
    target_chr = _r(record, 'rname')
    query_chr = _r(record, 'qname')
    query_seq = _r(record, 'qseq') or ''
    cigar     = _r(record, 'cigar')
    mapq      = _r(record, 'mapq')
    ref       = _r(record, 'rstart')
    query_pos = 0

    cigar_numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    cigar_codes   = re.findall(r'[A-Z]', cigar)

    for code, length in zip(cigar_codes, cigar_numbers):
        if code == 'M':
            ref += length
            query_pos += length
        elif code == 'D':
            if min_len <= length <= max_len:
                results.append(
                    f'{target_chr}\t{query_chr}\t{ref}\t{ref + length - 1}'
                    f'\t{length}\t{mapq}'
                    f'\t{target_chr}:{ref}-{ref + length - 1}_DEL={length}'
                    f'\tDEL\t"*"\n')
            ref += length
        elif code == 'I':
            if min_len <= length <= max_len:
                ins_seq = (query_seq[query_pos:query_pos + length]
                           if query_seq else '*')
                results.append(
                    f'{target_chr}\t{query_chr}\t{ref}\t{ref + 1}'
                    f'\t{length}\t{mapq}'
                    f'\t{target_chr}:{ref}-{ref + 1}_INS={length}'
                    f'\tINS\t{ins_seq}\n')
            query_pos += length
        elif code in ('S', 'H'):
            query_pos += length

    return results


# ---------------------------------------------------------------------------
# Layer 2 — Segment helpers
# ---------------------------------------------------------------------------

def _build_segment(record):
    """Convert a record (pysam or dict) -> internal segment list.

    Format: [readname, flag, chrom, start, end,
             [leftclip, mapped_len, rightclip], mapq, sequence]
    """
    clipinfo = [0, 0, 0]
    flag = _r(record, 'flag')
    cigars = _r(record, 'cigars')
    if flag != 4 and cigars:
        if cigars[0][0] in (4, 5):
            clipinfo[0] = cigars[0][1]
        if cigars[-1][0] in (4, 5):
            clipinfo[2] = cigars[-1][1]
        clipinfo[1] = _r(record, 'qalen')
    return [
        _r(record, 'qname'),
        flag,
        _r(record, 'rname'),
        _r(record, 'rstart'),
        _r(record, 'rend'),
        clipinfo,
        _r(record, 'mapq'),
        _r(record, 'qseq') or '*',
    ]


def _strand_from_flag(flag):
    return (flag % 32) > 15


def _make_svcall(chrom, readname, start, end, size, mapq, svid, sv_type,
                 sequence='*'):
    return (f'{chrom}\t{readname}\t{start}\t{end}\t{size}\t{mapq}'
            f'\t{svid}\t{sv_type}\t{sequence}')


# ---------------------------------------------------------------------------
# Consolidated segment-pair SV detector (pairwise, symmetric)
# ---------------------------------------------------------------------------

def call_segment_svs(seg_list, min_size, max_size, chromosome_list,
                     minimaq, chrom_size=None, dtype='cr'):
    """Detect INS/DEL/DUP/INV/TRA from alignment segments of one query.

    Two detection strategies depending on data type:

      Assembly mode (dtype='cr'/'sr'):
        Every segment is an equally valid chunk of the same contig.
        Tests all unique unordered pairs — no primary/supplementary hierarchy.

      Long-read mode (dtype='pb'/'ont'/'hifi'):
        The primary alignment (flag ≤ 16) is the true read mapping; SA-tag
        segments are alternative/chimeric alignments.  Only primary↔supp
        pairs are tested to avoid comparing two alternative mappings.

    Args:
        seg_list: list of [readname, flag, chrom, start, end,
                   [leftclip, maplen, rightclip], mapq, sequence]
        min_size, max_size: SV size bounds
        chromosome_list: ordered chromosome names (for TRA naming)
        minimaq: minimum MAPQ
        chrom_size: optional int — chromosome size (INV filter)

    Returns list of SV record strings.
    """
    n = len(seg_list)
    if n <= 1:
        return []

    # Filter out segments with very short mapped length
    segs = [s for s in seg_list if s[5][1] >= 250]
    n = len(segs)
    if n <= 1:
        return []

    readname = segs[0][0]
    is_asm = dtype in ('cr', 'sr')

    if is_asm:
        return _call_segment_svs_all_pairs(
            segs, n, readname, min_size, max_size, chromosome_list,
            minimaq, chrom_size)
    else:
        return _call_segment_svs_primary_anchored(
            segs, n, readname, min_size, max_size, chromosome_list,
            minimaq, chrom_size)


# ---------------------------------------------------------------------------
# Assembly-mode: all unique unordered pairs (no primary concept)
# ---------------------------------------------------------------------------

def _call_segment_svs_all_pairs(segs, n, readname, min_size, max_size,
                                chromosome_list, minimaq, chrom_size):
    """Test every unordered pair of segments. For genome-vs-genome assembly."""
    results = []

    for i in range(n):
        seg_a = segs[i]
        chrom_a   = seg_a[2]
        strand_a  = _strand_from_flag(seg_a[1])

        for j in range(i + 1, n):
            seg_b = segs[j]
            chrom_b  = seg_b[2]
            strand_b = _strand_from_flag(seg_b[1])

            maq = (int(seg_a[6]) + int(seg_b[6])) // 2
            if maq < minimaq:
                continue

            if chrom_a != chrom_b:
                _detect_tra(results, seg_a, seg_b, readname, maq, chromosome_list)
            elif strand_a != strand_b:
                _detect_inv(results, seg_a, seg_b, readname, maq,
                            min_size, max_size, chrom_size)
            else:
                _detect_samedir(results, seg_a, seg_b, readname, maq,
                                min_size, max_size)

    return results


# ---------------------------------------------------------------------------
# Long-read mode: primary-anchored (primary ↔ supplementary only)
# ---------------------------------------------------------------------------

def _call_segment_svs_primary_anchored(segs, n, readname, min_size, max_size,
                                       chromosome_list, minimaq, chrom_size):
    """Only test pairs where one segment is the true primary (flag ≤ 16).

    For long-read data, the primary is the read's real mapping; SA-tag
    segments are alternative/chimeric alignments.  Comparing two alternative
    mappings (supp↔supp) would create false positives.
    """
    # Find the true primary
    primary = None
    prim_idx = None
    for i, s in enumerate(segs):
        flag = int(s[1])
        if (flag & 0x800) == 0 and (flag & 0x100) == 0 and flag <= 16:
            primary = s
            prim_idx = i
            break

    if primary is None:
        return []   # no valid primary to anchor from

    pri_chrom  = primary[2]
    pri_strand = _strand_from_flag(primary[1])
    results = []

    for j in range(n):
        if j == prim_idx:
            continue
        supp = segs[j]
        supp_chrom  = supp[2]
        supp_strand = _strand_from_flag(supp[1])

        maq = (int(primary[6]) + int(supp[6])) // 2
        if maq < minimaq:
            continue

        if supp_chrom != pri_chrom:
            _detect_tra(results, primary, supp, readname, maq, chromosome_list)
        elif supp_strand != pri_strand:
            _detect_inv(results, primary, supp, readname, maq,
                        min_size, max_size, chrom_size)
        else:
            _detect_samedir(results, primary, supp, readname, maq,
                            min_size, max_size)

    return results


# ---------------------------------------------------------------------------
# Pair-level detectors (shared by both modes)
# ---------------------------------------------------------------------------

def _detect_samedir(results, seg_a, seg_b, readname, maq, min_size, max_size):
    """INS / DEL / DUP from a same-strand pair. Sorts by reference position."""
    left, right = (seg_a, seg_b) if seg_a[3] < seg_b[3] else (seg_b, seg_a)
    chrom = left[2]
    sh1, len1, sh2 = left[5]
    sh3, len2, sh4 = right[5]

    # --- INS: close on reference, gap in read ---
    if abs(right[3] - left[4]) <= 300:
        overlapmap = right[3] - left[4]
        ins_size = sh3 - len1 - sh1 - overlapmap
        if min_size <= ins_size <= max_size:
            sv_start = min(right[3], left[4])
            svid = f'{chrom}:{sv_start}-{sv_start + 1}_INSLEN={ins_size}'
            results.append(_make_svcall(
                chrom, readname, sv_start, sv_start + 1,
                ins_size, maq, svid, 'INS'))

    # --- DEL: reference gap bridged by read overlap ---
    overlapmap = sh1 + len1 - sh3
    if -200 < overlapmap < 1500:
        del_size = right[3] - left[4] + overlapmap
        if min_size <= del_size <= max_size:
            sv_start = max(0, left[4] - max(0, overlapmap))
            sv_end   = sv_start + del_size - 1
            svid = f'{chrom}:{sv_start}-{sv_end}_DELLEN={del_size}'
            results.append(_make_svcall(
                chrom, readname, sv_start, sv_end, del_size,
                maq, svid, 'DEL'))

    # --- DUP lap1: left 3' end overlaps right 5' start on query ---
    lap1 = sh1 + len1 - sh3
    if -200 < lap1 < 500 and (left[4] - right[3]) >= max(50, lap1):
        dup_size = left[4] - right[3] - max(lap1, 0)
        if min_size <= dup_size <= max_size:
            sv_start, sv_end = right[3], right[3] + dup_size
            svid = f'{chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
            results.append(_make_svcall(
                chrom, readname, sv_start, sv_end, dup_size,
                maq, svid, 'DUP'))

    # --- DUP lap2: right 3' end overlaps left 5' start ---
    lap2 = sh3 + len2 - sh1
    if -200 < lap2 < 500 and (right[4] - left[3]) >= max(50, lap2):
        dup_size = right[4] - left[3] - max(lap2, 0)
        if min_size <= dup_size <= max_size:
            sv_start, sv_end = left[3], left[3] + dup_size
            svid = f'{chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
            results.append(_make_svcall(
                chrom, readname, sv_start, sv_end, dup_size,
                maq, svid, 'DUP'))


def _detect_inv(results, seg_a, seg_b, readname, maq,
                min_size, max_size, chrom_size):
    """INV from an opposite-strand pair."""
    chrom = seg_a[2]  # same chrom for both

    # Sort by reference position
    left, right = (seg_a, seg_b) if seg_a[3] < seg_b[3] else (seg_b, seg_a)

    # Query-position consistency: when reference order flips,
    # query positions should stay consistent (both ends shift together)
    if not ((right[3] > left[3] and (right[4] - left[4]) > -200) or
            (right[3] < left[3] and (left[4] - right[4]) > -200)):
        return

    sh1, len1, sh2 = left[5]
    sh3, len2, sh4 = right[5]

    # INV lap1
    lap1 = sh3 + len2 - sh2
    if -200 < lap1 < 500 and (right[4] - left[4]) > max(100, lap1):
        inv_size = right[4] - left[4] - lap1
        if min_size <= inv_size <= max_size:
            if chrom_size is None or inv_size < chrom_size / 4:
                sv_start, sv_end = left[4], left[4] + inv_size
                svid = f'{chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                results.append(_make_svcall(
                    chrom, readname, sv_start, sv_end, inv_size,
                    maq, svid, 'INV'))

    # INV lap2
    lap2 = sh4 + len2 - sh1
    if -200 < lap2 < 500 and (right[3] - left[3]) >= max(100, lap2):
        inv_size = right[3] - left[3] - lap2
        if min_size <= inv_size <= max_size:
            if chrom_size is None or inv_size < chrom_size / 4:
                sv_start, sv_end = left[3], left[3] + inv_size
                svid = f'{chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                results.append(_make_svcall(
                    chrom, readname, sv_start, sv_end, inv_size,
                    maq, svid, 'INV'))


def _detect_tra(results, seg_a, seg_b, readname, maq, chromosome_list):
    """TRA from a different-chromosome pair."""
    chrom_a, chrom_b = seg_a[2], seg_b[2]
    sh1, len1, sh2 = seg_a[5]
    sh3, len2, sh4 = seg_b[5]

    # Find query-adjacent edges between the two segments
    bp_a = ''
    if abs(sh1 - (sh3 + len2)) <= 500:
        bp_a = seg_a[3]                        # A.start ≈ B.end on query
    elif abs((sh1 + len1) - sh3) <= 500:
        bp_a = seg_a[4]                        # A.end ≈ B.start on query

    bp_b = ''
    if abs(sh3 - (sh1 + len1)) <= 500:
        bp_b = seg_b[3]                        # B.start ≈ A.end on query
    elif abs((sh3 + len2) - sh1) <= 500:
        bp_b = seg_b[4]                        # B.end ≈ A.start on query

    if bp_a == '' or bp_b == '':
        return

    if chromosome_list.index(chrom_a) < chromosome_list.index(chrom_b):
        svid = f'{chrom_b}:{bp_b}_{chrom_a}:{bp_a}'
        results.append(_make_svcall(
            chrom_a, readname, bp_a, bp_b, 0, maq, svid, 'TRA'))
    else:
        svid = f'{chrom_a}:{bp_a}_{chrom_b}:{bp_b}'
        results.append(_make_svcall(
            chrom_b, readname, bp_b, bp_a, 0, maq, svid, 'TRA'))


# ---------------------------------------------------------------------------
# Cross-read segment SV caller (cr mode, second pass)
# ---------------------------------------------------------------------------

def call_segment_svs_cross_read(records, min_size, max_size, chromosome_list):
    """Group all records by query_name and run segment-pair analysis per query.

    Accepts pysam records or picklable dicts. Uses minimaq=60 for cross-read
    (consistent with DeBreak-inspired filtering for assembly mode).
    """
    supp_dict = defaultdict(list)
    for rec in records:
        if _r(rec, 'flag') == 4:
            continue
        supp_dict[_r(rec, 'qname')].append(_build_segment(rec))

    all_svs = []
    for seg_list in supp_dict.values():
        if len(seg_list) < 2:
            continue
        all_svs.extend(call_segment_svs(
            seg_list, min_size, max_size, chromosome_list, minimaq=60))
    return all_svs


# ---------------------------------------------------------------------------
# Per-record processor
# ---------------------------------------------------------------------------

def process_record(record, min_len, max_len, minimaq, msv, dtype,
                   chromosome_list, chrom_size=None):
    """Process one record for CIGAR SVs, coverage, breakpoints, and SA-based SVs.

    dtype in {'cr','sr'} uses assembly-mode CIGAR detection.
    dtype in {'pb','ont','hifi'} uses long-read-mode CIGAR + breakpoint detection.

    Returns (cigar_svs, cov_line, supp_svs, breakpoints).
    """
    cigar_svs  = []
    supp_svs   = []
    breakpoints = []
    cov_line   = None

    flag  = _r(record, 'flag')
    if flag == 4:
        return cigar_svs, cov_line, supp_svs, breakpoints

    mapq   = _r(record, 'mapq')
    rname  = _r(record, 'rname')
    rstart = _r(record, 'rstart')
    qname  = _r(record, 'qname')
    cigar  = _r(record, 'cigar')

    is_asm = dtype in ('cr', 'sr')
    is_lr  = dtype in ('pb', 'ont', 'hifi')

    # --- CIGAR-based DEL/INS + coverage ---
    if flag in (0, 16, 2048, 2064) and mapq >= minimaq:
        cigar_svs = call_cigar_del_ins(record, min_len, max_len)

        rend = _r(record, 'rend')
        cov_line = (f'{qname}\t{flag}\t{rname}'
                    f'\t{rstart}\t{rend}'
                    f'\t{mapq}\t{cigar}\n')

        # Long-read mode: breakpoints from large clips
        if is_lr:
            cnums  = [int(x) for x in re.findall(r'\d+', cigar)]
            ccodes = re.findall(r'[A-Z]', cigar)
            ref = rstart
            for code, length in zip(ccodes, cnums):
                if code in ('S', 'H') and length >= 500:
                    breakpoints.append(f'{qname}\t{flag}\t{rname}\t{ref}\n')
                if code in ('M', 'D', 'N', '=', 'X'):
                    ref += length

    # --- SA-based segment SV detection ---
    if msv == 'yes' and flag in (0, 16):
        sa_tag = _r(record, 'sa_tag')
        if sa_tag:
            seg = _build_segment(record)
            seg_list = [seg]

            for sa_entry in sa_tag.split(';'):
                if not sa_entry:
                    continue
                parts = sa_entry.split(',')
                if len(parts) < 5:
                    continue
                sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq = (
                    parts[0], int(parts[1]), parts[2], parts[3], int(parts[4]))

                if sa_chr == rname and sa_pos == rstart:
                    continue

                sflag   = 2048 if sa_strand == '+' else 2064
                sa_clip = parse_cigar2clipinfo(sa_cigar)
                sa_end  = sa_pos + sa_clip[1]

                seg_list.append([
                    qname, sflag, sa_chr, sa_pos, sa_end,
                    sa_clip, sa_mapq, '*',
                ])

            if 2 <= len(seg_list) <= 20:
                supp_svs = call_segment_svs(
                    seg_list, min_len, max_len, chromosome_list, minimaq,
                    chrom_size, dtype=dtype)

    return cigar_svs, cov_line, supp_svs, breakpoints


# ---------------------------------------------------------------------------
# Depth helper
# ---------------------------------------------------------------------------

def compute_depth(total_mapped_bases, chrom_size):
    return round(total_mapped_bases / chrom_size, 2) if chrom_size else 0.0
