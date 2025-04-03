import re
import subprocess
from collections import defaultdict
from math import ceil
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types
def sam_parser2Breaks(region_sam, min_maq):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    total_map_reads = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = int(row.reference_start)
        chr = row.reference_name
        cigar = row.cigarstring
        if (flag & 0x4) or maq < min_maq:
        #if (flag & 0x4):
            continue
        total_map_reads += 1
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start
        align_end = row.reference_end
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:  ## use clip length to filter ???
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chr, breakpoints_to_update)
    return breakpoints, total_map_reads

def sam_parser2Breaks_Del(region_sam, min_maq, sv_size, breakpoint_pos,region, span_bp):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    effective_spans = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        if (flag & 0x4) or maq < min_maq:
        #if (flag & 0x4):
            continue
        total_map_reads += 1
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start 
        align_end = row.reference_end
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if sv_size >= 100:
                if (align_start + span_bp < breakpoint_pos) and (align_end - span_bp > breakpoint_pos):
                    effective_spans += 1
            else:## small sv full cover
                if align_start - 10 < region[0] and align_end - 10 > region[1]:
                    effective_spans += 1
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':  # Deletion
                if  sv_size - 5 <= length <= sv_size + 5 :  # Only count size equally 
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    effective_spans -= 1
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 1
                    deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    return breakpoints, deletions, total_map_reads, effective_spans

def sam_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, span_bp):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    effective_span = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        if (flag & 0x4) or maq < min_maq:
        #if (flag & 0x4):
            continue
        total_map_reads += 1
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start 
        align_end = row.reference_end
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if (align_start + span_bp < sv_s)  and (align_end - span_bp > sv_s):
                effective_span += 1
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I':  # Insertion
                if sv_size - 5 <= length <= sv_size + 5 :  # Only count size equally ## to get close del points
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    insertions[insert_key] += 1  # Count the insertion span
                    effective_span -= 1
    return breakpoints, insertions, total_map_reads, effective_span

def get_spans(chromosome, start, cigar_types, cigar_numbers):
    """Generate spans based on CIGAR information."""
    current_position = start
    spans = {}
    for ctype, cnum in zip(cigar_types, cigar_numbers):
        if ctype == 'M':
            # For matches, calculate the span and update current position
            end_position = current_position + cnum
            span_key = f"{chromosome}:{current_position}-{end_position}"
            spans[span_key] = (current_position, end_position)
            current_position = end_position  # Move current position to end of match
        elif ctype == 'D':
            # For deletions, move the current position forward
            current_position += cnum  # Skip over the deleted region
    return spans

def sam_parser2SVInDel_CutSpan(region_sam, min_maq):
    """Process the SAM file to extract spans and count overlaps with a given range."""
    spans = defaultdict(int)
    for row in  region_sam:
        chrom, start, maq, cigar = row.reference_name, row.reference_start, row.mapping_quality, row.cigarstring
        if maq >= min_maq:
            cigar_numbers, cigar_types = parse_cigar(cigar)
            new_spans = get_spans(chrom, start, cigar_types, cigar_numbers)
            for span_key in new_spans:
                spans[span_key] += 1 
    return spans

def calculate_coverage(cut_span, sv_point):
    cov = 0
    for span in cut_span:
        cov_s = int(span.split(":")[1].split("-")[0])
        cov_e = int(span.split(":")[1].split("-")[1])
        # Check coverage around sv_s
        if sv_point - 75 > cov_s and sv_point + 75 < cov_e:
            cov += 1
    return cov

def determine_genotype(entry_ratio, homo_rate, ref_rate):
    """
    ONT easy lead to 0/1 and FP, here we try modify.
    """
    if entry_ratio  >= homo_rate:
        return "1/1"
    elif entry_ratio <  ref_rate:
        return "0/0"
    else:
        return "0/1"

def breaksCallGT(break_l_ratio, break_r_ratio):
    """
    INV, TRA, DUP, 
    Take two breakpoints for Genotype
    """
    if max(break_l_ratio, break_r_ratio) > 0.75 and min(break_l_ratio, break_r_ratio) > 0.4:
        return "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.2 or min(break_l_ratio, break_r_ratio) < 0.1:
        return "0/0"
    else:
        return "0/1"



def determine_dupGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.25 and break_r_ratio >= 0.25:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.3 and min(break_l_ratio, break_r_ratio) >= 0.1:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype

def determine_invGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.6 and break_r_ratio >= 0.6:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.75  and min(break_l_ratio, break_r_ratio) >= 0.3:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype


def determine_traGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.6 and break_r_ratio >= 0.6:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.75  and min(break_l_ratio, break_r_ratio) >= 0.3:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype


def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift=100, span_bp=50):
    info_return = []
    genotype = "0/0"  # Default genotype
    #sv_start_shift = set(range(sv_s - shift, sv_s + shift+100))
    #sv_end_shift   = set(range(sv_e - shift, sv_e + shift+100))
    sv_start_shift = set(range(sv_s - shift, sv_s + 120))
    sv_end_shift   = set(range(sv_e - shift, sv_e + 120))
    sv_size = abs(sv_size)
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, effective_spans = sam_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, span_bp)
    if total_map_reads == 0:
        info_return.append('./.')
        info_return.append(f"total_map_reads={total_map_reads}")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        cut_span = sam_parser2SVInDel_CutSpan(region_sam, min_maq)
        count_break_and_Ins = 0
        covIns = calculate_coverage(cut_span, sv_s)
        if inserts:  # Check if there are insertion entries
            for pos in inserts.keys():
                ins_s, ins_e = map(int, pos.split(":")[1].split("-"))
                # Check if insertions are within the shifted range
                if ins_s in sv_start_shift and ins_e in sv_end_shift:
                    count_break_and_Ins += inserts[pos]
        if breakpoints:  # No ins recording, check breakpoints
            for breakpoint in breakpoints.get(chrome, {}).keys():
                if breakpoint in sv_start_shift and breakpoint in sv_end_shift: ### and ?  or ?
                    count_break_and_Ins += breakpoints[chrome][breakpoint]
        ins_ratio = round(count_break_and_Ins / total_map_reads, 3)
        covIns_ratio = round(covIns / total_map_reads, 3)
        
        genotype = determine_genotype(ins_ratio, homo_rate, ref_rate)
        if genotype == "0/1" and effective_spans <= 0.01*total_map_reads+1:
            genotype = "1/1"
        elif genotype == "1/1" and effective_spans >= 0.05*total_map_reads+1:
            genotype = "0/1"
        elif genotype == "0/0":
            if effective_spans == 0:
                genotype = "1/1"
            #else:
            #    genotype = "0/1"
        print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\tIns_points_covered_ratio:{covIns_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads};effective_spans={effective_spans}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return
def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift=100, span_bp=50):
    ############ SVDel Case ##############
    info_return = []
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, effective_spans_l = sam_parser2Breaks_Del(left_sam, min_maq, sv_size,  sv_s, [sv_s,sv_e], span_bp)
    breakpoints_r, deles_r, total_map_reads_r, effective_spans_r = sam_parser2Breaks_Del(right_sam, min_maq, sv_size, sv_e, [sv_s,sv_e], span_bp)
    cut_span_l = sam_parser2SVInDel_CutSpan(left_sam, min_maq)
    cut_span_r = sam_parser2SVInDel_CutSpan(right_sam, min_maq)
    covDel_l = calculate_coverage(cut_span_l, sv_s)
    covDel_r = calculate_coverage(cut_span_l, sv_e)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("1/1") ## no reads properly dele
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0")
        info_return.append(f"deles_l_ratio=0;deles_r_ratio=0")
        return info_return
    if deles_l:  # Check if there are deletion entries ####### if there are deles but not the target deles
        for pos in deles_l.keys():
            dele_s, dele_e = map(int, pos.split(":")[1].split("-"))
            # Check if deletions are within the shifted range left right will have the same results #
            if dele_s in sv_start_shift:
                count_break_and_deles_l += deles_l[pos]
    if deles_r:
        for pos in deles_r.keys():
            dele_s, dele_e = map(int, pos.split(":")[1].split("-"))
            if dele_e in sv_end_shift:
                count_break_and_deles_r += deles_r[pos]
    if breakpoints_l:
        for breakpoint in breakpoints_l.get(chrome, {}).keys():
            if breakpoint in sv_start_shift:
                count_break_and_deles_l += breakpoints_l[chrome][breakpoint]
    if breakpoints_r:
        for breakpoint in breakpoints_r.get(chrome, {}).keys():
            if breakpoint in sv_end_shift: 
                count_break_and_deles_r += breakpoints_r[chrome][breakpoint]
    if total_map_reads_l:
        deles_l_ratio = round(count_break_and_deles_l / total_map_reads_l, 3)
        covDel_l_ratio =round( covDel_l / total_map_reads_l, 3)
    else:
        deles_l_ratio = 0
        covDel_l_ratio = 0
    if total_map_reads_r:
        deles_r_ratio = round(count_break_and_deles_r / total_map_reads_r, 3)
        covDel_r_ratio = round( covDel_r / total_map_reads_r, 3)
    else:
        deles_r_ratio = 0
        covDel_r_ratio = 0
    deles_ratio = max(deles_l_ratio, deles_r_ratio)
    genotype = determine_genotype(deles_ratio, homo_rate, ref_rate)
    if genotype == "0/1":
        if effective_spans_l + effective_spans_r == 0:
            genotype = "1/1"
    elif genotype == "1/1":
        if (effective_spans_l + effective_spans_r) >= (0.3 * total_map_reads + 2):
            genotype = "0/1"
    elif genotype == "0/0":
        if effective_spans_l + effective_spans_r == 0:
            genotype = '1/1'
    print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\tdele_l_covered_ratio:{covDel_l_ratio}\tdele_r_covered_ratio:{covDel_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l},span_l={effective_spans_l};total_map_reads_r={total_map_reads_r},span_r={effective_spans_r}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2GT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=50):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ INV Case ##############
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2 = sam_parser2Breaks(bp2_sam, min_maq)
    cut_span_bp1 = sam_parser2SVInDel_CutSpan(bp1_sam, min_maq)
    cut_span_bp2 = sam_parser2SVInDel_CutSpan(bp2_sam, min_maq)
    cov_break_bp1 = calculate_coverage(cut_span_bp1, bp1)
    cov_break_bp2 = calculate_coverage(cut_span_bp2, bp2)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.") 
        info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")

    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]
    if breakpoints_bp2:
        for breakpoint in breakpoints_bp2.get(chrome2, {}).keys():
            if breakpoint in bp2_shift: 
                count_break_bp2 += breakpoints_bp2[chrome2][breakpoint]
    if total_map_reads_bp1:
        break1_ratio =   round(count_break_bp1 / total_map_reads_bp1, 3)
        cov_break1_ratio =round( cov_break_bp1 / total_map_reads_bp1, 3)
    else:
        break1_ratio = 0
        cov_break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 3)
        cov_break2_ratio = round( cov_break_bp2 / total_map_reads_bp2, 3)
    else:
        break2_ratio = 0
        cov_break2_ratio = 0
    max_break_ratio = max(break1_ratio, break2_ratio)
    if sv_type == "DUP":
        genotype = determine_dupGT(break1_ratio, break2_ratio)
    elif sv_type == "INV":
        genotype = determine_invGT(break1_ratio, break2_ratio)
    elif sv_type == "TRA":
        genotype = determine_traGT(break1_ratio, break2_ratio)
    print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}_covered_ratio={cov_break1_ratio}\t{breakpoint2}_covered_ratio={cov_break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
    return info_return

def sam2readsID(region_sam):
    readsID = []
    maqs = 0
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        readID = row.query_name
        readsID.append(readID)
    if not readsID:
        return [], 0
    else:
        return readsID, ceil(maqs / len(readsID))

def traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=50):
    genotype = "0/0"
    info_return = []
    breaks_dict = {}
    bp1_readsID, maq1 = sam2readsID(bp1_sam)
    bp2_readsID, maq2 = sam2readsID(bp2_sam)

    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    if max(bp1_tra, bp2_tra) > 0.90:
        genotype = "1/1"
    elif bp1_tra >= 0.8 and bp2_tra >= 0.8:
        genotype = "1/1"
    elif bp1_tra + bp2_tra < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    info_return.append(genotype)
    print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return

