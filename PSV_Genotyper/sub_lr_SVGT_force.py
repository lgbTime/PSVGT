import re
import subprocess
from collections import defaultdict
from math import ceil,floor
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types

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

def sam_parser2Breaks(region_sam, min_maq):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq = row.mapping_quality
        maqs  += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (flag & 0x4) or maq < min_maq:
        if flag & 0x4:
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
            update_breakpoints(chrom, breakpoints_to_update)
    if total_map_reads >0:
        #print( breakpoints, total_map_reads, ceil(maqs / total_map_reads))
        return breakpoints, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, 0, 0


def sam_primary_parser2Breaks_Del(region_sam, min_maq, sv_size):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = int(row.reference_start)
        chr = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
            update_breakpoints(chr, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                if 0.5 * sv_size < length < 2 * sv_size:
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 1
                    else:
                        deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    if total_map_reads >0:
        return breakpoints, deletions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_start):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    spans = 0
    maqs = 0
    effective_span = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chr = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
            update_breakpoints(chr, breakpoints_to_update)
        else:
            if (align_start + 1000 < sv_start)  and (align_end - 1000 > sv_start):
                effective_span += 1

        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I' and (0.5 * sv_size < length < 2 * sv_size):
                insertion_start = current_start - 1
                insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                if insert_key not in insertions:
                    insertions[insert_key] = 1
                else:
                    insertions[insert_key] += 1
                effective_span -= 1

    if total_map_reads >0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads), effective_span
    else:
        return {},{},0,0,effective_span
def sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size):
    ## some dup signal may hide in I cigar ##
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
            continue
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start
        align_end = row.reference_end
        total_map_reads += 1
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  
            elif ctype == 'I':  # Insertion
                if  0.7 * sv_size < length < 1.5 * sv_size  :  # Only count size equally 
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    else:
                        insertions[insert_key] += 1  # Count the insertion span
    if total_map_reads > 0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def determine_genotype(breaks, depth, homo_rate=0.75,ref_rate = 0.05):
    """
    depth and ratio base genotype
    """
    if depth == 0:
        return "./."
    if depth <= 5:
        if floor(homo_rate*depth) + 1 <= breaks:
            return "1/1"
        elif breaks / depth <= 0.2: ## 2 reads --> 5X
            return "0/0"
        else:
            return "0/1"

    if 5 < depth <= 10:
        if homo_rate * depth + 1 <=  breaks:
            return "1/1"
        elif 0.05 * depth + 1 > breaks:
            return "0/0"
        else:
            return "0/1"
    
    if depth > 10:
        if breaks / depth >= homo_rate: ## or 0.65   or 0.625
            return "1/1"
        elif breaks / depth < ref_rate:
            return "0/0"
        else:
            return "0/1"

def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift):
    info_return = []
    genotype = "0/0"  # Default genotype
    shift = min(sv_size, 500)
    sv_start_shift = set(range(sv_s - shift, sv_s + shift ))
    dup_shift = set(range(sv_s - 1000, sv_s - 1000)) ## samll dup capture as ins
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, maq, effective_span = sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s)
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads},maq=0")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        count_break_and_Ins = 0
        if inserts:  
            for pos in inserts.keys():
                ins_s, ins_e = map(int, pos.split(":")[1].split("-"))
                if ins_s in sv_start_shift or ins_s in dup_shift:
                    count_break_and_Ins += inserts[pos]
        if breakpoints:   #check breakpoints
            for breakpoint in breakpoints.get(chrome, {}).keys():
                if breakpoint in sv_start_shift: 
                    count_break_and_Ins += breakpoints[chrome][breakpoint]
        ins_ratio = round(count_break_and_Ins / total_map_reads, 3)
        genotype = determine_genotype(count_break_and_Ins, total_map_reads, homo_rate, ref_rate)
        if genotype == "1/1":
            if floor(0.1*total_map_reads) + 1 <=  effective_span:
                genotype = "0/1"
                #print(f"***************Correting SVINS {chrome}:{sv_s}-{sv_e} genotype to 0/1 since it has {effective_span} span reads*****************")
        #print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads},maq={maq}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return

def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift):
    ############ SVDel Case ##############
    info_return = []
    breaks_dict = {}
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, maq_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq, sv_size)
    breakpoints_r, deles_r, total_map_reads_r, maq_r = sam_primary_parser2Breaks_Del(right_sam, min_maq, sv_size)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0,maq=0")
        info_return.append(f"deles_l_ratio=0,deles_r_ratio=0;DEL")
        return  info_return
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
    else:
        deles_l_ratio = 0
    if total_map_reads_r:
        deles_r_ratio = round(count_break_and_deles_r / total_map_reads_r, 3)
    else:
        deles_r_ratio = 0
    deles_ratio = max(deles_l_ratio, deles_r_ratio)
    breaks_dict[count_break_and_deles_l ] = total_map_reads_l
    breaks_dict[count_break_and_deles_r ] = total_map_reads_r
    max_breaks = max(count_break_and_deles_l,count_break_and_deles_r)
    genotype = determine_genotype(max_breaks,breaks_dict[max_breaks], homo_rate, ref_rate)
    #print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l};total_map_reads_r={total_map_reads_r};maq={max(maq_l,maq_r)}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
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
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 3)
    else:
        break2_ratio = 0
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.8:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return

def little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, sv_type, shift=500):
    """
    Taking total local mapping  to genotyping the small(<7k) dup SV;
    1st Scan all insertion cigar; 2nd Capture all breakpoints;
    """
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome}:{bp1}", f"{chrome}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1  + shift))
    bp2_shift = set(range(bp2 - shift, bp2  + shift))
    ins_shift1 = set(range(bp1 - sv_size, bp1 + sv_size))
    ins_shift2 = set(range(bp2 - sv_size, bp2 + sv_size))

    breakpoints_bp, inserts, total_map_reads, maq = sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size)

    count_break_bp = 0
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts:
        for pos in inserts.keys():
            ins_s, _ = map(int, pos.split(":")[1].split("-"))
            if ins_s in ins_shift1 or ins_s in ins_shift2:
                count_break_bp += inserts[pos]
    if breakpoints_bp:
        for breakpoint in breakpoints_bp.get(chrome, {}).keys():
            if breakpoint in bp1_shift or breakpoint in bp2_shift:
                count_break_bp += breakpoints_bp[chrome][breakpoint]
    if total_map_reads:
        break_ratio =   round(count_break_bp / total_map_reads, 3)
    else:
        break_ratio = 0
    
    if break_ratio > 0.75:
        genotype = "1/1"
    elif break_ratio < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads={total_map_reads};\tbreakpoint_ratio={break_ratio}\tbp1={breakpoint1},bp2={breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads={total_map_reads};maq={maq}")
    info_return.append(f"bp1={breakpoint1},bp2={breakpoint2},bp_ratio={break_ratio};{sv_type}")
    return info_return


def dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints local mapping info to genotyping the big SV;
    
    """
    breaks_dict = {}
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift = set(range(bp2 - shift, bp2 + shift))
    ins_shift1 = set(range(bp1 - shift, bp1 + ceil(shift+ sv_size/2)))
    ins_shift2 = set(range(bp2 - shift - ceil(sv_size / 2), bp2 + shift))

    breakpoints_bp1, inserts_bp1, total_map_reads_bp1, maq1 = sam_primary_parser2Breaks_dup(bp1_sam, min_maq, sv_size)
    breakpoints_bp2, inserts_bp2, total_map_reads_bp2, maq2 = sam_primary_parser2Breaks_dup(bp2_sam, min_maq, sv_size)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
    
    if inserts_bp1:  # Check if there are insertion entries
        for pos in inserts_bp1.keys():
            if pos in ins_shift1:
                count_break_bp1 += inserts_bp1[pos]
    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]

    if inserts_bp2:  # Check if there are insertion entries
        for pos in inserts_bp2.keys():
            if pos in ins_shift2:
                count_break_bp2 += inserts_bp2[pos]
    if breakpoints_bp2:
        for breakpoint in breakpoints_bp2.get(chrome2, {}).keys():
            if breakpoint in bp2_shift:
                count_break_bp2 += breakpoints_bp2[chrome2][breakpoint]

    if total_map_reads_bp1:
        break1_ratio =   round(count_break_bp1 / total_map_reads_bp1, 2)
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 2)
    else:
        break2_ratio = 0

    max_break_ratio = max(break1_ratio, break2_ratio)
    breaks_dict[count_break_bp1] = total_map_reads_bp1
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_breaks = max(count_break_bp1, count_break_bp2)
    if break1_ratio >=0.4 and break2_ratio >= 0.4:
        genotype = "1/1"
    elif break1_ratio + break2_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio + break2_ratio < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
    return info_return

def parse_cigar2clipinfo(cigarstring):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigarstring)]
    cigars = re.findall(r'[MIDNSHP=X]', cigarstring)
    leftclip = 0
    rightclip = 0
    read_len = sum(length for length, ctype in zip(numbers, cigars) if ctype in 'MNP=XI')
    if cigars[0] in "SH":
        leftclip = numbers[0]
    if cigars[-1] in "SH":
        rightclip = numbers[-1]
    return [leftclip, read_len, rightclip]

def traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=3000):
    """"
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints, 
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal  
    """
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
    #print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return

def breaks2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
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
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 3)
    else:
        break2_ratio = 0
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio >= 0.85 and break2_ratio >= 0.85:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.10:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"***breakpoint to genotype TRA********* {sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2} ****************")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return


def supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
        end = line.reference_end
        if (bp1 - shift) <= end <= (bp1 + shift):
            direction = "-" if int(flag) & 0x10 else "+"
            if line.has_tag("SA"):
                supps = line.get_tag("SA").split(";")[:-1]
                for supp in supps:
                    chrom, cigars,start,maq1 =supp.split(",")[0], supp.split(",")[3], int(supp.split(",")[1]), int(supp.split(',')[4])
                    if supp.split(",")[2] != direction and ref_chr == chrom:
                        if (bp2-shift) <= start <= (bp2 + shift):
                            genotype += 1
    
    if bp1_map != 0:
        bp1_rate = genotype / bp1_map
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate, maq_bp1 = 0,0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
        if (bp2 - shift) <= start <= (bp2 + shift):
            direction = "-" if int(flag) & 0x10 else "+"
            if line.has_tag("SA"):
                supps = line.get_tag("SA").split(";")[:-1]
                for supp in supps:
                    chrom, cigars,start,maq1 =supp.split(",")[0], supp.split(",")[3], int(supp.split(",")[1]), int(supp.split(',')[4])
                    if supp.split(",")[2] != direction and ref_chr == chrom:
                        clipinfo = parse_cigar2clipinfo(cigars)
                        end = start + clipinfo[1]
                        if (bp1 - shift) <= end <= (bp1 + shift):
                            genotype += 1
    if bp2_map != 0:
        bp2_rate = genotype / bp2_map
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
        maq_bp2 = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    rate, bp1_map, bp2_map, maq = supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    info_return = []
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.05:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

def supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0

    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
        end = line.reference_end
        if (bp1 - shift) <= end <= (bp1 + shift):
            direction = "-" if int(flag) & 0x10 else "+"
            if line.has_tag("SA"):
                supps = line.get_tag("SA").split(";")[:-1]
                for supp in supps:
                    chrom, cigars,start,maq1 =supp.split(",")[0], supp.split(",")[3], int(supp.split(",")[1]), int(supp.split(',')[4])
                    if supp.split(",")[2] == direction and ref_chr == chrom:
                        if (bp2-shift) <= start <= (bp2 + shift):
                            genotype += 1
    if bp1_map != 0:
        bp1_rate = genotype / bp1_map
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
        if (bp2 - shift) <= start <= (bp2 + shift):
            direction = "-" if int(flag) & 0x10 else "+"
            if line.has_tag("SA"):
                supps = line.get_tag("SA").split(";")[:-1]
                for supp in supps:
                    chrom, cigars,start,maq1 =supp.split(",")[0], supp.split(",")[3], int(supp.split(",")[1]), int(supp.split(',')[4])
                    if supp.split(",")[2] == direction and ref_chr == chrom:
                        clipinfo = parse_cigar2clipinfo(cigars)
                        end = start + clipinfo[1]
                        if (bp1 - shift) <= end <= (bp1 + shift):
                            genotype += 1
    if bp2_map != 0:
        bp2_rate = genotype / bp2_map
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    info_return = []
    rate, bp1_map, bp2_map, maq = supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.625:
        genotype = '1/1'
    elif rate < 0.1:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return
