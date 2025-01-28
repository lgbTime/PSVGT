import re
import subprocess
from collections import defaultdict
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types

def sam2readsID(region_sam):
    readsID = []
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        readID = row.query_name
        readsID.append(readID)
    return readsID

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


def sam_primary_parser2Breaks_Del(region_sam, min_maq, sv_size):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
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
        if (row.flag & 0x4) or (row.mapping_quality < min_maq):
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
            elif ctype == 'D':  # Deletion
            #elif ctype == 'D' and (0.5 * sv_size < length < 2 * sv_size):
                if  length >= 35 :  # Only count size equally 
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 0
                    deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    return breakpoints, deletions, total_map_reads

def sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = row.reference_start
        chr = row.reference_name
        cigar = row.cigarstring
        if (row.flag & 0x4) or (row.mapping_quality < min_maq):
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
                current_start += length  # Increment position past deletion
            elif ctype == 'I' and (0.5 * sv_size < length < 2 * sv_size):
                if  length > 35 :  # Only count size equally ## to get close del points
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 0
                    insertions[insert_key] += 1  # Count the insertion span
    return breakpoints, insertions, total_map_reads
def sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size):
    ## some dup signal may hide in I cigar ##
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = row.reference_start
        chr = row.reference_name
        cigar = row.cigarstring
        if (row.flag & 0x4) or (row.mapping_quality < min_maq):
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
                current_start += length  # Increment position past deletion
            elif ctype == 'I':  # Insertion
                if   0.5 * sv_size < length < 2 * sv_size  :  # Only count size equally 
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 0
                    insertions[insert_key] += 1  # Count the insertion span
    return breakpoints, insertions, total_map_reads

def determine_genotype(entry_ratio):
    """
    ONT easy lead to 0/1 and FP, here we try modify.
    """
    if entry_ratio >= 0.650: ## or 0.65   or 0.625
        return "1/1"
    elif entry_ratio <  0.10:
        return "0/0"
    else:
        return "0/1"

def breaksCallGT(break_l_ratio, break_r_ratio, sv_type):
    """
    INV, DUP, 
    Take two breakpoints for Genotype, as the shift parameter setting to a larger number the region clip reads ratio will be reduced, 
    here we must ensure the shifting and ratio should have a good adjust
    """
    if sv_type == "DUP":
        if max(break_l_ratio, break_r_ratio)   > 0.7:
            return "1/1"
        elif max(break_l_ratio, break_r_ratio) < 0.2:
            return "0/0"
        else:
            return "0/1"

    if sv_type == "INV":
        if max(break_l_ratio, break_r_ratio) > 0.7 and min(break_l_ratio, break_r_ratio) >= 0.55:
            return "1/1"
        elif max(break_l_ratio, break_r_ratio) < 0.3 or  min(break_l_ratio, break_r_ratio) < 0.15:
            return "0/0"
        else:
            return "0/1"


def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, shift):
    info_return = []
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    sv_size = abs(sv_size)
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads = sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size)
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads}")
        info_return.append(f"INS_rate=0;INS")
    else:
        count_break_and_Ins = 0
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
        genotype = determine_genotype(ins_ratio)
        print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return


def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, shift):
    ############ SVDel Case ##############
    info_return = []
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq, sv_size)
    breakpoints_r, deles_r, total_map_reads_r = sam_primary_parser2Breaks_Del(right_sam, min_maq, sv_size)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads}")
        info_return.append(f"DEL_rate=0;DEL")
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
    genotype = determine_genotype(deles_ratio)
    print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l};total_map_reads_r={total_map_reads_r}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2GT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0")
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
    max_break_ratio = max(break1_ratio, break2_ratio)
    genotype = breaksCallGT(break1_ratio, break2_ratio, sv_type)
    print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return
def dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    <<< how about using supp align to genotyoe >>> 
    """
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, inserts_bp1, total_map_reads_bp1 = sam_primary_parser2Breaks_dup(bp1_sam, min_maq, sv_size)
    breakpoints_bp2, inserts_bp2, total_map_reads_bp2 = sam_primary_parser2Breaks_dup(bp2_sam, min_maq, sv_size)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
    
    if inserts_bp1:  # Check if there are insertion entries
        for pos in inserts_bp1.keys():
            ins_s, ins_e = map(int, pos.split(":")[1].split("-"))
            # Check if insertions are within the shifted range
            if ins_s in bp1_shift and ins_e in bp1_shift:
                count_break_bp1 += inserts_bp1[pos]
    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]

    if inserts_bp2:  # Check if there are insertion entries
        for pos in inserts_bp2.keys():
            ins_s, ins_e = map(int, pos.split(":")[1].split("-"))
            # Check if insertions are within the shifted range
            if ins_s in bp2_shift and ins_e in bp2_shift:
                count_break_bp2 += inserts_bp2[pos]
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
    max_break_ratio = max(break1_ratio, break2_ratio)
    genotype = breaksCallGT(break1_ratio, break2_ratio, sv_type)
    print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
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
    genotype = '0/0'
    info_return = []
    bp1_readsID, bp2_readsID = sam2readsID(bp1_sam), sam2readsID(bp2_sam)
    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    genotype = determine_genotype(max(bp1_tra, bp2_tra))
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)}')
    info_return.append(f"{chrome1}:{bp1}_supp_ratio={bp1_tra},{chrome2}:{bp2}_supp_ratio={bp2_tra};{'TRA'}")
    return info_return

def supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=700):
    bp1_map = 0
    bp2_map = 0
    genotype = 0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
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
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
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
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map

def invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    rate, bp1_map, bp2_map = supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    info_return = []
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.15:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map}')
    info_return.append(f"{chrome1}:{bp1}_rate={rate},{chrome2}:{bp2}_rate={rate};segment_captured;{sv_type}")
    return info_return

def supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map = 0
    bp2_map = 0
    genotype = 0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
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
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
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
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map

def supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    info_return = []
    rate, bp1_map, bp2_map = supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.2:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map}')
    info_return.append(f"{chrome1}:{bp1}_rate={rate},{chrome2}:{bp2}_rate={rate};segment_captured;{sv_type}")
    return info_return
