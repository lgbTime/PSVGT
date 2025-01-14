import pysam
import re

def process_chromosome(chromosome,chromosome_list, bamfile_path, minLen, maxLen, min_maq, SVsignal_out_path, cov_out_path,dtype,msv):
    samfile = pysam.AlignmentFile(bamfile_path, 'rb')
    with open(f"{SVsignal_out_path}_{chromosome}.record.txt", 'w') as indel_out, open(f"{cov_out_path}_{chromosome}", 'w') as cov_out, open(f"{SVsignal_out_path}_{chromosome}.record.txt.suppAlign", 'w') as supp_sv_out:
        # Fetch all reads from the chromosome
        lines = samfile.fetch(chromosome)
        for line in lines:
            if dtype in ['sr', 'cr']:
                result = svInDel4asm(line, minLen, min_maq)
                if result is None:
                    continue
                results, cov_line = result
                if results:
                    indel_out.writelines(results)
                cov_out.writelines([cov_line])
            elif dtype in ['pb', 'ont', 'hifi']:
                svInDels, cov_line, supp_svsignal = svInDel4lr(line, minLen, min_maq, maxLen, msv, chromosome_list)
                if svInDels:
                    indel_out.writelines(svInDels)
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")
                if cov_line:
                    cov_out.writelines(cov_line)
        indel_out.close()
        cov_out.close()
        supp_sv_out.close()

def svInDel4asm(line, minLen, min_maq):
    if line.flag == 4 or line.mapping_quality < min_maq:
        return None  
    else:
        query_chr = line.query_name
        query_seq = line.query_sequence
        target_chr = line.reference_name
        target_start = line.reference_start
        target_end = line.reference_end
        maq = line.mapping_quality
        flag = line.flag
        strand = "-" if flag & 0x10 else "+"
        cigar = line.cigarstring
        cigar_numbers = list(map(int, re.findall(r'\d+', cigar)))
        cigar_codes = re.findall(r'[A-Z]', cigar)
        # Initialize variables for tracking reference and query positions
        ref = target_start
        query = 0
        results = []
        for code, length in zip(cigar_codes, cigar_numbers):
            if code == "M":  # Match
                ref += length
                query += length
            elif code == "D":  # Deletion
                if length >= minLen:
                    results.append(f'{target_chr}\t{query_chr}\t{ref}\t{ref + length - 1}\t{length}\t{maq}\t{target_chr}:{ref}-{ref+length-1}_DEL={length}\tDEL\t"*"\n')
                ref += length
            elif code == "I":  # Insertion
                if length >= minLen:
                    ins_seq = query_seq[query : query + length ]
                    results.append(f"{target_chr}\t{query_chr}\t{ref}\t{ref + 1}\t{length}\t{maq}\t{target_chr}:{ref}-{ref+1}_INS={length}\tINS\t{ins_seq}\n")
                query += length
            #elif code in {"N", "S", "H"}:
            #    ref += length if code == "N" else 0
            #    query += length
        return results, f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n'  # Return results and coverage line


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

def segmentsv4lr(supp_list,min_size, max_size, chromosome_list):
    """
    Collecting the supplementary alignment to identify the INS, DEL, DUP, INV, TRA signals
    LongRead align to reference that has "SA" will be recorded as dict as follow:
    { 'm64144': [['m64144',2048,'Db-Chr4',8949283,8949610,[0, 327, 11157],'43'],
                 ['m64144',0,   'Db-Chr4',8949649,8957962,[3162, 8322, 0],60,'AAXXXCTAATT'],
                 ['m64144',2064,'Db-Chr3',12497386,12502117,[5171, 4731, 1582],'60']]
                 }
    these code is refenced on DeBreak, the INS and DEL may be discarded in the future
    """
    if not any(int(supp[1]) <= 16 for supp in supp_list) or len(supp_list) <= 1:
        return []
    else:
        svsignal_supp = []
        def add_svcall(chrom,readname, start, end, size, maq, svid,sv_type, sequence=None):
            if maq == 60:
                sv_record = f"{chrom}\t{readname}\t{start}\t{end}\t{size}\t{maq}\t{svid}\t{sv_type}"
                if sequence:
                    sv_record += f"\t{sequence}"
                svsignal_supp.append(sv_record)
        primary_map = next(supp for supp in supp_list if int(supp[1]) <= 16)
        supps = [maplist for maplist in supp_list if maplist != primary_map]
        readseq = primary_map[7]
        pri_chrom = primary_map[2]
        pri_flag = int(primary_map[1]) % 32 > 15
        samedir, invdir, diffchr = [], [], []
        for supp in supps:
            chrom, supp_flag, supp_maplen = supp[2], int(supp[1]) % 32 > 15 , supp[5][1]
            if supp_maplen < 300: ### hifi may be setting to 500 also ok
                continue
            if chrom != pri_chrom:
                diffchr.append(supp)
            elif supp_flag != pri_flag:
                invdir.append(supp)
            else:
                samedir.append(supp)
        for supp_map in samedir:
            leftmap,rightmap = (primary_map, supp_map) if supp_map[3] > primary_map[3] else (supp_map, primary_map)
            sh1,len1,sh2 = leftmap[5][0],  leftmap[5][1], leftmap[5][2]
            sh3,len2,sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
            maq1,maq2 = leftmap[6], rightmap[6]
            readname = primary_map[0]
            maq = (int(maq1) + int(maq2)) // 2
            ## DEL ##
            overlapmap = sh1 + len1 - sh3
            if -200 < overlapmap < 1500:   ## check or search for better
                del_size = rightmap[3] - leftmap[4] + overlapmap
                if min_size <= del_size <= max_size:
                    sv_start = max(0, leftmap[4] - max(0, overlapmap))
                    sv_end = sv_start + del_size -1
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_DELLEN={del_size}'
                    add_svcall(pri_chrom, readname, sv_start, 
                               sv_end, del_size, maq, svid,"DEL",'*')
            ## INS ##
            if abs(rightmap[3] - leftmap[4]) <= 300: ### how about 200 ?
                overlapmap = rightmap[3] - leftmap[4]
                ins_size = sh3 - len1 - sh1 - overlapmap
                if min_size <= ins_size <= max_size:
                    sv_start = min(rightmap[3],leftmap[4])
                    sv_end = sv_start + 1
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_INSLEN={ins_size}'
                    seq = readseq[sh1 + len1: sh3 - overlapmap]
                    add_svcall(pri_chrom, readname, sv_start,
                               sv_end, ins_size, maq, svid,"INS",seq)
            ## DUP ##
            lap1 = sh1 + len1 - sh3
            if -200 < lap1 < 500 and (leftmap[4] - rightmap[3]) >= max(50, lap1):
                dup_size = leftmap[4] - rightmap[3] - max(lap1, 0)
                if min_size <= dup_size <= max_size:
                    sv_start, sv_end = rightmap[3],rightmap[3] + dup_size
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                    add_svcall(pri_chrom, readname, sv_start,
                               sv_end, dup_size, maq, svid,"DUP","*")
            lap2 = sh3 + len2 - sh1
            if -200 < lap2 < 500 and (rightmap[4] - leftmap[3]) >= max(1000, lap2):
                dup_size = rightmap[4] - leftmap[3] - lap2
                if min_size <= dup_size <= max_size:
                    sv_start,sv_end = leftmap[3], leftmap[3] + dup_size
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                    add_svcall(pri_chrom, readname, sv_start,
                            sv_end, dup_size, maq, svid, 'DUP', "*")
        ## INV ##
        for supp in invdir:
            if (supp[3] > primary_map[3] and (supp[4] - primary_map[4]) > -200) or \
           (supp[3] < primary_map[3] and (primary_map[4] - supp[4]) > -200):
                leftmap = primary_map if supp[3] > primary_map[3] else supp
                rightmap = supp if supp[3] > primary_map[3] else primary_map
            
                sh1,len1,sh2 = leftmap[5][0],  leftmap[5][1], leftmap[5][2]
                sh3,len2,sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
                maq1,maq2 = leftmap[6], rightmap[6]
                readname = primary_map[0]
                maq = (int(maq1) + int(maq2)) // 2
                
                lap1 = sh3 + len2 - sh2
                if -200 < lap1 < 500 and (rightmap[4] - leftmap[4]) > max(100, lap1):
                    inv_size = rightmap[4] - leftmap[4] - lap1
                    if min_size <= inv_size <= max_size:
                        sv_start,sv_end = leftmap[4],leftmap[4] + inv_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                            sv_end, inv_size, maq, svid,'INV', "*")
            
                lap2 = sh4 + len2 - sh1
                if -200 < lap2 < 500 and  (rightmap[3]-leftmap[3])>=max(100,lap2):
                    inv_size = rightmap[3] - leftmap[3] - lap2
                    if min_size <= inv_size <= max_size:
                        sv_start,sv_end = leftmap[3],leftmap[3] + inv_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                            sv_end, inv_size, maq, svid, "INV", "*")
        ## TRA ##
        ### win size 500 may fit into the TRA ###
        for supp in diffchr:
            maq1,maq2 = supp[6], primary_map[6]
            readname = primary_map[0]
            maq = (int(maq1) + int(maq2)) // 2
            sh1,len1,sh2 = primary_map[5][0],  primary_map[5][1], primary_map[5][2]
            sh3,len2,sh4 = supp[5][0], supp[5][1], supp[5][2]
            breakpoint1, breakpoint2 = '', ''
            if abs(sh1-sh3-len2)<= 500 or abs(sh1-len2-sh4) <= 500:
                chrom1=primary_map[2]
                breakpoint1=primary_map[3]
            elif abs(sh2-sh3-len2)<= 500 or abs(sh2-len2-sh4)<= 500:
                chrom1=primary_map[2]
                breakpoint1=primary_map[4]
            if abs(sh3-sh1-len1)<= 500 or abs(sh3-len1-sh2) <= 500:
                chrom2=supp[2]
                breakpoint2=supp[3]
            elif abs(sh4-sh1-len1) <= 500 or abs(sh4-len1-sh2)<= 500:
                chrom2=supp[2]
                breakpoint2=supp[4]
                ## Chr4	ERR3415829.505585   15689499    9991420	0	60	Chr4:15689499_Chr1:9991420_TRA TRA
            if breakpoint1!='' and breakpoint2!='' and maq == 60 and (chromosome_list.index(chrom1) < chromosome_list.index(chrom2)):
                svid = chrom2+':'+str(breakpoint2)+'_'+chrom1+':'+str(breakpoint1) 
                svsignal_supp +=[chrom1+'\t'+primary_map[0]+'\t'+str(breakpoint1)+'\t'+str(breakpoint2)+'\t0\t'+str(maq)+"\t"+svid+'\tTRA' + "\t*"]
            elif breakpoint1!='' and breakpoint2!='' and maq == 60 and (chromosome_list.index(chrom1) > chromosome_list.index(chrom2)):
                svid = chrom1+':'+str(breakpoint1)+'_'+chrom2+':'+str(breakpoint2) 
                svsignal_supp +=[chrom2+'\t'+primary_map[0]+'\t'+str(breakpoint2)+'\t'+str(breakpoint1)+'\t0\t'+str(maq)+"\t"+svid+'\tTRA'+ "\t*"]
        return svsignal_supp

def svInDel4lr(line, minLen, min_maq, maxLen, msv, chromosome_list):
    supp_dict = {}
    svInDels = []
    supp_svsignal = []
    covinfo = ''
    readname = line.query_name
    maq = line.mapping_quality
    target_start = line.reference_start
    refend = line.reference_end
    clipinfo = [0,0,0]
    if line.flag != 4:
        if line.cigar[0][0] in [4,5]:  ## left most cigar
            clipinfo[0] = line.cigar[0][1]
        if line.cigar[-1][0] in [4,5]:  ## right most cigar
            clipinfo[2] = line.cigar[-1][1]
        clipinfo[1] = line.query_alignment_length
    if msv == "yes":
        if line.is_supplementary:
            supp = [line.query_name, line.flag, line.reference_name,target_start, refend, clipinfo,maq,'']
            #print(supp)
            if readname in supp_dict.keys():
                supp_dict[readname] += [supp]
            else:
                supp_dict[readname]  = [supp]

    if line.flag  in [0,16] and line.mapping_quality >= min_maq:
        query_chr = line.query_name
        query_seq = line.query_sequence
        target_chr = line.reference_name
        target_start = line.reference_start
        target_end   = line.reference_end
        maq = line.mapping_quality
        flag = line.flag
        strand = "-" if flag & 0x10 else "+"
        cigar = line.cigarstring
        cigar_numbers = list(map(int, re.findall(r'\d+', cigar)))
        cigar_codes = re.findall(r'[A-Z]', cigar)
        # Initialize variables for tracking reference and query positions
        ref = target_start
        query = 0
        for code, length in zip(cigar_codes, cigar_numbers):
            if code == "M":  # Match
                ref += length
                query += length
            elif code == "D":  # Deletion
                if length >= minLen:
                    svInDels.append(f'{target_chr}\t{query_chr}\t{ref}\t{ref+length-1}\t{length}\t{maq}\t{target_chr}:{ref}-{ref+length-1}_DEL={length}\tDEL\t"*"\n')
                ref += length
            elif code == "I":  # Insertion  base in strandness and position to get Insertions
                if length >= minLen:
                    ins_seq = query_seq[query : query + length ]
                    #if len(ins_seq) < 50:
                    #    print("!!!!!!!!! ins_seq extract error !!!!!!!!!")
                    svInDels.append(f"{target_chr}\t{query_chr}\t{ref}\t{ref + 1}\t{length}\t{maq}\t{target_chr}:{ref}-{ref+1}_INS={length}\tINS\t{ins_seq}\n")
                query += length
            #elif code in {"N", "S", "H"}:
            #    ref += length if code == "N" else 0
            #    query += length
        covinfo = f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n'
        if msv == "yes":
            if line.has_tag("SA"):
                primary = [readname, line.flag, line.reference_name,target_start,refend,clipinfo, maq, line.query_sequence ]
                if readname in supp_dict.keys():
                    supp_dict[readname] += [primary]
                else:
                    supp_dict[readname] = [primary]
                #print(primary)
                supps = line.get_tag("SA").split(";")[:-1]
                diffchr = []
                for supp in supps:
                    chrom, cigars,start,maq =supp.split(",")[0], supp.split(",")[3], int(supp.split(",")[1]), int(supp.split(',')[4])
                    if supp.split(",") == line.reference_name:
                        continue
                    if supp.split(",")[2] == "+":
                        sflag = 2048
                    else:
                        sflag = 2064
                    clipinfo = parse_cigar2clipinfo(cigars)
                    end = start + clipinfo[1]
                    suppinfo = [readname, sflag, chrom, start, end, clipinfo, maq]
                    #print(suppinfo)
                    supp_dict[readname] += [suppinfo]
        if supp_dict:
            if 2<= len(supp_dict[readname]) <= 20:
                supp_svsignal = segmentsv4lr(supp_dict[readname],minLen, maxLen, chromosome_list)
    return svInDels, covinfo , supp_svsignal

