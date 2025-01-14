from os.path import basename
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from tqdm import tqdm
import pandas as pd
import pysam
from math import ceil
def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def get_cov(opened_sam, reference, start, end):
    num_maps = opened_sam.count(reference, start, end)
    #total_depth = [sum(depths) for depths in zip(*coverage)]
    #average_depth = sum(total_depth) / len(total_depth) if total_depth else 0
    return num_maps

def determine_genotype(entry_ratio):
    """
    ONT easy lead to 0/1 and FP, here we try modify.
    """
    if entry_ratio > 0.6:
        return "1/1"
    elif entry_ratio <  0.1:
        return "0/0"
    else:
        return "0/1"

def process_sv(sv_line, mapfile, name, min_maq, shift=100):
    opened_sam = pysam.AlignmentFile(mapfile, "rb")
    sampleID = name
    info = sv_line.strip().split("\t")
    sv_size, signal_num = int(info[3]),int(info[7])

    if "TRA" not in info[5]:
        chrome, bp1, bp2= info[0], int(info[1]), int(info[2])
        left_most = max(bp1 - 2*shift,0) ## aviod the region start at position value less than shift
        left_cov  = get_cov(opened_sam, reference=chrome, start=left_most, end= max(bp1 - shift, left_most))
        right_cov = get_cov(opened_sam, reference=chrome, start=bp2 + shift, end=bp2 + 2*shift)

        if "DEL" in info[5]:
            if left_cov == 0:
                del_l =0
            else:
                del_l = signal_num / left_cov
            if right_cov == 0:
                del_r = 0
            else:
                del_r = signal_num / right_cov
            del_rate =max(del_l,del_r)
            del_gt =  determine_genotype(del_rate)
            out = [chrome, bp1, bp2, sv_size, -sv_size, del_gt,f'total_map_reads_l={left_cov};total_map_reads_r={right_cov}', f'deles_l_ratio={round(del_r,2)},deles_r_ratio={round(del_r,2)};DEL' ]
            return "\t".join(map(str,out))
        
        elif "INS" in info[5]:
            cov = ceil((left_cov + right_cov) / 2)
            if cov == 0:
                ins_rate = 1
            else:
                ins_rate = signal_num / cov
            ins_gt = determine_genotype(ins_rate)
            out = [chrome, bp1, bp2, -sv_size, sv_size, ins_gt, f'total_map_reads={cov}', f'ins_ratio={round(ins_rate,2)};INS']
            return "\t".join(map(str, out))

        elif "DUP" in info[5]:
            if left_cov ==0:
                bp1_rate = 0
            else:
                bp1_rate = round(signal_num / left_cov, 2)
            if right_cov ==0:
                bp2_rate = 0
            else:
                bp2_rate = round(signal_num / right_cov, 2)
            dup_gt = determine_genotype(max(bp1_rate,bp2_rate))
            out = [chrome, bp1, bp2, -sv_size, sv_size, dup_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'{chrome}:{bp1}_ratio={bp1_rate},{chrome}:{bp2}_ratio={bp2_rate};DUP']
            return "\t".join(map(str, out))

        elif "INV" in info[5]:
            if left_cov ==0:
                bp1_rate = 0
            else:
                bp1_rate = round(signal_num / left_cov, 2)
            if right_cov ==0:
                bp2_rate = 0
            else:
                bp2_rate = round(signal_num / right_cov,2)
            inv_gt = determine_genotype(max(bp1_rate,bp2_rate))
            out = [chrome, bp1, bp2, -sv_size, sv_size, inv_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'{chrome}:{bp1}_ratio={bp1_rate},{chrome}:{bp2}_ratio={bp2_rate};INV']
            return "\t".join(map(str, out))
    elif "TRA" in info[5]:
        svid = info[4]
        ## gemome chr should not has _,since we take _ to extract chr in svid
        chr2_s2, chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1, bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2, bp2 = chr2_s2[0], int(chr2_s2[1])
        bp1_left = max(bp1 - 2*shift ,0) ## aviod the region start at position value less than shift
        bp2_left = max(bp2 - 2*shift ,0) ## aviod the region start at position value less than shift
        bp1_cov  = get_cov(opened_sam, chrome1, bp1_left,  max(bp1 - shift, bp1_left))
        bp2_cov =  get_cov(opened_sam, chrome2, bp2_left,  max(bp2 - shift, bp2_left))
        if bp1_cov == 0:
            tra_rate1 = 0
        else:
            tra_rate1 = round(signal_num / bp1_cov, 2)
        if bp2_cov ==0:
            tra_rate2 = 0
        else:
            tra_rate2 = round(signal_num / bp2_cov, 2)
        tra_gt = determine_genotype(max(tra_rate1, tra_rate2))
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size, tra_gt, f'total_map_reads_bp1={bp1_cov};total_map_reads_bp2={bp2_cov}',f'{chrome2}:{bp2}_ratio={tra_rate2},{chrome1}:{bp2}_ratio={tra_rate2};TRA']
        return "\t".join(map(str, out))
def svGenotyper(supp_sv_table, mapf, name, outdir, min_maq, shift):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(supp_sv_table, "r") as svf:
        sv_lines = svf.readlines()

    output_lines = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        future_to_svgt = {executor.submit(process_sv, sv_line, mapf, name, min_maq, shift): sv_line for sv_line in sv_lines[1:]}
        for future in tqdm(as_completed(future_to_svgt), total=len(future_to_svgt)):
            result = future.result()
            if result:
                output_lines.append(result)

    output_file_path = f"{outdir}/2_tmp_{name}_genotype.txt"
    with open(output_file_path, "w") as output_file:  # Open file in write mode
        # Write header first
        output_file.write(header_line + "\n")
        for line in output_lines:
            output_file.write(line + "\n")

if __name__ == "__main__":
    import argparse
    from time import time
    parser = argparse.ArgumentParser("the SVGT for long reads samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-i", dest="sv_info",required=True, help="sv call set")
    IN.add_argument("-m", dest="maq", default=50,type=int, help="the mini mapping quality of the contigs, range from 45-60, default is 50")
    IN.add_argument("-s", dest="shift",default=50,type=int, help="breakpoints may shifting, here we shift 100bp as default")
    IN.add_argument("-mapf", dest="mapf",required = True, help="the mapping file of sample in sam or bam format")
    IN.add_argument("-n", dest="ACC",required = True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir",required = True, help="the output dir")
    args = parser.parse_args()
    start_t = time()
    svGenotyper(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.shift)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")
