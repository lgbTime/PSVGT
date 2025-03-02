from os.path import basename
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from tqdm import tqdm
import pandas as pd
#import pysam
from math import ceil
def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def determine_genotype(entry_ratio):
    """
    ONT easy lead to 0/1 and FP, here we try modify.
    """
    if entry_ratio > 0.65:
        return "1/1"
    elif entry_ratio <  0.1:
        return "0/0"
    else:
        return "0/1"

def process_sv(sv_line, name):
    sampleID = name
    info = sv_line.strip().split("\t")
    sv_size, signal_num, signal_rate = int(info[3]),int(info[7]), float(info[8])
    left_cov,right_cov, cov = ceil(signal_num / signal_rate),  ceil(signal_num / signal_rate),  ceil(signal_num / signal_rate),
    if "TRA" not in info[5]:
        chrome, bp1, bp2= info[0], int(info[1]), int(info[2])
        if "DEL" in info[5]:
            del_gt =  determine_genotype(signal_rate)
            out = [chrome, bp1, bp2, sv_size, -sv_size, del_gt,f'total_map_reads_l={left_cov};total_map_reads_r={right_cov}', f'deles_l_ratio={round(signal_rate,2)},deles_r_ratio={round(signal_rate,2)};DEL' ]
            return "\t".join(map(str,out))
        
        elif "INS" in info[5]:
            ins_gt = determine_genotype(signal_rate)
            out = [chrome, bp1, bp2, -sv_size, sv_size, ins_gt, f'total_map_reads={cov}', f'ins_ratio={round(signal_rate,2)};INS']
            return "\t".join(map(str, out))

        elif "DUP" in info[5]:
            dup_gt = determine_genotype(signal_rate)
            out = [chrome, bp1, bp2, -sv_size, sv_size, dup_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'bp1={chrome}:{bp1},bp1_ratio={signal_rate},bp2={chrome}:{bp2},bp2_ratio={signal_rate};DUP']
            return "\t".join(map(str, out))

        elif "INV" in info[5]:
            inv_gt = determine_genotype(signal_rate)
            out = [chrome, bp1, bp2, -sv_size, sv_size, inv_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'bp1={chrome}:{bp1},bp1_ratio={signal_rate},bp2={chrome}:{bp2},bp2_ratio={signal_rate};INV']
            return "\t".join(map(str, out))
    elif "TRA" in info[5]:
        svid = info[4]
        ## gemome chr should not has _,since we take _ to extract chr in svid
        chr2_s2, chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1, bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2, bp2 = chr2_s2[0], int(chr2_s2[1])
        bp1_cov  = cov
        bp2_cov =  cov
        tra_gt = determine_genotype(signal_rate)
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size, tra_gt, f'total_map_reads_bp1={bp1_cov};total_map_reads_bp2={bp2_cov}',f'bp2={chrome2}:{bp2},bp2_ratio={signal_rate},bp1={chrome1}:{bp2},bp1_ratio={signal_rate};TRA']
        return "\t".join(map(str, out))
def svGenotyper(supp_sv_table, name, outdir):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(supp_sv_table, "r") as svf:
        sv_lines = svf.readlines()

    output_lines = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        future_to_svgt = {executor.submit(process_sv, sv_line, name): sv_line for sv_line in sv_lines[1:]}
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
    IN.add_argument("-n", dest="ACC",required = True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir",required = True, help="the output dir")
    args = parser.parse_args()
    start_t = time()
    svGenotyper(args.sv_info, args.ACC, args.dir)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")
