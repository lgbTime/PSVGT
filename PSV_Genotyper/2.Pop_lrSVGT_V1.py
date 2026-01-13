from os.path import basename
import subprocess
from multiprocessing import Pool
import os
from tqdm import tqdm
from sub_lr_SVGT import delGT, insGT, dupGT, supp_dupGT, traGT,breaks2traGT, invGT, breaks2invGT, little_dupGT
import pandas as pd
import pysam
from math import ceil
import time

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def process_sv(sv_line, opened_sam, name, min_maq, homo_rate, ref_rate, shift=200,minfilt=35):
    sampleID = name
    info = sv_line.strip().split("\t")
    if "DEL" in info[5]:
        chrome, sv_s, sv_e = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift, 0)
        left_sam = opened_sam.fetch(reference=chrome, start=left_most, end=sv_s + shift)
        right_sam = opened_sam.fetch(reference=chrome, start=sv_e - shift, end=sv_e + shift)
        genotype = delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift,minfilt)
        out = [chrome, sv_s, sv_e, sv_size, -sv_size]
        return "\t".join(map(str, out + genotype))

    if "INS" in info[5]:
        chrome = info[0]
        sv_s = int(info[1])
        sv_e = int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift, 0)
        region_sam = opened_sam.fetch(reference=chrome, start=left_most, end=sv_e + shift)
        genotype = insGT(sampleID, region_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift, minfilt)
        out = [chrome, sv_s, sv_e, -sv_size, sv_size]
        return "\t".join(map(str, out + genotype))

    if "DUP" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]), int(info[2])
        bp1_left = max(bp1 - shift, 0)
        bp1_sam = opened_sam.fetch(reference=chrome, start=bp1_left, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start=bp2 - shift, end=bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, -sv_size, sv_size]
        if sv_size <= 3000:
            shift = 500
            left_most = max(bp1 - shift, 0)
            region_sam = opened_sam.fetch(reference=chrome, start=left_most, end=bp2 + shift)
            genotype = little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift = shift)
        else:
            genotype = dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift = shift)
        return "\t".join(map(str, out + genotype))

    if "INV" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        if sv_size < 1000:
            shift = ceil(sv_size / 2)
        else:
            shift = 500
        left_most = max(bp1 - shift, 0)
        bp1_sam = opened_sam.fetch(reference=chrome, start=left_most, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start=bp2 - shift, end=bp2 + shift)
        out = [chrome, bp1, bp2, sv_size, sv_size]
        genotype = breaks2invGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "INV", shift=shift)
        return "\t".join(map(str, out + genotype))

    if "TRA" in info[5]:
        svid = info[4]
        #Lsat_v11_chr4	127272099	8722760	0	Lsat_v11_chr7:8722760_Lsat_v11_chr4:127272099	TRA
        chrome1, bp1 = info[0], int(info[1])
        chrome2, bp2 = svid.split(":")[0], int(info[2])
        shift = 500
        bp1_left = max(bp1 - shift, 0)
        bp2_left = max(bp2 - shift, 0)
        bp1_sam = opened_sam.fetch(reference=chrome1, start=bp1_left, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2_left, end=bp2 + shift)
        sv_size = 0
        genotype = traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, "TRA", shift=500)
        ### breakpoint to genotype TRA
        if genotype[0] == "0/1":
            bp1_sam = opened_sam.fetch(reference=chrome1, start=bp1_left, end=bp1 + shift)
            bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2_left, end=bp2 + shift)
            genotype = breaks2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, "TRA", shift=200) ### breakpoint calling
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size]
        return "\t".join(map(str, out + genotype))


def process_chromosome(args):
    chromosome_lines, mapf, name, min_maq, homo_rate, ref_rate, shift,minsv = args
    opened_sam = pysam.AlignmentFile(mapf, "rb")
    chromosome_output = []
    for line in chromosome_lines:
        result = process_sv(line, opened_sam, name, min_maq, homo_rate, ref_rate, shift, minsv)
        if result:
            chromosome_output.append(result)
    opened_sam.close()
    return chromosome_output


def svGenotyper(supp_sv_table, mapf, name, outdir, min_maq, homo_rate, ref_rate, shift, workers,minsv):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(supp_sv_table, "r") as svf:
        sv_lines = svf.readlines()

    chromosome_groups = {}
    for line in sv_lines[1:]:
        info = line.strip().split("\t")
        chromosome = info[0]
        if chromosome not in chromosome_groups:
            chromosome_groups[chromosome] = []
        chromosome_groups[chromosome].append(line)

    output_lines = []
    with Pool(processes=workers) as pool:
        args_list = [(lines, mapf, name, min_maq, homo_rate, ref_rate, shift,minsv) for _, lines in chromosome_groups.items()]
        results = list(tqdm(pool.imap(process_chromosome, args_list), total=len(args_list)))
        for result in results:
            output_lines.extend(result)

    output_file_path = f"{outdir}/2_tmp_{name}_genotype.txt"
    with open(output_file_path, "w") as output_file:
        output_file.write(header_line + "\n")
        for line in output_lines:
            output_file.write(line + "\n")


if __name__ == "__main__":
    import argparse
    from time import time
    parser = argparse.ArgumentParser("the SVGT for long reads samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-i", dest="sv_info", required=True, help="sv call set")
    IN.add_argument("-m", dest="maq", default=20, type=int, help="the mini mapping quality of the contigs, range from 20 - 60, default is 20")
    IN.add_argument("-s", dest="shift", default=150, type=int, help="breakpoints may shifting, here we shift 150bp as default")
    IN.add_argument("-mapf", dest="mapf", required=True, help="the mapping file of sample in sam or bam format")
    IN.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1")
    IN.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.05, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 5, the genotype will be 0/0")
    IN.add_argument("-n", dest="ACC", required=True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir", required=True, help="the output dir")
    IN.add_argument("-w", dest="workers", help="Number of worker processes", default=5, type=int)
    IN.add_argument("-minsv",dest="minsv",help="the minimum size of sv filter in GT", default=35,type=int)
    args = parser.parse_args()
    start_t = time()
    svGenotyper(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.lr_homo_rate, args.lr_ref_rate, args.shift, args.workers, args.minsv)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")
