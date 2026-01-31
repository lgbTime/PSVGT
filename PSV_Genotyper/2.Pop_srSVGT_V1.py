from os import system
from os.path import basename
import subprocess
import os
from tqdm import tqdm
from sub_srSVGT import delGT, insGT, traGT, breaks2GT
import pandas as pd
import pysam
import multiprocessing


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)


def split_region(start, end):
    step = (end - start) / 7
    points = [start + i * step for i in range(8)]
    del points[0]
    del points[-1]
    return points


def usepysam_count_localreads(opened_sam, chrom, start, end, maq):
    reads = opened_sam.fetch(reference=chrom, start=start, end=end)
    high_quality_reads = [read for read in reads if read.mapping_quality >= maq]
    count = len(high_quality_reads)
    return count


def five_local_reads_count4SVDel(chrom, start, end, min_maq, opened_sam):
    reads_counts = []
    points = split_region(start, end)
    for point in points:
        start, end = point, point + 7
        count = usepysam_count_localreads(opened_sam, chrom, start, end, min_maq)
        reads_counts.append(count)
    return reads_counts


def process_sv(sv_line, opened_sam, name, min_maq, homo_rate, ref_rate, shift, span):
    sampleID = name
    info = sv_line.strip().split("\t")
    if "DEL" in info[5]:
        chrome, sv_s, sv_e = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        left_sam = opened_sam.fetch(reference=chrome, start=sv_s - shift, end=sv_s + shift)
        right_sam = opened_sam.fetch(reference=chrome, start=sv_e - shift, end=sv_e + shift)
        genotype = delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift, span)
        if genotype[0] != "1/1":
            local_5counts = five_local_reads_count4SVDel(chrome, sv_s, sv_e, min_maq, opened_sam)
            if local_5counts.count(0) >= 2:
                genotype[0] = "1/1"
                print(f"********* DEL {chrome}:{sv_s}-{sv_e} genotype to 1/1, since at least two local region count 0 reads **********")
        out = [chrome, sv_s, sv_e, sv_size, -sv_size]
        return "\t".join(map(str, out + genotype))
    if "INS" in info[5]:
        chrome = info[0]
        sv_s = int(info[1])
        sv_e = int(info[2])
        sv_size = int(info[3])
        region_sam = opened_sam.fetch(reference=chrome, start=sv_s - shift, end=sv_e + shift)
        genotype = insGT(sampleID, region_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift, span)
        out = [chrome, sv_s, sv_e, -sv_size, sv_size]
        return "\t".join(map(str, out + genotype))
    if "TRA" in info[5]:
        svid = info[4]
        chrome1, bp1 = info[0], int(info[1])
        chrome2, bp2 = svid.split(":")[0], int(info[2])
        bp1_left = max(bp1 - shift, 0)
        bp2_left = max(bp2 - shift, 0)
        bp1_sam = opened_sam.fetch(reference=chrome1, start=bp1_left, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2_left, end=bp2 + shift)
        sv_size = 0
        genotype = traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, "TRA", shift)
        if genotype[0] != "1/1":
            bp1_sam = opened_sam.fetch(reference=chrome1, start=bp1 - shift, end=bp1 + shift)
            bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2 - shift, end=bp2 + shift)
            genotype = breaks2GT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, "TRA", shift)
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size]
        return "\t".join(map(str, out + genotype))
    if "DUP" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]), int(info[2])
        bp1_sam = opened_sam.fetch(reference=chrome, start=bp1 - shift, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start=bp2 - shift, end=bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, -sv_size, sv_size]
        genotype = breaks2GT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift)
        return "\t".join(map(str, out + genotype))
    if "INV" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]), int(info[2])
        bp1_sam = opened_sam.fetch(reference=chrome, start=bp1 - shift, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start=bp2 - shift, end=bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, sv_size, sv_size]
        genotype = breaks2GT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "INV", shift)
        return "\t".join(map(str, out + genotype))


def process_chromosome(chromosome_lines, mapf_path, name, min_maq, homo_rate, ref_rate, shift, span):
    output_lines = []
    try:
        opened_sam = pysam.AlignmentFile(mapf_path)
        for sv_line in chromosome_lines:
            result = process_sv(sv_line, opened_sam, name, min_maq, homo_rate, ref_rate, shift, span)
            output_lines.append(result)
    except Exception as e:
        print(f"Error processing chromosome: {e}")
    finally:
        if 'opened_sam' in locals():
            opened_sam.close()
    return output_lines


def srSVGT(sv_table, mapf, name, outdir, min_maq, homo_rate, ref_rate, shift, span, cpu):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(sv_table, "r") as svf:
        sv_lines = svf.readlines()

    chromosome_groups = {}
    for sv_line in sv_lines[1:]:
        chrom = sv_line.strip().split("\t")[0]
        if chrom not in chromosome_groups:
            chromosome_groups[chrom] = []
        chromosome_groups[chrom].append(sv_line)
    tasks = [(chromosome_lines, mapf, name, min_maq, homo_rate, ref_rate, shift, span) for chromosome_lines in chromosome_groups.values()]

    nproc = min(cpu,len(tasks))
    with multiprocessing.Pool(processes=nproc) as pool:
        results = list(tqdm(pool.starmap(process_chromosome, tasks), total=len(tasks)))

    all_output_lines = [line for sublist in results for line in sublist]
    output_file_path = f"{outdir}/2_tmp_{name}_bpgenotype.txt"
    with open(output_file_path, "w") as output_file:
        output_file.write(header_line + "\n")
        for line in set(all_output_lines):
            output_file.write(line + "\n")


if __name__ == "__main__":
    import argparse
    from time import time

    parser = argparse.ArgumentParser("srGT Genotyping", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-i", dest="sv_info", required=True,
                    help="SV record table can be generated by PSVGT step 1")
    IN.add_argument("-maq", dest="maq", default=1, type=int,
                    help="the mini mapping quality of the contigs, range from 0 to 60, default is 45")
    IN.add_argument("-s", dest="shift", default=30, type=int,
                    help="breakpoints may shifting, if the shifting range is expand, the clip reads ratio to call genotype should be smaller, please check sub_svGT ratio to determine the genotype")
    IN.add_argument("-span", dest="span", default=50, type=int,
                    help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")
    IN.add_argument("-ref_rate", dest="ref_rate", default=0.05, type=float,
                    help="if local breakpoints/maps <= 0.05, the genotype should will 0/0, you can increaing the float to filter putative heterzygous SV")
    IN.add_argument("-homo_rate", dest="homo_rate", default=0.65, type=float,                    
                    help="if local breakpoints/maps >= 0.65, the genotype should will 1/1, base on the known heterozygous rate of the species, you can lower down the value")

    IN.add_argument("-mapf", dest="mapf", required=True, help="the mapping file in bam format")
    IN.add_argument("-n", dest="ACC", required=True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir", required=True, help="the output dir")
    IN.add_argument("-srcpu", dest="cpu", type=int, default=10,help="Max CPUs to use in genotyping")
    args = parser.parse_args()
    start_t = time()
    srSVGT(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.homo_rate, args.ref_rate, args.shift, args.span, args.cpu)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")
