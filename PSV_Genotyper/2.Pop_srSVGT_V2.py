from os import system
from os.path import basename
import subprocess
import os
from tqdm import tqdm
from sub_srSVGT import delGT, insGT, breaks2GT
import pandas as pd
import pysam
import multiprocessing

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def split_region(start, end):
    # Calculate the step size for dividing the region into 5 equal parts
    step = (end - start) / 5
    # Generate the five split points
    points = [start + i * step for i in range(6)]  # 6 points (start + 5 splits)
    del points[0]
    del points[-1]
    return points

def usepysam_count_localreads(opened_sam, chrom, start, end, maq):
    reads = opened_sam.fetch(reference=chrom, start=start, end=end)
    high_quality_reads = [read for read in reads if read.mapping_quality >= maq]
    count = len(high_quality_reads)
    return count

def four_local_reads_count4SVDel(chrom, start, end, min_maq, opened_sam):
    """
    to fix the big shift of DEL points,
    here five space of the region will be cutted,
    reads will be count within the five position.
    """
    reads_counts = []
    points = split_region(start, end)
    for point in points:
        start, end = point, point + 9
        count = usepysam_count_localreads(opened_sam, chrom, start, end, min_maq)
        reads_counts.append(count)
    reads_counts.sort()
    return reads_counts

def process_sv(sv_line, opened_sam, name, min_maq, shift):
    sampleID = name
    info = sv_line.strip().split("\t")
    if "DEL" in info[5]:
        chrome, sv_s, sv_e = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        left_sam = opened_sam.fetch(reference=chrome, start=sv_s - shift, end=sv_s + shift)
        right_sam = opened_sam.fetch(reference=chrome, start=sv_e - shift, end=sv_e + shift)
        genotype = delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, shift)
        if genotype[0] != "1/1":
            local_4counts = four_local_reads_count4SVDel(chrome, sv_s, sv_e, min_maq, opened_sam)
            if local_4counts[1] == 0:  ##  3 inter  regions count 0 reads
                genotype[0] = "1/1"
        out = [chrome, sv_s, sv_e, sv_size, -sv_size]
        return "\t".join(map(str, out + genotype))
    if "INS" in info[5]:
        chrome = info[0]
        sv_s = int(info[1])
        sv_e = int(info[2])
        sv_size = int(info[3])
        region_sam = opened_sam.fetch(reference=chrome, start=sv_s - shift, end=sv_e + shift)
        genotype = insGT(sampleID, region_sam, chrome, sv_s, sv_e, sv_size, min_maq, shift)
        out = [chrome, sv_s, sv_e, -sv_size, sv_size]
        return "\t".join(map(str, out + genotype))

    if "TRA" in info[5]:
        svid = info[4]
        ## gemome chr should not has _,since we take _ to extract chr in svid
        chr2_s2, chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1, bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2, bp2 = chr2_s2[0], int(chr2_s2[1])
        bp1_sam = opened_sam.fetch(reference=chrome1, start=bp1 - shift, end=bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome2, start=bp2 - shift, end=bp2 + shift)
        sv_size = 0
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

def process_chromosome(chromosome_lines, opened_sam, name, min_maq, shift):
    output_lines = []
    for sv_line in chromosome_lines:
        result = process_sv(sv_line, opened_sam, name, min_maq, shift)
        output_lines.append(result)
    return output_lines

def srSVGT(sv_table, mapf, name, outdir, min_maq, shift):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(sv_table, "r") as svf:
        sv_lines = svf.readlines()

    # Group SV lines by chromosome
    chromosome_groups = {}
    for sv_line in sv_lines[1:]:
        info = sv_line.strip().split("\t")
        chromosome = info[0]
        if chromosome not in chromosome_groups:
            chromosome_groups[chromosome] = []
        chromosome_groups[chromosome].append(sv_line)

    # Create a pool of processes
    pool = multiprocessing.Pool()
    results = []
    for chromosome, lines in chromosome_groups.items():
        result = pool.apply_async(process_chromosome, args=(lines, mapf, name, min_maq, shift))
        results.append(result)

    # Collect the results
    all_output_lines = []
    for result in results:
        all_output_lines.extend(result.get())

    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()

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
                    help="InDel record table can be generated by PSVGT step 1")
    IN.add_argument("-m", dest="maq", default=45, type=int,
                    help="the mini mapping quality of the contigs, range from 0-60, default is 45")
    IN.add_argument("-s", dest="shift", default=30, type=int,
                    help="breakpoints may shifting, if the shifting range is expand, the clip reads ratio to call genotype should be smaller, please check sub_svGT ratio to determine the genotype")
    IN.add_argument("-mapf", dest="mapf", required=True, help="the mapping file in bam format")
    IN.add_argument("-n", dest="ACC", required=True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir", required=True, help="the output dir")
    args = parser.parse_args()
    start_t = time()
    opened_sam = pysam.AlignmentFile(args.mapf)
    srSVGT(args.sv_info, opened_sam, args.ACC, args.dir, args.maq, args.shift)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")

