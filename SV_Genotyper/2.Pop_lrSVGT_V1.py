from os.path import basename
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from tqdm import tqdm
from sub_lr_SVGT import delGT, insGT,dupGT, supp_dupGT, traGT,invGT, breaks2GT
import pandas as pd
import pysam
def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def fetch_sam_data_with_retry(opened_sam, reference, start, end, max_retries=25, delay=0.1):
    retries = 0
    while retries < max_retries:
        try:
            print(f"Attempting to fetch: {reference} {start}-{end}")
            return opened_sam.fetch(reference=reference, start=start, end=end)
        except BlockingIOError:
            print(f"BlockingIOError, retrying {retries + 1}/{max_retries}")
            time.sleep(delay)  # Increase the delay to 1 second
            retries += 1
        except OSError as e:
            print(f"OSError encountered: {e}")
            raise RuntimeError(f"OSError: Failed to fetch from BAM file due to corruption.")
    raise RuntimeError(f"Failed to fetch data after {max_retries} retries.")

def process_sv(sv_line, mapfile, name, min_maq, shift=100):
    opened_sam = pysam.AlignmentFile(mapfile, "rb")
    sampleID = name
    info = sv_line.strip().split("\t")
    if "DEL" in info[5]:
        chrome,sv_s, sv_e = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift,0) ## aviod the region start at position value less than shift 
        left_sam  = opened_sam.fetch(reference=chrome, start=left_most, end=sv_s + shift)
        right_sam = opened_sam.fetch(reference=chrome, start=sv_e - shift, end=sv_e + shift)
        genotype = delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, shift=200)
        out = [chrome, sv_s, sv_e, sv_size, -sv_size]
        return "\t".join(map(str, out + genotype))
    
    if "INS" in info[5]:
        chrome = info[0]
        sv_s = int(info[1])
        sv_e = int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift,0) ## aviod the region start at position value less than shift
        region_sam  = opened_sam.fetch(reference=chrome, start=left_most, end=sv_e + shift)
        genotype = insGT(sampleID, region_sam, chrome, sv_s, sv_e, sv_size, min_maq, shift=200)
        out = [chrome, sv_s, sv_e, -sv_size, sv_size]
        return "\t".join(map(str, out + genotype))

    if "DUP" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]),int(info[2])
        bp1_left = max(bp1 - shift,0) ## aviod the region start at position value less than shift
        bp1_sam = opened_sam.fetch(reference=chrome, start= bp1_left, end= bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start= bp2 - shift, end= bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, -sv_size, sv_size]
        if abs(sv_size) < 10000:
            genotype = dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP",shift=500)
        else:
            genotype = supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP",shift=500)
        return "\t".join(map(str, out + genotype))

    if "INV" in info[5]:
        chrome,bp1,bp2 = info[0], int(info[1]),int(info[2])
        sv_size = int(info[3])
        if sv_size <= 10000:
            shift = sv_size * 0.2 ## take sive into  ??
        else:
            shift = 2000
        left_most = max(bp1 - shift, 0)  ## aviod the region start at position value less than shift
        bp1_sam = opened_sam.fetch(reference=chrome, start= left_most, end= bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start= bp2 - shift, end= bp2 + shift)
        out = [chrome, bp1, bp2, sv_size, sv_size]
        genotype = invGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "INV",shift=shift)
        return "\t".join(map(str, out + genotype))

    if "TRA" in info[5]:
        svid = info[4]
        ## gemome chr should not has _,since we take _ to extract chr in svid
        chr2_s2,chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1,bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2,bp2 = chr2_s2[0], int(chr2_s2[1])
        bp1_left = max(bp1 - shift ,0) ## aviod the region start at position value less than shift
        bp2_left = max(bp2 - shift ,0) ## aviod the region start at position value less than shift
        bp1_sam = opened_sam.fetch(reference=chrome1, start= bp1_left, end= bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome2, start= bp2_left, end= bp2 + shift)
        sv_size = 0
        genotype = traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq,"TRA",shift=2000)
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size]
        return "\t".join(map(str, out + genotype))

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
    IN.add_argument("-s", dest="shift",default=100,type=int, help="breakpoints may shifting, here we shift 100bp as default")
    IN.add_argument("-mapf", dest="mapf",required = True, help="the mapping file of sample in sam or bam format")
    IN.add_argument("-n", dest="ACC",required = True, help="Accession name of the Individual")
    IN.add_argument("-o", dest="dir",required = True, help="the output dir")
    args = parser.parse_args()
    start_t = time()
    svGenotyper(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.shift)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")
