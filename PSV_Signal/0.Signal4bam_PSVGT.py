import argparse
from time import time
import pandas as pd
from multiprocessing import Pool
import sub_Signal4bam_PSVGT
def parse_arguments():
    parser = argparse.ArgumentParser("SV signal extract from sam file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Signal Capture")
    IN.add_argument("-b", dest="bam", required=True, help="sorted bam and index bam only")
    IN.add_argument("-o", dest="out", required=True, help="Output SV signal info file chromsomely")
    IN.add_argument("-m", dest="min", help="SV min length", default=40, type=int)
    IN.add_argument("-maq", dest="maq", help="The min mapping quality", default=50, type=int)
    IN.add_argument("-dtype", dest="dtype", required=True, help="The sequence type (ont, hifi, pb, cr, sr)")
    IN.add_argument("-M", dest="max", help="SV max length", default=10000000, type=int)
    IN.add_argument("-fai", dest="fai", help="Chromosome fai index file", type=str)
    IN.add_argument("-msv", dest="msv", help="Detecting complex SVs (INV, DP, TRA, INS, DEL) from supplementary alignment", default="no", type=str)
    return parser.parse_args()

def main():
    args = parse_arguments()
    # Load chromosome list from the fai file
    fai = pd.read_csv(args.fai, header=None, sep="\t",dtype=str, index_col=None)
    chromosome_list =  fai[0].tolist()
    chrom_size_dict = dict(zip(fai[0], fai[1].astype(int)))
    if args.out[-4:] == ".bam":
        args.out = args.out.replace('.bam', '')
    SVsignal_out_path = f"{args.out}"
    # Use multiprocessing to parallelize the chromosome processing
    with Pool() as pool:
        pool.starmap(sub_Signal4bam_PSVGT.process_chromosome, 
                     [(chromosome,chrom_size_dict[chromosome],chromosome_list, args.bam, args.min,args.max, args.maq, SVsignal_out_path,args.dtype,args.msv) 
                      for chromosome in chromosome_list])
    exe_end = time()
    print(f"{'*' * 40} done SV searching {'*' * 40}\ncost time: {exe_end - exe_start}")

if __name__ == "__main__":
    exe_start = time()
    main()
