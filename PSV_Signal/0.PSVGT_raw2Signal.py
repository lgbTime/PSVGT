import os
import subprocess
import argparse
from os.path import basename
from concurrent.futures import ThreadPoolExecutor
from time import time
sub_scripts_folder = os.path.dirname(os.path.abspath(__file__))

def run_command(cmd):
    """Run a shell command."""
    subprocess.run(cmd, shell=True, check=True)

def check_dir(directory):
    """Check if the output directory exists; if not, create it."""
    if not os.path.isdir(directory):
        os.mkdir(directory)
def construct_minimap2_cmd(ref, sample, out, dtype, cpu, secondary=False):
    secondary_option = "--secondary=no" if not secondary else ""
    if dtype in ["hifi", "pb", "ont"]:
         mapping_type = f"map-{dtype}"  # Long-read mapping types
    elif dtype in ['sr']:
        mapping_type = f"map-hifi"
    elif dtype in ["cr"]: ## since sr data will  assemble to contigs
        mapping_type = "map-hifi"
        #mapping_type = "minimap2 -a -xasm5 --cs -t 10 "
    else:
        raise ValueError("Invalid dtype provided. Choose from ['hifi', 'pb', 'ont', 'sr', 'cr'].")
    if dtype in ['cr']:
        cmd = f"minimap2 -ax {mapping_type} -t {cpu} {secondary_option} -o {out}.sam {ref} {sample}"
    else:
        cmd = f"minimap2 -ax {mapping_type} -t {cpu} {secondary_option} -o {out}.sam {ref} {sample}"
    return cmd

def svsignal(args):
    sample = args.sample
    check_dir(args.outdir)
    outSV = "0_tmp_" + basename(sample)
    ref = args.ref
    minimapCMD = construct_minimap2_cmd(args.ref, args.sample, f"{args.outdir}/{outSV}", args.dtype, args.minimapCPU, secondary=False)
    sort_cmd = f"samtools sort {args.outdir}/{outSV}.sam -@ {args.minimapCPU} -o {args.outdir}/{outSV}.bam && samtools index {args.outdir}/{outSV}.bam -@ {args.minimapCPU} && rm {args.outdir}/{outSV}.sam"
    cmd1 = minimapCMD + "&&" + sort_cmd
    print(cmd1)
    run_command(cmd1)
    start_t = time()
    if args.msv == "yes" and args.dtype in ["hifi", "ont", "pb"]:
        cmd2 = f"python {sub_scripts_folder}/0.Signal4bam_PSVGT.py -b {args.outdir}/{outSV}.bam -o {args.outdir}/{outSV} -m {args.min} -maq {args.maq} -dtype {args.dtype} -msv yes -fai {args.ref}.fai"
    else:
        cmd2 = f"python {sub_scripts_folder}/0.Signal4bam_PSVGT.py -b {args.outdir}/{outSV}.bam -o {args.outdir}/{outSV} -m {args.min} -maq {args.maq} -dtype {args.dtype} -msv no -fai {args.ref}.fai "
    print(f'********************* runing cmd ************************\n{cmd2}')
    run_command(cmd2)
    end_t = time()
    return outSV, f"{args.outdir}/{outSV}.bam"

if __name__ == "__main__":
    from time import time
    parser = argparse.ArgumentParser("Individual SVInDel Searching", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file")
    IN.add_argument("-i", dest="sample",required=True, help="Individual should be a contig fasta file or long sequencing data")
    IN.add_argument("-dtype", dest="dtype",required=True, help="sample type in sr/hifi/ont/pb/cr")
    IN.add_argument("-r", dest="ref", required=True, help="The reference genome fasta file")
    IN.add_argument("-m", dest="min", help="The min length of an SVInDel", default=40, type=int)
    IN.add_argument("-maq", dest='maq', help="The min mapping quality", default=50,type=int)
    IN.add_argument("-o", dest="outdir", help="The output directory", default="./")
    IN.add_argument("-minimapCPU", dest="minimapCPU",default=20,type=int, help="The CPU use in minimap mapping")
    IN.add_argument("-msv", dest="msv",default="no",type=str, help="collect supplementary alignment reads to dectect complex sv(TRA, INV, DUP, INS, DEL) for hifi/ont/CLR data")
    
    args = parser.parse_args()
    start = time()
    outInDel, ibam = svsignal(args)
    end = time()
    print(f"########################## Time Cost in SVInDel Indentification: {end - start} #####################################")
