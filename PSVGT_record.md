## ./
```shell
total 316
-rwxr-xr-x 1 lgb xinwang 38580 Jul 13 17:30 0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py
-rwxr-xr-x 1 lgb xinwang 37677 May 25 15:36 backup_20250525_PSVGT1.0.py
-rwxr-xr-x 1 lgb xinwang 21122 May 22 16:55 back_v1_PSVGT1.0.py
-rwxr-xr-x 1 lgb xinwang 22654 May 22 16:55 back_v2_PSVGT1.0.py
drwxr-xr-x 2 lgb xinwang  4096 May 22 16:55 CapsPop
drwxr-xr-x 2 lgb xinwang  4096 May 30 16:14 img
-rwxr-xr-x 1 lgb xinwang  2707 Jul 12 17:13 merge_close_SV.py
drwxr-xr-x 3 lgb xinwang  4096 Sep 12 15:09 PSV_Genotyper
-rwxr-xr-x 1 lgb xinwang 33891 Jul 23 09:11 PSVGT1.0_Dynamic_Window.py
-rwxr-xr-x 1 lgb xinwang 33739 Jul 23 09:05 PSVGT1.0_Fix_Window.py
-rwxr-xr-x 1 lgb xinwang 35062 Sep 15 10:59 PSVGT1.0.py
-rw-r--r-- 1 lgb xinwang     0 Sep 16 07:33 PSVGT_record.md
drwxr-xr-x 4 lgb xinwang  4096 Sep 15 09:51 PSV_Signal
drwxr-xr-x 2 lgb xinwang  4096 May 22 16:55 __pycache__
-rw-r--r-- 1 lgb xinwang 15859 Jul  4 22:39 README.md
-rwxr-xr-x 1 lgb xinwang 21318 May 22 16:55 SSVGT1.0.py
drwxr-xr-x 2 lgb xinwang  4096 Aug 26 12:27 SVInDel_Anno
drwxr-xr-x 2 lgb xinwang  4096 May 22 17:37 SVInDel_Imputation
drwxr-xr-x 2 lgb xinwang  4096 Jul  4 11:50 SVInDel_Primer
drwxr-xr-x 2 lgb xinwang  4096 Aug 14 17:26 utils

```
### ./merge_close_SV.py
```python
import pandas as pd
import sys
import argparse


def merge_svs(input_file, output_file, distant=800):
    try:
        df = pd.read_csv(input_file, sep='\t')
        df['SVType'] = df['SV_support'].str.split(';').str[-1]
        df = df.sort_values(by=['SVType', 'Target_start'])

        merged_rows = []
        i = 0
        while i < len(df):
            current_row = df.iloc[i]
            current_chrom = current_row['#Target_name']
            current_start = current_row['Target_start']
            current_size = current_row['Target_size']
            current_svtype = current_row['SVType']

            j = i + 1
            while j < len(df):
                next_row = df.iloc[j]
                next_chrom = next_row['#Target_name']
                next_start = next_row['Target_start']
                next_size = next_row['Target_size']
                next_svtype = next_row['SVType']
                # Use the 'distant' parameter instead of hardcoded 800
                if (next_chrom == current_chrom and 
                    next_svtype == current_svtype and 
                    abs(next_start - current_start) < distant):
                    current_size += next_size
                    j += 1
                else:
                    break

            current_row['Target_size'] = current_size
            merged_rows.append(current_row)
            i = j

        merged_df = pd.DataFrame(merged_rows)
        merged_df.drop(columns=['SVType'], inplace=True)
        merged_df.to_csv(output_file, sep='\t', na_rep='nan', index=False)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # Initialize argument parser with a description
    parser = argparse.ArgumentParser(
        description='Merge structural variants (SVs) based on proximity and type.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # Shows default values in help
    )
    
    # Required positional arguments
    parser.add_argument('input_file', 
                        help='Path to the input TSV file containing SV data')
    parser.add_argument('output_file', 
                        help='Path to save the merged output TSV file')
    
    # Optional argument with default value
    parser.add_argument('--distant', type=int, default=800,
                        help='Maximum distance (bp) between SVs to be merged')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the merge function with parsed arguments
    merge_svs(args.input_file, args.output_file, args.distant)

```

### ./PSVGT1.0.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
import multiprocessing
import concurrent
from functools import partial
from os.path import basename
from concurrent.futures import ThreadPoolExecutor, as_completed

PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from Sub_readfa2Dict import readfa2Dict

def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("PSVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-fix_window", "--fix_window", default="no", help="PSVGT allow clustering signal in both flexible(size adaptive) and fix window sliding to cluster signal, for case of haplotype genomes analysis, we suggest setting to yes, while long reads sequencing datasets we suggest no")
    parser.add_argument("-win", "--window",default=500,type=int,help="Window size to parse signal, this parameter is to merge fragment SV due to aligment, default fix window size is 500bp to merge cluster the SV with a sv_start less than window size")
    parser.add_argument("-diploid", "--diploid", help="for diploid resolved assembly, to get a phased genotype please provide table list each line in format hap1\thap2\tSampleName in cr directory")
    parser.add_argument("-polyploid", "--polyploid", help="for polyploid haplotype resolved assembly like potato(4 haplotype assemblies available), to get a merge genotype of samples please provide table list each line in format hap1\thap2\thap3\thapn\tSampleName")


    parser.add_argument("-w", "--max_workers", default = 6,type=int, help="the max workers thread pool excutor, 6 means run 6 samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="no", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=30,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.10, type=float, help= "the percent of reads that support a candidate SV (0.10 means at a depth 20X region, a SV signal should have at least 2 reads support, this parameter is for the variaty depth of hifi/ont/pb samples")
    parser.add_argument("-nreads",  "--nreads", type=int, help= "the number of reads to support a candidate SV (SV signal should have at least numbers reads support, this parameter is for the various depth of hifi/ont/pb samples")
    parser.add_argument("--num_hap", "--num_hap", default=2, type=int, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    parser.add_argument("-msv","--msv_mode",default="no", help= "Setting msv to `yes` signals of INS,DEL,INV,DUP,TRA will captured, while setting to `no` Only insertions and deletions will be capture from alignment ")
    parser.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1, for species like polyploid potato we suggest 0.8")
    parser.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 0.10, the genotype will be 0/0")
    parser.add_argument("-sr_homo_rate", dest="sr_homo_rate",default=0.65, type=float, help="to determine a homozygous site, if 0.65 of the local mapping signal suport the sv the genotyoe will be 1/1, you can lower down the value if your specise have a low heterozygous rate")
    parser.add_argument("-sr_ref_rate", dest="sr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 10 percent, the genotype will be 0/0, you can increase the value to filter putative heterozygous SV")

    parser.add_argument("-span", dest="span", default=50, type=int, help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")

    args = parser.parse_args()
    if args.fix_window == "yes":
        KLOOK = "0.KLookCluster_LocalDepthAdaptive.py"
    else:
        KLOOK = "0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py"
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    

    done_analysor = file_capture(args.outdir, ".record.txt")
    def run_clu2fil_cmd(clu2fil_cmd):
        try:
            subprocess.run(clu2fil_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {clu2fil_cmd}\n{e}")
    def process_single_contig(contig, dtype, args, PSVGT, fai):
        try:
            ## final chromsome result ##
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"
            ## skip samples ##
            if os.path.exists(done_name):
                print(f"Skipping processed contig: {basename(contig)}")
                return

            ## minimaping and signal detect ##
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py '
                f'-i {contig} -dtype {dtype} -r {args.refGenome} '
                f'-m {args.min} -maq {maq} -o {args.outdir} '
                f'-minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            ## chromosome pool run ##
            with multiprocessing.Pool(processes=args.max_workers) as pool:
                if dtype in ['hifi','ont','pb'] and args.csv and args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/{KLOOK} '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} --nreads {args.nreads}'
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                
                elif dtype in ['hifi','ont','pb'] and args.csv and not args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/{KLOOK} '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                elif dtype in ['hifi','ont','pb'] and args.nreads and not args.csv:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/{KLOOK} '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--nreads {args.nreads} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                else:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/{KLOOK} '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]

                pool.map(run_clu2fil_cmd, clu2fil_cmds)
            ## ACC SV ##
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/0_tmp_{basename(contig)} '
                f'-fai {args.refGenome}.fai -M {args.max} --nrate {args.csv}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(contig)}: {str(e)}")
    def add_commands4fq(files, dtype):
        max_workers = min(len(files), args.max_workers)
        if max_workers >0:
            print(f'{max_workers} workers')
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        process_single_contig,
                        contig, dtype, args, PSVGT, fai
                    ) for contig in files
                ]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Thread error: {str(e)}")

    def process_single_bam(bam, dtype, args, PSVGT, fai):
        try:## last chromosome result ##
            done_flag = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if os.path.exists(done_flag):
                print(f"Skipping processed BAM: {basename(bam)}")
                return

            # Phase 1: signal
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py '
                f'-b {bam} -o {args.outdir}/{basename(bam)} '
                f'-m {args.min} -maq {maq} -dtype {dtype} '
                f'-msv {args.msv_mode} -fai {args.refGenome}.fai'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            # Phase 2: parafly chrs
            base_prefix = basename(bam).replace(".bam", "")
            if args.nreads and args.csv and dtype in ['hifi', 'ont', 'pb']:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/{KLOOK} '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--nreads {args.nreads} --rate_depth {args.csv} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.csv and dtype in ['hifi', 'ont', 'pb'] and not args.nreads:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/{KLOOK} '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--rate_depth {args.csv} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.nreads and dtype in ['hifi', 'ont', 'pb'] and not args.csv:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/{KLOOK} '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--nreads {args.nreads}  '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            else:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/{KLOOK} '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]


            # dynamic cpu
            with multiprocessing.Pool(processes=min(len(fai[0]), os.cpu_count()//5)) as pool:
                pool.map(run_clu2fil_cmd, chrom_commands)

            # Phase 3: merge signal
            nhap ={1:0.99, 2:0.49, 3:0.33, 4:0.24, 5:0.19, 6:0.16, 7:0.14, 8:0.12}
            ## for haplotype genomes
            if dtype in ['cr']:
                if args.diploid or args.polyploid: 
                    nrate = nhap[1] ## for haplotype genome
                else:
                    nrate = nhap[args.num_hap] ## for all haplotypes in one fasta genome
                ACC_SV_cmd = (
                    f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                    f'-preffix {args.outdir}/{base_prefix} '
                    f'-fai {args.refGenome}.fai --nrate {nrate} --minimaq {args.maq}'
                )
            elif dtype in ['sr']: ## currently sr assemble can only capture one haplotype
                nrate = nhap[1]
                ACC_SV_cmd = (
                    f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                    f'-preffix {args.outdir}/{base_prefix} '
                    f'-fai {args.refGenome}.fai --nrate {nrate} --minimaq {args.maq}'
                )
            else:
                ACC_SV_cmd = (
                    f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                    f'-preffix {args.outdir}/{base_prefix} '
                    f'-fai {args.refGenome}.fai --nrate {args.csv} --minimaq {args.maq}'
                )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(bam)}: {str(e)}")

    def add_commands4bam(files, dtype):
        if len(files)>0:
            max_threads = min(len(files),args.max_workers)
            with ThreadPoolExecutor(max_workers=max_threads) as executor:
                task_func = partial(
                    process_single_bam,
                    dtype=dtype,
                    args=args,
                    PSVGT=PSVGT,
                    fai=fai
                )
                futures = {executor.submit(task_func, bam): bam for bam in files}
                for future in concurrent.futures.as_completed(futures):
                    bam_file = futures[future]
                    try:
                        future.result()
                        print(f"Success: {basename(bam_file)}")
                    except Exception as e:
                        print(f"Failed processing {basename(bam_file)}: {str(e)}")

    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        add_commands4fq(file_capture(args.crdir, ".gz"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    ## step1 no-redundant population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file} -m {args.maq} -lr_homo_rate {args.lr_homo_rate} -lr_ref_rate {args.lr_ref_rate}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf {args.refGenome}.fai"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")

    ###################### haplotype resoved assembly genotype phased ###################
    if args.diploid:
        print("*************** diploid calling **************")
        with open(args.diploid, 'r') as f:
            lines = f.readlines()
        for hh in lines:
            h1 = hh.strip().split("\t")[0]
            h1 = "2_tmp_" + h1.replace(".bam", "") + "_genotype.txt" if  h1[-4:] == ".bam" else "2_tmp_" + h1 + "_genotype.txt"
            h2 = hh.strip().split("\t")[1]
            h2 = "2_tmp_" + h2.replace(".bam", "") + "_genotype.txt" if  h2[-4:] == ".bam" else "2_tmp_" + h2 + "_genotype.txt"
            samplename =  hh.strip().split("\t")[2]
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_diploid_asm.py {args.outdir}/{h1} {args.outdir}/{h2} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f"***************** try to phased hap1: {args.outdir}/{h1} and hap2: {args.outdir}/{h2} to {args.outdir}/2_tmp_{samplename}_genotype.txt ******************")
            tab2vcf =    f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(tab2vcf)
    if args.polyploid:
        print("*************** polyploid genotype merging **************")
        with open(args.polyploid, 'r') as f:
            lines = f.readlines()
        path_haps = ""
        for hh in lines:
            path_haps = ""
            haps = hh.strip().split("\t")
            total_haps = len(haps) - 1
            samplename = haps[-1]
            for i in range(total_haps):
                hap = haps[i] 
                hap_file = "2_tmp_" + hap.replace(".bam", "") + "_genotype.txt" if  hap[-4:] == ".bam" else "2_tmp_" + hap + "_genotype.txt"
                path_haps += f'{args.outdir}/{hap_file} '
            print(f"*************************** phased polyploid {samplename} ***************************** ")
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_polyploid_genome_gt.py {path_haps} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f'********************** phased polyploid ************************\n{phased_cmd}')
            vcf_cmd = f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(vcf_cmd)


    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {bam} -maq {args.maq} -span {args.span} -s 50 -n {sampleID} -o {args.outdir} -homo_rate {args.sr_homo_rate} -ref_rate {args.sr_ref_rate} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf {args.refGenome}.fai"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), args.max_workers)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 80 --max 600 --frank 400 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 80 600 400 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py
```python
import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil,floor

def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cost time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        try:
            num_maps = covinfo.count(str(reference), start, end)
        except AttributeError:
            print(f"Error: covinfo does not have a 'count' method.")
            num_maps = 0
    return num_maps

def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)


def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_mode = subset.mode()
        if not subset_mode.empty:
            return subset_mode.iloc[0]  
        else:
            return ceil(subset.mean())
    else:
        # if subset empty return original mode
        original_mode = series.mode()
        if not original_mode.empty:
            return original_mode.iloc[0]
        else:
            # no mode, return mean
            return ceil(series.mean())

def candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate=0.1,add=1):
    """
    Solve clusters dataframe
    Get each cluster's sv breakpoints, SV length
    For low depth clus, the cluster must be vary strict and intense, else it will generate a lot of false positive sv. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts() ## default reversed count
    min_clusters = min(num_hap, len(cluster_counts))
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.05:
        print("proportion is too low, return []")
        return []
    else: 
        print(f"top1 cluster percent is {proportion}")
    
    top_clusters = cluster_counts.head(min_clusters).index
    clusdf = clusdf[clusdf[cluster_col].isin(top_clusters)]
    clusters_to_process =  clusdf[cluster_col].unique()

    msv = []
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = mode_or_median(clu_df['Target_end'])
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        
        if sv_eye < 2:
            continue
        if svtype in ["INS", "DUP"]:
            break_cov = local_cov(opened_bam, sv_chrom, sv_start-70, sv_start-50) + 1
        elif svtype in ["DEL", "INV"]:
            midpoint = int(sv_start + sv_len*0.5)
            break_cov = local_cov(opened_bam, sv_chrom, midpoint, midpoint+20) + 1
        else:
            print("!!!!!!!!!!!!!!!!!!!!! SV type error !!!!!!!!!!!!!!!!!!!!!!!!")

        SV_rate = round(sv_eye / break_cov, 2)
        if sv_eye >= break_cov*support_rate + add:
            print(f"{svid} sv signal propertion support")
            clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
            msv.append(svid)
        if sv_len > 3000 and svtype=="INS" and sv_eye >= 2 and sv_eye >= break_cov*support_rate*0.6: ## to capture big INS
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} big INS,breakpoint coverage is {break_cov} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        if sv_len > 500 and svtype=="DUP" and sv_eye >= 2 and sv_eye >= break_cov*support_rate*0.5: ## to capture more DUP
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} DUP,breakpoint coverage is {break_cov} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        
        if sv_eye >= nreads -1: ## filter 
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid},but local coverage is {break_cov}, breakdepth not support this sv")
    if len(set(msv)) >= 2:
        outmsv = "\t".join(list(set(msv)))
        print(f'msv:\t{outmsv}')
    return clus

def klook_clusters(clusdf, max_diff_func, len_fold=0.8):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    num_signals = len(clusdf)
    if num_signals == 1:
        clusdf['shift_cluster'] = -1
        return clusdf
    
    # Precompute the 'Target_start' and 'Target_end' columns for quicker access
    start_values = clusdf['Target_start'].values
    end_values = clusdf['Target_end'].values
    svlen_values = clusdf['SVlen'].values
    
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    
    for i in range(1, len(clusdf)):
        current_svlen = svlen_values[i]
        current_start = start_values[i]
        current_end = end_values[i]
        found_cluster = False
        # Adjust the previous cluster range dynamically
        start_index = max(0, i - ceil(num_signals * 0.8))  # previous 80% record

        # Vectorized approach: avoid nested for loops
        for j in range(i - 1, start_index - 1, -1):
            old_len = svlen_values[j]
            relate_size = min(current_svlen / old_len, old_len/current_svlen)
            max_diff = max_diff_func(old_len)
            len_condition = relate_size > len_fold
            pos_condition = (abs(current_start - start_values[j]) <= max_diff or
                             abs(current_end - end_values[j]) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf

def max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 50 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 100 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 450 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 1000

def max_diff_func4INS(old_len):
    if old_len <= 100:
        return 40 + (old_len-50) * 0.8
    elif 100 < old_len <= 500:
        return 80 + (old_len-100) * 0.8
    elif 500 < old_len <= 5000:
        return 400 + (old_len-500) * 0.1
    elif old_len > 5000:
        return 1000


def lowdepth_max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 25 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 50 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 100 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 800

def lowdepth_max_diff_func4INS(old_len):
    if old_len <= 100:
        return 25 + (old_len-50) * 0.4
    elif 100 < old_len <= 500:
        return 50 + (old_len-100) * 0.4
    elif 500 < old_len <= 5000:
        return 200 + (old_len-500) * 0.05
    elif old_len > 5000:
        return 800

def onedepth_all_clus(all_signal,svtype, opened_bam, msvfold, nrate=0.25):
    """
    Process the contig or genome level SV signal.
    Strict condition for merge, single chromosome mode.
    """
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    if all_signal.shape[0] > 1:
        all_signal = merge_and_sort_cr(all_signal)
        print(f"******************signal after merging ****************\n {all_signal}")
        if all_signal.shape[0] > 1:
            all_clus = klook_clusters(all_signal, max_diff_func, msvfold)
        else:
            all_signal['shift_cluster'] = -1
            all_clus = all_signal
    else:
        all_signal['shift_cluster'] = -1
        all_clus = all_signal
    print(all_clus[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    sv_chrom = all_clus["#Target_name"].unique()[0]
    svtype = all_clus["SVType"].unique()[0]
    cluster_col = all_clus.columns[-1]
    clus = []
    msv = []
    for clu in all_clus[cluster_col].unique():
        clu_df = all_clus[all_clus[cluster_col] == clu] 
        print(clu_df[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])   
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len =   mode_or_median(clu_df['SVlen'])
        sv_end =   mode_or_median(clu_df['Target_end'])
        maq =      mode_or_median(clu_df['maq'])
        readsname = clu_df['Query_name'].tolist()
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 250), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 250)
        depth = max([start_local_map,end_local_map])
        if depth > 0:
            SV_rate = round(sv_eye / depth, 2)
        else:
            SV_rate = 1
        if SV_rate < nrate: ## 0.25 to meet Tetraploid heterozygous
            continue
        print(f"************************cluster done***************************\n{[sv_chrom, sv_start, sv_end, sv_len, svid, svtype, '*', sv_eye, SV_rate, maq, readsname]}")
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
        msv.append(svid)
    if len(msv) >=2:
        outmsv = "\t".join(msv)
        print(f'msv:\t{outmsv}')
    return clus

def lowdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, msvfold, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= lowdepth_max_diff_func4INS
    else:
        max_diff_func= lowdepth_max_diff_func4DEL    

    clusdf = klook_clusters(clusdf, max_diff_func, msvfold)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def highdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, msvfold, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    clusdf = klook_clusters(clusdf, max_diff_func, msvfold)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])
    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def filter_overlaps(df):
    """
    to filter overlap SV within same queryname
    """
    filtered_groups = []
    for query, group in df.groupby('Query_name'):
        if len(group) <= 1:
            filtered_groups.append(group)
            continue
        sorted_group = group.sort_values(
            by=['Target_start', 'SVlen'],
            ascending=[True, False]
        )
        kept = []
        last_end = -1
        for _, row in sorted_group.iterrows():
            start = row['Target_start']
            end = row['Target_end']
            if start > last_end:  
                kept.append(row)
                last_end = end
            else:
                 if end > last_end:
                     kept[-1] = row
                     last_end = end
        filtered_groups.append(pd.DataFrame(kept))    
    return pd.concat(filtered_groups, ignore_index=True)

def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates( subset=["Query_name", "Target_start", "Target_end"], keep="first",inplace=True) ## we dont keep the same reads multiple segment
    print(clusdf[["#Target_name","Target_start","Target_end","SVlen","SVType"]])
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf


def merge_and_sort_cr(clusdf):
    """
    filter overlap SV on same query name;
    merge sv based on query name;
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates(
        subset=["Query_name", "Target_start", "Target_end"],
        keep="first",
        inplace=True
    )
    clusdf = filter_overlaps(clusdf)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'min',       
        'Target_end': 'max',    
        'SVlen': 'sum',              
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv[['#Target_name', 'Target_start', 'Target_end', 'SVlen', 'SVType', 'maq', 'seq', 'Query_name']]
    clusdf = clusdf.astype({
        'Target_start': np.int32,
        'Target_end': np.int32,
        'SVlen': np.int32,
        'maq': np.int16
    })
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf

def windows_slide4asm(dfs, svtype, window_size=500):
    """
    fix window gap for assemblies
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        j = i  
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom:
            last_in_window_start = dfs.at[j-1, 'Target_start'] if j > i else current_start
            if svtype == 'DEL':
                dynamic_window_end = last_in_window_start + window_size + 250
            else:
                dynamic_window_end = last_in_window_start + window_size
    
            if dfs.at[j, 'Target_start'] < dynamic_window_end:
                j += 1  
            else:
                break  
        window_df = dfs[i:j].copy()
        if 1 <= len(window_df['Query_name'].unique()):
            window_df.index = range(len(window_df))  
            windows[current_start] = window_df  
        i = j  
    return windows


def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if svtype == 'DEL':
            window_end = current_start + window_size + 250 ## avoid small fragment deletions
        else:
            window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 100:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows

##def windows_slide(dfs, depth, svtype, nreads=2,window_size=500, bigINS=False):
def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if current_svlen < 100:
            window_size = 200 + current_svlen * 0.5
        elif 100 < current_svlen <= 500:
            window_size = 400 + current_svlen * 0.5
        elif 500 < current_svlen <= 1000:
            window_size = 600 + current_svlen * 0.5
        else:
            window_size = 1500
        window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        ## allow big sv in one redas support ##
        if svtype == "INS" and window_df['SVlen'].max() >=5000:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
            print("!!!!!!!!!!!! big SV in window !!!!!!!!!!!!")
        
        if nreads <= len(window_df['Query_name'].unique()) < depth * 10:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows


def windows_slide4tra(df, shift=2000000):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df = df.copy()
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs = dfs.sort_values(by=['Target_start1', 'Target_start2'])
        dfs.index = range(len(dfs))
        i = 0
        while i < len(dfs):
            current_chrom = str(dfs.at[i, '#Target_name1'])
            current_start = dfs.at[i, 'Target_start1']
            window_size = shift
            window_end = current_start + window_size
            j = i
            while j < len(dfs) and dfs.at[j, 'Target_start1'] < window_end:
                j += 1
            window_df = dfs[i:j].copy()
            window_df = window_df.reset_index(drop=True)
            key = f'{chrom1}_{current_start}_{chrom2}'
            windows[key] = window_df
            i = j
    return windows

def klook_clu_tra(win):
    """
    tra clus, the clusdf should be chromosome pair mode,
    each clusdf only allow two chromosome.
    """
    if len(win) == 1:
        win['shift_cluster'] = -1
        return win
    cluster_id = 0
    win.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(win)):
        current_start1 = win.loc[i, 'Target_start1']
        current_start2 = win.loc[i, 'Target_start2']
        found_cluster = False
        start_index = max(0, i - 40)
        for j in range(i - 1, start_index - 1, -1):
            max_diff = 1000
            ## already ensure chrom1 != chrom2
            pos_condition = (abs(current_start1 - win.loc[j, 'Target_start1']) <= max_diff and
                             abs(current_start2 - win.loc[j, 'Target_start2']) <= max_diff)
            if  pos_condition:
                win.loc[i, 'shift_cluster'] = win.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            win.loc[i, 'shift_cluster'] = cluster_id
    return win

def candidate_tra(window, opened_bam, dtype):
    """
    Process translocation data and find candidate translocations.
    """
    win_clus = klook_clu_tra(window)
    tras = []
    for clu in win_clus['shift_cluster'].unique():
        win = win_clus[win_clus['shift_cluster']==clu]
        chr1 = win['#Target_name1'].iloc[0]
        chr1_start = mode_or_median(win['Target_start1'])
        sv_len = 0
        chr2 = win['#Target_name2'].iloc[0]
        chr2_start = mode_or_median(win['Target_start2'])
        maq = mode_or_median(win['maq'])
        readsname = set(win['Query_name'])
        svid = f'{chr2}:{chr2_start}_{chr1}:{chr1_start}'
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, chr1, max(0, chr1_start - 250), max(chr1_start - 150, 0))
        end_local_map = local_cov(opened_bam, chr2, chr2_start + 150, chr2_start + 250)
        local_depth = np.mean([start_local_map, end_local_map])
        if local_depth >0:
            SV_rate = round(sv_eye / local_depth, 2)
        else:
            SV_rate = 1
        if dtype in ['ont','hifi', 'pb']:
            ## to one depth ##
            if sv_eye >= (local_depth * 0.1 + 0.7 ):
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
        elif dtype in ['cr', 'sr']:
            if sv_eye > local_depth * 0.25:
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tras

def load_and_process_sv_data(args):
    from math import floor
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    except pd.errors.EmptyDataError:
        print(f"Warning: File {args.raw_signal} is empty.")
        sv_indel_data = pd.DataFrame()
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else ceil(float(depth_stat.iloc[0, 3])+0.3)
        if args.nreads:
            nreads_fil = args.nreads
        else:
            nreads_fil = floor(depth / 10)
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth, nreads_fil
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data = sv_indel_data.copy()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_indel_data['Target_end'] = sv_indel_data['Target_end'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"***************** empty {supp_align_file} file ******************")
        else:
            if not msv.empty:
                print(f'*************** {args.raw_signal}.suppAlign has {msv.shape[0]} rows ********************')
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                #msv = msv.drop_duplicates()
                print(f'*************** {args.raw_signal}.suppAlign after duplicates drop has {msv.shape[0]} rows ********************')
                #msv = msv.reset_index(drop=True)
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['Target_end'] = sv_data[svtype]['Target_end'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    print(f'******************** all chromosomes list {chroms} ****************************')
    return sv_data, chroms, depth, nreads_fil

def process_svtype(args, sv_data, chroms, svtype, depth, nreads, minLen):
    """
        cov info parse by covfile or bam file
    """
    print(f"start klook for {args.raw_signal}  SV type: {svtype}") 
    if args.dtype in ['cr', 'sr']:
        try:
            covinfo = pd.read_csv(args.covfile, sep="\t", index_col=None, header=None)
            covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
            covinfo['target_chr'] = covinfo['target_chr'].astype(str)
            bam_path = covinfo
        except FileNotFoundError:
            print(f"Error: Coverage file {args.covfile} not found.")
            return [], [], [], [], []
    else:
        bam_path = args.bam

    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    if args.dtype in ['pb', 'ont', 'hifi']:
        print(f'data type is {args.dtype}')
        try:
            with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
                for chrom in chroms:
                    chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                    def process_sv_type(svtype):
                        if svtype == "TRA":
                            log = open("log_tra", 'w')
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            windows = windows_slide4tra(chrom_data[svtype], 2000000)
                            print(f'The TRA windows by 2M: \n{windows}', file=log)
                            tra_list = []
                            for win in windows.values():
                                win_clus = klook_clu_tra(win)
                                print(win_clus.iloc[:,[0,2,3,4,5,-1]], file=log)
                                tra = candidate_tra(win_clus, opened_bam, args.dtype)
                                print(tra, file=log)

                                if tra:
                                    tra_list.extend(tra)
                            log.close()
                            return tra_list
                        else:
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                            sv_windows = windows_slide(sv_dfs, depth, svtype, nreads,args.window_size)
                            candidate_svs = []
                            if args.dtype in ['hifi', 'ont']:
                                hdepth = 5
                            else:
                                hdepth = 10
                            for sv_window in sv_windows.values():
                                if len(sv_window['Query_name']) > hdepth: 
                                    candisv = highdepth_clu(sv_window,args.num_hap, svtype, opened_bam, nreads, args.msvfold, args.rate_depth, 1)
                                else:
                                    candisv = lowdepth_clu(sv_window, args.num_hap, svtype, opened_bam, nreads, args.msvfold, args.rate_depth, 1)
                                candidate_svs.extend(candisv)
                            return candidate_svs

                    if svtype == "DEL":
                        del_clus.extend(process_sv_type("DEL"))
                    elif svtype == "INS":
                        ins_clus.extend(process_sv_type("INS"))
                    elif svtype == "INV":
                        inv_clus.extend(process_sv_type("INV"))
                    elif svtype == "DUP":
                        dup_clus.extend(process_sv_type("DUP"))
                    elif svtype == "TRA":
                        tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f"Error processing BAM file: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    elif args.dtype in ['sr', 'cr']:
        print(f'data type is {args.dtype}')
        try:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                def process_sv_type(svtype):
                    if svtype != "TRA":
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                        if sv_dfs.empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        sv_dfs = sv_dfs.sort_values(by=['Target_start', 'SVlen'])
                        sv_dfs.index = range(len(sv_dfs))
                        print(sv_dfs.head(10))
                        print("**************************** Calling one depth all clustering ***********************")
                        candidate_svs = []
                        win_dfs = windows_slide4asm(sv_dfs, svtype, args.window_size)
                        for win_df in win_dfs.values():
                            print(f'**************************window signals*********************\n{win_df}')
                            candidate_sv = onedepth_all_clus(win_df, svtype, bam_path, args.msvfold, max(args.rate_depth,0.25))
                            candidate_svs.extend(candidate_sv)
                        return candidate_svs
                    else:
                        if chrom_data[svtype].empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        windows = windows_slide4tra(chrom_data[svtype], 1000000)
                        tra_list = []
                        for value in windows.values():
                            win = value
                            tra = candidate_tra(win, bam_path,args.dtype)
                            if tra:
                                tra_list.extend(tra)
                        print(f'****************************** the tra condidate ***************************')
                        print(tra_list)
                        return tra_list
                if svtype == "DEL":
                    del_clus.extend(process_sv_type("DEL"))
                elif svtype == "INS":
                    ins_clus.extend(process_sv_type("INS"))
                elif svtype == "INV":
                    inv_clus.extend(process_sv_type("INV"))
                elif svtype == "DUP":
                    dup_clus.extend(process_sv_type("DUP"))
                elif svtype == "TRA":
                    tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f'Error in processing {args.dtype}: {args.raw_signal}')
            print(f"Error processing {svtype} for chromosome {chrom}: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        print(f"Error: sequence data dtype error: {args.dtype} is not in [pb,hifi,ont,sr,cr], please check parameter -dtype")
        return [],[],[],[],[]

def candidateSV(args):
    sv_data, chroms, depth, nreads = load_and_process_sv_data(args)
    print(chroms, depth, nreads)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth, nreads, args.min) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            if result is not None:
                tra_clus += result[0]
                del_clus += result[1]
                ins_clus += result[2]
                inv_clus += result[3]
                dup_clus += result[4]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        return [], [], [], [], []


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True,
                    help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=800, type=int,
                    help="the distance shift of breakpoint to cluster the TRA/big_INV/big_DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=10000000, help="the max SV length")
    IN.add_argument("-m", dest="min", type=int, default=45, help="the minimum SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    IN.add_argument("--nreads", dest="nreads", type=int, help="the minimum numbers of reads to support SV, if not provided, we use average_depth / 10 as threshold" )
    IN.add_argument("--rate_depth", dest="rate_depth", type=float, default=0.1, help="the sv supports of local depth ratio to support sv, 0.1 means the percent of local reads shoule support sv")
    IN.add_argument("--window", dest="window_size", type=int, default=500, help="the window size of signal to parse in klook cluster, 500bp suggested")
    IN.add_argument("--num_hap", dest="num_hap", type=int, default=2, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    IN.add_argument("--msvfold", dest="msvfold", type=float, default=0.8, help="the SV size fold change of multi sv alleles")

    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

```

### ./PSVGT1.0_Dynamic_Window.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
import multiprocessing
import concurrent
from functools import partial
from os.path import basename
from concurrent.futures import ThreadPoolExecutor, as_completed

PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from Sub_readfa2Dict import readfa2Dict

def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("PSVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-win", "--window",default=500,type=int,help="Window size to parse signal, this parameter is to merge fragment SV due to aligment, default window size is 500bp to merge cluster the SV with a sv_start less than window size")
    parser.add_argument("-diploid", "--diploid", help="for diploid resolved assembly, to get a phased genotype please provide table list each line in format hap1\thap2\tSampleName in cr directory")
    parser.add_argument("-polyploid", "--polyploid", help="for polyploid haplotype resolved assembly like potato(4 haplotype assemblies available), to get a merge genotype of samples please provide table list each line in format hap1\thap2\thap3\thapn\tSampleName")


    parser.add_argument("-w", "--max_workers", default = 6,type=int, help="the max workers thread pool excutor, 6 means run 6 samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="no", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=30,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.10, type=float, help= "the percent of reads that support a candidate SV (0.10 means at a depth 20X region, a SV signal should have at least 2 reads support, this parameter is for the variaty depth of hifi/ont/pb samples")
    parser.add_argument("-nreads",  "--nreads", type=int, help= "the number of reads to support a candidate SV (SV signal should have at least numbers reads support, this parameter is for the various depth of hifi/ont/pb samples")
    parser.add_argument("--num_hap", "--num_hap", default=2, type=int, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    parser.add_argument("-msv","--msv_mode",default="no", help= "Setting msv to `yes` signals of INS,DEL,INV,DUP,TRA will captured, while setting to `no` Only insertions and deletions will be capture from alignment ")
    parser.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1, for species like polyploid potato we suggest 0.8")
    parser.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 0.10, the genotype will be 0/0")
    parser.add_argument("-sr_homo_rate", dest="sr_homo_rate",default=0.65, type=float, help="to determine a homozygous site, if 0.65 of the local mapping signal suport the sv the genotyoe will be 1/1, you can lower down the value if your specise have a low heterozygous rate")
    parser.add_argument("-sr_ref_rate", dest="sr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 10 percent, the genotype will be 0/0, you can increase the value to filter putative heterozygous SV")

    parser.add_argument("-span", dest="span", default=50, type=int, help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    

    done_analysor = file_capture(args.outdir, ".record.txt")
    def run_clu2fil_cmd(clu2fil_cmd):
        try:
            subprocess.run(clu2fil_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {clu2fil_cmd}\n{e}")
    def process_single_contig(contig, dtype, args, PSVGT, fai):
        try:
            ## final chromsome result ##
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"
            ## skip samples ##
            if os.path.exists(done_name):
                print(f"Skipping processed contig: {basename(contig)}")
                return

            ## minimaping and signal detect ##
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py '
                f'-i {contig} -dtype {dtype} -r {args.refGenome} '
                f'-m {args.min} -maq {maq} -o {args.outdir} '
                f'-minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            ## chromosome pool run ##
            with multiprocessing.Pool(processes=args.max_workers) as pool:
                if dtype in ['hifi','ont','pb'] and args.csv and args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} --nreads {args.nreads}'
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                
                elif dtype in ['hifi','ont','pb'] and args.csv and not args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                elif dtype in ['hifi','ont','pb'] and args.nreads and not args.csv:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--nreads {args.nreads} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                else:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]

                pool.map(run_clu2fil_cmd, clu2fil_cmds)
            ## ACC SV ##
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/0_tmp_{basename(contig)} '
                f'-fai {args.refGenome}.fai -M {args.max} --nrate {args.csv}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(contig)}: {str(e)}")
    def add_commands4fq(files, dtype):
        max_workers = min(len(files), args.max_workers)
        if max_workers >0:
            print(f'{max_workers} workers')
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        process_single_contig,
                        contig, dtype, args, PSVGT, fai
                    ) for contig in files
                ]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Thread error: {str(e)}")

    def process_single_bam(bam, dtype, args, PSVGT, fai):
        try:## last chromosome result ##
            done_flag = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if os.path.exists(done_flag):
                print(f"Skipping processed BAM: {basename(bam)}")
                return

            # Phase 1: signal
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py '
                f'-b {bam} -o {args.outdir}/{basename(bam)} '
                f'-m {args.min} -maq {maq} -dtype {dtype} '
                f'-msv {args.msv_mode} -fai {args.refGenome}.fai'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            # Phase 2: parafly chrs
            base_prefix = basename(bam).replace(".bam", "")
            if args.nreads and args.csv and dtype in ['hifi', 'ont', 'pb']:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--nreads {args.nreads} --rate_depth {args.csv} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.csv and dtype in ['hifi', 'ont', 'pb'] and not args.nreads:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--rate_depth {args.csv} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.nreads and dtype in ['hifi', 'ont', 'pb'] and not args.csv:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--nreads {args.nreads}  '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            else:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]


            # dynamic cpu
            with multiprocessing.Pool(processes=min(len(fai[0]), os.cpu_count()//5)) as pool:
                pool.map(run_clu2fil_cmd, chrom_commands)

            # Phase 3: merge signal
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/{base_prefix} '
                f'-fai {args.refGenome}.fai --nrate {args.csv} --minimaq {args.maq}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(bam)}: {str(e)}")

    def add_commands4bam(files, dtype):
        if len(files)>0:
            max_threads = min(len(files),args.max_workers)
            with ThreadPoolExecutor(max_workers=max_threads) as executor:
                task_func = partial(
                    process_single_bam,
                    dtype=dtype,
                    args=args,
                    PSVGT=PSVGT,
                    fai=fai
                )
                futures = {executor.submit(task_func, bam): bam for bam in files}
                for future in concurrent.futures.as_completed(futures):
                    bam_file = futures[future]
                    try:
                        future.result()
                        print(f"Success: {basename(bam_file)}")
                    except Exception as e:
                        print(f"Failed processing {basename(bam_file)}: {str(e)}")

    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        add_commands4fq(file_capture(args.crdir, ".gz"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file} -m {args.maq} -lr_homo_rate {args.lr_homo_rate} -lr_ref_rate {args.lr_ref_rate}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf {args.refGenome}.fai"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")

    ###################### haplotype resoved assembly genotype phased ###################
    if args.diploid:
        print("*************** diploid calling **************")
        with open(args.diploid, 'r') as f:
            lines = f.readlines()
        for hh in lines:
            h1 = hh.strip().split("\t")[0]
            h1 = "2_tmp_" + h1.replace(".bam", "") + "_genotype.txt" if  h1[-4:] == ".bam" else "2_tmp_" + h1 + "_genotype.txt"
            h2 = hh.strip().split("\t")[1]
            h2 = "2_tmp_" + h2.replace(".bam", "") + "_genotype.txt" if  h2[-4:] == ".bam" else "2_tmp_" + h2 + "_genotype.txt"
            samplename =  hh.strip().split("\t")[2]
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_diploid_asm.py {args.outdir}/{h1} {args.outdir}/{h2} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f"***************** try to phased hap1: {args.outdir}/{h1} and hap2: {args.outdir}/{h2} to {args.outdir}/2_tmp_{samplename}_genotype.txt ******************")
            tab2vcf =    f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(tab2vcf)
    if args.polyploid:
        print("*************** polyploid genotype merging **************")
        with open(args.polyploid, 'r') as f:
            lines = f.readlines()
        path_haps = ""
        for hh in lines:
            path_haps = ""
            haps = hh.strip().split("\t")
            total_haps = len(haps) - 1
            samplename = haps[-1]
            for i in range(total_haps):
                hap = haps[i] 
                hap_file = "2_tmp_" + hap.replace(".bam", "") + "_genotype.txt" if  hap[-4:] == ".bam" else "2_tmp_" + hap + "_genotype.txt"
                path_haps += f'{args.outdir}/{hap_file} '
            print(f"*************************** phased polyploid {samplename} ***************************** ")
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_polyploid_genome_gt.py {path_haps} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f'********************** phased polyploid ************************\n{phased_cmd}')
            vcf_cmd = f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(vcf_cmd)


    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {bam} -maq {args.maq} -span {args.span} -s 50 -n {sampleID} -o {args.outdir} -homo_rate {args.sr_homo_rate} -ref_rate {args.sr_ref_rate} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf {args.refGenome}.fai"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), args.max_workers)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 80 --max 600 --frank 400 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 80 600 400 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./back_v1_PSVGT1.0.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from concurrent.futures import ThreadPoolExecutor, as_completed
from os.path import basename
from Sub_readfa2Dict import readfa2Dict
def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("Pop SVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-w", "--max_workers", default = 4,type=int, help="the max workers thread pool excutor, 4 means run for samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="yes", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=30,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.25, type=float, help= "the percent of reads that support a candidate SV (0.25 means at a depth 20X region, a SV signal should have at least 5 reads support, this parameter is for the variaty depth of HIFI/ONT/PB samples")
    parser.add_argument("-msv","--msv_mode",default="no", help= "In msv mode signals of INS,DEL,INV,DUP,TRA will captured from ont/hifi/pb, while for assemble contig from short reads or genome lelve samples we detect SVInDel Only. If no hifi or ont or pacbio data is provided, please setting -msv no, PSVGT will detect SVInDel Only")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    
    ############# SVInDel population mode step0 #############
    Pop_SV_Analysor_cmds = []
    done_analysor = file_capture(args.outdir, ".record.txt")
    def add_commands4fq(files, dtype):
        for contig in files:
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt" ## last chrom
            if done_name not in done_analysor:
                maq = min(args.maq, 60) 
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py -i {contig} -dtype {dtype} -r {args.refGenome} -m {args.min}  -maq {maq} -o {args.outdir} -minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
                clu2fil_cmds = ''
                for chrom in fai[0]:
                    clu2fil_cmd = f'&& python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthPASS.py -f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max} --b {args.outdir}/0_tmp_{basename(contig)}.bam --cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov'
                    clu2fil_cmds += clu2fil_cmd
                ACC_SV_cmd = f'&& python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/0_tmp_{basename(contig)} -fai  {args.refGenome}.fai'
                final_cmd = signal_cmd + clu2fil_cmds + ACC_SV_cmd
                Pop_SV_Analysor_cmds.append(final_cmd)

    def add_commands4bam(files, dtype):
        for bam in files:
            done_name = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if done_name not in done_analysor:
                maq = min(args.maq + 5, 60) ## biger than illumina breakpoints quality
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py -b {bam} -o {args.outdir}/{basename(bam)} -m {args.min} -maq {args.maq} -dtype {dtype} -msv {args.msv_mode} -fai {args.refGenome}.fai' 
                clu2fil_cmds =''
                for chrom in fai[0]:
                    clu2fil_cmd = f'&& python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthPASS.py -f {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max}  --b {bam} --cov {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt.cov'
                    clu2fil_cmds += clu2fil_cmd
                ACC_SV_cmd = f' &&  python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/{basename(bam).replace(".bam", "")} -fai  {args.refGenome}.fai'
                final_cmd = signal_cmd + clu2fil_cmds + ACC_SV_cmd
                Pop_SV_Analysor_cmds.append(final_cmd)

    # Capture files for different types
    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    # Execute commands using ThreadPool
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures1 = [executor.submit(execute_commands, cmd) for cmd in Pop_SV_Analysor_cmds]
        results1 =[]
        for future in as_completed(futures1):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=all_log)
            else:
                print(stdout,file=all_log)
            results1.append((stdout,stderr,returncode,cmd))
    
    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")
    
    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_clustered_Record.txt -mapf {bam} -s 50 -n {sampleID} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), 40)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 50 --max 500 --frank 500 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 50 500 500 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./PSVGT1.0_Fix_Window.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
import multiprocessing
import concurrent
from functools import partial
from os.path import basename
from concurrent.futures import ThreadPoolExecutor, as_completed

PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from Sub_readfa2Dict import readfa2Dict

def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("PSVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-win", "--window",default=500,type=int,help="Window size to parse signal, this parameter is to merge fragment SV due to aligment, default window size is 500bp to merge cluster the SV with a sv_start less than window size")
    parser.add_argument("-diploid", "--diploid", help="for diploid resolved assembly, to get a phased genotype please provide table list each line in format hap1\thap2\tSampleName in cr directory")
    parser.add_argument("-polyploid", "--polyploid", help="for polyploid haplotype resolved assembly like potato(4 haplotype assemblies available), to get a merge genotype of samples please provide table list each line in format hap1\thap2\thap3\thapn\tSampleName")


    parser.add_argument("-w", "--max_workers", default = 6,type=int, help="the max workers thread pool excutor, 6 means run 6 samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="no", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=30,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.10, type=float, help= "the percent of reads that support a candidate SV (0.10 means at a depth 20X region, a SV signal should have at least 2 reads support, this parameter is for the variaty depth of hifi/ont/pb samples")
    parser.add_argument("-nreads",  "--nreads", type=int, help= "the number of reads to support a candidate SV (SV signal should have at least numbers reads support, this parameter is for the various depth of hifi/ont/pb samples")
    parser.add_argument("--num_hap", "--num_hap", default=2, type=int, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    parser.add_argument("-msv","--msv_mode",default="no", help= "Setting msv to `yes` signals of INS,DEL,INV,DUP,TRA will captured, while setting to `no` Only insertions and deletions will be capture from alignment ")
    parser.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1, for species like polyploid potato we suggest 0.8")
    parser.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 0.10, the genotype will be 0/0")
    parser.add_argument("-sr_homo_rate", dest="sr_homo_rate",default=0.65, type=float, help="to determine a homozygous site, if 0.65 of the local mapping signal suport the sv the genotyoe will be 1/1, you can lower down the value if your specise have a low heterozygous rate")
    parser.add_argument("-sr_ref_rate", dest="sr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 10 percent, the genotype will be 0/0, you can increase the value to filter putative heterozygous SV")

    parser.add_argument("-span", dest="span", default=50, type=int, help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    

    done_analysor = file_capture(args.outdir, ".record.txt")
    def run_clu2fil_cmd(clu2fil_cmd):
        try:
            subprocess.run(clu2fil_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {clu2fil_cmd}\n{e}")
    def process_single_contig(contig, dtype, args, PSVGT, fai):
        try:
            ## final chromsome result ##
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"
            ## skip samples ##
            if os.path.exists(done_name):
                print(f"Skipping processed contig: {basename(contig)}")
                return

            ## minimaping and signal detect ##
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py '
                f'-i {contig} -dtype {dtype} -r {args.refGenome} '
                f'-m {args.min} -maq {maq} -o {args.outdir} '
                f'-minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            ## chromosome pool run ##
            with multiprocessing.Pool(processes=args.max_workers) as pool:
                if dtype in ['hifi','ont','pb'] and args.csv and args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} --nreads {args.nreads}'
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                
                elif dtype in ['hifi','ont','pb'] and args.csv and not args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                elif dtype in ['hifi','ont','pb'] and args.nreads and not args.csv:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--nreads {args.nreads} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                else:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]

                pool.map(run_clu2fil_cmd, clu2fil_cmds)
            ## ACC SV ##
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/0_tmp_{basename(contig)} '
                f'-fai {args.refGenome}.fai -M {args.max} --nrate {args.csv}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(contig)}: {str(e)}")
    def add_commands4fq(files, dtype):
        max_workers = min(len(files), args.max_workers)
        if max_workers >0:
            print(f'{max_workers} workers')
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        process_single_contig,
                        contig, dtype, args, PSVGT, fai
                    ) for contig in files
                ]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Thread error: {str(e)}")

    def process_single_bam(bam, dtype, args, PSVGT, fai):
        try:## last chromosome result ##
            done_flag = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if os.path.exists(done_flag):
                print(f"Skipping processed BAM: {basename(bam)}")
                return

            # Phase 1: signal
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py '
                f'-b {bam} -o {args.outdir}/{basename(bam)} '
                f'-m {args.min} -maq {maq} -dtype {dtype} '
                f'-msv {args.msv_mode} -fai {args.refGenome}.fai'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            # Phase 2: parafly chrs
            base_prefix = basename(bam).replace(".bam", "")
            if args.nreads and args.csv and dtype in ['hifi', 'ont', 'pb']:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--nreads {args.nreads} --rate_depth {args.csv} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.csv and dtype in ['hifi', 'ont', 'pb'] and not args.nreads:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--rate_depth {args.csv} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.nreads and dtype in ['hifi', 'ont', 'pb'] and not args.csv:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--nreads {args.nreads}  '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            else:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]


            # dynamic cpu
            with multiprocessing.Pool(processes=min(len(fai[0]), os.cpu_count()//5)) as pool:
                pool.map(run_clu2fil_cmd, chrom_commands)

            # Phase 3: merge signal
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/{base_prefix} '
                f'-fai {args.refGenome}.fai --nrate {args.csv} --minimaq {args.maq}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(bam)}: {str(e)}")

    def add_commands4bam(files, dtype):
        if len(files)>0:
            max_threads = min(len(files),args.max_workers)
            with ThreadPoolExecutor(max_workers=max_threads) as executor:
                task_func = partial(
                    process_single_bam,
                    dtype=dtype,
                    args=args,
                    PSVGT=PSVGT,
                    fai=fai
                )
                futures = {executor.submit(task_func, bam): bam for bam in files}
                for future in concurrent.futures.as_completed(futures):
                    bam_file = futures[future]
                    try:
                        future.result()
                        print(f"Success: {basename(bam_file)}")
                    except Exception as e:
                        print(f"Failed processing {basename(bam_file)}: {str(e)}")

    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        add_commands4fq(file_capture(args.crdir, ".gz"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file} -m {args.maq} -lr_homo_rate {args.lr_homo_rate} -lr_ref_rate {args.lr_ref_rate}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf {args.refGenome}.fai"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")

    ###################### haplotype resoved assembly genotype phased ###################
    if args.diploid:
        print("*************** diploid calling **************")
        with open(args.diploid, 'r') as f:
            lines = f.readlines()
        for hh in lines:
            h1 = hh.strip().split("\t")[0]
            h1 = "2_tmp_" + h1.replace(".bam", "") + "_genotype.txt" if  h1[-4:] == ".bam" else "2_tmp_" + h1 + "_genotype.txt"
            h2 = hh.strip().split("\t")[1]
            h2 = "2_tmp_" + h2.replace(".bam", "") + "_genotype.txt" if  h2[-4:] == ".bam" else "2_tmp_" + h2 + "_genotype.txt"
            samplename =  hh.strip().split("\t")[2]
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_diploid_asm.py {args.outdir}/{h1} {args.outdir}/{h2} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f"***************** try to phased hap1: {args.outdir}/{h1} and hap2: {args.outdir}/{h2} to {args.outdir}/2_tmp_{samplename}_genotype.txt ******************")
            tab2vcf =    f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(tab2vcf)
    if args.polyploid:
        print("*************** polyploid genotype merging **************")
        with open(args.polyploid, 'r') as f:
            lines = f.readlines()
        path_haps = ""
        for hh in lines:
            path_haps = ""
            haps = hh.strip().split("\t")
            total_haps = len(haps) - 1
            samplename = haps[-1]
            for i in range(total_haps):
                hap = haps[i] 
                hap_file = "2_tmp_" + hap.replace(".bam", "") + "_genotype.txt" if  hap[-4:] == ".bam" else "2_tmp_" + hap + "_genotype.txt"
                path_haps += f'{args.outdir}/{hap_file} '
            print(f"*************************** phased polyploid {samplename} ***************************** ")
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_polyploid_genome_gt.py {path_haps} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f'********************** phased polyploid ************************\n{phased_cmd}')
            vcf_cmd = f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(vcf_cmd)


    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {bam} -maq {args.maq} -span {args.span} -s 50 -n {sampleID} -o {args.outdir} -homo_rate {args.sr_homo_rate} -ref_rate {args.sr_ref_rate} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf {args.refGenome}.fai"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), args.max_workers)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 80 --max 600 --frank 400 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 80 600 400 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./SSVGT1.0.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
PSVGT_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT_folder}/SV_Genotyper")
from concurrent.futures import ThreadPoolExecutor, as_completed
from os.path import basename
from Sub_readfa2Dict import readfa2Dict
def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            captures.append(os.path.join(dir, file))
    return captures 


def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        check_dir(f"00_bwa_mem_out")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("SSVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-w", "--max_workers", default = 4,type=int, help="the max workers thread pool excutor, 4 means run for samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=50, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="yes", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=45,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-support",  "--support",default=0.2, type=float, help= "the percent of reads that support a candidate SVInDel(0.15 means 20X data should have at least 3 reads support SVINDel), this parameter is for the variaty depth of HIFI/ONT/PB samples")
    parser.add_argument("-msv","--msv_mode",default="no", help= "In msv mode signals of INS,DEL,INV,DUP,TRA will captured from ont/hifi/pb, while for assemble contig from short reads or genome lelve samples we detect SVInDel Only. If no hifi or ont or pacbio data is provided, please setting -msv no, PSVGT will detect SVInDel Only")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai',sep="\t",index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    
    ############# SVInDel population mode step0 #############
    Pop_SV_Analysor_cmds = []
    done_analysor = file_capture(args.outdir, ".record.txt")
    def add_commands4fq(files, dtype):
        for contig in files:
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[0]}.record.txt" ## last chrome
            if done_name not in done_analysor:
                maq = min(args.maq + 5, 60) ## biger than illumina breakpoints quality
                signal_cmd = f'python {PSVGT_folder}/SV_Genotyper/0.PSVGT_raw2Signal.py -i {contig} -dtype {dtype} -r {args.refGenome} -m {args.min}  -maq {maq} -o {args.outdir} -minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
                clu2fil_cmds = ''
                for chrom in fai[0]:
                    clu2fil_cmd = f'&& python {PSVGT_folder}/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt -dtype {dtype} --b {args.outdir}/0_tmp_{basename(contig)}.bam --cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov -s 1000 -M {args.max} -csv 0.2'
                    clu2fil_cmds += clu2fil_cmd
                final_cmd = signal_cmd + clu2fil_cmds + f"&& python {PSVGT_folder}/SV_Genotyper/1.SSV_signal_cluster.py  -preffix {args.outdir}/0_tmp_{basename(contig)} -s 50 -fai {args.refGenome}.fai"
                Pop_SV_Analysor_cmds.append(final_cmd)

    def add_commands4bam(files, dtype):
        for bam in files:
            done_name = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[0]}.record.txt"
            if done_name not in done_analysor:
                maq = min(args.maq + 5, 60) ## biger than illumina breakpoints quality
                signal_cmd = f'python {PSVGT_folder}/SV_Genotyper/0.Signal4bam_PSVGT.py -b {bam} -o {args.outdir}/{basename(bam)} -m {args.min} -maq {args.maq} -dtype {dtype} -msv {args.msv_mode} -fai {args.refGenome}.fai'
                clu2fil_cmds =''
                for chrom in fai[0]:
                    clu2fil_cmd = f'&& python {PSVGT_folder}/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt -dtype {dtype} --b {bam} -s 1000 -M {args.max} -csv 0.2 --cov {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt.cov'
                    clu2fil_cmds += clu2fil_cmd
                final_cmd = signal_cmd + clu2fil_cmds + f"&& python {PSVGT_folder}/SV_Genotyper/1.SSV_signal_cluster.py  -preffix {args.outdir}/{basename(bam).replace('.bam', '')} -s 50 -fai {args.refGenome}.fai"
                Pop_SV_Analysor_cmds.append(final_cmd)
    # Capture files for different types
    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use ont mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.hifidir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.hifidir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        already_maps += file_capture(args.crdir, ".bam")

    # Execute commands using ThreadPool
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures1 = [executor.submit(execute_commands, cmd) for cmd in Pop_SV_Analysor_cmds]
        results1 =[]
        for future in as_completed(futures1):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=all_log)
            else:
                print(stdout,file=all_log)
            results1.append((stdout,stderr,returncode,cmd))
    
    
    ## step2 genotyping by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(f'*************************************** Gentypeing the out file is {if_done_name} **************************************************')
            acc_name = basename(mapinfo_file).replace('0_tmp_', '').replace(".bam", '')
            cluster_signal_file   =  f"{args.outdir}/{basename(mapinfo_file).replace('.bam', '')}_clustered_Record.txt" ## only capture target signal file
            print(cluster_signal_file)
            cmd = f"python {PSVGT_folder}/SV_Genotyper/2.Single_lrSVGT_V1.py -i {cluster_signal_file} -mapf {mapinfo_file}  -n {acc_name} -o {args.outdir} && python {PSVGT_folder}/SV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")
    
    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        with open(f"{args.outdir}/bwa_bams.txt", "r") as bwas:
            bwa_bams_lines = bwas.readlines()
            for line in bwa_bams_lines:
                bam = line.strip()
                sampleID = basename(bam)[:-4]
                cluster_signal = f'{args.outdir}/0_tmp_{basename(bam).replace("_sorted.bam", "")}.contigs.fa.record.txt_Clustered_Record.txt'
                cluster_signal_files  = file_capture(f'{args.outdir}', '_clustered_Record.txt')
                cluster_signal_file   =  [f for f in cluster_signal_files if sampleID.replace("_sorted", "") in f] ## only capture target signal file
                print(cluster_signal_files, cluster_signal_file)

                bpgt_cmd =  f"python {PSVGT_folder}/SV_Genotyper/2.Pop_srSVGT_V1.py -i {cluster_signal_file[0]} -mapf {bam} -s 25 -n {sampleID} -o {args.outdir} && python {PSVGT_folder}/SV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf"
                print(bpgt_cmd)
                breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), 40)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT_folder}/SV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT_folder}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 50 --max 500 --frank 500 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT_folder}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 50 500 500 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT_folder}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT_folder}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel && rm {final_gt}_tmp.tab &")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT_folder}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT_folder}/CapsPop/pop_maf0.05_caps.py {PSVGT_folder}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./back_v2_PSVGT1.0.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
import multiprocessing
PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from concurrent.futures import ThreadPoolExecutor, as_completed
from os.path import basename
from Sub_readfa2Dict import readfa2Dict
def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("Pop SVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-w", "--max_workers", default = 4,type=int, help="the max workers thread pool excutor, 4 means run for samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="yes", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=1,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.25, type=float, help= "the percent of reads that support a candidate SV (0.25 means at a depth 20X region, a SV signal should have at least 5 reads support, this parameter is for the variaty depth of HIFI/ONT/PB samples")
    parser.add_argument("-msv","--msv_mode",default="no", help= "In msv mode signals of INS,DEL,INV,DUP,TRA will captured from ont/hifi/pb, while for assemble contig from short reads or genome lelve samples we detect SVInDel Only. If no hifi or ont or pacbio data is provided, please setting -msv no, PSVGT will detect SVInDel Only")
    parser.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1")
    parser.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.05, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 0.05, the genotype will be 0/0")
    parser.add_argument("-span", dest="span", default=50, type=int, help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    
    done_analysor = file_capture(args.outdir, ".record.txt")
    def run_clu2fil_cmd(clu2fil_cmd):
        try:
            subprocess.run(clu2fil_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {clu2fil_cmd}\n{e}")

    def add_commands4fq(files, dtype):
        for contig in files:
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"  ## last chrom
            if done_name not in done_analysor:
                maq = min(args.maq, 60)
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py -i {contig} -dtype {dtype} -r {args.refGenome} -m {args.min}  -maq {maq} -o {args.outdir} -minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
                try:
                    subprocess.run(signal_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running signal command: {signal_cmd}\n{e}")
                    continue

                clu2fil_cmds = []
                for chrom in fai[0]:
                    clu2fil_cmd = f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthPASS.py -f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max} --b {args.outdir}/0_tmp_{basename(contig)}.bam --cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov'
                    clu2fil_cmds.append(clu2fil_cmd)
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                pool.map(run_clu2fil_cmd, clu2fil_cmds)
                pool.close()
                pool.join()

                ACC_SV_cmd = f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/0_tmp_{basename(contig)} -fai  {args.refGenome}.fai'
                try:
                    subprocess.run(ACC_SV_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running ACC_SV command: {ACC_SV_cmd}\n{e}")

    def add_commands4bam(files, dtype):
        for bam in files:
            done_name = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if done_name not in done_analysor:
                maq = min(args.maq, 60)  ## biger than illumina breakpoints quality
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py -b {bam} -o {args.outdir}/{basename(bam)} -m {args.min} -maq {args.maq} -dtype {dtype} -msv {args.msv_mode} -fai {args.refGenome}.fai'
                try:
                    subprocess.run(signal_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running signal command: {signal_cmd}\n{e}")
                    continue
                clu2fil_cmds = []
                for chrom in fai[0]:
                    clu2fil_cmd = f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthPASS.py -f {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max}  --b {bam} --cov {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt.cov'
                    clu2fil_cmds.append(clu2fil_cmd)
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                pool.map(run_clu2fil_cmd, clu2fil_cmds)
                pool.close()
                pool.join()

                ACC_SV_cmd = f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/{basename(bam).replace(".bam", "")} -fai  {args.refGenome}.fai'
                try:
                    subprocess.run(ACC_SV_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running ACC_SV command: {ACC_SV_cmd}\n{e}")
    # Capture files for different types
    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        add_commands4fq(file_capture(args.crdir, ".gz"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file} -m {args.maq} -lr_homo_rate {args.lr_homo_rate} -lr_ref_rate {args.lr_ref_rate}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")
    
    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {bam} -m {args.maq} -span {args.span} -s 50 -n {sampleID} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), 40)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 50 --max 500 --frank 500 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 50 500 500 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

### ./backup_20250525_PSVGT1.0.py
```python
import pandas as pd
import os
import subprocess
import re
import sys
import shutil
import multiprocessing
import concurrent
from functools import partial
from os.path import basename
from concurrent.futures import ThreadPoolExecutor, as_completed

PSVGT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{PSVGT}/PSV_Genotyper")
sys.path.append(f"{PSVGT}/PSV_Signal")
from Sub_readfa2Dict import readfa2Dict

def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        stdout = stdout.decode('utf-8')
    except UnicodeDecodeError:
        stdout = stdout.decode('ISO-8859-1')
    try:
        stderr = stderr.decode('utf-8')
    except UnicodeDecodeError:
        stderr = stderr.decode('ISO-8859-1')
    return stdout, stderr, process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

def check_dir(dir, log=None):
    if os.path.isdir(dir) == False:
        os.mkdir(dir)
def rm_dir(dir, log=None):
    if os.path.isdir(dir) == True:
        shutil.rmtree(dir)
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)
def recursive_listdir(path):
    fqs = []
    files = os.listdir(path)
    for file in files:
        file_path = os.path.join(path, file)
        if os.path.isfile(file_path):
            fqs.append(file_path)
        elif os.path.isdir(file_path):
          recursive_listdir(file_path)
    fqs.sort()
    return fqs
def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            print(f"************ capture {suffix} as suffix file in {dir} *************")
            captures.append(os.path.join(dir, file))
    if captures:
        print(f'************ captured file ************** \n{captures}')
    return captures

def pairend2contig(path, threads, ref):
    """
    require BWA, picard megahit program,please pre-install them
    """
    megahit_sample_dir = []
    contigs = []
    bwa_bams = []
    asm_cmds = []
    map_sort_cmds = []
    check_dir("./00_megahit")
    check_dir(f"00_bwa_mem_out")
    check_dir("./00_megahit_log")
    fqs = recursive_listdir(path)
    for i  in range(0,len(fqs)-1,2):
        fq1,fq2 = fqs[i], fqs[i+1]
        fq1_name = basename(fq1)
        # Split by "r1" or "R1"
        sample_parts = re.split(r'_?R1.fastq|_?R1.gz|_?r1.gz|_?1.fastq.gz|_?1.fq.gz|_?R1.fq.gz|_?R1.fastq.gz|_?r1.fq.gz|_?r1.fastq.gz|_?R1.clean.fastq.gz|_?R1.clean.fq.gz', fq1_name)
        sample = sample_parts[0]
        print(f"please ensure the paired end data, R1: {fq1} ; R2: {fq2} sample name extract is {sample}")
        assembly_cmd = f"megahit -t {threads} -1 {fq1} -2 {fq2} -o 00_megahit/{sample} --out-prefix {sample} 1>00_megahit_log/{sample}.err 2>00_megahit_log/{sample}.log ; rm -r 00_megahit/{sample}/intermediate_contigs" 
        asm_cmds.append(assembly_cmd)
        megahit_sample_dir.append(f"00_megahit/{sample}")
        contigs.append(f"00_megahit/{sample}/{sample}.contigs.fa")
        bwa_cmd = f"bwa mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:{sample}\\tPU:run' -t {threads} {ref} {fq1} {fq2} | samtools view -buhS | samtools sort -@ 10 -o 00_bwa_mem_out/{sample}_sorted.bam && samtools index 00_bwa_mem_out/{sample}_sorted.bam -@ 10"
        picard_cmd = f"picard MarkDuplicates -I 00_bwa_mem_out/{sample}_sorted.bam -O 00_bwa_mem_out/{sample}.dedup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M 00_bwa_mem_out/{sample}.dedup.metrics && rm {sample}_sorted.bam"
        #map_sort_cmd = bwa_cmd + " && " + picard_cmd
        #map_sort_cmds.append(map_sort_cmd)
        map_sort_cmds.append(bwa_cmd)
        bwa_bams.append(f'00_bwa_mem_out/{sample}_sorted.bam')
        #bwa_bams.append(f'00_bwa_mem_out/{sample}_dedup.bam')
    return asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir

if __name__ == "__main__":
    import pyfiglet
    def print_large_SVInDel():
        ascii_art = pyfiglet.figlet_format("Pop SVGT", font="slant")
        print(ascii_art)
    print_large_SVInDel()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="This is SVGT in population working flow", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file or indexed bam files of hifi")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file or indexed bam files of ont")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads files or indexed bam files of pb")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly contig/chromosome level or indexed bam files of samples ")
    parser.add_argument("-win", "--window",default=500,type=int,help="Window size to parse signal, this parameter is to merge fragment SV due to aligment, default window size is 500bp to merge cluster the SV with a sv_start less than window size")
    parser.add_argument("-diploid", "--diploid", help="for diploid resolved assembly, to get a phased genotype please provide table list each line in format hap1\thap2\tSampleName in cr directory")
    parser.add_argument("-polyploid", "--polyploid", help="for polyploid haplotype resolved assembly like potato(4 haplotype assemblies available), to get a merge genotype of samples please provide table list each line in format hap1\thap2\thap3\thapn\tSampleName")


    parser.add_argument("-w", "--max_workers", default = 4,type=int, help="the max workers thread pool excutor, 4 means run 4 samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome ")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=40, help= "The min length of SV ")
    parser.add_argument("-M",  "--max", default=10000000, help= "The max length of  SV ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="no", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=1,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-csv",  "--csv",default=0.15, type=float, help= "the percent of reads that support a candidate SV (0.15 means at a depth 20X region, a SV signal should have at least 3 reads support, this parameter is for the variaty depth of hifi/ont/pb samples")
    parser.add_argument("-nreads",  "--nreads", type=int, help= "the number of reads to support a candidate SV (SV signal should have at least numbers reads support, this parameter is for the various depth of hifi/ont/pb samples")
    parser.add_argument("--num_hap", "--num_hap", default=2, type=int, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    parser.add_argument("-msv","--msv_mode",default="no", help= "In msv mode signals of INS,DEL,INV,DUP,TRA will captured from ont/hifi/pb, while for assemble contig from short reads or genome lelve samples we detect SVInDel Only. If no hifi or ont or pacbio data is provided, please setting -msv no, PSVGT will detect SVInDel Only")
    parser.add_argument("-lr_homo_rate", dest="lr_homo_rate",default=0.75, type=float, help="to determine a homozygous site, if 0.75 of the local mapping signal suport the sv the genotyoe will be 1/1, for species like polyploid potato we suggest 0.8")
    parser.add_argument("-lr_ref_rate", dest="lr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 0.10, the genotype will be 0/0")
    parser.add_argument("-sr_homo_rate", dest="sr_homo_rate",default=0.65, type=float, help="to determine a homozygous site, if 0.65 of the local mapping signal suport the sv the genotyoe will be 1/1, you can lower down the value if your specise have a low heterozygous rate")
    parser.add_argument("-sr_ref_rate", dest="sr_ref_rate",default=0.10, type=float, help="to determine reference allele, in a 100X data, if suport of local signal less than 10 percent, the genotype will be 0/0, you can increase the value to filter putative heterozygous SV")

    parser.add_argument("-span", dest="span", default=50, type=int, help="heterzygous evdence, a read (maping start - 50) < breakpoint < (mapping end - 50) will be taken as span the breakpoint, for 150 bp reads we suggest 50, for 125bp may be 45 will be better")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
    fai = pd.read_csv(f'{args.refGenome}.fai', sep='\t',index_col=None,header=None)
    genome_size = 0
    for key,seq in fa.items():
        genome_size += len(seq)
    if args.srdir:
        asm_cmds, map_sort_cmds, contigs, bwa_bams, megahit_sample_dir  = pairend2contig(args.srdir, args.threads, args.refGenome)
        ## store bwa bam to list file
        with open(f"{args.outdir}/bwa_bams.txt", 'w') as bamlst:
            for bam in bwa_bams:
                bamlst.write(f'{bam}\n')
        ## remove the done jobs ##
        for i in range(len(asm_cmds) - 1, -1, -1):  # Iterate backwards
            if os.path.isfile(contigs[i]):  # Replace with your condition
                del asm_cmds[i]
            else:
                rm_dir(megahit_sample_dir[i]) ## to avoid re assmbely and err
        for i in range(len(bwa_bams)-1,-1,-1):
            if os.path.isfile(bwa_bams[i]):
                del map_sort_cmds[i]
        if args.breaker=="yes":
            bwa_index_suffix = ["ann", "pac", "amb", "bwt", "sa"]
            should_index = 0
            for suffix in bwa_index_suffix:
                if not os.path.isfile(f"{args.refGenome}.{suffix}"):
                    should_index += 1
            if should_index !=0:
                run_command(f"bwa index {args.refGenome}")
            run_cmds = asm_cmds + map_sort_cmds
        else:
            run_cmds = asm_cmds
        if len(run_cmds) >0:
            print(f"*************** running minimap and bwa cmd are {run_cmds} *******************")
            for cmd in run_cmds:
                print(cmd, file=all_log)
        # Using ThreadPoolExecutor to run command1 first
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures0 = [executor.submit(execute_commands, cmd) for cmd in run_cmds]
            results0 =[]
            for future in as_completed(futures0):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results0.append((stdout,stderr,returncode,cmd))
    
    done_analysor = file_capture(args.outdir, ".record.txt")
    def run_clu2fil_cmd(clu2fil_cmd):
        try:
            subprocess.run(clu2fil_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {clu2fil_cmd}\n{e}")
    
    def process_single_contig(contig, dtype, args, PSVGT, fai):
        try:
            ## final chromsome result ##
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"
            ## skip samples ##
            if os.path.exists(done_name):
                print(f"Skipping processed contig: {basename(contig)}")
                return

            ## minimaping and signal detect ##
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py '
                f'-i {contig} -dtype {dtype} -r {args.refGenome} '
                f'-m {args.min} -maq {maq} -o {args.outdir} '
                f'-minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            ## chromosome pool run ##
            with multiprocessing.Pool(processes=args.max_workers) as pool:
                if dtype in ['hifi','ont','pb'] and args.csv and args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} --nreads {args.nreads}'
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                
                elif dtype in ['hifi','ont','pb'] and args.csv and not args.nreads:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--rate_depth {args.csv} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                elif dtype in ['hifi','ont','pb'] and args.nreads and not args.csv:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--nreads {args.nreads} '
                        f'--window {args.window} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov > {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]
                else:
                    clu2fil_cmds = [
                        f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                        f'-f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt '
                        f'-dtype {dtype} -s 800 -M {args.max} '
                        f'--b {args.outdir}/0_tmp_{basename(contig)}.bam '
                        f'--window {args.window} '
                        f'--num_hap {args.num_hap} '
                        f'--cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov >{args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                        for chrom in fai[0]
                    ]

                pool.map(run_clu2fil_cmd, clu2fil_cmds)
            ## ACC SV ##
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/0_tmp_{basename(contig)} '
                f'-fai {args.refGenome}.fai -M {args.max} --nrate {args.csv}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(contig)}: {str(e)}")

    def add_commands4fq(files, dtype):
        for contig in files:
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}_{fai[0].iloc[-1]}.record.txt"  ## last chrom
            if done_name not in done_analysor:
                maq = min(args.maq, 60)
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.PSVGT_raw2Signal.py -i {contig} -dtype {dtype} -r {args.refGenome} -m {args.min}  -maq {maq} -o {args.outdir} -minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
                try:
                    subprocess.run(signal_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running signal command: {signal_cmd}\n{e}")
                    continue

                clu2fil_cmds = []
                for chrom in fai[0]:
                    clu2fil_cmd = f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py -f {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max} --b {args.outdir}/0_tmp_{basename(contig)}.bam --cov {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.cov --window {args.window} --num_hap {args.num_hap} {args.outdir}/0_tmp_{basename(contig)}_{chrom}.record.txt.log'
                    clu2fil_cmds.append(clu2fil_cmd)
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                pool.map(run_clu2fil_cmd, clu2fil_cmds)
                pool.close()
                pool.join()

                ACC_SV_cmd = f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/0_tmp_{basename(contig)} -fai  {args.refGenome}.fai --nrate {args.csv}'
                try:
                    subprocess.run(ACC_SV_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running ACC_SV command: {ACC_SV_cmd}\n{e}")
    
    def add_commands4fq(files, dtype):
        max_workers = min(len(files), args.max_workers)
        if max_workers >0:
            print(f'{max_workers} workers')
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        process_single_contig,
                        contig, dtype, args, PSVGT, fai
                    ) for contig in files
                ]

        
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Thread error: {str(e)}")

    def add_commands4bam(files, dtype):
        for bam in files:
            done_name = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if done_name not in done_analysor:
                maq = min(args.maq, 60)  ## biger than illumina breakpoints quality
                signal_cmd = f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py -b {bam} -o {args.outdir}/{basename(bam)} -m {args.min} -maq {args.maq} -dtype {dtype} -msv {args.msv_mode} -fai {args.refGenome}.fai'
                try:
                    subprocess.run(signal_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running signal command: {signal_cmd}\n{e}")
                    continue
                clu2fil_cmds = []
                for chrom in fai[0]:
                    clu2fil_cmd = f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py -f {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt -dtype {dtype} -s 800 -M {args.max}  --b {bam} --cov {args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt.cov --window {args.window} --num_hap {args.num_hap} >{args.outdir}/{basename(bam).replace(".bam", "")}_{chrom}.record.txt.log'
                    clu2fil_cmds.append(clu2fil_cmd)
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                pool.map(run_clu2fil_cmd, clu2fil_cmds)
                pool.close()
                pool.join()

                ACC_SV_cmd = f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix {args.outdir}/{basename(bam).replace(".bam", "")} -fai  {args.refGenome}.fai --nrate {args.csv}'
                try:
                    subprocess.run(ACC_SV_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running ACC_SV command: {ACC_SV_cmd}\n{e}")
    
    def process_single_bam(bam, dtype, args, PSVGT, fai):
        try:## last chromosome result ##
            done_flag = f"{args.outdir}/0_tmp_{basename(bam)}_{fai[0].iloc[-1]}.record.txt"
            if os.path.exists(done_flag):
                print(f"Skipping processed BAM: {basename(bam)}")
                return

            # Phase 1: signal
            maq = min(args.maq, 60)
            signal_cmd = (
                f'python {PSVGT}/PSV_Signal/0.Signal4bam_PSVGT.py '
                f'-b {bam} -o {args.outdir}/{basename(bam)} '
                f'-m {args.min} -maq {maq} -dtype {dtype} '
                f'-msv {args.msv_mode} -fai {args.refGenome}.fai'
            )
            subprocess.run(signal_cmd, shell=True, check=True)

            # Phase 2: parafly chrs
            base_prefix = basename(bam).replace(".bam", "")
            if args.nreads and args.csv and dtype in ['hifi', 'ont', 'pb']:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--nreads {args.nreads} --rate_depth {args.csv} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.csv and dtype in ['hifi', 'ont', 'pb'] and not args.nreads:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--rate_depth {args.csv} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            elif args.nreads and dtype in ['hifi', 'ont', 'pb'] and not args.csv:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--nreads {args.nreads}  '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]
            else:
                chrom_commands = [
                    f'python {PSVGT}/PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py '
                    f'-f {args.outdir}/{base_prefix}_{chrom}.record.txt '
                    f'-dtype {dtype} -s 800 -M {args.max} '
                    f'--window {args.window} '
                    f'--num_hap {args.num_hap} '
                    f'--b {bam} --cov {args.outdir}/{base_prefix}_{chrom}.record.txt.cov >{args.outdir}/{base_prefix}_{chrom}.record.txt.log'
                for chrom in fai[0]
                ]


            # dynamic cpu
            with multiprocessing.Pool(processes=min(len(fai[0]), os.cpu_count()//5)) as pool:
                pool.map(run_clu2fil_cmd, chrom_commands)

            # Phase 3: merge signal
            ACC_SV_cmd = (
                f'python {PSVGT}/PSV_Signal/1.ACCSV_Signal_Cluster.py '
                f'-preffix {args.outdir}/{base_prefix} '
                f'-fai {args.refGenome}.fai --nrate {args.csv}'
            )
            subprocess.run(ACC_SV_cmd, shell=True, check=True)

        except subprocess.CalledProcessError as e:
            print(f"Error processing {basename(bam)}: {str(e)}")

    def add_commands4bam(files, dtype):
        if len(files)>0:
            max_threads = min(len(files),args.max_workers)
            with ThreadPoolExecutor(max_workers=max_threads) as executor:
                task_func = partial(
                    process_single_bam,
                    dtype=dtype,
                    args=args,
                    PSVGT=PSVGT,
                    fai=fai
                )
                futures = {executor.submit(task_func, bam): bam for bam in files}
                for future in concurrent.futures.as_completed(futures):
                    bam_file = futures[future]
                    try:
                        future.result()
                        print(f"Success: {basename(bam_file)}")
                    except Exception as e:
                        print(f"Failed processing {basename(bam_file)}: {str(e)}")

    already_maps = []
    if args.srdir:
        add_commands4fq(contigs, "sr") #### the short reads assembly reads use hifi mode to mapping ####
    if args.hifidir:
        add_commands4fq(file_capture(args.hifidir, ".gz"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fastq"), "hifi")
        add_commands4fq(file_capture(args.hifidir, ".fq"), "hifi")
        add_commands4bam(file_capture(args.hifidir, ".bam"), "hifi")
        already_maps += file_capture(args.hifidir, ".bam")
    if args.ontdir:
        add_commands4fq(file_capture(args.ontdir, ".gz"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fastq"), "ont")
        add_commands4fq(file_capture(args.ontdir, ".fq"), "ont")
        add_commands4bam(file_capture(args.ontdir, ".bam"), "ont")
        already_maps += file_capture(args.ontdir, ".bam")

    if args.pbdir:
        add_commands4fq(file_capture(args.pbdir, ".gz"), 'pb')
        add_commands4fq(file_capture(args.pbdir, ".fastq"), "pb")
        add_commands4fq(file_capture(args.pbdir, ".fq"), "pb")
        add_commands4bam(file_capture(args.pbdir, ".bam"),'pb')
        already_maps += file_capture(args.pbdir, ".bam")
    if args.crdir:
        add_commands4fq(file_capture(args.crdir, ".fasta"), "cr")
        add_commands4fq(file_capture(args.crdir, ".fa"), "cr")
        add_commands4bam(file_capture(args.crdir, ".bam"), "cr")
        add_commands4fq(file_capture(args.crdir, ".gz"), "cr")
        already_maps += file_capture(args.crdir, ".bam")


    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT}/PSV_Signal/1.PSV_signal_cluster.py -d {args.outdir} -s 50")
    
    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam") + already_maps
    mapinfo_files.sort()
    gt_cmds = []
    #### not repeat the genotyping #####
    doneGT  = file_capture(f"{args.outdir}", "_genotype.txt")
    if len(doneGT) > 0:
        print(f"The sample file {doneGT} has been genotype before, if you want to update genotyping results, please remove the file in the list")
    for mapinfo_file in mapinfo_files:
        if_done_name = f"{args.outdir}/2_tmp_{basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')}_genotype.txt"
        if if_done_name not in doneGT:
            print(if_done_name)
            acc_name = basename(mapinfo_file).replace('.bam', '').replace('0_tmp_', '')

            cmd = f"python {PSVGT}/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {mapinfo_file} -m {args.maq} -lr_homo_rate {args.lr_homo_rate} -lr_ref_rate {args.lr_ref_rate}  -n {acc_name} -o {args.outdir} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf {args.refGenome}.fai"
            print(cmd)
            gt_cmds.append(cmd)
    if len(gt_cmds) >0:
        with open ("gt_by_longseq_log.txt", 'w') as longseq_gt_log:
            max_workers = min(len(gt_cmds), 5)
            with ThreadPoolExecutor(max_workers= max_workers) as executor: #### in the way cpu 104 will reach 104 * 2 ####
                futures2 = [executor.submit(execute_commands, cmd) for cmd in gt_cmds ]
                results2 =[]
                for future in as_completed(futures2):
                    try:
                        stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                        if returncode != 0:
                            print(f"Error executing command: {cmd}" , file=longseq_gt_log)
                            print(stderr, file=longseq_gt_log)
                        else:
                            print(f"{cmd}", file=longseq_gt_log)
                            print(stdout,file=longseq_gt_log)
                        results2.append((stdout,stderr,returncode,cmd))
                    except Exception as e:
                        print(f'An error occurred: {e}', file=longseq_gt_log)
    else:
        print(f"all samples has been genotype before, if you want to repeat genotype please remove the 2_tmp_XXX_genotype.txt files in the {args.outdir}")

    ###################### haplotype resoved assembly genotype phased ###################
    if args.diploid:
        print("*************** diploid calling **************")
        with open(args.diploid, 'r') as f:
            lines = f.readlines()
        for hh in lines:
            h1 = hh.strip().split("\t")[0]
            h1 = "2_tmp_" + h1.replace(".bam", "") + "_genotype.txt" if  h1[-4:] == ".bam" else "2_tmp_" + h1 + "_genotype.txt"
            h2 = hh.strip().split("\t")[1]
            h2 = "2_tmp_" + h2.replace(".bam", "") + "_genotype.txt" if  h2[-4:] == ".bam" else "2_tmp_" + h2 + "_genotype.txt"
            samplename =  hh.strip().split("\t")[2]
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_diploid_asm.py {args.outdir}/{h1} {args.outdir}/{h2} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f"***************** try to phased hap1: {args.outdir}/{h1} and hap2: {args.outdir}/{h2} to {args.outdir}/2_tmp_{samplename}_genotype.txt ******************")
            tab2vcf =    f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(tab2vcf)
    if args.polyploid:
        print("*************** polyploid genotype merging **************")
        with open(args.polyploid, 'r') as f:
            lines = f.readlines()
        path_haps = ""
        for hh in lines:
            path_haps = ""
            haps = hh.strip().split("\t")
            total_haps = len(haps) - 1
            samplename = haps[-1]
            for i in range(total_haps):
                hap = haps[i] 
                hap_file = "2_tmp_" + hap.replace(".bam", "") + "_genotype.txt" if  hap[-4:] == ".bam" else "2_tmp_" + hap + "_genotype.txt"
                path_haps += f'{args.outdir}/{hap_file} '
            print(f"*************************** phased polyploid {samplename} ***************************** ")
            phased_cmd = f'python {PSVGT}/PSV_Genotyper/phased_polyploid_genome_gt.py {path_haps} {args.outdir}/2_tmp_{samplename}_genotype.txt'
            print(f'********************** phased polyploid ************************\n{phased_cmd}')
            vcf_cmd = f'python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{samplename}_genotype.txt {args.outdir}/2_tmp_{samplename}_genotype.vcf {args.refGenome}.fai'
            run_command(phased_cmd)
            run_command(vcf_cmd)


    ###################### calling illumina breaker #########################
    if args.breaker == "yes":
        breaker_gt_cmds = []
        bams = file_capture(f"00_bwa_mem_out", ".bam")
        for bam in bams:
            sampleID = basename(bam)[:-4]
            bpgt_cmd =  f"python {PSVGT}/PSV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_Candidate_Record.txt -mapf {bam} -maq {args.maq} -span {args.span} -s 50 -n {sampleID} -o {args.outdir} -homo_rate {args.sr_homo_rate} -ref_rate {args.sr_ref_rate} && python {PSVGT}/PSV_Genotyper/SVGT_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf {args.refGenome}.fai"
            print(bpgt_cmd)
            breaker_gt_cmds.append(bpgt_cmd)
        with open("gt_sv_by_bwa_bam_log.txt", 'w') as sr_bpgt_log:
            max_workers =  min(len(breaker_gt_cmds), args.max_workers)
            sr_svgt_results = threading_cmd(breaker_gt_cmds, sr_bpgt_log, worker=max_workers)
    
    
    ###################### merging the vcf files in vcf.list #######################
    print(f'###################### merging the vcf files  #######################')
    merge_cmd = f"python {PSVGT}/PSV_Genotyper/merge_vcf_by_pandas.py -d {args.outdir} -o {args.outdir}/PSVGT_all.vcf2"
    print(merge_cmd)
    run_command(merge_cmd)
    if args.popInDel == "yes":
        run_command(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} --min 80 --max 600 --frank 400 --maf 0.01 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
        print(f"python {PSVGT}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 80 600 400 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
    ###################### Annotaion SVInDel For  Population #########################
    final_gt = f"{args.outdir}/PSVGT_all.vcf2.SVInDel"
    if args.gff:
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Annotation.py -g {args.gff} -s  {args.outdir}/PSVGT_all.vcf2.SVInDel -m ID -c Parent -o {args.outdir}/SVInDels_Lead_Gene_Variant.txt &")
        run_command(f"python {PSVGT}/SVInDel_Anno/SV_Features_Position.py {args.gff} {final_gt}_tmp.tab {args.outdir}/PSVInDel")
    if args.popcaps == "yes":
        popcaps_cmds = []
        out_lst = open("bam_lst.txt", "w")
        contig_bams = file_capture(args.outdir, ".bam")
        for bam in contig_bams:   ################################# use bwa bam or minimap bam   ########################################
            print(bam, file=out_lst)
        out_lst.close()
        for chrom in fa.keys():
            popcaps_cmd  = f"samtools mpileup -b bam_lst.txt -q 55 -Q 30 -r {chrom} |python {PSVGT}/CapsPop/mpileup_stdin4popcasp.py > PopCaps_{chrom}_input.txt && python {PSVGT}/CapsPop/pop_maf0.05_caps.py {PSVGT}/CapsPop/common_enzyme.list {args.refGenome} PopCaps_{chrom}_input.txt Out_PopCaps_{chrom}_maf0.05.txt 300 && rm PopCaps_{chrom}_input.txt"
            popcaps_cmds.append(popcaps_cmd)
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures4 = [executor.submit(execute_commands, cmd) for cmd in popcaps_cmds ]
            results4 =[]
            for future in as_completed(futures4):
                stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
                if returncode != 0:
                    print(f"Error executing command: {cmd}" , file=all_log)
                    print(stderr, file=all_log)
                else:
                    print(f"Output of command '{cmd}':", file=all_log)
                    print(stdout,file=all_log)
                results4.append((stdout,stderr,returncode,cmd))
    end_t = time()
    print(f'{"*" * 20} Total Time Cost In SVGT Program: {end_t - start_t}s\t{"*" * 20}')

```

## ./SVInDel_Primer
```shell
total 24
-rw-r--r-- 1 lgb xinwang  604 May 22 16:55 input.txt
-rwxr-xr-x 1 lgb xinwang 8084 May 22 16:55 SVInDel_GenotypeTab2primer.py
-rwxr-xr-x 1 lgb xinwang 9035 Jul  4 11:50 vcf2primer.py

```
### ./SVInDel_Primer/SVInDel_GenotypeTab2primer.py
```python
import argparse
import primer3
import subprocess
import gzip
from collections import deque
from gzip import BadGzipFile
import re

def reverse_complement(primer):
    """Return the reverse complement of a DNA primer."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(primer))

def find_primer_positions(sequence, primer, strand):
    """
    Find the start position of the forward primer and the end position of the reverse primer in a sequence.
    Parameters:
    - sequence (str): The DNA sequence to search within.
    - primer_f (str): The forward primer.
    - primer_r (str): The reverse primer.
    Returns:
    - tuple: (start position of F, end position of R) or (None, None) if not found.
    """
    if strand == "+":
        start_f = sequence.find(primer)
        return start_f
    else:
        primer_r_rev_comp = reverse_complement(primer)
        start_r = sequence.find(primer_r_rev_comp)
        end_r = start_r + len(primer_r_rev_comp) - 1  
        if not end_r:
            print(primer, sequence)
        return end_r

def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

def my_primer(info, seq, seqid, seq_args_left, seq_args_right, global_args, sv_size):
    results_left = primer3.bindings.design_primers(seq_args=seq_args_left,global_args=global_args)
    primers_L=[]
    for i in range(results_left['PRIMER_RIGHT_NUM_RETURNED']):
        left_seq = results_left[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right_seq = results_left[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        product_size = results_left[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        primers_L.append((left_seq))

    # Run the primer3-py primer design algorithm with the updated global and sequence-specific arguments for the right subregion
    results_right = primer3.bindings.design_primers(seq_args=seq_args_right,global_args=global_args)
    primers_R=[]
    for i in range(results_right['PRIMER_LEFT_NUM_RETURNED']):
        left_seq = results_right[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right_seq = results_right[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        product_size = results_right[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        primers_R.append((right_seq))
    if "N" not in seq: 
        if len(primers_L) >0  and len(primers_R) > 0:
            pos_f = []
            pos_f_list = []
            for f in primers_L:
                pos = find_primer_positions(seq, f, "+")
                pos_f.append(f"{f}_start_at_{str(pos)}")
                pos_f_list.append(pos)
            best_f_index = pos_f_list.index(max(pos_f_list))
            best_forward_primer = pos_f[best_f_index]
            
            pos_r = []
            pos_r_list = []
            for r in primers_R:
                pos = find_primer_positions(seq, r, "-")
                pos_r.append(f"{r}_star_at_{str(pos)}")
                pos_r_list.append(pos)
            best_r_index = pos_r_list.index(min(pos_r_list))
            best_reverse_primer = pos_r[best_r_index]
            miniPCR1 = min(pos_r_list) - max(pos_f_list)
            miniPCR2 = miniPCR1 - sv_size
            outinfo = "\t".join(info)
            print(f"{outinfo}\t{seqid}\t{seq}\tForward primers: {set(pos_f)}\tReverse primers: {set(pos_r)}\tMiniPCR_Pair_Primer:{(best_forward_primer,best_reverse_primer)}\tRefPCRsize:{miniPCR1}\tQueryPCRsize:{miniPCR2}")
# Define the argparse parameters
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="the SVInDel record first col could be anything, 2nd col is chr, 3rd is start, 4th col is end, 5th col is target sv size")
parser.add_argument("ref", help="the reference genome to extract sequence design primers")
args = parser.parse_args()

global_args={
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER': 10,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_RIGHT_PRIMER': 10,
    'PRIMER_OPT_SIZE': 21,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 53.0,
    'PRIMER_MAX_TM': 65.0,
    'PRIMER_MIN_GC': 40.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
}

def is_valid_dna_sequence(sequence):
    """
    Check if the DNA sequence is valid.    
    Parameters:
    sequence (str): The DNA sequence to check.
    Returns:
    bool: True if the sequence is valid, False otherwise.
    """
    # Define valid nucleotides
    valid_nucleotides = {'A', 'G', 'C', 'T', 'N'}
    # Use set comparison to check for any invalid characters
    return all(char in valid_nucleotides for char in sequence)

# Read the input file and design primers for each sequence
with open(args.input_file) as f:
    lines = f.readlines()
## get seq fields ##
fa = readfa2Dict(args.ref)
header = ""
for line in lines:
    if "#" in line:
        print(line.strip(),"Seqid\tSeq.\tForward primers\tReverse primers\tMiniPCR Product Primer\tRef PCR size\tQuery PCR size",sep="\t") 
        continue
    info = line.strip().split("\t")
    chrom, startI, endI, sv_size =info[1], int(info[2]), int(info[3]), int(info[4])
    if startI+20 < 600:
        continue
    InDelLocal = f'{chrom}:{startI - 600}-{endI + 600}'
    seq = fa[chrom][startI - 600 : endI + 600].upper()
    pos = 600
    seqid = f'{chrom}:{startI - 600}-{endI + 600}'
    left_subregion_start = 0
    left_subregion_end = pos - 1
    right_subregion_start = len(seq) - pos + 1
    right_subregion_end = len(seq)

    # Define the seq_args dictionary to include the left subregion
    seq_args_left={
            'SEQUENCE_ID': seqid+'_left',
            'SEQUENCE_TEMPLATE': seq[left_subregion_start:left_subregion_end],
            'SEQUENCE_INCLUDED_REGION': [0, left_subregion_end - left_subregion_start],
        }

    # Define the seq_args dictionary to include the right subregion
    seq_args_right={
            'SEQUENCE_ID': seqid+'_right',
            'SEQUENCE_TEMPLATE': seq[right_subregion_start:right_subregion_end],
            'SEQUENCE_INCLUDED_REGION': [0, right_subregion_end - right_subregion_start],
        }
    #print("*********************************************\n",seq_args_right)
        # Design primers for the left and right subregions and print the results
    if is_valid_dna_sequence(seq):
        my_primer(info, seq, seqid, seq_args_left, seq_args_right, global_args, sv_size)

```

### ./SVInDel_Primer/vcf2primer.py
```python
import argparse
import primer3
import subprocess
import gzip
from collections import deque
from gzip import BadGzipFile
import re

def reverse_complement(primer):
    """Return the reverse complement of a DNA primer."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(primer))

def find_primer_positions(sequence, primer, strand):
    """
    Find the start position of the forward primer and the end position of the reverse primer in a sequence.
    Parameters:
    - sequence (str): The DNA sequence to search within.
    - primer_f (str): The forward primer.
    - primer_r (str): The reverse primer.
    Returns:
    - tuple: (start position of F, end position of R) or (None, None) if not found.
    """
    if strand == "+":
        start_f = sequence.find(primer)
        return start_f
    else:
        primer_r_rev_comp = reverse_complement(primer)
        start_r = sequence.find(primer_r_rev_comp)
        end_r = start_r + len(primer_r_rev_comp) - 1  
        if not end_r:
            print(primer, sequence)
        return end_r

def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

def my_primer(info, seq, seqid, seq_args_left, seq_args_right, global_args, sv_size):
    results_left = primer3.bindings.design_primers(seq_args=seq_args_left,global_args=global_args)
    primers_L=[]
    for i in range(results_left['PRIMER_RIGHT_NUM_RETURNED']):
        left_seq = results_left[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right_seq = results_left[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        product_size = results_left[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        primers_L.append((left_seq))

    # Run the primer3-py primer design algorithm with the updated global and sequence-specific arguments for the right subregion
    results_right = primer3.bindings.design_primers(seq_args=seq_args_right,global_args=global_args)
    primers_R=[]
    for i in range(results_right['PRIMER_LEFT_NUM_RETURNED']):
        left_seq = results_right[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right_seq = results_right[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        product_size = results_right[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        primers_R.append((right_seq))
    if "N" not in seq: 
        if len(primers_L) >0  and len(primers_R) > 0:
            pos_f = []
            pos_f_list = []
            for f in primers_L:
                pos = find_primer_positions(seq, f, "+")
                pos_f.append(f"{f}_start_at_{str(pos)}")
                pos_f_list.append(pos)
            best_f_index = pos_f_list.index(max(pos_f_list))
            best_forward_primer = pos_f[best_f_index]
            
            pos_r = []
            pos_r_list = []
            for r in primers_R:
                pos = find_primer_positions(seq, r, "-")
                pos_r.append(f"{r}_star_at_{str(pos)}")
                pos_r_list.append(pos)
            best_r_index = pos_r_list.index(min(pos_r_list))
            best_reverse_primer = pos_r[best_r_index]
            miniPCR1 = min(pos_r_list) - max(pos_f_list)
            miniPCR2 = miniPCR1 - sv_size
            outinfo = "\t".join(info)
            print(f"{outinfo}\t{seqid}\t{seq}\tForward primers: {set(pos_f)}\tReverse primers: {set(pos_r)}\tMiniPCR_Pair_Primer:{(best_forward_primer,best_reverse_primer)}\tRefPCRsize:{miniPCR1}\tQueryPCRsize:{miniPCR2}")
# Define the argparse parameters
parser = argparse.ArgumentParser("MAF InDel Marker Analysis For PSVGT InDel VCF", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("vcf_file", help="the vcf file from PSVGT")
parser.add_argument("ref", help="the reference genome to extract sequence design primers")
parser.add_argument("--min", help="the min size of InDel", default=50,type=int)
parser.add_argument("--max", help="the max size of InDel", default=600, type=int)
parser.add_argument("--frank", help="the exten size from the target region", default=300, type=int)
parser.add_argument("--maf", help="the minor allele frequency setting for the InDel Primer", default=0.01, type=float)

args = parser.parse_args()

global_args={
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER': 10,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_RIGHT_PRIMER': 10,
    'PRIMER_OPT_SIZE': 21,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 53.0,
    'PRIMER_MAX_TM': 65.0,
    'PRIMER_MIN_GC': 40.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
}

def is_valid_dna_sequence(sequence):
    """
    Check if the DNA sequence is valid.    
    Parameters:
    sequence (str): The DNA sequence to check.
    Returns:
    bool: True if the sequence is valid, False otherwise.
    """
    # Define valid nucleotides
    valid_nucleotides = {'A', 'G', 'C', 'T', 'N'}
    # Use set comparison to check for any invalid characters
    return all(char in valid_nucleotides for char in sequence)

# Read the input file and design primers for each sequence
with open(args.vcf_file) as f:
    lines = f.readlines()
## get seq fields ##
fa = readfa2Dict(args.ref)
header = ""
for line in lines:
    if "#" in line:
        print(line.strip(),"Seqid\tSeq.\tForward primers\tReverse primers\tMiniPCR Product Primer\tRef PCR size\tQuery PCR size",sep="\t") 
        continue
    info = line.strip().split("\t")
    SVID = info[2] ## Chr1:22023-22242_220
    chrom  = SVID.split(":")[0]
    startI = int(SVID.split(":")[1].split("-")[0])
    endI   = int(SVID.split(":")[1].split("-")[1].split("_")[0])
    sv_size = int(SVID.split(":")[1].split("-")[1].split("_")[1])
    
    alt_maf = round(info[9:].count("1/1") / len(info[9:]), 2)
    if len(info) > 10: ## more than one samples
        if alt_maf <= args.maf:
            continue
        if alt_maf > 1 - args.maf:
            continue
    if sv_size >= args.max:
        continue
    if sv_size <= args.min:
        continue
    if startI+20 < 1000:
        continue
    InDelLocal = f'{chrom}:{startI - args.frank}-{endI + args.frank}'
    seq = fa[chrom][startI - args.frank : endI + args.frank].upper()
    pos = args.frank
    seqid = f'{chrom}:{startI - args.frank}-{endI + args.frank}'
    left_subregion_start = 0
    left_subregion_end = pos - 1
    right_subregion_start = len(seq) - pos + 1
    right_subregion_end = len(seq)

    # Define the seq_args dictionary to include the left subregion
    seq_args_left={
            'SEQUENCE_ID': seqid+'_left',
            'SEQUENCE_TEMPLATE': seq[left_subregion_start:left_subregion_end],
            'SEQUENCE_INCLUDED_REGION': [0, left_subregion_end - left_subregion_start],
        }

    # Define the seq_args dictionary to include the right subregion
    seq_args_right={
            'SEQUENCE_ID': seqid+'_right',
            'SEQUENCE_TEMPLATE': seq[right_subregion_start:right_subregion_end],
            'SEQUENCE_INCLUDED_REGION': [0, right_subregion_end - right_subregion_start],
        }
    #print("*********************************************\n",seq_args_right)
        # Design primers for the left and right subregions and print the results
    if is_valid_dna_sequence(seq):
        my_primer(info[0:5]+info[9:]+[str(alt_maf)], seq, seqid, seq_args_left, seq_args_right, global_args, sv_size)

```

## ./PSV_Signal
```shell
total 220
-rwxr-xr-x 1 lgb xinwang 38464 Jul 27 11:12 0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py
-rwxr-xr-x 1 lgb xinwang 35047 Jul  5 15:30 0.KLookCluster_LocalDepthAdaptive.py
-rwxr-xr-x 1 lgb xinwang 34478 May 23 10:51 0.KLookCluster_SumAlleles4LocalDepthAdaptive.py
-rwxr-xr-x 1 lgb xinwang  3824 May 25 16:57 0.PSVGT_raw2Signal.py
-rwxr-xr-x 1 lgb xinwang  2198 Jun 11 18:12 0.Signal4bam_PSVGT.py
-rw-r--r-- 1 lgb xinwang  2205 May 22 16:55 0.Signal4bam_PSVGT_v2.py
-rwxr-xr-x 1 lgb xinwang 10950 Jul 10 22:42 1.ACCSV_Signal_Cluster.py
-rwxr-xr-x 1 lgb xinwang  8735 May 22 16:55 1.PSV_signal_cluster.py
drwxr-xr-x 3 lgb xinwang  4096 Sep  4 21:11 parse_sam_signal
drwxr-xr-x 2 lgb xinwang  4096 Jul 10 10:35 __pycache__
-rwxr-xr-x 1 lgb xinwang 18167 May 22 16:55 sub_PSVGT_Signaling.py
-rwxr-xr-x 1 lgb xinwang  1683 May 22 16:55 Sub_readfa2Dict.py
-rwxr-xr-x 1 lgb xinwang 38978 Jul 10 10:30 sub_Signal4bam_PSVGT.py

```
### ./PSV_Signal/0.PSVGT_raw2Signal.py
```python
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
    if args.msv == "yes":
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

```

### ./PSV_Signal/sub_PSVGT_Signaling.py
```python
import pandas as pd
import re
from tqdm import tqdm
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from os import system
def inOut(args):
    samfile = pysam.AlignmentFile(args.sam)
    lines = samfile.fetch()
    if args.out[-4:] == ".bam":
        args.out = args.out.replace('.bam', '')
    indel_out = open(args.out, 'w')
    cov_out = open(args.out + ".cov", 'w')
    supp_sv_out = open(f"{args.out}.suppAlign", 'w')
    return indel_out, lines, cov_out, supp_sv_out

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
                ## should i lower down the maq to capture more signal ???
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
                    ins_seq = query_seq[query + 1: query + length -1]
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
            ## 1 or 2
            if 1<= len(supp_dict[readname]) <= 20:
                supp_svsignal = segmentsv4lr(supp_dict[readname],minLen, maxLen, chromosome_list)
    return svInDels, covinfo , supp_svsignal
def svInDel4asm(line, minLen, min_maq):
    if line.flag == 4 or line.mapping_quality < min_maq:
        return None  # Skip header lines
    else:
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
                    ins_seq = query_seq[query + 1: query + length -1]
                    results.append(f"{target_chr}\t{query_chr}\t{ref}\t{ref + 1}\t{length}\t{maq}\t{target_chr}:{ref}-{ref+1}_INS={length}\tINS\t{ins_seq}\n")
                query += length
            #elif code in {"N", "S", "H"}:
            #    ref += length if code == "N" else 0
            #    query += length
        return results, f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n'  # Return results and coverage line

if __name__ == "__main__":
    import argparse
    from time import time
    parser = argparse.ArgumentParser("SV signal extract from sam file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Sigal Capture")
    IN.add_argument("-s", dest="sam",required=True, help="sam file from genome aln (minimap2 -ax asm5 -t 20 --secondary=no -o out.sam target.fa query.fa)")
    IN.add_argument("-o", dest="out",required=True, help="output InDel info file")
    IN.add_argument("-m", dest="min", help="InDel min length", default=50, type=int)
    IN.add_argument("-maq", dest="maq", help="the min mapping quality", default = 50, type=int)
    IN.add_argument("-dtype", dest="seqtype",required=True, help="the sequence type(ont,hifi,cr)")
    IN.add_argument("-M", dest="max", help="SV max length", default=10000000, type=int)
    IN.add_argument("-fai", dest="fai", help="chromsome fai index file", type=str)
    IN.add_argument("-msv", dest="msv", help="detecting complex sv(INV,DP,TRA,INS,DEL) from supplementary alignment", default="no", type=str)
    args = parser.parse_args()
    exe_start = time()
    chromosome_list = pd.read_csv(args.fai,header=None,sep="\t",index_col=None)[0].tolist()
    indel_out, lines, cov_out, supp_sv_out = inOut(args)
    if args.seqtype in ["ont","hifi",'pb']:
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(svInDel4lr, line, args.min, args.maq, args.max, args.msv, chromosome_list): line for line in lines}
            cov_lines = []  # To collect coverage lines
            for future in tqdm(as_completed(futures), total=len(futures)):
                result = future.result()
                if result is None:
                    continue
                svInDels, cov_line, supp_svsignal = result
                if svInDels:
                    for svInDel in svInDels:
                        indel_out.writelines(f'{svInDel}')
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")
                cov_lines.append(cov_line)
        # Write coverage data after processing all lines
        cov_out.writelines(cov_lines)
        indel_out.close()
        cov_out.close()
        supp_sv_out.close()
        exe_end = time()
    elif args.seqtype in ["sr", "cr"]:
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(svInDel4asm, line, args.min, args.maq): line for line in lines}
            cov_lines = []  # To collect coverage lines
            for future in tqdm(as_completed(futures), total=len(futures)):
                result = future.result()
                if result is None:
                    continue
                results, cov_line = result
                if results:
                    indel_out.writelines(results)
                cov_lines.append(cov_line)
        # Write coverage data after processing all lines
        cov_out.writelines(cov_lines)
        indel_out.close()
        cov_out.close()
        supp_sv_out.close()
        exe_end = time()
    print(f"{'*' * 40} done SV searching {'*' * 40}\ncost time: {exe_end - exe_start}s")


```

### ./PSV_Signal/sub_Signal4bam_PSVGT.py
```python
import pysam
import re
def process_chromosome(chromosome,chrom_size,chromosome_list, bamfile_path, minLen, maxLen, min_maq, SVsignal_out_path,dtype,msv):
    samfile = pysam.AlignmentFile(bamfile_path, 'rb')
    with open(f"{SVsignal_out_path}_{chromosome}.record.txt", 'w') as indel_out, open(f"{SVsignal_out_path}_{chromosome}.record.txt.cov", 'w') as cov_out, open(f"{SVsignal_out_path}_{chromosome}.record.txt.suppAlign", 'w') as supp_sv_out, open(f'{SVsignal_out_path}_{chromosome}.record.txt.depth','w') as depth_out:
        # Fetch all reads from the chromosome
        lines = samfile.fetch(chromosome)
        total_map_lens = 0
        for line in lines:
            if dtype in ['sr', 'cr']:
                result = svInDel4asm(line, minLen, min_maq,msv, maxLen, chromosome_list)
                if result is None:
                    continue
                results, cov_line, supp_svsignal = result
                if results:
                    indel_out.writelines(results)
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")   
                cov_out.writelines([cov_line])
            elif dtype in ['pb', 'ont', 'hifi']:
                svInDels, cov_line, supp_svsignal = svInDel4lr(line, minLen, min_maq, maxLen, msv, chromosome_list)
                if cov_line:
                    #print(cov_line)
                    total_map_lens += int(cov_line.split("\t")[4]) - int(cov_line.split("\t")[3])
                if svInDels:
                    indel_out.writelines(svInDels)
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")
                if cov_line:
                    cov_out.writelines(cov_line)
        depth = round( total_map_lens / chrom_size, 2 )
        depth_out.write(f'{chromosome}\t{chrom_size}\t{total_map_lens}\t{depth}\n')
        depth_out.close()
        indel_out.close()
        cov_out.close()
        supp_sv_out.close()

        if dtype == 'cr' and msv == 'yes':
            with open(f"{SVsignal_out_path}_{chromosome}.record.txt.suppAlign", 'a') as out:
                lines = samfile.fetch(chromosome)
                svs = segment_segment_sv(lines, minLen, maxLen, chromosome_list)
                for sv in svs:
                    out.writelines(f"{sv}\n")
            out.close()

def segmentsv4asm(supp_list, min_size, max_size, chromosome_list):
    """
    Collecting the supplementary alignment to identify the INS, DEL, DUP, INV, TRA signals
    LongRead align to reference that has "SA" will be recorded as dict as follow:
    { 'm64144': [['m64144',2048,'Db-Chr4',8949283,8949610,[0, 327, 11157],'43'],
                 ['m64144',0,   'Db-Chr4',8949649,8957962,[3162, 8322, 0],60,'AAXXXCTAATT'],
                 ['m64144',2064,'Db-Chr3',12497386,12502117,[5171, 4731, 1582],'60']]
                 }
    these code is referenced on DeBreak, the INS and DEL may be discarded in the future
    """
    if not any(int(supp[1]) <= 16 for supp in supp_list) or len(supp_list) <= 1:
        return []

    # Sort the supplementary list by chromosome and then by the start position
    supp_list_sorted = sorted(supp_list, key=lambda x: (x[2], int(x[3])))
    svsignal_supp = []

    def add_svcall(chrom, readname, start, end, size, maq, svid, sv_type, sequence=None):
        if maq == 60:
            sv_record = f"{chrom}\t{readname}\t{start}\t{end}\t{size}\t{maq}\t{svid}\t{sv_type}"
            if sequence:
                sv_record += f"\t{sequence}"
            svsignal_supp.append(sv_record)

    # Iterate through the sorted list, each smaller start position is treated as the primary segment
    for idx, primary_map in enumerate(supp_list_sorted):
        readname = primary_map[0]
        pri_chrom = primary_map[2]
        pri_start = int(primary_map[3])
        readseq = primary_map[7] if len(primary_map)==8 else "*"
        #print(f"********************************* flag: {primary_map[1]}*************************************")
        pri_flag = int(primary_map[1]) % 32 > 15
        
        # Separate the supplementary segments
        supps = [supp for supp in supp_list_sorted if supp != primary_map]
        
        samedir, invdir, diffchr = [], [], []

        # Classify supplementary segments into samedir, invdir, and diffchr categories
        for supp in supps:
            chrom, supp_flag, supp_maplen = supp[2], int(supp[1]) % 32 > 15, supp[5][1]
            if supp_maplen < 250:  # Ignore if the map length is less than 250
                continue
            if chrom != pri_chrom:
                diffchr.append(supp)
            elif supp_flag != pri_flag:
                invdir.append(supp)
            else:
                samedir.append(supp)

        # Process samedir: same direction alignments
        for supp_map in samedir:
            leftmap, rightmap = (primary_map, supp_map) if supp_map[3] > primary_map[3] else (supp_map, primary_map)
            sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
            sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
            maq1, maq2 = leftmap[6], rightmap[6]
            maq = (int(maq1) + int(maq2)) // 2
            
            # INS (Insertion) Detection
            if abs(rightmap[3] - leftmap[4]) <= 300:
                overlapmap = rightmap[3] - leftmap[4]
                ins_size = sh3 - len1 - sh1 - overlapmap
                if min_size <= ins_size <= max_size:
                    sv_start = min(rightmap[3], leftmap[4])
                    sv_end = sv_start + 1
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_INSLEN={ins_size}'
                    seq = readseq[sh1 + len1: sh3 - overlapmap]
                    add_svcall(pri_chrom, readname, sv_start, sv_end, ins_size, maq, svid, "INS", "*")
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

def segmentsv4lr(supp_list,min_size, max_size, chromosome_list, minimaq):
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
            if maq >= minimaq:
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
            if maq < minimaq:
                continue
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
        #if maq < min_maq:
        #    continue
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
            if 2<= len(supp_dict[readname])<=20:
                supp_svsignal = segmentsv4lr(supp_dict[readname],minLen, maxLen, chromosome_list,min_maq)
    return svInDels, covinfo , supp_svsignal


def svInDel4asm(line, minLen, min_maq, msv, maxLen, chromosome_list):
    if line.flag == 4 or line.mapping_quality < min_maq:
        return None  
    else:
        clipinfo = [0,0,0]
        supp_svsignal = []
        if line.cigar[0][0] in [4,5]:  ## left most cigar
            clipinfo[0] = line.cigar[0][1]
        if line.cigar[-1][0] in [4,5]:  ## right most cigar
            clipinfo[2] = line.cigar[-1][1]
        clipinfo[1] = line.query_alignment_length
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
        
        if msv == "no":
            return results, f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n', []  
        
        elif msv == "yes":
            supp_dict = {}
            readname, refend = query_chr,target_end
            if line.has_tag("SA"):
                primary = [readname, line.flag, line.reference_name,target_start,refend,clipinfo, maq, line.query_sequence ]
                if readname in supp_dict.keys():
                    supp_dict[readname] += [primary]
                else:
                    supp_dict[readname] = [primary]
                #print(primary[0])
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
                    suppinfo = [readname, sflag, chrom, start, end, clipinfo, maq, "*"]
                    #print(suppinfo)
                    supp_dict[readname] += [suppinfo]
            if supp_dict:
                if 2<= len(supp_dict[readname]):
                    supp_svsignal = segmentsv4asm(supp_dict[readname],minLen, maxLen, chromosome_list)
            return results,  f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n', supp_svsignal



def segment_segment_sv(lines, min_size, max_size, chromosome_list):
    """
    Collecting the supplementary alignment to identify the INS, DEL, DUP, INV, TRA signals
    LongRead align to reference that has "SA" will be recorded as dict as follow:
    { 'm64144': [['m64144',2048,'Db-Chr4',8949283,8949610,[0, 327, 11157],'43'],
                 ['m64144',0,   'Db-Chr4',8949649,8957962,[3162, 8322, 0],60,'AAXXXCTAATT'],
                 ['m64144',2064,'Db-Chr3',12497386,12502117,[5171, 4731, 1582],'60']]
                 }
    """
    supp_dict = {}
    for line in lines:
        clipinfo = [0, 0, 0]
        if line.flag != 4:
            if line.cigar[0][0] in [4, 5]:  ## left most cigar
                clipinfo[0] = line.cigar[0][1]
            if line.cigar[-1][0] in [4, 5]:  ## right most cigar
                clipinfo[2] = line.cigar[-1][1]
            clipinfo[1] = line.query_alignment_length
            query_chr = line.query_name
            query_seq = line.query_sequence
            target_chr = line.reference_name
            target_start = line.reference_start
            target_end = line.reference_end
            maq = line.mapping_quality
            flag = line.flag
            supp = [query_chr, flag, target_chr, target_start, target_end, clipinfo, maq, query_seq]
            if query_chr in supp_dict.keys():
                supp_dict[query_chr] += [supp]
            else:
                supp_dict[query_chr] = [supp]

    svsignal_supp = []

    def add_svcall(chrom, readname, start, end, size, maq, svid, sv_type, sequence=None):
        if maq == 60:
            sv_record = f"{chrom}\t{readname}\t{start}\t{end}\t{size}\t{maq}\t{svid}\t{sv_type}"
            if sequence:
                sv_record += f"\t{sequence}"
            svsignal_supp.append(sv_record)
            
    print("******************** Calling segment segment SV *************************")

    for query_chr, query_supp_list in supp_dict.items():
        supp_list_sorted = sorted(query_supp_list, key=lambda x: (x[2], int(x[3])))

        for idx, primary_map in enumerate(supp_list_sorted):
            readname = primary_map[0]
            pri_chrom = primary_map[2]
            pri_start = int(primary_map[3])
            readseq = primary_map[7] if len(primary_map) == 8 else "*"
            #print(f"********************************* flag: {primary_map[1]}*************************************")
            try:
                pri_flag = int(primary_map[1]) % 32 > 15
            except ValueError:
                print(f"Invalid flag value: {primary_map[1]}. Skipping this entry.")
                continue

            # Separate the supplementary segments
            supps = [supp for supp in supp_list_sorted if supp != primary_map]

            samedir, invdir, diffchr = [], [], []

            # Classify supplementary segments into samedir, invdir, and diffchr categories
            for supp in supps:
                #print(f'{primary_map[0:5]}\t{supp[0:5]}')
                try:
                    chrom, supp_flag, supp_maplen = supp[2], int(supp[1]) % 32 > 15, supp[5][1]
                except ValueError:
                    print(f"Invalid flag value: {supp[1]}. Skipping this supp entry.")
                    continue
                if supp_maplen < 250:  # Ignore if the map length is less than 250
                    continue
                if chrom != pri_chrom:
                    diffchr.append(supp)
                elif supp_flag != pri_flag:
                    invdir.append(supp)
                else:
                    samedir.append(supp)

            # Process samedir: same direction alignments
            for supp_map in samedir:
                leftmap, rightmap = (primary_map, supp_map) if supp_map[3] > primary_map[3] else (supp_map, primary_map)
                sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
                sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
                maq1, maq2 = leftmap[6], rightmap[6]
                maq = (int(maq1) + int(maq2)) // 2

                ## INS ##
                if abs(rightmap[3] - leftmap[4]) <= 300:
                    overlapmap = rightmap[3] - leftmap[4]
                    ins_size = sh3 - len1 - sh1 - overlapmap
                    if min_size <= ins_size <= max_size:
                        sv_start = min(rightmap[3], leftmap[4])
                        sv_end = sv_start + 1
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_INSLEN={ins_size}'
                        if readseq:
                            seq = readseq[sh1 + len1: sh3 - overlapmap]
                        else:
                            seq = "*"
                        add_svcall(pri_chrom, readname, sv_start, sv_end, ins_size, maq, svid, "INS", "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                ## DEL ##
                overlapmap = sh1 + len1 - sh3
                if -200 < overlapmap < 1500:  ## check or search for better
                    del_size = rightmap[3] - leftmap[4] + overlapmap
                    if min_size <= del_size <= max_size:
                        sv_start = max(0, leftmap[4] - max(0, overlapmap))
                        sv_end = sv_start + del_size - 1
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DELLEN={del_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, del_size, maq, svid, "DEL", '*')
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                ## DUP ##
                lap1 = sh1 + len1 - sh3
                if -200 < lap1 < 500 and (leftmap[4] - rightmap[3]) >= max(50, lap1):
                    dup_size = leftmap[4] - rightmap[3] - max(lap1, 0)
                    if min_size <= dup_size <= max_size:
                        sv_start, sv_end = rightmap[3], rightmap[3] + dup_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, dup_size, maq, svid, "DUP", "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                lap2 = sh3 + len2 - sh1
                if -200 < lap2 < 500 and (rightmap[4] - leftmap[3]) >= max(1000, lap2):
                    dup_size = rightmap[4] - leftmap[3] - lap2
                    if min_size <= dup_size <= max_size:
                        sv_start, sv_end = leftmap[3], leftmap[3] + dup_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, dup_size, maq, svid, 'DUP', "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
            ## INV ##
            for supp in invdir:
                if (supp[3] > primary_map[3] and (supp[4] - primary_map[4]) > -200) or \
                        (supp[3] < primary_map[3] and (primary_map[4] - supp[4]) > -200):
                    leftmap = primary_map if supp[3] > primary_map[3] else supp
                    rightmap = supp if supp[3] > primary_map[3] else primary_map

                    sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
                    sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
                    maq1, maq2 = leftmap[6], rightmap[6]
                    readname = primary_map[0]
                    maq = (int(maq1) + int(maq2)) // 2

                    lap1 = sh3 + len2 - sh2
                    if -200 < lap1 < 500 and (rightmap[4] - leftmap[4]) > max(100, lap1):
                        inv_size = rightmap[4] - leftmap[4] - lap1
                        if min_size <= inv_size <= max_size:
                            sv_start, sv_end = leftmap[4], leftmap[4] + inv_size
                            svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                            add_svcall(pri_chrom, readname, sv_start,
                                       sv_end, inv_size, maq, svid, 'INV', "*")
                            print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')

                    lap2 = sh4 + len2 - sh1
                    if -200 < lap2 < 500 and (rightmap[3] - leftmap[3]) >= max(100, lap2):
                        inv_size = rightmap[3] - leftmap[3] - lap2
                        if min_size <= inv_size <= max_size:
                            sv_start, sv_end = leftmap[3], leftmap[3] + inv_size
                            svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                            add_svcall(pri_chrom, readname, sv_start,
                                       sv_end, inv_size, maq, svid, "INV", "*")
                            print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
            ## TRA ##
            ### win size 500 may fit into the TRA ###
            for supp in diffchr:
                maq1, maq2 = supp[6], primary_map[6]
                readname = primary_map[0]
                maq = (int(maq1) + int(maq2)) // 2
                sh1, len1, sh2 = primary_map[5][0], primary_map[5][1], primary_map[5][2]
                sh3, len2, sh4 = supp[5][0], supp[5][1], supp[5][2]
                breakpoint1, breakpoint2 = '', ''
                if abs(sh1 - sh3 - len2) <= 500 or abs(sh1 - len2 - sh4) <= 500:
                    chrom1 = primary_map[2]
                    breakpoint1 = primary_map[3]
                elif abs(sh2 - sh3 - len2) <= 500 or abs(sh2 - len2 - sh4) <= 500:
                    chrom1 = primary_map[2]
                    breakpoint1 = primary_map[4]
                if abs(sh3 - sh1 - len1) <= 500 or abs(sh3 - len1 - sh2) <= 500:
                    chrom2 = supp[2]
                    breakpoint2 = supp[3]
                elif abs(sh4 - sh1 - len1) <= 500 or abs(sh4 - len1 - sh2) <= 500:
                    chrom2 = supp[2]
                    breakpoint2 = supp[4]
                    ## Chr4	ERR3415829.505585   15689499    9991420	0	60	Chr4:15689499_Chr1:9991420_TRA TRA
                if breakpoint1 != '' and breakpoint2 != '' and maq == 60 and (
                        chromosome_list.index(chrom1) < chromosome_list.index(chrom2)):
                    svid = chrom2 + ':' + str(breakpoint2) + '_' + chrom1 + ':' + str(breakpoint1)
                    svsignal_supp += [
                        chrom1 + '\t' + primary_map[0] + '\t' + str(breakpoint1) + '\t' + str(breakpoint2) + '\t0\t' + str(
                            maq) + "\t" + svid + '\tTRA' + "\t*"]
                elif breakpoint1 != '' and breakpoint2 != '' and maq == 60 and (
                        chromosome_list.index(chrom1) > chromosome_list.index(chrom2)):
                    svid = chrom1 + ':' + str(breakpoint1) + '_' + chrom2 + ':' + str(breakpoint2)
                    svsignal_supp += [
                        chrom2 + '\t' + primary_map[0] + '\t' + str(breakpoint2) + '\t' + str(breakpoint1) + '\t0\t' + str(
                            maq) + "\t" + svid + '\tTRA' + "\t*"]

    return svsignal_supp

```

### ./PSV_Signal/0.KLookCluster_SumAlleles4LocalDepthAdaptive.py
```python
import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil,floor

def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cost time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        try:
            num_maps = covinfo.count(str(reference), start, end)
        except AttributeError:
            print(f"Error: covinfo does not have a 'count' method.")
            num_maps = 0
    return num_maps

def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)

def candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate=0.1,add=1):
    """
    Solve clusters dataframe
    Get each cluster's sv breakpoints, SV length
    For low depth clus, the cluster must be vary strict and intense, else it will generate a lot of false positive sv. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts() ## default reversed count
    min_clusters = min(num_hap, len(cluster_counts))
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.1:
        print("proportion is too low, return []")
        return []
    else: 
        print(f"top1 cluster percent is {proportion}")
    
    top_clusters = cluster_counts.head(min_clusters).index
    clusdf = clusdf[clusdf[cluster_col].isin(top_clusters)]
    Total_signal = clusdf.shape[0]
    clusters_to_process =  clusdf[cluster_col].unique()

    msv = []
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = mode_or_median(clu_df['Target_end'])
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        if sv_eye < 2:
            continue
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 150), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 150)
        
        SV_rate = round(sv_eye / min(start_local_map+1,end_local_map+1), 2)
        if sv_eye >= min([start_local_map,end_local_map])*support_rate + add:
            print(f"{svid} sv signal propertion support")
            clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
            msv.append(svid)
        if sv_len > 2000 and svtype=="INS" and sv_eye >= 2 and sv_eye >= min([start_local_map,end_local_map])*support_rate*0.5: ## to capture big INS
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} big INS,local_map is {min([start_local_map,end_local_map])} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        if sv_eye >= nreads -1: ## filter 
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid},local coverage is {min([start_local_map,end_local_map])} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
            else:
                print(f"{svid}: number sv reads lower than required")
    if len(set(msv)) >= 2:
        outmsv = "\t".join(list(set(msv)))
        print(f'msv:\t{outmsv}')
    return clus

def klook_clusters(clusdf, max_diff_func, len_fold=0.8):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    num_signals = len(clusdf)
    if num_signals == 1:
        clusdf['shift_cluster'] = -1
        return clusdf
    
    # Precompute the 'Target_start' and 'Target_end' columns for quicker access
    start_values = clusdf['Target_start'].values
    end_values = clusdf['Target_end'].values
    svlen_values = clusdf['SVlen'].values
    
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    
    for i in range(1, len(clusdf)):
        current_svlen = svlen_values[i]
        current_start = start_values[i]
        current_end = end_values[i]
        found_cluster = False
        # Adjust the previous cluster range dynamically
        start_index = max(0, i - ceil(num_signals * 0.8))  # previous 80% record

        # Vectorized approach: avoid nested for loops
        for j in range(i - 1, start_index - 1, -1):
            old_len = svlen_values[j]
            relate_size = min(current_svlen / old_len, old_len/current_svlen)
            max_diff = max_diff_func(old_len)
            len_condition = relate_size > len_fold
            pos_condition = (abs(current_start - start_values[j]) <= max_diff or
                             abs(current_end - end_values[j]) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf

def max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 50 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 100 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 450 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 1000

def max_diff_func4INS(old_len):
    if old_len <= 100:
        return 50 + (old_len-50) * 0.5
    elif 100 < old_len <= 500:
        return 100 + (old_len-100) * 0.5
    elif 500 < old_len <= 5000:
        return 300 + (old_len-500) * 0.1
    elif old_len > 5000:
        return 1000

def onedepth_all_clus(all_signal,svtype, opened_bam, nrate=0.25):
    """
    Process the contig or genome level SV signal.
    Strict condition for merge, single chromosome mode.
    """
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    if all_signal.shape[0] > 1:
        all_signal = merge_and_sort_cr(all_signal)
        print(f"******************signal after merging ****************\n {all_signal}")
        if all_signal.shape[0] > 1:
            all_clus = klook_clusters(all_signal, max_diff_func, 0.8)
        else:
            all_signal['shift_cluster'] = -1
            all_clus = all_signal
    else:
        all_signal['shift_cluster'] = -1
        all_clus = all_signal
    print(all_clus[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    sv_chrom = all_clus["#Target_name"].unique()[0]
    svtype = all_clus["SVType"].unique()[0]
    cluster_col = all_clus.columns[-1]
    clus = []
    msv = []
    for clu in all_clus[cluster_col].unique():
        clu_df = all_clus[all_clus[cluster_col] == clu] 
        print(clu_df[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])   
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len =   mode_or_median(clu_df['SVlen'])
        sv_end =   mode_or_median(clu_df['Target_end'])
        maq =      mode_or_median(clu_df['maq'])
        readsname = clu_df['Query_name'].tolist()
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 250), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 250)
        depth = max([start_local_map,end_local_map])
        if depth > 0:
            SV_rate = round(sv_eye / depth, 2)
        else:
            SV_rate = 1
        if SV_rate < nrate: ## 0.25 to meet Tetraploid heterozygous
            continue
        print(f"************************cluster done***************************\n{[sv_chrom, sv_start, sv_end, sv_len, svid, svtype, '*', sv_eye, SV_rate, maq, readsname]}")
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
        msv.append(svid)
    if len(msv) >=2:
        outmsv = "\t".join(msv)
        print(f'msv:\t{outmsv}')
    return clus

def lowdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL    

    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def highdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])
    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates( subset=["Query_name", "Target_start", "Target_end"], keep="first",inplace=True) ## we dont keep the same reads multiple segment
    print(clusdf[["#Target_name","Target_start","Target_end","SVlen","SVType"]])
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf

def merge_and_sort_cr(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    Optimized version for better performance with large datasets.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates(subset=["Query_name", "Target_start", "Target_end"], keep="first", inplace=True)
    clusdf.sort_values(['Query_name', 'Target_start'], inplace=True)
    clusdf['same_query'] = clusdf['Query_name'] == clusdf['Query_name'].shift(1)
    clusdf['overlap'] = (clusdf['Target_start'] < clusdf['Target_end'].shift(1)) & clusdf['same_query']
    clusdf['group_key'] = (clusdf['Query_name'] != clusdf['Query_name'].shift(1)) | clusdf['overlap']
    
    clusdf['group_key'] = clusdf['group_key'].cumsum()
    merged_sv = clusdf.groupby('group_key').agg({
        'Query_name': 'first',
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index(drop=True)
    merged_sv = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    merged_sv.sort_values(by=['Target_start'], inplace=True)
    merged_sv.reset_index(drop=True, inplace=True)
    merged_sv['shift_cluster'] = -1
    return merged_sv


def windows_slide4asm(dfs, svtype, window_size=500):
    """
    fix window gap for assemblies
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        j = i  
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom:
            last_in_window_start = dfs.at[j-1, 'Target_start'] if j > i else current_start
            if svtype == 'DEL':
                dynamic_window_end = last_in_window_start + window_size + 250
            else:
                dynamic_window_end = last_in_window_start + window_size
    
            if dfs.at[j, 'Target_start'] < dynamic_window_end:
                j += 1  
            else:
                break  
        window_df = dfs[i:j].copy()
        if 1 <= len(window_df['Query_name'].unique()):
            window_df.index = range(len(window_df))  
            windows[current_start] = window_df  
        i = j  
    return windows

def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if svtype == 'DEL':
            window_end = current_start + window_size + 250 ## avoid small fragment deletions
        else:
            window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 100:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows

def windows_slide4tra(df, shift=2000000):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df = df.copy()
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs = dfs.sort_values(by=['Target_start1', 'Target_start2'])
        dfs.index = range(len(dfs))
        i = 0
        while i < len(dfs):
            current_chrom = str(dfs.at[i, '#Target_name1'])
            current_start = dfs.at[i, 'Target_start1']
            window_size = shift
            window_end = current_start + window_size
            j = i
            while j < len(dfs) and dfs.at[j, 'Target_start1'] < window_end:
                j += 1
            window_df = dfs[i:j].copy()
            window_df = window_df.reset_index(drop=True)
            key = f'{chrom1}_{current_start}_{chrom2}'
            windows[key] = window_df
            i = j
    return windows

def klook_clu_tra(win):
    """
    tra clus, the clusdf should be chromosome pair mode,
    each clusdf only allow two chromosome.
    """
    if len(win) == 1:
        win['shift_cluster'] = -1
        return win
    cluster_id = 0
    win.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(win)):
        current_start1 = win.loc[i, 'Target_start1']
        current_start2 = win.loc[i, 'Target_start2']
        found_cluster = False
        start_index = max(0, i - 40)
        for j in range(i - 1, start_index - 1, -1):
            max_diff = 1000
            ## already ensure chrom1 != chrom2
            pos_condition = (abs(current_start1 - win.loc[j, 'Target_start1']) <= max_diff and
                             abs(current_start2 - win.loc[j, 'Target_start2']) <= max_diff)
            if  pos_condition:
                win.loc[i, 'shift_cluster'] = win.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            win.loc[i, 'shift_cluster'] = cluster_id
    return win

def candidate_tra(window, opened_bam, dtype):
    """
    Process translocation data and find candidate translocations.
    """
    win_clus = klook_clu_tra(window)
    tras = []
    for clu in win_clus['shift_cluster'].unique():
        win = win_clus[win_clus['shift_cluster']==clu]
        chr1 = win['#Target_name1'].iloc[0]
        chr1_start = mode_or_median(win['Target_start1'])
        sv_len = 0
        chr2 = win['#Target_name2'].iloc[0]
        chr2_start = mode_or_median(win['Target_start2'])
        maq = mode_or_median(win['maq'])
        readsname = set(win['Query_name'])
        svid = f'{chr2}:{chr2_start}_{chr1}:{chr1_start}'
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, chr1, max(0, chr1_start - 250), max(chr1_start - 150, 0))
        end_local_map = local_cov(opened_bam, chr2, chr2_start + 150, chr2_start + 250)
        local_depth = np.mean([start_local_map, end_local_map])
        if local_depth >0:
            SV_rate = round(sv_eye / local_depth, 2)
        else:
            SV_rate = 1
        if dtype in ['ont','hifi', 'pb']:
            ## to one depth ##
            if sv_eye >= (local_depth * 0.1 + 0.7 ):
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
        elif dtype in ['cr', 'sr']:
            if sv_eye > local_depth * 0.25:
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tras

def load_and_process_sv_data(args):
    from math import floor
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    except pd.errors.EmptyDataError:
        print(f"Warning: File {args.raw_signal} is empty.")
        sv_indel_data = pd.DataFrame()
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else ceil(float(depth_stat.iloc[0, 3])+0.3)
        if args.nreads:
            nreads_fil = args.nreads
        else:
            nreads_fil = floor(depth / 10)
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth, nreads_fil
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data = sv_indel_data.copy()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_indel_data['Target_end'] = sv_indel_data['Target_end'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"***************** empty {supp_align_file} file ******************")
        else:
            if not msv.empty:
                print(f'*************** {args.raw_signal}.suppAlign has {msv.shape[0]} rows ********************')
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                #msv = msv.drop_duplicates()
                print(f'*************** {args.raw_signal}.suppAlign after duplicates drop has {msv.shape[0]} rows ********************')
                #msv = msv.reset_index(drop=True)
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['Target_end'] = sv_data[svtype]['Target_end'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    print(f'******************** all chromosomes list {chroms} ****************************')
    return sv_data, chroms, depth, nreads_fil

def process_svtype(args, sv_data, chroms, svtype, depth, nreads, minLen):
    """
        cov info parse by covfile or bam file
    """
    print(f"start klook for {args.raw_signal}  SV type: {svtype}") 
    if args.dtype in ['cr', 'sr']:
        try:
            covinfo = pd.read_csv(args.covfile, sep="\t", index_col=None, header=None)
            covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
            covinfo['target_chr'] = covinfo['target_chr'].astype(str)
            bam_path = covinfo
        except FileNotFoundError:
            print(f"Error: Coverage file {args.covfile} not found.")
            return [], [], [], [], []
    else:
        bam_path = args.bam

    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    if args.dtype in ['pb', 'ont', 'hifi']:
        print(f'data type is {args.dtype}')
        try:
            with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
                for chrom in chroms:
                    chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                    def process_sv_type(svtype):
                        if svtype == "TRA":
                            log = open("log_tra", 'w')
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            windows = windows_slide4tra(chrom_data[svtype], 2000000)
                            print(f'The TRA windows by 2M: \n{windows}', file=log)
                            tra_list = []
                            for win in windows.values():
                                win_clus = klook_clu_tra(win)
                                print(win_clus.iloc[:,[0,2,3,4,5,-1]], file=log)
                                tra = candidate_tra(win_clus, opened_bam, args.dtype)
                                print(tra, file=log)

                                if tra:
                                    tra_list.extend(tra)
                            log.close()
                            return tra_list
                        else:
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                            sv_windows = windows_slide(sv_dfs, depth, svtype, nreads,args.window_size)
                            candidate_svs = []
                            if args.dtype == 'hifi':
                                hdepth = 5
                            else:
                                hdepth = 10
                            for sv_window in sv_windows.values():
                                if len(sv_window['Query_name']) > hdepth: 
                                    candisv = highdepth_clu(sv_window,args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
                                else:
                                    candisv = lowdepth_clu(sv_window, args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
                                candidate_svs.extend(candisv)
                            return candidate_svs

                    if svtype == "DEL":
                        del_clus.extend(process_sv_type("DEL"))
                    elif svtype == "INS":
                        ins_clus.extend(process_sv_type("INS"))
                    elif svtype == "INV":
                        inv_clus.extend(process_sv_type("INV"))
                    elif svtype == "DUP":
                        dup_clus.extend(process_sv_type("DUP"))
                    elif svtype == "TRA":
                        tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f"Error processing BAM file: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    elif args.dtype in ['sr', 'cr']:
        print(f'data type is {args.dtype}')
        try:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                def process_sv_type(svtype):
                    if svtype != "TRA":
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                        if sv_dfs.empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        sv_dfs = sv_dfs.sort_values(by=['Target_start', 'SVlen'])
                        sv_dfs.index = range(len(sv_dfs))
                        print(sv_dfs.head(10))
                        print("**************************** Calling one depth all clustering ***********************")
                        candidate_svs = []
                        win_dfs = windows_slide4asm(sv_dfs, svtype, args.window_size)
                        for win_df in win_dfs.values():
                            print(f'**************************window signals*********************\n{win_df}')
                            candidate_sv = onedepth_all_clus(win_df, svtype, bam_path,max(args.rate_depth,0.25))
                            candidate_svs.extend(candidate_sv)
                        return candidate_svs
                    else:
                        if chrom_data[svtype].empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        windows = windows_slide4tra(chrom_data[svtype], 1000000)
                        tra_list = []
                        for value in windows.values():
                            win = value
                            tra = candidate_tra(win, bam_path,args.dtype)
                            if tra:
                                tra_list.extend(tra)
                        print(f'****************************** the tra condidate ***************************')
                        print(tra_list)
                        return tra_list
                if svtype == "DEL":
                    del_clus.extend(process_sv_type("DEL"))
                elif svtype == "INS":
                    ins_clus.extend(process_sv_type("INS"))
                elif svtype == "INV":
                    inv_clus.extend(process_sv_type("INV"))
                elif svtype == "DUP":
                    dup_clus.extend(process_sv_type("DUP"))
                elif svtype == "TRA":
                    tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f'Error in processing {args.dtype}: {args.raw_signal}')
            print(f"Error processing {svtype} for chromosome {chrom}: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        print(f"Error: sequence data dtype error: {args.dtype} is not in [pb,hifi,ont,sr,cr], please check parameter -dtype")
        return [],[],[],[],[]

def candidateSV(args):
    sv_data, chroms, depth, nreads = load_and_process_sv_data(args)
    print(chroms, depth, nreads)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth, nreads, args.min) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            if result is not None:
                tra_clus += result[0]
                del_clus += result[1]
                ins_clus += result[2]
                inv_clus += result[3]
                dup_clus += result[4]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        return [], [], [], [], []


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True,
                    help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=800, type=int,
                    help="the distance shift of breakpoint to cluster the TRA/big_INV/big_DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=10000000, help="the max SV length")
    IN.add_argument("-m", dest="min", type=int, default=45, help="the minimum SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    IN.add_argument("--nreads", dest="nreads", type=int, help="the minimum numbers of reads to support SV, if not provided, we use average_depth / 10 as threshold" )
    IN.add_argument("--rate_depth", dest="rate_depth", type=float, default=0.1, help="the sv supports of local depth ratio to support sv, 0.1 means the percent of local reads shoule support sv")
    IN.add_argument("--window", dest="window_size", type=int, default=500, help="the window size of signal to parse in klook cluster, 500bp suggested")
    IN.add_argument("--num_hap", dest="num_hap", type=int, default=2, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

```

### ./PSV_Signal/Sub_readfa2Dict.py
```python
import gzip
from collections import deque
from gzip import BadGzipFile
import re
def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

```

### ./PSV_Signal/1.PSV_signal_cluster.py
```python
import pandas as pd
import argparse
from time import time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
def mode_or_median(series):
    """
    Returns the most common value from a series. If there are ties, it will return the first mode.
    """
    return series.mode().iloc[0] if not series.mode().empty else math.floor(series.median())

def most_common(series):
    """
    Returns the most common value from a series. If there are ties, it will return the first mode.
    """
    return series.mode().iloc[0] if not series.mode().empty else series.iloc[0]

def svs_clu(chrsvdf, svtype, chrom, max_diff=50):
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start'])
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1,len(df)):
        if df.loc[i, '#Target_name'] == df.loc[i-1, '#Target_name']:
            old_len = df.loc[i-1, 'SVlen']
            now_len = df.loc[i, 'SVlen']
            relate_size = old_len / now_len
            if abs(df.loc[i, 'Target_start'] - df.loc[i-1, 'Target_start']) <= max_diff and abs(df.loc[i, 'Target_end'] - df.loc[i-1, 'Target_end']) <= max_diff and (0.8 < relate_size < 1.2):
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id
    print(f'The {chrom} {svtype} data {df.shape} after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        maq = int(cs['maq'].mean())
        sv_rate = most_common(cs['sv_rate'])
        Target_name = most_common(cs['#Target_name'])
        Target_start = most_common(cs['Target_start'])
        Target_end = most_common(cs['Target_end'])
        cluster_size = most_common(cs['cluster_size'])
        SVlen =  most_common(cs['SVlen'])
        SVID =  most_common(cs['SVID'])
        SVType =   most_common(cs['SVType'])
        seq = most_common(cs['seq'])
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start,
            'Target_end': Target_end,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq':maq,
            'cluster_size_prevalent': cluster_size,
            'sv_rate_prevalent': sv_rate
                })
    return clu

def tra_clu(tradf, chrom, max_diff=100):
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','maq','cluster_size','sv_rate']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(int)
    df["Target_start2"] = df["Target_start2"].astype(int)
    df.sort_values(by=['#Target_name1', '#Target_name2','Target_start1', 'Target_start2'],inplace=True)
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1, len(df)):
        if df.loc[i, '#Target_name1'] == df.loc[i-1, '#Target_name1'] and df.loc[i, '#Target_name2'] == df.loc[i-1, '#Target_name2']:
            if abs(df.loc[i, 'Target_start1'] - df.loc[i-1, 'Target_start1']) <= max_diff and abs(df.loc[i, 'Target_start2'] - df.loc[i-1, 'Target_start2']) <= max_diff:
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id

    print(f'The {chrom} TRA signal data {tradf.shape} after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        Target_name = most_common(cs['#Target_name1'])
        Target_start1 = most_common(cs['Target_start1'])
        Target_start2 = most_common(cs['Target_start2'])
        cluster_size =  most_common(cs['cluster_size'])
        sv_rate = most_common(cs['sv_rate'])
        SVlen =  0
        maq = cs['maq'].mean()
        SVID =  most_common(cs['SVID'])
        SVType =  "TRA"
        seq = "*"
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq':int(maq),
            'cluster_size_prevalent': cluster_size,
            'sv_rate_prevalent':sv_rate
            })
    return clu

def file_capture(dir, suffix):
    import os
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            captures.append(os.path.join(dir, file))
    return captures
def read_file(file_name):
    import os
    # Check if the file exists and is not empty
    if os.path.exists(file_name) and os.path.getsize(file_name) > 0:
        try:
            svi = pd.read_csv(file_name, header=0, index_col=None, sep="\t")
            return svi
        except pd.errors.EmptyDataError:
            print(f"The file {file_name} is empty or malformed.")
            return None
    else:
        print(f"The file {file_name} does not exist or is empty.")
        return None

def candidateSV(args):
    file_lists = file_capture(args.sv_dir, "_Clustered_Record.txt")
    file_lists.sort()
    pop_num = len(file_lists)
    print(f'************************** The total population number found is {pop_num}\n{file_lists} ***************************')
    sv = pd.read_csv(file_lists[0],header=0,index_col=None,sep="\t")
    for file_name in file_lists[1:]:
        svi = read_file(file_name)
        sv = pd.concat([sv,svi],axis=0)
    ori = sv.shape
    sv_types = sv["SVType"].unique()
    sv.sort_values(by=["#Target_name","Target_start","SVlen","SVID"], inplace=True)
    def process_chromosome(chrom):
        inschr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INS")]
        delchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DEL")]
        trachr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="TRA")]
        invchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INV")]
        dupchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DUP")]
        if  delchr.empty:
            dels = []
        else:
            dels = svs_clu(delchr,'DEL', chrom, args.shift)

        if inschr.empty:
            ins = []
        else:
            ins = svs_clu(inschr,'INS', chrom, args.shift)

        if invchr.empty:
            inv = []
        else:
            inv = svs_clu(invchr, 'INV', chrom, args.shift)

        if dupchr.empty:
            dup =[]
        else:
            dup = svs_clu(dupchr,'DUP', chrom, args.shift)

        if  trachr.empty:
            tra = []
        else:
            tra = tra_clu(trachr, chrom, args.shift*2)
        return tra + dels + ins + inv + dup
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_chromosome, sv['#Target_name'].unique())
    out = []
    for result in results:
        out += result
    sv_out = pd.DataFrame(out)
    print(sv_out.head())
    sv_out["SVlen"] = sv_out["SVlen"].astype(int)
    sv_out = sv_out.sort_values(by=['SVID'],inplace=False)
    sv_fil = sv_out.drop_duplicates(subset='SVID', keep='last', inplace=False) 
    sv_fil[sv_fil["SVlen"] <= args.max].to_csv(f"{args.sv_dir}/PopSV_Candidate_Record.txt",header=True,index=None,sep="\t")
    now = sv_fil.shape
    print(f'Original data shape: {ori}, after clustering: {now}\nOutput file is {args.sv_dir}/PopSV_Candidate_Record.txt')


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-d", dest="sv_dir", required=True, help="the PSVGT output directory")
    IN.add_argument("-s", dest="shift", default=30, type=int, help="the distance of shifting the breakpoints ")
    IN.add_argument("-M", dest="max", default=6868886, type=int, help="the max SV length ")
    args = parser.parse_args()
    start_t = time()
    candidateSV(args)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")


```

### ./PSV_Signal/1.ACCSV_Signal_Cluster.py
```python
import pandas as pd
import argparse
from time import time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
def svs_clu(chrsvdf, svtype,max_diff=30):
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start','SVlen'])
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1,len(df)):
        if df.loc[i, '#Target_name'] == df.loc[i-1, '#Target_name']:
            old_len = df.loc[i-1, 'SVlen']
            now_len = df.loc[i, 'SVlen']
            relate_size = old_len / now_len
            if abs(df.loc[i, 'Target_start'] - df.loc[i-1, 'Target_start']) <= max_diff and abs(df.loc[i, 'Target_end'] - df.loc[i-1, 'Target_end']) <= max_diff and (0.8 < relate_size < 1.2):
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id
    print(f'The {svtype} data shape: {df.shape}\nThe {svtype} data shape after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        max_clu_row = cs[cs['cluster_size'] == cs['cluster_size'].max()].iloc[0]
        Target_name = max_clu_row['#Target_name']
        Target_start1 = max_clu_row['Target_start']
        Target_start2 = max_clu_row['Target_end']
        cluster_size = cs['cluster_size'].sum()
        SVlen =  cs['SVlen'].max()
        maq = int(cs['maq'].mean())
        SVID = max_clu_row['SVID']
        SVType = max_clu_row['SVType']
        seq = "*"
        sv_rate = cs['sv_rate'].sum()
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq':maq,
            'cluster_size': cluster_size,
            'sv_rate': sv_rate
         })
    return clu


def tra_clu(tradf,  max_diff=800):
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','cluster_size','sv_rate', 'maq', 'readsID']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(int)
    df["Target_start2"] = df["Target_start2"].astype(int)
    df.sort_values(by=['#Target_name1', '#Target_name2','Target_start1', 'Target_start2'],inplace=True)
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1, len(df)):
        if df.loc[i, '#Target_name1'] == df.loc[i-1, '#Target_name1'] and df.loc[i, '#Target_name2'] == df.loc[i-1, '#Target_name2']:
            if abs(df.loc[i, 'Target_start1'] - df.loc[i-1, 'Target_start1']) <= max_diff and abs(df.loc[i, 'Target_start2'] - df.loc[i-1, 'Target_start2']) <= max_diff:
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id

    print(f'The TRA signal shape: {df.shape}\nThe TRA signal data shape after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        max_clu_row = cs[cs['cluster_size'] == cs['cluster_size'].max()].iloc[0]
        Target_name = max_clu_row['#Target_name1']
        Target_start1 = max_clu_row['Target_start1']
        Target_start2 = max_clu_row['Target_start2']
        cluster_size = cs['cluster_size'].sum()
        maq = int(cs['maq'].mean())
        SVlen = 0
        SVID = max_clu_row['SVID']
        SVType = "TRA"
        seq = "*"
        sv_rate = cs['sv_rate'].sum()
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq': maq,
            'cluster_size': cluster_size,
            'sv_rate': sv_rate
         })
    return clu


def drop_ins_from_dup(df):
    out=open("log", 'w')
    ins = df[df['SVType'].isin(['INS', 'DUP'])]
    ins = ins.sort_values(by=['#Target_name', 'Target_start'])
    ins.index = range(len(ins))
    dup_rows = ins[ins['SVType'] == 'DUP']
    if dup_rows.empty:
        return df

    dropINS = []
    for index, dup_row in dup_rows.iterrows():
        start_index = max(0, index - 10)
        end_index = min(len(df), index + 10)
        nearby_ins_rows = ins[(ins.index >= start_index) & (ins.index < end_index)]
        print(nearby_ins_rows, file=out)
        for ins_index, ins_row in nearby_ins_rows.iterrows():
            if (dup_row['Target_start'] < ins_row['Target_start'] - 50 < dup_row['Target_end'] and
                    dup_row['#Target_name'] == ins_row['#Target_name']):
                ratio = min(ins_row['SVlen'] / dup_row['SVlen'], dup_row['SVlen']/ins_row['SVlen'])
                if ratio >= 0.8:
                    dropINS.append(ins_row['SVID'])
    filtered_df = df[~df['SVID'].isin(dropINS)]
    print(f'******************************* drop {len(dropINS)} INS SV since they are DUP **************************************')
    print(f'******************************* drop {dropINS} INS SV since they are DUP **************************************')
    out.close()
    return filtered_df

def read_file(file_name):
    import os
    if os.path.exists(file_name) and os.path.getsize(file_name) > 0:
        try:
            svi = pd.read_csv(file_name, header=0, index_col=None, sep="\t")
            return svi
        except pd.errors.EmptyDataError:
            print(f"The file {file_name} is empty or malformed.")
            return None
    else:
        print(f"The file {file_name} does not exist or is empty.")
        return None

def candidateSV(args):
    from os.path import basename, exists
    import pandas as pd
    chrs = pd.read_csv(args.fai,sep="\t",index_col=None,header=None)[0].tolist()
    file_lists = []
    for sv in ['INS', 'DEL', 'TRA', 'DUP','INV']:
        for chrom in chrs:
            file_name = f"{args.preffix}_{chrom}.record.txt_{sv}.signal"
            if exists(file_name):
                file_lists.append(file_name)
    print(file_lists)

    sv = read_file(file_lists[0])
    if args.nreads_nrate:
        sv = sv[(sv["cluster_size"] >= args.nreads) & (sv["sv_rate"] >= args.nrate)]
    else:
        if args.nrate and not args.nreads:
            sv = sv[sv["sv_rate"] >= args.nrate]
        elif args.nreads and not args.nrate:
            sv = sv[sv["cluster_size"] >= args.nreads]
        elif args.nreads and args.nrate:
            sv = sv[(sv["sv_rate"] >= args.nrate) | (sv["cluster_size"] >= args.nreads)]

    for file_name in file_lists[1:]:
        svi = read_file(file_name)
        sv = pd.concat([sv,svi],axis=0)
    if args.nreads_nrate:
        sv = sv[(sv["cluster_size"] >= args.nreads) & (sv["sv_rate"] >= args.nrate)]
    else:
        if args.nrate and not args.nreads:
            sv = sv[sv["sv_rate"] >= args.nrate]
        elif args.nreads and not args.nrate:
            sv = sv[sv["cluster_size"] >= args.nreads]
        elif args.nreads and args.nrate:
            sv = sv[(sv["sv_rate"] >= args.nrate) | (sv["cluster_size"] >= args.nreads)]

    ori = sv.shape
    sv_types = sv["SVType"].unique()
    ori = sv.shape
    sv_types = sv["SVType"].unique()
    sv.sort_values(by=["#Target_name","Target_start","SVlen","SVID"], inplace=True)
    def process_chromosome(chrom):
        inschr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INS")]
        delchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DEL")]
        trachr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="TRA")]
        invchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INV")]
        dupchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DUP")]
        if  delchr.empty:
            dels = []
        else:
            dels = svs_clu(delchr,'DEL', args.shift)

        if inschr.empty:
            ins = []
        else:
            ins = svs_clu(inschr,'INS', args.shift)

        if invchr.empty:
            inv = []
        else:
            inv = svs_clu(invchr, 'INV', args.shift)

        if dupchr.empty:
            dup =[]
        else:
            dup = svs_clu(dupchr,'DUP', args.shift)

        if  trachr.empty:
            tra = []
        else:
            tra = tra_clu(trachr, args.shift*2)
        return tra + dels + ins + inv + dup
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_chromosome, sv['#Target_name'].unique())
    out = []
    for result in results:
        out += result
    svs_df = pd.DataFrame(out)
    sv_out = drop_ins_from_dup(svs_df)
    sv_out.loc[:,"SVlen"] = sv_out["SVlen"].astype(int)
    sv_out[(sv_out["SVlen"] <= args.max) & (sv_out['maq']>=args.minimaq)].to_csv(f"{args.preffix}_Clustered_Record.txt",header=True,index=None,sep="\t")
    now = sv_out.shape
    print(f'Original data shape: {ori}, after clustering: {now}\nOutput file is {args.preffix}_Clustered_Record.txt')


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-preffix", dest="preffix", required=True, help="the sample preffix name of SV raw sginal")
    IN.add_argument("-fai", dest="fai",required=True, help="the reference faidx file")
    IN.add_argument("-s", dest="shift", default=100, type=int, help="the distance of shifting the breakpoints ")
    IN.add_argument("-M", dest="max", default=6868886, type=int, help="the max SV length ")
    IN.add_argument("--nreads", dest="nreads", type=int, help="filter signal by minimum support reads")
    IN.add_argument("--nrate", dest="nrate",   type=float, help="filter signal by minimum ratio of signal")
    IN.add_argument("--nn", dest="nreads_nrate",action="store_true", help="strict filter for signal use both meet threshold of nreads and nrate, if not --nn, keep SV meet nreads or nrate")
    IN.add_argument("--minimaq", dest="minimaq", type=int, default=50,help="minimum mapping quality of SV")
    args = parser.parse_args()
    start_t = time()
    candidateSV(args)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

```

### ./PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py
```python
import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil,floor

def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cost time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        try:
            num_maps = covinfo.count(str(reference), start, end)
        except AttributeError:
            print(f"Error: covinfo does not have a 'count' method.")
            num_maps = 0
    return num_maps

def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)


def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_mode = subset.mode()
        if not subset_mode.empty:
            return subset_mode.iloc[0]  
        else:
            return ceil(subset.mean())
    else:
        # if subset empty return original mode
        original_mode = series.mode()
        if not original_mode.empty:
            return original_mode.iloc[0]
        else:
            # no mode, return mean
            return ceil(series.mean())

def candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate=0.1,add=1):
    """
    Solve clusters dataframe
    Get each cluster's sv breakpoints, SV length
    For low depth clus, the cluster must be vary strict and intense, else it will generate a lot of false positive sv. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts() ## default reversed count
    min_clusters = min(num_hap, len(cluster_counts))
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.05:
        print("proportion is too low, return []")
        return []
    else: 
        print(f"top1 cluster percent is {proportion}")
    
    top_clusters = cluster_counts.head(min_clusters).index
    clusdf = clusdf[clusdf[cluster_col].isin(top_clusters)]
    clusters_to_process =  clusdf[cluster_col].unique()

    msv = []
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = clu_df['Target_end'].max()
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        
        if sv_eye < 2:
            continue
        if svtype in ["INS", "DUP"]:
            break_cov = local_cov(opened_bam, sv_chrom, sv_start-70, sv_start-50) + 1
        elif svtype in ["DEL", "INV"]:
            midpoint = int(sv_start + sv_len*0.5)
            break_cov = local_cov(opened_bam, sv_chrom, midpoint, midpoint+20) + 1
        else:
            print("!!!!!!!!!!!!!!!!!!!!! SV type error !!!!!!!!!!!!!!!!!!!!!!!!")

        SV_rate = round(sv_eye / break_cov, 2)
        if sv_eye >= break_cov*support_rate + add:
            print(f"{svid} sv signal propertion support")
            clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
            msv.append(svid)
        if sv_len > 3000 and svtype=="INS" and sv_eye >= 2 and sv_eye >= break_cov*support_rate*0.6: ## to capture big INS
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} big INS,breakpoint coverage is {break_cov} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        if sv_len > 500 and svtype=="DUP" and sv_eye >= 2 and sv_eye >= break_cov*support_rate*0.5: ## to capture more DUP
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} DUP,breakpoint coverage is {break_cov} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        
        if sv_eye >= nreads -1: ## filter 
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid},but local coverage is {break_cov}, breakdepth not support this sv")
    if len(set(msv)) >= 2:
        outmsv = "\t".join(list(set(msv)))
        print(f'msv:\t{outmsv}')
    return clus

def klook_clusters(clusdf, max_diff_func, len_fold=0.8):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    num_signals = len(clusdf)
    if num_signals == 1:
        clusdf['shift_cluster'] = -1
        return clusdf
    
    # Precompute the 'Target_start' and 'Target_end' columns for quicker access
    start_values = clusdf['Target_start'].values
    end_values = clusdf['Target_end'].values
    svlen_values = clusdf['SVlen'].values
    
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    
    for i in range(1, len(clusdf)):
        current_svlen = svlen_values[i]
        current_start = start_values[i]
        current_end = end_values[i]
        found_cluster = False
        # Adjust the previous cluster range dynamically
        start_index = max(0, i - ceil(num_signals * 0.8))  # previous 80% record

        # Vectorized approach: avoid nested for loops
        for j in range(i - 1, start_index - 1, -1):
            old_len = svlen_values[j]
            relate_size = min(current_svlen / old_len, old_len/current_svlen)
            max_diff = max_diff_func(old_len)
            len_condition = relate_size > len_fold
            pos_condition = (abs(current_start - start_values[j]) <= max_diff or
                             abs(current_end - end_values[j]) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf

def max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 50 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 100 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 450 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 1000

def max_diff_func4INS(old_len):
    if old_len <= 100:
        return 40 + (old_len-50) * 0.8
    elif 100 < old_len <= 500:
        return 80 + (old_len-100) * 0.8
    elif 500 < old_len <= 5000:
        return 400 + (old_len-500) * 0.1
    elif old_len > 5000:
        return 1000


def lowdepth_max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 25 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 50 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 100 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 800

def lowdepth_max_diff_func4INS(old_len):
    if old_len <= 100:
        return 25 + (old_len-50) * 0.4
    elif 100 < old_len <= 500:
        return 50 + (old_len-100) * 0.4
    elif 500 < old_len <= 5000:
        return 200 + (old_len-500) * 0.05
    elif old_len > 5000:
        return 800

def onedepth_all_clus(all_signal,svtype, opened_bam, msvfold, nrate=0.25):
    """
    Process the contig or genome level SV signal.
    Strict condition for merge, single chromosome mode.
    """
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    if all_signal.shape[0] > 1:
        all_signal = merge_and_sort_cr(all_signal)
        print(f"******************signal after merging ****************\n {all_signal}")
        if all_signal.shape[0] > 1:
            all_clus = klook_clusters(all_signal, max_diff_func, msvfold)
        else:
            all_signal['shift_cluster'] = -1
            all_clus = all_signal
    else:
        all_signal['shift_cluster'] = -1
        all_clus = all_signal
    print(all_clus[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    sv_chrom = all_clus["#Target_name"].unique()[0]
    svtype = all_clus["SVType"].unique()[0]
    cluster_col = all_clus.columns[-1]
    clus = []
    msv = []
    for clu in all_clus[cluster_col].unique():
        clu_df = all_clus[all_clus[cluster_col] == clu] 
        print(clu_df[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])   
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len =   mode_or_median(clu_df['SVlen'])
        sv_end =   mode_or_median(clu_df['Target_end'])
        maq =      mode_or_median(clu_df['maq'])
        readsname = clu_df['Query_name'].tolist()
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 250), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 250)
        depth = max([start_local_map,end_local_map])
        if depth > 0:
            SV_rate = round(sv_eye / depth, 2)
        else:
            SV_rate = 1
        if SV_rate < nrate: ## 0.25 to meet Tetraploid heterozygous
            continue
        print(f"************************cluster done***************************\n{[sv_chrom, sv_start, sv_end, sv_len, svid, svtype, '*', sv_eye, SV_rate, maq, readsname]}")
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
        msv.append(svid)
    if len(msv) >=2:
        outmsv = "\t".join(msv)
        print(f'msv:\t{outmsv}')
    return clus

def lowdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, msvfold, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= lowdepth_max_diff_func4INS
    else:
        max_diff_func= lowdepth_max_diff_func4DEL    

    clusdf = klook_clusters(clusdf, max_diff_func, msvfold)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def highdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, msvfold, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    clusdf = klook_clusters(clusdf, max_diff_func, msvfold)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])
    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def filter_overlaps(df):
    """
    to filter overlap SV within same queryname
    """
    filtered_groups = []
    for query, group in df.groupby('Query_name'):
        if len(group) <= 1:
            filtered_groups.append(group)
            continue
        sorted_group = group.sort_values(
            by=['Target_start', 'SVlen'],
            ascending=[True, False]
        )
        kept = []
        last_end = -1
        for _, row in sorted_group.iterrows():
            start = row['Target_start']
            end = row['Target_end']
            if start > last_end:  
                kept.append(row)
                last_end = end
            else:
                 if end > last_end:
                     kept[-1] = row
                     last_end = end
        filtered_groups.append(pd.DataFrame(kept))    
    return pd.concat(filtered_groups, ignore_index=True)

def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates( subset=["Query_name", "Target_start", "Target_end"], keep="first",inplace=True) ## we dont keep the same reads multiple segment
    print(clusdf[["#Target_name","Target_start","Target_end","SVlen","SVType"]])
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf


def merge_and_sort_cr(clusdf):
    """
    filter overlap SV on same query name;
    merge sv based on query name;
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates(
        subset=["Query_name", "Target_start", "Target_end"],
        keep="first",
        inplace=True
    )
    clusdf = filter_overlaps(clusdf)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'min',       
        'Target_end': 'max',    
        'SVlen': 'sum',              
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv[['#Target_name', 'Target_start', 'Target_end', 'SVlen', 'SVType', 'maq', 'seq', 'Query_name']]
    clusdf = clusdf.astype({
        'Target_start': np.int32,
        'Target_end': np.int32,
        'SVlen': np.int32,
        'maq': np.int16
    })
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf

def windows_slide4asm(dfs, svtype, window_size=500):
    """
    fix window gap for assemblies
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        j = i  
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom:
            last_in_window_start = dfs.at[j-1, 'Target_start'] if j > i else current_start
            if svtype == 'DEL':
                dynamic_window_end = last_in_window_start + window_size + 250
            else:
                dynamic_window_end = last_in_window_start + window_size
    
            if dfs.at[j, 'Target_start'] < dynamic_window_end:
                j += 1  
            else:
                break  
        window_df = dfs[i:j].copy()
        if 1 <= len(window_df['Query_name'].unique()):
            window_df.index = range(len(window_df))  
            windows[current_start] = window_df  
        i = j  
    return windows


def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if svtype == 'DEL':
            window_end = current_start + window_size + 250 ## avoid small fragment deletions
        else:
            window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 100:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows


def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if current_svlen < 100:
            window_size = 200 + current_svlen * 0.5
        elif 100 < current_svlen <= 500:
            window_size = 400 + current_svlen * 0.5
        elif 500 < current_svlen <= 1000:
            window_size = 600 + current_svlen * 0.5
        else:
            window_size = 1500
        window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        ## allow big sv in one redas support ##
        if svtype == "INS" and window_df['SVlen'].max() >=5000:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
            print("!!!!!!!!!!!! big SV in window !!!!!!!!!!!!")
        if nreads >= 5:
            nfil = 3
        else:
            nfil = 2
        if  nfil <= len(window_df['Query_name'].unique()) < depth * 10:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows

def windows_slide4tra(df, shift=2000000):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df = df.copy()
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs = dfs.sort_values(by=['Target_start1', 'Target_start2'])
        dfs.index = range(len(dfs))
        i = 0
        while i < len(dfs):
            current_chrom = str(dfs.at[i, '#Target_name1'])
            current_start = dfs.at[i, 'Target_start1']
            window_size = shift
            window_end = current_start + window_size
            j = i
            while j < len(dfs) and dfs.at[j, 'Target_start1'] < window_end:
                j += 1
            window_df = dfs[i:j].copy()
            window_df = window_df.reset_index(drop=True)
            key = f'{chrom1}_{current_start}_{chrom2}'
            windows[key] = window_df
            i = j
    return windows

def klook_clu_tra(win):
    """
    tra clus, the clusdf should be chromosome pair mode,
    each clusdf only allow two chromosome.
    """
    if len(win) == 1:
        win['shift_cluster'] = -1
        return win
    cluster_id = 0
    win.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(win)):
        current_start1 = win.loc[i, 'Target_start1']
        current_start2 = win.loc[i, 'Target_start2']
        found_cluster = False
        start_index = max(0, i - 40)
        for j in range(i - 1, start_index - 1, -1):
            max_diff = 1000
            ## already ensure chrom1 != chrom2
            pos_condition = (abs(current_start1 - win.loc[j, 'Target_start1']) <= max_diff and
                             abs(current_start2 - win.loc[j, 'Target_start2']) <= max_diff)
            if  pos_condition:
                win.loc[i, 'shift_cluster'] = win.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            win.loc[i, 'shift_cluster'] = cluster_id
    return win

def candidate_tra(window, opened_bam, dtype):
    """
    Process translocation data and find candidate translocations.
    """
    win_clus = klook_clu_tra(window)
    tras = []
    for clu in win_clus['shift_cluster'].unique():
        win = win_clus[win_clus['shift_cluster']==clu]
        chr1 = win['#Target_name1'].iloc[0]
        chr1_start = mode_or_median(win['Target_start1'])
        sv_len = 0
        chr2 = win['#Target_name2'].iloc[0]
        chr2_start = mode_or_median(win['Target_start2'])
        maq = mode_or_median(win['maq'])
        readsname = set(win['Query_name'])
        svid = f'{chr2}:{chr2_start}_{chr1}:{chr1_start}'
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, chr1, max(0, chr1_start - 250), max(chr1_start - 150, 0))
        end_local_map = local_cov(opened_bam, chr2, chr2_start + 150, chr2_start + 250)
        local_depth = np.mean([start_local_map, end_local_map])
        if local_depth >0:
            SV_rate = round(sv_eye / local_depth, 2)
        else:
            SV_rate = 1
        if dtype in ['ont','hifi', 'pb']:
            ## to one depth ##
            if sv_eye >= (local_depth * 0.1 + 0.7 ):
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
        elif dtype in ['cr', 'sr']:
            if sv_eye > local_depth * 0.25:
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tras

def load_and_process_sv_data(args):
    from math import floor
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    except pd.errors.EmptyDataError:
        print(f"Warning: File {args.raw_signal} is empty.")
        sv_indel_data = pd.DataFrame()
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else ceil(float(depth_stat.iloc[0, 3])+0.3)
        if args.nreads:
            nreads_fil = args.nreads
        else:
            nreads_fil = floor(depth / 10)
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth, nreads_fil
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data = sv_indel_data.copy()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_indel_data['Target_end'] = sv_indel_data['Target_end'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"***************** empty {supp_align_file} file ******************")
        else:
            if not msv.empty:
                print(f'*************** {args.raw_signal}.suppAlign has {msv.shape[0]} rows ********************')
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                print(f'*************** {args.raw_signal}.suppAlign after duplicates drop has {msv.shape[0]} rows ********************')
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['Target_end'] = sv_data[svtype]['Target_end'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    print(f'******************** all chromosomes list {chroms} ****************************')
    return sv_data, chroms, depth, nreads_fil

def process_svtype(args, sv_data, chroms, svtype, depth, nreads, minLen):
    """
        cov info parse by covfile or bam file
    """
    print(f"start klook for {args.raw_signal}  SV type: {svtype}") 
    if args.dtype in ['cr', 'sr']:
        try:
            covinfo = pd.read_csv(args.covfile, sep="\t", index_col=None, header=None)
            covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
            covinfo['target_chr'] = covinfo['target_chr'].astype(str)
            bam_path = covinfo
        except FileNotFoundError:
            print(f"Error: Coverage file {args.covfile} not found.")
            return [], [], [], [], []
    else:
        bam_path = args.bam

    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    if args.dtype in ['pb', 'ont', 'hifi']:
        print(f'data type is {args.dtype}')
        try:
            with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
                for chrom in chroms:
                    chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                    def process_sv_type(svtype):
                        if svtype == "TRA":
                            log = open("log_tra", 'w')
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            windows = windows_slide4tra(chrom_data[svtype], 2000000)
                            print(f'The TRA windows by 2M: \n{windows}', file=log)
                            tra_list = []
                            for win in windows.values():
                                win_clus = klook_clu_tra(win)
                                print(win_clus.iloc[:,[0,2,3,4,5,-1]], file=log)
                                tra = candidate_tra(win_clus, opened_bam, args.dtype)
                                print(tra, file=log)

                                if tra:
                                    tra_list.extend(tra)
                            log.close()
                            return tra_list
                        else:
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                            sv_windows = windows_slide(sv_dfs, depth, svtype, nreads,args.window_size)
                            candidate_svs = []
                            if args.dtype in ['hifi', 'ont']:
                                hdepth = 5
                            else:
                                hdepth = 10
                            for sv_window in sv_windows.values():
                                if len(sv_window['Query_name']) > hdepth: 
                                    candisv = highdepth_clu(sv_window,args.num_hap, svtype, opened_bam, nreads, args.msvfold, args.rate_depth, 1)
                                else:
                                    candisv = lowdepth_clu(sv_window, args.num_hap, svtype, opened_bam, nreads, args.msvfold, args.rate_depth, 1)
                                candidate_svs.extend(candisv)
                            return candidate_svs

                    if svtype == "DEL":
                        del_clus.extend(process_sv_type("DEL"))
                    elif svtype == "INS":
                        ins_clus.extend(process_sv_type("INS"))
                    elif svtype == "INV":
                        inv_clus.extend(process_sv_type("INV"))
                    elif svtype == "DUP":
                        dup_clus.extend(process_sv_type("DUP"))
                    elif svtype == "TRA":
                        tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f"Error processing BAM file: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    elif args.dtype in ['sr', 'cr']:
        print(f'data type is {args.dtype}')
        try:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                def process_sv_type(svtype):
                    if svtype != "TRA":
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                        if sv_dfs.empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        sv_dfs = sv_dfs.sort_values(by=['Target_start', 'SVlen'])
                        sv_dfs.index = range(len(sv_dfs))
                        print(sv_dfs.head(10))
                        print("**************************** Calling one depth all clustering ***********************")
                        candidate_svs = []
                        win_dfs = windows_slide4asm(sv_dfs, svtype, args.window_size)
                        for win_df in win_dfs.values():
                            print(f'**************************window signals*********************\n{win_df}')
                            candidate_sv = onedepth_all_clus(win_df, svtype, bam_path, args.msvfold, max(args.rate_depth,0.25))
                            candidate_svs.extend(candidate_sv)
                        return candidate_svs
                    else:
                        if chrom_data[svtype].empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        windows = windows_slide4tra(chrom_data[svtype], 1000000)
                        tra_list = []
                        for value in windows.values():
                            win = value
                            tra = candidate_tra(win, bam_path,args.dtype)
                            if tra:
                                tra_list.extend(tra)
                        print(f'****************************** the tra condidate ***************************')
                        print(tra_list)
                        return tra_list
                if svtype == "DEL":
                    del_clus.extend(process_sv_type("DEL"))
                elif svtype == "INS":
                    ins_clus.extend(process_sv_type("INS"))
                elif svtype == "INV":
                    inv_clus.extend(process_sv_type("INV"))
                elif svtype == "DUP":
                    dup_clus.extend(process_sv_type("DUP"))
                elif svtype == "TRA":
                    tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f'Error in processing {args.dtype}: {args.raw_signal}')
            print(f"Error processing {svtype} for chromosome {chrom}: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        print(f"Error: sequence data dtype error: {args.dtype} is not in [pb,hifi,ont,sr,cr], please check parameter -dtype")
        return [],[],[],[],[]

def candidateSV(args):
    sv_data, chroms, depth, nreads = load_and_process_sv_data(args)
    print(chroms, depth, nreads)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth, nreads, args.min) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            if result is not None:
                tra_clus += result[0]
                del_clus += result[1]
                ins_clus += result[2]
                inv_clus += result[3]
                dup_clus += result[4]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        return [], [], [], [], []


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True,
                    help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=800, type=int,
                    help="the distance shift of breakpoint to cluster the TRA/big_INV/big_DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=10000000, help="the max SV length")
    IN.add_argument("-m", dest="min", type=int, default=45, help="the minimum SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    IN.add_argument("--nreads", dest="nreads", type=int, help="the minimum numbers of reads to support SV, if not provided, we use average_depth / 10 as threshold" )
    IN.add_argument("--rate_depth", dest="rate_depth", type=float, default=0.1, help="the sv supports of local depth ratio to support sv, 0.1 means the percent of local reads shoule support sv")
    IN.add_argument("--window", dest="window_size", type=int, default=500, help="the window size of signal to parse in klook cluster, 500bp suggested")
    IN.add_argument("--num_hap", dest="num_hap", type=int, default=2, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    IN.add_argument("--msvfold", dest="msvfold", type=float, default=0.8, help="the SV size fold change of multi sv alleles")

    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

```

### ./PSV_Signal/0.Signal4bam_PSVGT.py
```python
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
    IN.add_argument("-w", dest="workers", help="Number of worker processes", default=10, type=int)
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
    with Pool(processes=args.workers) as pool:
        pool.starmap(sub_Signal4bam_PSVGT.process_chromosome, 
                     [(chromosome,chrom_size_dict[chromosome],chromosome_list, args.bam, args.min,args.max, args.maq, SVsignal_out_path,args.dtype,args.msv) 
                      for chromosome in chromosome_list])
    exe_end = time()
    print(f"{'*' * 40} done SV searching {'*' * 40}\ncost time: {exe_end - exe_start}")

if __name__ == "__main__":
    exe_start = time()
    main()

```

### ./PSV_Signal/0.Signal4bam_PSVGT_v2.py
```python
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
    IN.add_argument("-w", dest="workers", help="Number of worker processes", default=10, type=int)
    return parser.parse_args()


def main():
    args = parse_arguments()
    # Load chromosome list from the fai file
    fai = pd.read_csv(args.fai, header=None, sep="\t", dtype=str, index_col=None)
    chromosome_list = fai[0].tolist()
    chrom_size_dict = dict(zip(fai[0], fai[1].astype(int)))
    if args.out[-4:] == ".bam":
        args.out = args.out.replace('.bam', '')
    SVsignal_out_path = f"{args.out}"
    # Use multiprocessing to parallelize the chromosome processing
    with Pool(processes=args.workers) as pool:
        pool.starmap(sub_Signal4bam_PSVGT.process_chromosome,
                     [(chromosome, chrom_size_dict[chromosome], chromosome_list, args.bam, args.min, args.max, args.maq, SVsignal_out_path, args.dtype, args.msv)
                      for chromosome in chromosome_list])
    exe_end = time()
    print(f"{'*' * 40} done SV searching {'*' * 40}\ncost time: {exe_end - exe_start}")


if __name__ == "__main__":
    exe_start = time()
    main()

```

### ./PSV_Signal/0.KLookCluster_LocalDepthAdaptive.py
```python
import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil,floor

def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cost time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        try:
            num_maps = covinfo.count(str(reference), start, end)
        except AttributeError:
            print(f"Error: covinfo does not have a 'count' method.")
            num_maps = 0
    return num_maps

def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)

def candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate=0.1,add=1):
    """
    Solve clusters dataframe
    Get each cluster's sv breakpoints, SV length
    For low depth clus, the cluster must be vary strict and intense, else it will generate a lot of false positive sv. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts() ## default reversed count
    min_clusters = min(num_hap, len(cluster_counts))
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.05:
        print("proportion is too low, return []")
        return []
    else: 
        print(f"top1 cluster percent is {proportion}")
    
    top_clusters = cluster_counts.head(min_clusters).index
    clusdf = clusdf[clusdf[cluster_col].isin(top_clusters)]
    clusters_to_process =  clusdf[cluster_col].unique()

    msv = []
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = mode_or_median(clu_df['Target_end'])
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        if sv_eye < 2:
            continue
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 150), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 150)

        
        SV_rate = round(sv_eye / min(start_local_map+1,end_local_map+1), 2)
        if sv_eye >= min([start_local_map,end_local_map])*support_rate + add:
            print(f"{svid} sv signal propertion support")
            clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
            msv.append(svid)
        if sv_len > 2000 and svtype=="INS" and sv_eye >= 2 and sv_eye >= min([start_local_map,end_local_map])*support_rate*0.5: ## to capture big INS
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} big INS,local_map is {min([start_local_map,end_local_map])} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        if sv_eye >= nreads -1: ## filter 
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid},local coverage is {min([start_local_map,end_local_map])} ")
                #clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                #msv.append(svid)
            else:
                print(f"{svid}: number sv reads lower than required")
    if len(set(msv)) >= 2:
        outmsv = "\t".join(list(set(msv)))
        print(f'msv:\t{outmsv}')
    return clus

def klook_clusters(clusdf, max_diff_func, len_fold=0.8):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    num_signals = len(clusdf)
    if num_signals == 1:
        clusdf['shift_cluster'] = -1
        return clusdf
    
    # Precompute the 'Target_start' and 'Target_end' columns for quicker access
    start_values = clusdf['Target_start'].values
    end_values = clusdf['Target_end'].values
    svlen_values = clusdf['SVlen'].values
    
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    
    for i in range(1, len(clusdf)):
        current_svlen = svlen_values[i]
        current_start = start_values[i]
        current_end = end_values[i]
        found_cluster = False
        # Adjust the previous cluster range dynamically
        start_index = max(0, i - ceil(num_signals * 0.8))  # previous 80% record

        # Vectorized approach: avoid nested for loops
        for j in range(i - 1, start_index - 1, -1):
            old_len = svlen_values[j]
            relate_size = min(current_svlen / old_len, old_len/current_svlen)
            max_diff = max_diff_func(old_len)
            len_condition = relate_size > len_fold
            pos_condition = (abs(current_start - start_values[j]) <= max_diff or
                             abs(current_end - end_values[j]) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf

def max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 50 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 100 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 450 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 1000

def max_diff_func4INS(old_len):
    if old_len <= 100:
        return 50 + (old_len-50) * 0.5
    elif 100 < old_len <= 500:
        return 100 + (old_len-100) * 0.5
    elif 500 < old_len <= 5000:
        return 300 + (old_len-500) * 0.1
    elif old_len > 5000:
        return 1000

def onedepth_all_clus(all_signal,svtype, opened_bam, nrate=0.25):
    """
    Process the contig or genome level SV signal.
    Strict condition for merge, single chromosome mode.
    """
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    if all_signal.shape[0] > 1:
        all_signal = merge_and_sort_cr(all_signal)
        print(f"******************signal after merging ****************\n {all_signal}")
        if all_signal.shape[0] > 1:
            all_clus = klook_clusters(all_signal, max_diff_func, 0.8)
        else:
            all_signal['shift_cluster'] = -1
            all_clus = all_signal
    else:
        all_signal['shift_cluster'] = -1
        all_clus = all_signal
    print(all_clus[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    sv_chrom = all_clus["#Target_name"].unique()[0]
    svtype = all_clus["SVType"].unique()[0]
    cluster_col = all_clus.columns[-1]
    clus = []
    msv = []
    for clu in all_clus[cluster_col].unique():
        clu_df = all_clus[all_clus[cluster_col] == clu] 
        print(clu_df[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])   
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len =   mode_or_median(clu_df['SVlen'])
        sv_end =   mode_or_median(clu_df['Target_end'])
        maq =      mode_or_median(clu_df['maq'])
        readsname = clu_df['Query_name'].tolist()
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 250), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 250)
        depth = max([start_local_map,end_local_map])
        if depth > 0:
            SV_rate = round(sv_eye / depth, 2)
        else:
            SV_rate = 1
        if SV_rate < nrate: ## 0.25 to meet Tetraploid heterozygous
            continue
        print(f"************************cluster done***************************\n{[sv_chrom, sv_start, sv_end, sv_len, svid, svtype, '*', sv_eye, SV_rate, maq, readsname]}")
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
        msv.append(svid)
    if len(msv) >=2:
        outmsv = "\t".join(msv)
        print(f'msv:\t{outmsv}')
    return clus

def lowdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL    

    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def highdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])
    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def filter_overlaps(df):
    """
    to filter overlap SV within same queryname
    """
    filtered_groups = []
    for query, group in df.groupby('Query_name'):
        if len(group) <= 1:
            filtered_groups.append(group)
            continue
        sorted_group = group.sort_values(
            by=['Target_start', 'SVlen'],
            ascending=[True, False]
        )
        kept = []
        last_end = -1
        for _, row in sorted_group.iterrows():
            start = row['Target_start']
            end = row['Target_end']
            if start > last_end:  
                kept.append(row)
                last_end = end
            else:
                 if end > last_end:
                     kept[-1] = row
                     last_end = end
        filtered_groups.append(pd.DataFrame(kept))    
    return pd.concat(filtered_groups, ignore_index=True)

def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates( subset=["Query_name", "Target_start", "Target_end"], keep="first",inplace=True) ## we dont keep the same reads multiple segment
    print(clusdf[["#Target_name","Target_start","Target_end","SVlen","SVType"]])
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf


def merge_and_sort_cr(clusdf):
    """
    filter overlap SV on same query name;
    merge sv based on query name;
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates(
        subset=["Query_name", "Target_start", "Target_end"],
        keep="first",
        inplace=True
    )
    clusdf = filter_overlaps(clusdf)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'min',       
        'Target_end': 'max',    
        'SVlen': 'sum',              
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv[['#Target_name', 'Target_start', 'Target_end', 'SVlen', 'SVType', 'maq', 'seq', 'Query_name']]
    clusdf = clusdf.astype({
        'Target_start': np.int32,
        'Target_end': np.int32,
        'SVlen': np.int32,
        'maq': np.int16
    })
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf

def windows_slide4asm(dfs, svtype, window_size=500):
    """
    fix window gap for assemblies
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        j = i  
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom:
            last_in_window_start = dfs.at[j-1, 'Target_start'] if j > i else current_start
            if svtype == 'DEL':
                dynamic_window_end = last_in_window_start + window_size + 250
            else:
                dynamic_window_end = last_in_window_start + window_size
    
            if dfs.at[j, 'Target_start'] < dynamic_window_end:
                j += 1  
            else:
                break  
        window_df = dfs[i:j].copy()
        if 1 <= len(window_df['Query_name'].unique()):
            window_df.index = range(len(window_df))  
            windows[current_start] = window_df  
        i = j  
    return windows

def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if svtype == 'DEL':
            window_end = current_start + window_size + 250 ## avoid small fragment deletions
        else:
            window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 100:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows

def windows_slide4tra(df, shift=2000000):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df = df.copy()
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs = dfs.sort_values(by=['Target_start1', 'Target_start2'])
        dfs.index = range(len(dfs))
        i = 0
        while i < len(dfs):
            current_chrom = str(dfs.at[i, '#Target_name1'])
            current_start = dfs.at[i, 'Target_start1']
            window_size = shift
            window_end = current_start + window_size
            j = i
            while j < len(dfs) and dfs.at[j, 'Target_start1'] < window_end:
                j += 1
            window_df = dfs[i:j].copy()
            window_df = window_df.reset_index(drop=True)
            key = f'{chrom1}_{current_start}_{chrom2}'
            windows[key] = window_df
            i = j
    return windows

def klook_clu_tra(win):
    """
    tra clus, the clusdf should be chromosome pair mode,
    each clusdf only allow two chromosome.
    """
    if len(win) == 1:
        win['shift_cluster'] = -1
        return win
    cluster_id = 0
    win.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(win)):
        current_start1 = win.loc[i, 'Target_start1']
        current_start2 = win.loc[i, 'Target_start2']
        found_cluster = False
        start_index = max(0, i - 40)
        for j in range(i - 1, start_index - 1, -1):
            max_diff = 1000
            ## already ensure chrom1 != chrom2
            pos_condition = (abs(current_start1 - win.loc[j, 'Target_start1']) <= max_diff and
                             abs(current_start2 - win.loc[j, 'Target_start2']) <= max_diff)
            if  pos_condition:
                win.loc[i, 'shift_cluster'] = win.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            win.loc[i, 'shift_cluster'] = cluster_id
    return win

def candidate_tra(window, opened_bam, dtype):
    """
    Process translocation data and find candidate translocations.
    """
    win_clus = klook_clu_tra(window)
    tras = []
    for clu in win_clus['shift_cluster'].unique():
        win = win_clus[win_clus['shift_cluster']==clu]
        chr1 = win['#Target_name1'].iloc[0]
        chr1_start = mode_or_median(win['Target_start1'])
        sv_len = 0
        chr2 = win['#Target_name2'].iloc[0]
        chr2_start = mode_or_median(win['Target_start2'])
        maq = mode_or_median(win['maq'])
        readsname = set(win['Query_name'])
        svid = f'{chr2}:{chr2_start}_{chr1}:{chr1_start}'
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, chr1, max(0, chr1_start - 250), max(chr1_start - 150, 0))
        end_local_map = local_cov(opened_bam, chr2, chr2_start + 150, chr2_start + 250)
        local_depth = np.mean([start_local_map, end_local_map])
        if local_depth >0:
            SV_rate = round(sv_eye / local_depth, 2)
        else:
            SV_rate = 1
        if dtype in ['ont','hifi', 'pb']:
            ## to one depth ##
            if sv_eye >= (local_depth * 0.1 + 0.7 ):
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
        elif dtype in ['cr', 'sr']:
            if sv_eye > local_depth * 0.25:
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tras

def load_and_process_sv_data(args):
    from math import floor
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    except pd.errors.EmptyDataError:
        print(f"Warning: File {args.raw_signal} is empty.")
        sv_indel_data = pd.DataFrame()
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else ceil(float(depth_stat.iloc[0, 3])+0.3)
        if args.nreads:
            nreads_fil = args.nreads
        else:
            nreads_fil = floor(depth / 10)
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth, nreads_fil
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data = sv_indel_data.copy()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_indel_data['Target_end'] = sv_indel_data['Target_end'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"***************** empty {supp_align_file} file ******************")
        else:
            if not msv.empty:
                print(f'*************** {args.raw_signal}.suppAlign has {msv.shape[0]} rows ********************')
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                #msv = msv.drop_duplicates()
                print(f'*************** {args.raw_signal}.suppAlign after duplicates drop has {msv.shape[0]} rows ********************')
                #msv = msv.reset_index(drop=True)
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['Target_end'] = sv_data[svtype]['Target_end'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    print(f'******************** all chromosomes list {chroms} ****************************')
    return sv_data, chroms, depth, nreads_fil

def process_svtype(args, sv_data, chroms, svtype, depth, nreads, minLen):
    """
        cov info parse by covfile or bam file
    """
    print(f"start klook for {args.raw_signal}  SV type: {svtype}") 
    if args.dtype in ['cr', 'sr']:
        try:
            covinfo = pd.read_csv(args.covfile, sep="\t", index_col=None, header=None)
            covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
            covinfo['target_chr'] = covinfo['target_chr'].astype(str)
            bam_path = covinfo
        except FileNotFoundError:
            print(f"Error: Coverage file {args.covfile} not found.")
            return [], [], [], [], []
    else:
        bam_path = args.bam

    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    if args.dtype in ['pb', 'ont', 'hifi']:
        print(f'data type is {args.dtype}')
        try:
            with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
                for chrom in chroms:
                    chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                    def process_sv_type(svtype):
                        if svtype == "TRA":
                            log = open("log_tra", 'w')
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            windows = windows_slide4tra(chrom_data[svtype], 2000000)
                            print(f'The TRA windows by 2M: \n{windows}', file=log)
                            tra_list = []
                            for win in windows.values():
                                win_clus = klook_clu_tra(win)
                                print(win_clus.iloc[:,[0,2,3,4,5,-1]], file=log)
                                tra = candidate_tra(win_clus, opened_bam, args.dtype)
                                print(tra, file=log)

                                if tra:
                                    tra_list.extend(tra)
                            log.close()
                            return tra_list
                        else:
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                            sv_windows = windows_slide(sv_dfs, depth, svtype, nreads,args.window_size)
                            candidate_svs = []
                            if args.dtype == 'hifi':
                                hdepth = 5
                            else:
                                hdepth = 10
                            for sv_window in sv_windows.values():
                                if len(sv_window['Query_name']) > hdepth: 
                                    candisv = highdepth_clu(sv_window,args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
                                else:
                                    candisv = lowdepth_clu(sv_window, args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
                                candidate_svs.extend(candisv)
                            return candidate_svs

                    if svtype == "DEL":
                        del_clus.extend(process_sv_type("DEL"))
                    elif svtype == "INS":
                        ins_clus.extend(process_sv_type("INS"))
                    elif svtype == "INV":
                        inv_clus.extend(process_sv_type("INV"))
                    elif svtype == "DUP":
                        dup_clus.extend(process_sv_type("DUP"))
                    elif svtype == "TRA":
                        tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f"Error processing BAM file: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    elif args.dtype in ['sr', 'cr']:
        print(f'data type is {args.dtype}')
        try:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                def process_sv_type(svtype):
                    if svtype != "TRA":
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                        if sv_dfs.empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        sv_dfs = sv_dfs.sort_values(by=['Target_start', 'SVlen'])
                        sv_dfs.index = range(len(sv_dfs))
                        print(sv_dfs.head(10))
                        print("**************************** Calling one depth all clustering ***********************")
                        candidate_svs = []
                        win_dfs = windows_slide4asm(sv_dfs, svtype, args.window_size)
                        for win_df in win_dfs.values():
                            print(f'**************************window signals*********************\n{win_df}')
                            candidate_sv = onedepth_all_clus(win_df, svtype, bam_path,max(args.rate_depth,0.25))
                            candidate_svs.extend(candidate_sv)
                        return candidate_svs
                    else:
                        if chrom_data[svtype].empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        windows = windows_slide4tra(chrom_data[svtype], 1000000)
                        tra_list = []
                        for value in windows.values():
                            win = value
                            tra = candidate_tra(win, bam_path,args.dtype)
                            if tra:
                                tra_list.extend(tra)
                        print(f'****************************** the tra condidate ***************************')
                        print(tra_list)
                        return tra_list
                if svtype == "DEL":
                    del_clus.extend(process_sv_type("DEL"))
                elif svtype == "INS":
                    ins_clus.extend(process_sv_type("INS"))
                elif svtype == "INV":
                    inv_clus.extend(process_sv_type("INV"))
                elif svtype == "DUP":
                    dup_clus.extend(process_sv_type("DUP"))
                elif svtype == "TRA":
                    tra_clus.extend(process_sv_type("TRA"))
        except Exception as e:
            print(f'Error in processing {args.dtype}: {args.raw_signal}')
            print(f"Error processing {svtype} for chromosome {chrom}: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        print(f"Error: sequence data dtype error: {args.dtype} is not in [pb,hifi,ont,sr,cr], please check parameter -dtype")
        return [],[],[],[],[]

def candidateSV(args):
    sv_data, chroms, depth, nreads = load_and_process_sv_data(args)
    print(chroms, depth, nreads)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth, nreads, args.min) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            if result is not None:
                tra_clus += result[0]
                del_clus += result[1]
                ins_clus += result[2]
                inv_clus += result[3]
                dup_clus += result[4]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        return [], [], [], [], []


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True,
                    help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=800, type=int,
                    help="the distance shift of breakpoint to cluster the TRA/big_INV/big_DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=10000000, help="the max SV length")
    IN.add_argument("-m", dest="min", type=int, default=45, help="the minimum SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    IN.add_argument("--nreads", dest="nreads", type=int, help="the minimum numbers of reads to support SV, if not provided, we use average_depth / 10 as threshold" )
    IN.add_argument("--rate_depth", dest="rate_depth", type=float, default=0.1, help="the sv supports of local depth ratio to support sv, 0.1 means the percent of local reads shoule support sv")
    IN.add_argument("--window", dest="window_size", type=int, default=500, help="the window size of signal to parse in klook cluster, 500bp suggested")
    IN.add_argument("--num_hap", dest="num_hap", type=int, default=2, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

```

## ./PSV_Signal/parse_sam_signal
```shell
total 53036
-rwxr-xr-x 1 lgb xinwang     2178 Aug 14 10:47 0.Signal4bam_PSVGT.py
-rw-r--r-- 1 lgb xinwang 26676270 Aug 14 16:03 cigar_error.log
-rwxr-xr-x 1 lgb xinwang      633 Aug 14 17:11 filter_long_cigar_sam.py
-rw-r--r-- 1 lgb xinwang 25335711 Aug 14 14:41 InDel
-rw-r--r-- 1 lgb xinwang  2216270 Aug 14 14:41 InDel.cov
-rw-r--r-- 1 lgb xinwang      680 Aug 14 17:15 longreads_265M_cigar_out.py
drwxr-xr-x 2 lgb xinwang     4096 Aug 14 16:02 __pycache__
-rwxr-xr-x 1 lgb xinwang    44059 Aug 14 16:02 sub_Signal4bam_PSVGT.py
-rw-r--r-- 1 lgb xinwang        0 Aug 14 16:02 test.record.txt
-rw-r--r-- 1 lgb xinwang        0 Aug 14 16:02 test.record.txt.cov
-rw-r--r-- 1 lgb xinwang        0 Aug 14 16:02 test.record.txt.depth
-rw-r--r-- 1 lgb xinwang        0 Aug 14 16:02 test.record.txt.suppAlign
-rw-r--r-- 1 lgb xinwang    11520 Sep  4 21:11 test.txt

```
### ./PSV_Signal/parse_sam_signal/longreads_265M_cigar_out.py
```python
import re
import sys
insam = sys.argv[1]
outsam = sys.argv[2]
MAX_HS_LEN = 265000000
cigar_pattern = re.compile(r'(\d+)([MIDNSHPX=])')
with open(insam, 'r') as fin, open(outsam, 'w') as fout:
    for line in fin:
        if line.startswith('@'):
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        cigar = parts[5]
        valid = True
        for length_str, op in cigar_pattern.findall(cigar):
            if op in ['H', 'S']:
                if int(length_str) > MAX_HS_LEN:
                    valid = False
                    break
        if valid:
            fout.write(line)

```

### ./PSV_Signal/parse_sam_signal/sub_Signal4bam_PSVGT.py
```python
import pysam
import re

def process_chromosome(chromosome,chrom_size,chromosome_list, bamfile_path, minLen, maxLen, min_maq, SVsignal_out_path,dtype,msv):
    samfile = pysam.AlignmentFile(bamfile_path, 'rb')
    with open(f"{SVsignal_out_path}.record.txt", 'w') as indel_out, open(f"{SVsignal_out_path}.record.txt.cov", 'w') as cov_out, open(f"{SVsignal_out_path}.record.txt.suppAlign", 'w') as supp_sv_out, open(f'{SVsignal_out_path}.record.txt.depth','w') as depth_out:
        # Fetch all reads from the sam file
        lines = samfile.fetch()
        total_map_lens = 0
        for line in lines:
            if dtype in ['sr', 'cr']:
                result = svInDel4asm(line, minLen, min_maq,msv, maxLen, chromosome_list)
                if result is None:
                    continue
                results, cov_line, supp_svsignal = result
                if results:
                    indel_out.writelines(results)
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")   
                cov_out.writelines([cov_line])
            elif dtype in ['pb', 'ont', 'hifi']:
                svInDels, cov_line, supp_svsignal = svInDel4lr(line, minLen, min_maq, maxLen, msv, chromosome_list)
                if cov_line:
                    #print(cov_line)
                    total_map_lens += int(cov_line.split("\t")[4]) - int(cov_line.split("\t")[3])
                if svInDels:
                    indel_out.writelines(svInDels)
                if supp_svsignal:
                    for sv in supp_svsignal:
                        supp_sv_out.writelines(f"{sv}\n")
                if cov_line:
                    cov_out.writelines(cov_line)
        depth = round( total_map_lens / chrom_size, 2 )
        depth_out.write(f'{chromosome}\t{chrom_size}\t{total_map_lens}\t{depth}\n')
        depth_out.close()
        indel_out.close()
        cov_out.close()
        supp_sv_out.close()

        if dtype == 'cr' and msv == 'yes':
            with open(f"{SVsignal_out_path}_{chromosome}.record.txt.suppAlign", 'a') as out:
                lines = samfile.fetch()
                svs = segment_segment_sv(lines, minLen, maxLen, chromosome_list)
                for sv in svs:
                    out.writelines(f"{sv}\n")
            out.close()

import pysam
import re
import logging
logging.basicConfig(
    filename='cigar_error.log',
    level=logging.WARNING,
    format='%(asctime)s - %(message)s'
)

def is_valid_cigar(cigar, max_operation_length=1000000):
    if not cigar:
        return False  
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    for length, op in cigar_ops:
        try:
            length = int(length)
            if length > max_operation_length:
                logging.warning(f"Invalid CIGAR: {length}{op} (exceeds max length {max_operation_length})")
                return False
        except ValueError:
            logging.warning(f"Invalid CIGAR number: {length}")
            return False
    return True

def process_chromosome(chromosome, chrom_size, chromosome_list, bamfile_path, minLen, maxLen, min_maq, SVsignal_out_path, dtype, msv):
    samfile = pysam.AlignmentFile(bamfile_path, 'rb', check_sq=False)
    try:
        with open(f"{SVsignal_out_path}.record.txt", 'w') as indel_out, \
             open(f"{SVsignal_out_path}.record.txt.cov", 'w') as cov_out, \
             open(f"{SVsignal_out_path}.record.txt.suppAlign", 'w') as supp_sv_out, \
             open(f'{SVsignal_out_path}.record.txt.depth', 'w') as depth_out:

            total_map_lens = 0
            # 遍历reads，过滤异常CIGAR
            for read in samfile.fetch():
                # 跳过未映射的read
                if read.is_unmapped:
                    continue
                # 检查CIGAR合法性
                cigar = read.cigarstring
                if not is_valid_cigar(cigar):
                    continue  # 跳过CIGAR异常的read
                # 检查 Mapping Quality
                if read.mapping_quality < min_maq:
                    continue  # 跳过低质量映射

                # 按dtype处理
                if dtype in ['sr', 'cr']:
                    try:
                        result = svInDel4asm(read, minLen, min_maq, msv, maxLen, chromosome_list)
                    except Exception as e:
                        logging.warning(f"Error processing read {read.query_name}: {e}")
                        continue
                    if result is None:
                        continue
                    results, cov_line, supp_svsignal = result
                    if results:
                        indel_out.writelines(results)
                    if supp_svsignal:
                        for sv in supp_svsignal:
                            supp_sv_out.write(f"{sv}\n")
                    if cov_line:
                        cov_out.write(cov_line)
                elif dtype in ['pb', 'ont', 'hifi']:
                    try:
                        svInDels, cov_line, supp_svsignal = svInDel4lr(read, minLen, min_maq, maxLen, msv, chromosome_list)
                    except Exception as e:
                        logging.warning(f"Error processing read {read.query_name}: {e}")
                        continue
                    if cov_line:
                        # 累加映射长度（假设cov_line第4-5列为start-end）
                        try:
                            start = int(cov_line.split("\t")[3])
                            end = int(cov_line.split("\t")[4])
                            total_map_lens += (end - start)
                        except (IndexError, ValueError):
                            logging.warning(f"Invalid cov_line format: {cov_line}")
                    if svInDels:
                        indel_out.writelines(svInDels)
                    if supp_svsignal:
                        for sv in supp_svsignal:
                            supp_sv_out.write(f"{sv}\n")
                    if cov_line:
                        cov_out.write(cov_line)

            # 计算深度并写入
            depth = round(total_map_lens / chrom_size, 2) if chrom_size != 0 else 0
            depth_out.write(f'{chromosome}\t{chrom_size}\t{total_map_lens}\t{depth}\n')

        # 处理cr类型且msv=yes的情况
        if dtype == 'cr' and msv == 'yes':
            with open(f"{SVsignal_out_path}_{chromosome}.record.txt.suppAlign", 'a') as out:
                # 重新遍历reads，再次过滤CIGAR
                samfile.reset()  # 重置文件指针到开头
                for read in samfile.fetch():
                    if read.is_unmapped:
                        continue
                    cigar = read.cigarstring
                    if not is_valid_cigar(cigar):
                        continue
                    if read.mapping_quality < min_maq:
                        continue
                    try:
                        svs = segment_segment_sv([read], minLen, maxLen, chromosome_list)  # 传入单个read
                        for sv in svs:
                            out.write(f"{sv}\n")
                    except Exception as e:
                        logging.warning(f"Error in segment_segment_sv for read {read.query_name}: {e}")
                        continue

    except Exception as e:
        logging.error(f"Fatal error in process_chromosome: {e}")
        raise  
    finally:
        samfile.close()  

def segmentsv4asm(supp_list, min_size, max_size, chromosome_list):
    """
    Collecting the supplementary alignment to identify the INS, DEL, DUP, INV, TRA signals
    LongRead align to reference that has "SA" will be recorded as dict as follow:
    { 'm64144': [['m64144',2048,'Db-Chr4',8949283,8949610,[0, 327, 11157],'43'],
                 ['m64144',0,   'Db-Chr4',8949649,8957962,[3162, 8322, 0],60,'AAXXXCTAATT'],
                 ['m64144',2064,'Db-Chr3',12497386,12502117,[5171, 4731, 1582],'60']]
                 }
    these code is referenced on DeBreak, the INS and DEL may be discarded in the future
    """
    if not any(int(supp[1]) <= 16 for supp in supp_list) or len(supp_list) <= 1:
        return []

    # Sort the supplementary list by chromosome and then by the start position
    supp_list_sorted = sorted(supp_list, key=lambda x: (x[2], int(x[3])))
    svsignal_supp = []

    def add_svcall(chrom, readname, start, end, size, maq, svid, sv_type, sequence=None):
        if maq == 60:
            sv_record = f"{chrom}\t{readname}\t{start}\t{end}\t{size}\t{maq}\t{svid}\t{sv_type}"
            if sequence:
                sv_record += f"\t{sequence}"
            svsignal_supp.append(sv_record)

    # Iterate through the sorted list, each smaller start position is treated as the primary segment
    for idx, primary_map in enumerate(supp_list_sorted):
        readname = primary_map[0]
        pri_chrom = primary_map[2]
        pri_start = int(primary_map[3])
        readseq = primary_map[7] if len(primary_map)==8 else "*"
        #print(f"********************************* flag: {primary_map[1]}*************************************")
        pri_flag = int(primary_map[1]) % 32 > 15
        
        # Separate the supplementary segments
        supps = [supp for supp in supp_list_sorted if supp != primary_map]
        
        samedir, invdir, diffchr = [], [], []

        # Classify supplementary segments into samedir, invdir, and diffchr categories
        for supp in supps:
            chrom, supp_flag, supp_maplen = supp[2], int(supp[1]) % 32 > 15, supp[5][1]
            if supp_maplen < 250:  # Ignore if the map length is less than 250
                continue
            if chrom != pri_chrom:
                diffchr.append(supp)
            elif supp_flag != pri_flag:
                invdir.append(supp)
            else:
                samedir.append(supp)

        # Process samedir: same direction alignments
        for supp_map in samedir:
            leftmap, rightmap = (primary_map, supp_map) if supp_map[3] > primary_map[3] else (supp_map, primary_map)
            sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
            sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
            maq1, maq2 = leftmap[6], rightmap[6]
            maq = (int(maq1) + int(maq2)) // 2
            
            # INS (Insertion) Detection
            if abs(rightmap[3] - leftmap[4]) <= 300:
                overlapmap = rightmap[3] - leftmap[4]
                ins_size = sh3 - len1 - sh1 - overlapmap
                if min_size <= ins_size <= max_size:
                    sv_start = min(rightmap[3], leftmap[4])
                    sv_end = sv_start + 1
                    svid = f'{pri_chrom}:{sv_start}-{sv_end}_INSLEN={ins_size}'
                    seq = readseq[sh1 + len1: sh3 - overlapmap]
                    add_svcall(pri_chrom, readname, sv_start, sv_end, ins_size, maq, svid, "INS", "*")
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

def segmentsv4lr(supp_list,min_size, max_size, chromosome_list, minimaq):
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
            if maq >= minimaq:
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
            if maq < minimaq:
                continue
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
        #if maq < min_maq:
        #    continue
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
            if 2<= len(supp_dict[readname])<=20:
                supp_svsignal = segmentsv4lr(supp_dict[readname],minLen, maxLen, chromosome_list,min_maq)
    return svInDels, covinfo , supp_svsignal


def svInDel4asm(line, minLen, min_maq, msv, maxLen, chromosome_list):
    if line.flag == 4 or line.mapping_quality < min_maq:
        return None  
    else:
        clipinfo = [0,0,0]
        supp_svsignal = []
        if line.cigar[0][0] in [4,5]:  ## left most cigar
            clipinfo[0] = line.cigar[0][1]
        if line.cigar[-1][0] in [4,5]:  ## right most cigar
            clipinfo[2] = line.cigar[-1][1]
        clipinfo[1] = line.query_alignment_length
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
        
        if msv == "no":
            return results, f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n', []  
        
        elif msv == "yes":
            supp_dict = {}
            readname, refend = query_chr,target_end
            if line.has_tag("SA"):
                primary = [readname, line.flag, line.reference_name,target_start,refend,clipinfo, maq, line.query_sequence ]
                if readname in supp_dict.keys():
                    supp_dict[readname] += [primary]
                else:
                    supp_dict[readname] = [primary]
                #print(primary[0])
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
                    suppinfo = [readname, sflag, chrom, start, end, clipinfo, maq, "*"]
                    #print(suppinfo)
                    supp_dict[readname] += [suppinfo]
            if supp_dict:
                if 2<= len(supp_dict[readname]):
                    supp_svsignal = segmentsv4asm(supp_dict[readname],minLen, maxLen, chromosome_list)
            return results,  f'{query_chr}\t{flag}\t{target_chr}\t{target_start}\t{target_end}\t{maq}\t{cigar}\n', supp_svsignal



def segment_segment_sv(lines, min_size, max_size, chromosome_list):
    """
    Collecting the supplementary alignment to identify the INS, DEL, DUP, INV, TRA signals
    LongRead align to reference that has "SA" will be recorded as dict as follow:
    { 'm64144': [['m64144',2048,'Db-Chr4',8949283,8949610,[0, 327, 11157],'43'],
                 ['m64144',0,   'Db-Chr4',8949649,8957962,[3162, 8322, 0],60,'AAXXXCTAATT'],
                 ['m64144',2064,'Db-Chr3',12497386,12502117,[5171, 4731, 1582],'60']]
                 }
    """
    supp_dict = {}
    for line in lines:
        clipinfo = [0, 0, 0]
        if line.flag != 4:
            if line.cigar[0][0] in [4, 5]:  ## left most cigar
                clipinfo[0] = line.cigar[0][1]
            if line.cigar[-1][0] in [4, 5]:  ## right most cigar
                clipinfo[2] = line.cigar[-1][1]
            clipinfo[1] = line.query_alignment_length
            query_chr = line.query_name
            query_seq = line.query_sequence
            target_chr = line.reference_name
            target_start = line.reference_start
            target_end = line.reference_end
            maq = line.mapping_quality
            flag = line.flag
            supp = [query_chr, flag, target_chr, target_start, target_end, clipinfo, maq, query_seq]
            if query_chr in supp_dict.keys():
                supp_dict[query_chr] += [supp]
            else:
                supp_dict[query_chr] = [supp]

    svsignal_supp = []

    def add_svcall(chrom, readname, start, end, size, maq, svid, sv_type, sequence=None):
        if maq == 60:
            sv_record = f"{chrom}\t{readname}\t{start}\t{end}\t{size}\t{maq}\t{svid}\t{sv_type}"
            if sequence:
                sv_record += f"\t{sequence}"
            svsignal_supp.append(sv_record)
            
    print("******************** Calling segment segment SV *************************")

    for query_chr, query_supp_list in supp_dict.items():
        supp_list_sorted = sorted(query_supp_list, key=lambda x: (x[2], int(x[3])))

        for idx, primary_map in enumerate(supp_list_sorted):
            readname = primary_map[0]
            pri_chrom = primary_map[2]
            pri_start = int(primary_map[3])
            readseq = primary_map[7] if len(primary_map) == 8 else "*"
            #print(f"********************************* flag: {primary_map[1]}*************************************")
            try:
                pri_flag = int(primary_map[1]) % 32 > 15
            except ValueError:
                print(f"Invalid flag value: {primary_map[1]}. Skipping this entry.")
                continue

            # Separate the supplementary segments
            supps = [supp for supp in supp_list_sorted if supp != primary_map]

            samedir, invdir, diffchr = [], [], []

            # Classify supplementary segments into samedir, invdir, and diffchr categories
            for supp in supps:
                #print(f'{primary_map[0:5]}\t{supp[0:5]}')
                try:
                    chrom, supp_flag, supp_maplen = supp[2], int(supp[1]) % 32 > 15, supp[5][1]
                except ValueError:
                    print(f"Invalid flag value: {supp[1]}. Skipping this supp entry.")
                    continue
                if supp_maplen < 250:  # Ignore if the map length is less than 250
                    continue
                if chrom != pri_chrom:
                    diffchr.append(supp)
                elif supp_flag != pri_flag:
                    invdir.append(supp)
                else:
                    samedir.append(supp)

            # Process samedir: same direction alignments
            for supp_map in samedir:
                leftmap, rightmap = (primary_map, supp_map) if supp_map[3] > primary_map[3] else (supp_map, primary_map)
                sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
                sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
                maq1, maq2 = leftmap[6], rightmap[6]
                maq = (int(maq1) + int(maq2)) // 2

                ## INS ##
                if abs(rightmap[3] - leftmap[4]) <= 300:
                    overlapmap = rightmap[3] - leftmap[4]
                    ins_size = sh3 - len1 - sh1 - overlapmap
                    if min_size <= ins_size <= max_size:
                        sv_start = min(rightmap[3], leftmap[4])
                        sv_end = sv_start + 1
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_INSLEN={ins_size}'
                        if readseq:
                            seq = readseq[sh1 + len1: sh3 - overlapmap]
                        else:
                            seq = "*"
                        add_svcall(pri_chrom, readname, sv_start, sv_end, ins_size, maq, svid, "INS", "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                ## DEL ##
                overlapmap = sh1 + len1 - sh3
                if -200 < overlapmap < 1500:  ## check or search for better
                    del_size = rightmap[3] - leftmap[4] + overlapmap
                    if min_size <= del_size <= max_size:
                        sv_start = max(0, leftmap[4] - max(0, overlapmap))
                        sv_end = sv_start + del_size - 1
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DELLEN={del_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, del_size, maq, svid, "DEL", '*')
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                ## DUP ##
                lap1 = sh1 + len1 - sh3
                if -200 < lap1 < 500 and (leftmap[4] - rightmap[3]) >= max(50, lap1):
                    dup_size = leftmap[4] - rightmap[3] - max(lap1, 0)
                    if min_size <= dup_size <= max_size:
                        sv_start, sv_end = rightmap[3], rightmap[3] + dup_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, dup_size, maq, svid, "DUP", "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
                lap2 = sh3 + len2 - sh1
                if -200 < lap2 < 500 and (rightmap[4] - leftmap[3]) >= max(1000, lap2):
                    dup_size = rightmap[4] - leftmap[3] - lap2
                    if min_size <= dup_size <= max_size:
                        sv_start, sv_end = leftmap[3], leftmap[3] + dup_size
                        svid = f'{pri_chrom}:{sv_start}-{sv_end}_DUPLEN={dup_size}'
                        add_svcall(pri_chrom, readname, sv_start,
                                   sv_end, dup_size, maq, svid, 'DUP', "*")
                        print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
            ## INV ##
            for supp in invdir:
                if (supp[3] > primary_map[3] and (supp[4] - primary_map[4]) > -200) or \
                        (supp[3] < primary_map[3] and (primary_map[4] - supp[4]) > -200):
                    leftmap = primary_map if supp[3] > primary_map[3] else supp
                    rightmap = supp if supp[3] > primary_map[3] else primary_map

                    sh1, len1, sh2 = leftmap[5][0], leftmap[5][1], leftmap[5][2]
                    sh3, len2, sh4 = rightmap[5][0], rightmap[5][1], rightmap[5][2]
                    maq1, maq2 = leftmap[6], rightmap[6]
                    readname = primary_map[0]
                    maq = (int(maq1) + int(maq2)) // 2

                    lap1 = sh3 + len2 - sh2
                    if -200 < lap1 < 500 and (rightmap[4] - leftmap[4]) > max(100, lap1):
                        inv_size = rightmap[4] - leftmap[4] - lap1
                        if min_size <= inv_size <= max_size:
                            sv_start, sv_end = leftmap[4], leftmap[4] + inv_size
                            svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                            add_svcall(pri_chrom, readname, sv_start,
                                       sv_end, inv_size, maq, svid, 'INV', "*")
                            print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')

                    lap2 = sh4 + len2 - sh1
                    if -200 < lap2 < 500 and (rightmap[3] - leftmap[3]) >= max(100, lap2):
                        inv_size = rightmap[3] - leftmap[3] - lap2
                        if min_size <= inv_size <= max_size:
                            sv_start, sv_end = leftmap[3], leftmap[3] + inv_size
                            svid = f'{pri_chrom}:{sv_start}-{sv_end}_INVLEN={inv_size}'
                            add_svcall(pri_chrom, readname, sv_start,
                                       sv_end, inv_size, maq, svid, "INV", "*")
                            print(f'{svid}: detect from: {leftmap[0:5]}\t{rightmap[0:5]}')
            ## TRA ##
            ### win size 500 may fit into the TRA ###
            for supp in diffchr:
                maq1, maq2 = supp[6], primary_map[6]
                readname = primary_map[0]
                maq = (int(maq1) + int(maq2)) // 2
                sh1, len1, sh2 = primary_map[5][0], primary_map[5][1], primary_map[5][2]
                sh3, len2, sh4 = supp[5][0], supp[5][1], supp[5][2]
                breakpoint1, breakpoint2 = '', ''
                if abs(sh1 - sh3 - len2) <= 500 or abs(sh1 - len2 - sh4) <= 500:
                    chrom1 = primary_map[2]
                    breakpoint1 = primary_map[3]
                elif abs(sh2 - sh3 - len2) <= 500 or abs(sh2 - len2 - sh4) <= 500:
                    chrom1 = primary_map[2]
                    breakpoint1 = primary_map[4]
                if abs(sh3 - sh1 - len1) <= 500 or abs(sh3 - len1 - sh2) <= 500:
                    chrom2 = supp[2]
                    breakpoint2 = supp[3]
                elif abs(sh4 - sh1 - len1) <= 500 or abs(sh4 - len1 - sh2) <= 500:
                    chrom2 = supp[2]
                    breakpoint2 = supp[4]
                    ## Chr4	ERR3415829.505585   15689499    9991420	0	60	Chr4:15689499_Chr1:9991420_TRA TRA
                if breakpoint1 != '' and breakpoint2 != '' and maq == 60 and (
                        chromosome_list.index(chrom1) < chromosome_list.index(chrom2)):
                    svid = chrom2 + ':' + str(breakpoint2) + '_' + chrom1 + ':' + str(breakpoint1)
                    svsignal_supp += [
                        chrom1 + '\t' + primary_map[0] + '\t' + str(breakpoint1) + '\t' + str(breakpoint2) + '\t0\t' + str(
                            maq) + "\t" + svid + '\tTRA' + "\t*"]
                elif breakpoint1 != '' and breakpoint2 != '' and maq == 60 and (
                        chromosome_list.index(chrom1) > chromosome_list.index(chrom2)):
                    svid = chrom1 + ':' + str(breakpoint1) + '_' + chrom2 + ':' + str(breakpoint2)
                    svsignal_supp += [
                        chrom2 + '\t' + primary_map[0] + '\t' + str(breakpoint2) + '\t' + str(breakpoint1) + '\t0\t' + str(
                            maq) + "\t" + svid + '\tTRA' + "\t*"]

    return svsignal_supp

```

### ./PSV_Signal/parse_sam_signal/0.Signal4bam_PSVGT.py
```python
import argparse
from time import time
import pandas as pd
from multiprocessing import Pool
import sub_Signal4bam_PSVGT
def parse_arguments():
    parser = argparse.ArgumentParser("SV signal extract from sam file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Signal Capture")
    IN.add_argument("-s", dest="sam", required=True, help="sam file ")
    IN.add_argument("-o", dest="out", required=True, help="Output SV signal info file chromsomely")
    IN.add_argument("-m", dest="min", help="SV min length", default=40, type=int)
    IN.add_argument("-maq", dest="maq", help="The min mapping quality", default=50, type=int)
    IN.add_argument("-dtype", dest="dtype", required=True, help="The sequence type (ont, hifi, pb, cr, sr)")
    IN.add_argument("-M", dest="max", help="SV max length", default=10000000, type=int)
    IN.add_argument("-fai", dest="fai", help="Chromosome fai index file", type=str)
    IN.add_argument("-msv", dest="msv", help="Detecting complex SVs (INV, DP, TRA, INS, DEL) from supplementary alignment", default="no", type=str)
    IN.add_argument("-w", dest="workers", help="Number of worker processes", default=10, type=int)
    return parser.parse_args()

def main():
    args = parse_arguments()
    # Load chromosome list from the fai file
    fai = pd.read_csv(args.fai, header=None, sep="\t",dtype=str, index_col=None)
    chromosome_list =  fai[0].tolist()
    chrom_size_dict = dict(zip(fai[0], fai[1].astype(int)))
    if args.out[-4:] == ".sam":
        args.out = args.out.replace('.sam', '')
    SVsignal_out_path = f"{args.out}"
    # Use multiprocessing to parallelize the chromosome processing
    with Pool(processes=args.workers) as pool:
        pool.starmap(sub_Signal4bam_PSVGT.process_chromosome, 
                     [(chromosome,chrom_size_dict[chromosome],chromosome_list, args.sam, args.min,args.max, args.maq, SVsignal_out_path,args.dtype,args.msv) 
                      for chromosome in chromosome_list])
    exe_end = time()
    print(f"{'*' * 40} done SV searching {'*' * 40}\ncost time: {exe_end - exe_start}")

if __name__ == "__main__":
    exe_start = time()
    main()

```

### ./PSV_Signal/parse_sam_signal/filter_long_cigar_sam.py
```python
import re
import sys
insam =sys.argv[1]
outsam = sys.argv[2]
MAX_LEN = 265000000
cigar_pattern = re.compile(r'(\d+)([MIDNSHPX=])')

with open(insam, 'r') as fin, open(outsam, 'w') as fout:
    for line in fin:
        if line.startswith('@'):  
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:  
            continue
        cigar = parts[5]
        valid = True
        for length_str, op in cigar_pattern.findall(cigar):
            if int(length_str) > MAX_LEN:
                valid = False
                break
        if valid:
            fout.write(line)

```

## ./PSV_Genotyper
```shell
total 244
-rwxr-xr-x 1 lgb xinwang  4591 May 22 16:55 2.ACC_lrSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang  6717 May 22 16:55 2.Pop_crSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang  7623 May 22 16:55 2.Pop_lrSVGT_Force.py
-rwxr-xr-x 1 lgb xinwang  7789 Aug 30 12:00 2.Pop_lrSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang  8444 Sep 12 15:09 2.Pop_srSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang  7101 May 22 16:55 2.Single_lrSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang  7781 May 22 16:55 back_2.Pop_lrSVGT_V1.py
-rwxr-xr-x 1 lgb xinwang 34837 May 27 17:40 back_sub_lr_SVGT.py
-rwxr-xr-x 1 lgb xinwang  2054 May 22 16:55 merge_vcf_by_pandas.py
-rwxr-xr-x 1 lgb xinwang  1653 May 22 16:55 phased_diploid_asm.py
-rwxr-xr-x 1 lgb xinwang  2368 May 22 16:55 phased_polyploid_genome_gt.py
drwxr-xr-x 2 lgb xinwang  4096 Jul 23 09:40 __pycache__
-rwxr-xr-x 1 lgb xinwang 25725 May 22 16:55 sub_cr_SVGT.py
-rwxr-xr-x 1 lgb xinwang 34484 May 22 16:55 sub_lr_SVGT_force.py
-rwxr-xr-x 1 lgb xinwang 34900 Jul 22 18:14 sub_lr_SVGT.py
-rwxr-xr-x 1 lgb xinwang  1683 May 22 16:55 Sub_readfa2Dict.py
-rwxr-xr-x 1 lgb xinwang 21262 May 22 16:55 sub_srSVGT.py
-rwxr-xr-x 1 lgb xinwang  3664 May 22 16:55 SVGT_tab2vcf.py

```
### ./PSV_Genotyper/phased_diploid_asm.py
```python
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process and combine haplotype data')
    parser.add_argument('hap1path', help='Path to the first haplotype file')
    parser.add_argument('hap2path', help='Path to the second haplotype file')
    parser.add_argument('out_diploid', help='output diploid phased genotype file')
    args = parser.parse_args()

    for path in [args.hap1path, args.hap2path]:
        try:
            with open(path, 'r'):
                pass
        except FileNotFoundError:
            raise SystemExit(f"错误：文件 {path} 不存在")

    hap1 = pd.read_csv(args.hap1path, sep='\t', header=0, index_col=None)
    hap2 = pd.read_csv(args.hap2path, sep='\t', header=0, index_col=None)
    gt_h1 = hap1.columns[5]  
    gt_h2 = hap2.columns[5]
    genotype_map = {
        "1/1": "1",
        "0/0": "0",
        "./.": ".",
        "0/1": "0/1"
    }

    hap1['GT'] = hap1[gt_h1].map(genotype_map).fillna(hap1[gt_h1])  # 未知基因型保留原值
    hap2['GT'] = hap2[gt_h2].map(genotype_map).fillna(hap2[gt_h2])
    # ---------------------- 合并相位基因型 ----------------------
    phased = hap1.copy()
    phased['GT'] = hap1['GT'] + "|" + hap2['GT']
    phase_standard = {
        '0/1|1': '0|1',  
        '0/1|0': '1|0',
        '0|0/1': '0|1',
        '1|0/1': '1|0',
	'0/1|0/1':"0/1"
    }
    phased['GT'] = phased['GT'].replace(phase_standard)
    phased[gt_h1] = phased['GT']
    phased.drop(columns=phased.columns[-1], inplace=True)
    phased.to_csv(args.out_diploid, sep='\t', index=False, header=True)

if __name__ == "__main__":
    main()

```

### ./PSV_Genotyper/2.Pop_lrSVGT_V1.py
```python
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

```

### ./PSV_Genotyper/2.Pop_crSVGT_V1.py
```python
from os.path import basename
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from tqdm import tqdm
from sub_cr_SVGT import delGT, insGT,dupGT,supp2traGT, invGT, breaks2GT
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
        genotype = supp2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq,"TRA",shift=2000)
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size]
        return "\t".join(map(str, out + genotype))

    if "DUP" in info[5]:
        chrome, bp1, bp2 = info[0], int(info[1]),int(info[2])
        bp1_left = max(bp1 - shift,0) ## aviod the region start at position value less than shift
        bp1_sam = opened_sam.fetch(reference=chrome, start= bp1_left, end= bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start= bp2 - shift, end= bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, -sv_size, sv_size]
        genotype = dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP",shift=200)
        return "\t".join(map(str, out + genotype))

    if "INV" in info[5]:
        chrome,bp1,bp2 = info[0], int(info[1]),int(info[2])
        left_most = max(bp1 - shift, 0)  ## aviod the region start at position value less than shift
        bp1_sam = opened_sam.fetch(reference=chrome, start= left_most, end= bp1 + shift)
        bp2_sam = opened_sam.fetch(reference=chrome, start= bp2 - shift, end= bp2 + shift)
        sv_size = int(info[3])
        out = [chrome, bp1, bp2, sv_size, sv_size]
        genotype = invGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "INV",shift=1000)
        return "\t".join(map(str, out + genotype))

def svGenotyper(supp_sv_table, mapf, name, outdir, min_maq, shift):
    header_line = f"#Target_name\tTarget_start\tTarget_end\tTarget_size\tQuery_size\t{name}\tTotal_Map_Reads\tSV_support"
    with open(supp_sv_table, "r") as svf:
        sv_lines = svf.readlines()

    output_lines = []
    with ThreadPoolExecutor(max_workers=25) as executor:
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

```

### ./PSV_Genotyper/merge_vcf_by_pandas.py
```python
import os
import pandas as pd
import gzip

def vcf2df(vcfpath, reset=False):
    is_gzipped = vcfpath.endswith('.gz')
    if is_gzipped:
        with gzip.open(vcfpath, 'rt') as f:
            header = next((i for i, line in enumerate(f) if not line.startswith('##')), 0)
    else:
        with open(vcfpath, 'r') as f:
            header = next((i for i, line in enumerate(f) if not line.startswith('##')), 0)
    vcf = pd.read_csv(vcfpath, header=header, sep="\t", dtype=str, compression='infer' if is_gzipped else None)
    return vcf


def file_capture(dir, suffix):
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            captures.append(os.path.join(dir, file))
    return captures

def merge(args):
    file_lists = file_capture(args.dir, '.vcf')
    print(file_lists)
    file_lists.sort()
    vcf = vcf2df(file_lists[0])
    vcf.sort_values(by=["#CHROM","POS"],inplace=True)
    vcf.index = vcf["ID"]
    vcf_gt = vcf[vcf.columns[-1]]
    
    for file_name in file_lists[1:]:
        vcfi = vcf2df(file_name)
        vcfi.sort_values(by=["#CHROM","POS"],inplace=True)
        vcfi.index = vcfi["ID"]
        vcfi_gt = vcfi[vcfi.columns[-1]]
        vcf_gt = pd.concat([vcf_gt, vcfi_gt],axis=1)
    vcf["INFO"] = "."
    vcf_out = pd.concat([vcf[vcf.columns[0:9]],vcf_gt], axis=1)
    print(vcf_out.head())
    vcf_out.to_csv(f"{args.out}",header=True,index=None,sep="\t")
    vcf_out[vcf_out["ALT"].isin(["<DEL>", "<INS>"])].to_csv(f"{args.out}.SVInDel",header=True,index=None,sep="\t")
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("merging the population scale SVGT vcf", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-d", dest="dir",required=True, help="a file list all the vcf file path which generated by SVGT")
    IN.add_argument("-o", dest="out",required=True, help="output file merged vcf")
    args  = parser.parse_args()
    merge(args)


```

### ./PSV_Genotyper/2.Pop_lrSVGT_Force.py
```python
from os.path import basename
import subprocess
from multiprocessing import Pool
import os
from tqdm import tqdm
from sub_lr_SVGT_force import delGT, insGT, dupGT, supp_dupGT, traGT,breaks2traGT, invGT, breaks2invGT, little_dupGT
import pandas as pd
import pysam
from math import ceil
import time

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def process_sv(sv_line, opened_sam, name, min_maq, homo_rate, ref_rate, shift=100):
    sampleID = name
    info = sv_line.strip().split("\t")
    if "DEL" in info[5]:
        chrome, sv_s, sv_e = info[0], int(info[1]), int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift, 0)
        left_sam = opened_sam.fetch(reference=chrome, start=left_most, end=sv_s + shift)
        right_sam = opened_sam.fetch(reference=chrome, start=sv_e - shift, end=sv_e + shift)
        genotype = delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift=200)
        out = [chrome, sv_s, sv_e, sv_size, -sv_size]
        return "\t".join(map(str, out + genotype))

    if "INS" in info[5]:
        chrome = info[0]
        sv_s = int(info[1])
        sv_e = int(info[2])
        sv_size = int(info[3])
        left_most = max(sv_s - shift, 0)
        region_sam = opened_sam.fetch(reference=chrome, start=left_most, end=sv_e + shift)
        genotype = insGT(sampleID, region_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift=200)
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
            genotype = little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift=shift)
        else:
            genotype = dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift=shift)
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
        chr2_s2, chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1, bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2, bp2 = chr2_s2[0], int(chr2_s2[1])
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
    chromosome_lines, mapf, name, min_maq, homo_rate, ref_rate, shift = args
    opened_sam = pysam.AlignmentFile(mapf, "rb")
    chromosome_output = []
    for line in chromosome_lines:
        result = process_sv(line, opened_sam, name, min_maq, homo_rate, ref_rate, shift)
        if result:
            chromosome_output.append(result)
    opened_sam.close()
    return chromosome_output


def svGenotyper(supp_sv_table, mapf, name, outdir, min_maq, homo_rate, ref_rate, shift, workers):
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
        args_list = [(lines, mapf, name, min_maq, homo_rate, ref_rate, shift) for _, lines in chromosome_groups.items()]
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
    args = parser.parse_args()
    start_t = time()
    svGenotyper(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.lr_homo_rate, args.lr_ref_rate, args.shift, args.workers)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")

```

### ./PSV_Genotyper/2.ACC_lrSVGT_V1.py
```python
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

```

### ./PSV_Genotyper/sub_cr_SVGT.py
```python
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


def sam_primary_parser2Breaks_Del(region_sam, min_maq):
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
                if  length >= 40 :  # Only count size equally 
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 0
                    deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    return breakpoints, deletions, total_map_reads

def sam_primary_parser2Breaks_Ins(region_sam, min_maq):
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
                if  length > 40 :  # Only count size equally ## to get close del points
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
                if   0.6 * sv_size < length < 1.5 * sv_size  :  # Only count size equally 
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
    if entry_ratio > 0.65:
        return "1/1"
    elif entry_ratio <  0.15:
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
        if max(break_l_ratio, break_r_ratio)   > 0.6:  
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
    breakpoints, inserts, total_map_reads = sam_primary_parser2Breaks_Ins(region_sam, min_maq)
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
    breakpoints_l, deles_l, total_map_reads_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq)
    breakpoints_r, deles_r, total_map_reads_r = sam_primary_parser2Breaks_Del(right_sam, min_maq)
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
    if bp1_map + bp2_map ==0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.2:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map}')
    info_return.append(f"{chrome1}:{bp1},{chrome2}:{bp2};segment_captured;{sv_type}")
    return info_return

def supp2tra(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=2000):
    bp1_map, bp2_map = 0, 0
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
                    if ref_chr != chrom:
                        clipinfo = parse_cigar2clipinfo(cigars)
                        end = start + clipinfo[1]
                        if (bp2-shift) <= start <= (bp2 + shift) or (bp2-shift) <= end <= (bp2 + shift) :
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
                    if ref_chr != chrom:
                        clipinfo = parse_cigar2clipinfo(cigars)
                        end = start + clipinfo[1]
                        if (bp1 - shift) <= start <= (bp1 + shift) or (bp1 - shift) <= end <= (bp1 + shift):
                            genotype += 1
    if bp2_map != 0:
        bp2_rate = genotype / bp2_map
    else:
        bp2_rate = 0
    return bp1_rate+bp2_rate, bp1_map, bp2_map

def supp2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=2000):
    """"
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    use in hifi, cr, genome
    """
    rate,bp1_map, bp2_map =  supp2tra(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=2000)
    if bp1_map + bp2_map == 0:
        genotype = './.'
    info_return = []
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.2:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map}')
    info_return.append(f"{chrome1}:{bp1},{chrome2}:{bp2};segment_captured;{sv_type}")
    return info_return

```

### ./PSV_Genotyper/phased_polyploid_genome_gt.py
```python
import pandas as pd
import argparse

def convert_gt(gt_str):
    """
    polyploid gt to stardard gt
    """
    alleles = set(gt_str.split('|'))
    if '0' in alleles and '1' in alleles:
        return "0/1"
    elif '1' in alleles and '.' not in alleles:
        return "1/1"
    elif '1' in alleles and '.' in alleles:
        return "1/1" ## in hifi reads it will be 1/1
    elif '0' in alleles and '.' in alleles:
        return "./0"
    elif '1' not in alleles and '.' not in alleles:
        return "0/0"
    else:
        return ".".join(alleles)

def main():
    parser = argparse.ArgumentParser(description='Process and combine multiple haplotypes to multi-allele genotype table')
    parser.add_argument('hap_paths', nargs='+', help='Paths to haplotype files (e.g., hap1.tsv hap2.tsv hap3.tsv ...)')
    parser.add_argument('out_multiploid', help='Output file for multi-haplotype phased genotype')
    args = parser.parse_args()
    for path in args.hap_paths:
        try:
            with open(path, 'r'):
                pass
        except FileNotFoundError:
            raise SystemExit(f"file {path} not found")
    haps = []
    for idx, path in enumerate(args.hap_paths, 1):
        df = pd.read_csv(path, sep='\t', header=0, index_col=None)
        df['hap_id'] = f"hap{idx}"
        haps.append(df)
    genotype_map = {
        "1/1": "1",
        "0/0": "0",
        "./.": ".",
        "0/1": "0"  # haplotype genome dont have heterozygous GT
    }
    for df in haps:
        gt_col = df.columns[5]
        df['GT'] = df[gt_col].map(genotype_map).fillna(df[gt_col])
    merged = haps[0].copy()
    gt_columns = [df['GT'] for df in haps]
    merged['GT_polyploid'] = merged.apply(lambda row: "|".join([str(gt.iloc[row.name]) for gt in gt_columns]), axis=1)
    merged['GT_standard'] = merged['GT_polyploid'].apply(convert_gt)
    merged[haps[0].columns[5]] = merged['GT_standard']
    merged.drop(columns=['hap_id', 'GT_polyploid', 'GT_standard'], inplace=False, errors='ignore').to_csv(args.out_multiploid, sep='\t', index=False, header=True)
    print(merged.head())
    merged[haps[0].columns[5]] = merged['GT_polyploid']
    merged.drop(columns=['hap_id', 'GT_polyploid', 'GT_standard'], inplace=False, errors='ignore').to_csv(f"{args.out_multiploid}.polyploid.gt", sep='\t', index=False, header=True)
if __name__ == "__main__":
    main()

```

### ./PSV_Genotyper/Sub_readfa2Dict.py
```python
import gzip
from collections import deque
from gzip import BadGzipFile
import re
def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

```

### ./PSV_Genotyper/sub_srSVGT.py
```python
import re
import subprocess
from collections import defaultdict
from math import ceil
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types
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
        #if (flag & 0x4):
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

def sam_parser2Breaks_Del(region_sam, min_maq, sv_size, breakpoint_pos,region, span_bp):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    effective_spans = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        if (flag & 0x4) or maq < min_maq:
        #if (flag & 0x4):
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
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if sv_size >= 100:
                if (align_start + span_bp < breakpoint_pos) and (align_end - span_bp > breakpoint_pos):
                    effective_spans += 1
            else:## small sv full cover
                if align_start - 10 < region[0] and align_end - 10 > region[1]:
                    effective_spans += 1
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':  # Deletion
                if  sv_size - 5 <= length <= sv_size + 5 :  # Only count size equally 
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    effective_spans -= 1
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 1
                    deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    return breakpoints, deletions, total_map_reads, effective_spans

def sam_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, span_bp):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    effective_span = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        if (flag & 0x4) or maq < min_maq:
        #if (flag & 0x4):
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
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if (align_start + span_bp < sv_s)  and (align_end - span_bp > sv_s):
                effective_span += 1
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I':  # Insertion
                if sv_size - 5 <= length <= sv_size + 5 :  # Only count size equally ## to get close del points
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    insertions[insert_key] += 1  # Count the insertion span
                    effective_span -= 1
    return breakpoints, insertions, total_map_reads, effective_span

def get_spans(chromosome, start, cigar_types, cigar_numbers):
    """Generate spans based on CIGAR information."""
    current_position = start
    spans = {}
    for ctype, cnum in zip(cigar_types, cigar_numbers):
        if ctype == 'M':
            # For matches, calculate the span and update current position
            end_position = current_position + cnum
            span_key = f"{chromosome}:{current_position}-{end_position}"
            spans[span_key] = (current_position, end_position)
            current_position = end_position  # Move current position to end of match
        elif ctype == 'D':
            # For deletions, move the current position forward
            current_position += cnum  # Skip over the deleted region
    return spans

def sam_parser2SVInDel_CutSpan(region_sam, min_maq):
    """Process the SAM file to extract spans and count overlaps with a given range."""
    spans = defaultdict(int)
    for row in  region_sam:
        chrom, start, maq, cigar = row.reference_name, row.reference_start, row.mapping_quality, row.cigarstring
        if maq >= min_maq:
            cigar_numbers, cigar_types = parse_cigar(cigar)
            new_spans = get_spans(chrom, start, cigar_types, cigar_numbers)
            for span_key in new_spans:
                spans[span_key] += 1 
    return spans

def calculate_coverage(cut_span, sv_point):
    cov = 0
    for span in cut_span:
        cov_s = int(span.split(":")[1].split("-")[0])
        cov_e = int(span.split(":")[1].split("-")[1])
        # Check coverage around sv_s
        if sv_point - 75 > cov_s and sv_point + 75 < cov_e:
            cov += 1
    return cov

def determine_genotype(entry_ratio, homo_rate, ref_rate):
    """
    ONT easy lead to 0/1 and FP, here we try modify.
    """
    if entry_ratio  >= homo_rate:
        return "1/1"
    elif entry_ratio <  ref_rate:
        return "0/0"
    else:
        return "0/1"

def breaksCallGT(break_l_ratio, break_r_ratio):
    """
    INV, TRA, DUP, 
    Take two breakpoints for Genotype
    """
    if max(break_l_ratio, break_r_ratio) > 0.75 and min(break_l_ratio, break_r_ratio) > 0.4:
        return "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.2 or min(break_l_ratio, break_r_ratio) < 0.1:
        return "0/0"
    else:
        return "0/1"



def determine_dupGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.25 and break_r_ratio >= 0.25:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.3 and min(break_l_ratio, break_r_ratio) >= 0.1:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype

def determine_invGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.6 and break_r_ratio >= 0.6:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.75  and min(break_l_ratio, break_r_ratio) >= 0.3:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype


def determine_traGT(break_l_ratio, break_r_ratio):
    """
    Take two breakpoints for Genotype
    """
    if break_l_ratio >= 0.6 and break_r_ratio >= 0.6:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio ) >= 0.75  and min(break_l_ratio, break_r_ratio) >= 0.3:
        genotype = "1/1"
    elif max(break_l_ratio, break_r_ratio) < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    return genotype


def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift=100, span_bp=50):
    info_return = []
    genotype = "0/0"  # Default genotype
    #sv_start_shift = set(range(sv_s - shift, sv_s + shift+100))
    #sv_end_shift   = set(range(sv_e - shift, sv_e + shift+100))
    sv_start_shift = set(range(sv_s - shift, sv_s + 120))
    sv_end_shift   = set(range(sv_e - shift, sv_e + 120))
    sv_size = abs(sv_size)
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, effective_spans = sam_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, span_bp)
    if total_map_reads == 0:
        info_return.append('./.')
        info_return.append(f"total_map_reads={total_map_reads}")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        cut_span = sam_parser2SVInDel_CutSpan(region_sam, min_maq)
        count_break_and_Ins = 0
        covIns = calculate_coverage(cut_span, sv_s)
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
        covIns_ratio = round(covIns / total_map_reads, 3)
        
        genotype = determine_genotype(ins_ratio, homo_rate, ref_rate)
        if genotype == "0/1" and effective_spans <= 0.01*total_map_reads+1:
            genotype = "1/1"
        elif genotype == "1/1" and effective_spans >= 0.05*total_map_reads+1:
            genotype = "0/1"
        #elif genotype == "0/0":
        #    if effective_spans == 0:
        #        genotype = "1/1"
            #else:
            #    genotype = "0/1"
        print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\tIns_points_covered_ratio:{covIns_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads};effective_spans={effective_spans}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return
def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift=100, span_bp=50):
    ############ SVDel Case ##############
    info_return = []
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, effective_spans_l = sam_parser2Breaks_Del(left_sam, min_maq, sv_size,  sv_s, [sv_s,sv_e], span_bp)
    breakpoints_r, deles_r, total_map_reads_r, effective_spans_r = sam_parser2Breaks_Del(right_sam, min_maq, sv_size, sv_e, [sv_s,sv_e], span_bp)
    cut_span_l = sam_parser2SVInDel_CutSpan(left_sam, min_maq)
    cut_span_r = sam_parser2SVInDel_CutSpan(right_sam, min_maq)
    covDel_l = calculate_coverage(cut_span_l, sv_s)
    covDel_r = calculate_coverage(cut_span_l, sv_e)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("1/1") ## no reads properly dele
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0")
        info_return.append(f"deles_l_ratio=0,deles_r_ratio=0;DEL")
        return info_return
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
        covDel_l_ratio =round( covDel_l / total_map_reads_l, 3)
    else:
        deles_l_ratio = 0
        covDel_l_ratio = 0
    if total_map_reads_r:
        deles_r_ratio = round(count_break_and_deles_r / total_map_reads_r, 3)
        covDel_r_ratio = round( covDel_r / total_map_reads_r, 3)
    else:
        deles_r_ratio = 0
        covDel_r_ratio = 0
    deles_ratio = max(deles_l_ratio, deles_r_ratio)
    genotype = determine_genotype(deles_ratio, homo_rate, ref_rate)
    if genotype == "0/1":
        if effective_spans_l + effective_spans_r == 0:
            genotype = "1/1"
    elif genotype == "1/1":
        if (effective_spans_l + effective_spans_r) >= (0.3 * total_map_reads + 2):
            genotype = "0/1"
    elif genotype == "0/0":
        if effective_spans_l + effective_spans_r == 0:
            genotype = '1/1'
    print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\tdele_l_covered_ratio:{covDel_l_ratio}\tdele_r_covered_ratio:{covDel_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l},span_l={effective_spans_l};total_map_reads_r={total_map_reads_r},span_r={effective_spans_r}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2GT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=50):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ INV Case ##############
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2 = sam_parser2Breaks(bp2_sam, min_maq)
    cut_span_bp1 = sam_parser2SVInDel_CutSpan(bp1_sam, min_maq)
    cut_span_bp2 = sam_parser2SVInDel_CutSpan(bp2_sam, min_maq)
    cov_break_bp1 = calculate_coverage(cut_span_bp1, bp1)
    cov_break_bp2 = calculate_coverage(cut_span_bp2, bp2)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.") 
        info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")

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
        cov_break1_ratio =round( cov_break_bp1 / total_map_reads_bp1, 3)
    else:
        break1_ratio = 0
        cov_break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 3)
        cov_break2_ratio = round( cov_break_bp2 / total_map_reads_bp2, 3)
    else:
        break2_ratio = 0
        cov_break2_ratio = 0
    max_break_ratio = max(break1_ratio, break2_ratio)
    if sv_type == "DUP":
        genotype = determine_dupGT(break1_ratio, break2_ratio)
    elif sv_type == "INV":
        genotype = determine_invGT(break1_ratio, break2_ratio)
    elif sv_type == "TRA":
        genotype = determine_traGT(break1_ratio, break2_ratio)
    print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}_covered_ratio={cov_break1_ratio}\t{breakpoint2}_covered_ratio={cov_break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
    return info_return

def sam2readsID(region_sam):
    readsID = []
    maqs = 0
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        readID = row.query_name
        readsID.append(readID)
    if not readsID:
        return [], 0
    else:
        return readsID, ceil(maqs / len(readsID))

def traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=50):
    genotype = "0/0"
    info_return = []
    breaks_dict = {}
    bp1_readsID, maq1 = sam2readsID(bp1_sam)
    bp2_readsID, maq2 = sam2readsID(bp2_sam)

    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    if max(bp1_tra, bp2_tra) > 0.90:
        genotype = "1/1"
    elif bp1_tra >= 0.8 and bp2_tra >= 0.8:
        genotype = "1/1"
    elif bp1_tra + bp2_tra < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    info_return.append(genotype)
    print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return


```

### ./PSV_Genotyper/back_sub_lr_SVGT.py
```python
import re
import subprocess
from collections import defaultdict
from math import ceil,floor
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types

def sam2readsID(region_sam):
    readsID = []
    maqs = 0
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        readID = row.query_name
        readsID.append(readID)
    if not readsID:
        return [], 0
    else:
        return readsID, ceil(maqs / len(readsID))

def sam_parser2Breaks(region_sam, min_maq):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq = row.mapping_quality
        maqs  += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (flag & 0x4) or maq < min_maq:
        if flag & 0x4:
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
            update_breakpoints(chrom, breakpoints_to_update)
    if total_map_reads >0:
        #print( breakpoints, total_map_reads, ceil(maqs / total_map_reads))
        return breakpoints, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, 0, 0


def sv_egde(current_svlen):
    if current_svlen < 100:
        window_size = 150 + current_svlen * 0.2
    elif 100 < current_svlen <= 500:
        window_size = 200 + current_svlen * 0.2
    elif 500 < current_svlen <= 1000:
        window_size = 300 + current_svlen * 0.2
    else:
        window_size = 500
    return window_size + 50



def sam_primary_parser2Breaks_Del(region_sam, min_maq, sv_size, sv_start, sv_end, minfilt):
    #sv_shift = sv_egde(sv_size)
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        del_size = 0
        readname = row.query_name
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):

            continue
        total_map_reads += 1
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start 
        align_end = row.reference_end
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS' and align_start > sv_end - 0.5*sv_size: ## right bp
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS' and align_end < sv_start + 0.5*sv_size: ## left breakpoint
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                if sv_start - 200 < current_start < sv_end + 200:
                    if sv_size*0.8 <= length <= sv_size*1.2:
                        deletions[readname] = length
                    if length >= minfilt:
                        del_size += length
                        if sv_size*0.8 <= del_size <= sv_size*1.2:
                            deletions[readname] = del_size
                current_start += length
    
        if 0.7*sv_size <= del_size <= 1.5*sv_size:
            deletions[readname] = del_size
    
    if total_map_reads >0:
        #print(f'{chrom}:{sv_start}\t{sv_size}:{deletions}')
        return breakpoints, deletions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_start, sv_end, minfilt):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    reads2ins_size = defaultdict(int) 
    total_map_reads = 0
    spans = 0
    maqs = 0
    effective_span = 0
    #sv_shift = sv_egde(sv_size)
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        ins_size = 0
        readname = row.query_name
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if (align_start + 1000 < sv_start)  and (align_end - 1000 > sv_end):
                effective_span += 1

        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I' and (sv_start-200 < current_start <= sv_end + 200) and length >=minfilt: ## since window size in cluster is 500
                if 0.8*sv_size <= length <= 1.2*sv_size:
                    reads2ins_size[readname] = length
                else:
                    ins_size += length
                    if 0.8*sv_size <= length <= 1.2*sv_size:
                        reads2ins_size[readname] = ins_size
        if 0.7*sv_size <= ins_size < 1.5*sv_size:
            reads2ins_size[readname] = ins_size
        
    effective_span -= len(reads2ins_size)
    if total_map_reads >0:
        #print(f'{chrom}:{sv_start}-{sv_end}\t{sv_size}:{reads2ins_size}')
        return breakpoints, reads2ins_size, total_map_reads, ceil(maqs / total_map_reads), effective_span
    else:
        return {},{},0,0,effective_span
def sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size):
    ## some dup signal may hide in I cigar ##
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
            continue
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start
        align_end = row.reference_end
        total_map_reads += 1
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  
            elif ctype == 'I':  # Insertion
                if  0.7 * sv_size < length < 1.5 * sv_size  :  # Only count size equally 
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    else:
                        insertions[insert_key] += 1  # Count the insertion span
    if total_map_reads > 0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def determine_genotype(breaks, depth, homo_rate=0.75,ref_rate = 0.05):
    """
    depth and ratio base genotype
    """
    if depth == 0:
        return "./."
    if depth <= 5:
        if floor(homo_rate*depth) + 1 <= breaks:
            return "1/1"
        elif breaks / depth <= 0.2: ## 2 reads --> 5X
            return "0/0"
        else:
            return "0/1"

    if 5 < depth <= 10:
        if homo_rate * depth + 1 <=  breaks:
            return "1/1"
        elif 0.05 * depth + 1 > breaks:
            return "0/0"
        else:
            return "0/1"
    
    if depth > 10:
        if breaks / depth >= homo_rate: ## or 0.65   or 0.625
            return "1/1"
        elif breaks / depth < ref_rate:
            return "0/0"
        else:
            return "0/1"

def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift, minfilt):
    info_return = []
    genotype = "0/0"  # Default genotype
    shift = min(sv_size, 500)
    sv_start_shift = set(range(sv_s - shift, sv_s + shift ))
    sv_end_shift=set(range(sv_e-shift, sv_e+shift))
    dup_shift = set(range(sv_s - 1000, sv_s - 1000)) ## samll dup capture as ins
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, maq, effective_span = sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, minfilt)
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads},maq=0")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        count_break_and_Ins = len(inserts)
        if breakpoints:   #check breakpoints
            for breakpoint in breakpoints.get(chrome, {}).keys():
                if breakpoint in sv_start_shift or breakpoint in sv_end_shift: 
                    count_break_and_Ins += breakpoints[chrome][breakpoint]
        ins_ratio = round(count_break_and_Ins / total_map_reads, 3)
        genotype = determine_genotype(count_break_and_Ins, total_map_reads, homo_rate, ref_rate)
        if genotype == "1/1":
            if floor(0.1*total_map_reads) + 1 <=  effective_span:
                genotype = "0/1"
                #print(f"***************Correting SVINS {chrome}:{sv_s}-{sv_e} genotype to 0/1 since it has {effective_span} span reads*****************")
        #print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads},maq={maq}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return

def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift,minfilt):
    ############ SVDel Case ##############
    info_return = []
    breaks_dict = {}
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, maq_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq, sv_size, sv_s, sv_e, minfilt)
    breakpoints_r, deles_r, total_map_reads_r, maq_r = sam_primary_parser2Breaks_Del(right_sam, min_maq, sv_size, sv_s, sv_e, minfilt)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0,maq=0")
        info_return.append(f"deles_l_ratio=0,deles_r_ratio=0;DEL")
        return  info_return
    
    count_break_and_deles_l += len(deles_l)
    count_break_and_deles_r += len(deles_r)
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
    breaks_dict[count_break_and_deles_l ] = total_map_reads_l
    breaks_dict[count_break_and_deles_r ] = total_map_reads_r
    max_breaks = max(count_break_and_deles_l,count_break_and_deles_r)
    genotype = determine_genotype(max_breaks,breaks_dict[max_breaks], homo_rate, ref_rate)
    #print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l};total_map_reads_r={total_map_reads_r};maq={max(maq_l,maq_r)}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
        return info_return
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.8:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return

def little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, sv_type, shift=500):
    """
    Taking total local mapping  to genotyping the small(<7k) dup SV;
    1st Scan all insertion cigar; 2nd Capture all breakpoints;
    """
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome}:{bp1}", f"{chrome}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1  + shift))
    bp2_shift = set(range(bp2 - shift, bp2  + shift))
    ins_shift1 = set(range(bp1 - sv_size, bp1 + sv_size))
    ins_shift2 = set(range(bp2 - sv_size, bp2 + sv_size))

    breakpoints_bp, inserts, total_map_reads, maq = sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size)

    count_break_bp = 0
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts:
        for pos in inserts.keys():
            ins_s, _ = map(int, pos.split(":")[1].split("-"))
            if ins_s in ins_shift1 or ins_s in ins_shift2:
                count_break_bp += inserts[pos]
    if breakpoints_bp:
        for breakpoint in breakpoints_bp.get(chrome, {}).keys():
            if breakpoint in bp1_shift or breakpoint in bp2_shift:
                count_break_bp += breakpoints_bp[chrome][breakpoint]
    if total_map_reads:
        break_ratio =   round(count_break_bp / total_map_reads, 3)
    else:
        break_ratio = 0
    
    if break_ratio > 0.75:
        genotype = "1/1"
    elif break_ratio < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads={total_map_reads};\tbreakpoint_ratio={break_ratio}\tbp1={breakpoint1},bp2={breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads={total_map_reads};maq={maq}")
    info_return.append(f"bp1={breakpoint1},bp2={breakpoint2},bp_ratio={break_ratio};{sv_type}")
    return info_return


def dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints local mapping info to genotyping the big SV;
    
    """
    breaks_dict = {}
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift = set(range(bp2 - shift, bp2 + shift))
    ins_shift1 = set(range(bp1 - shift, bp1 + ceil(shift+ sv_size/2)))
    ins_shift2 = set(range(bp2 - shift - ceil(sv_size / 2), bp2 + shift))

    breakpoints_bp1, inserts_bp1, total_map_reads_bp1, maq1 = sam_primary_parser2Breaks_dup(bp1_sam, min_maq, sv_size)
    breakpoints_bp2, inserts_bp2, total_map_reads_bp2, maq2 = sam_primary_parser2Breaks_dup(bp2_sam, min_maq, sv_size)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts_bp1:  # Check if there are insertion entries
        for pos in inserts_bp1.keys():
            if pos in ins_shift1:
                count_break_bp1 += inserts_bp1[pos]
    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]

    if inserts_bp2:  # Check if there are insertion entries
        for pos in inserts_bp2.keys():
            if pos in ins_shift2:
                count_break_bp2 += inserts_bp2[pos]
    if breakpoints_bp2:
        for breakpoint in breakpoints_bp2.get(chrome2, {}).keys():
            if breakpoint in bp2_shift:
                count_break_bp2 += breakpoints_bp2[chrome2][breakpoint]

    if total_map_reads_bp1:
        break1_ratio =   round(count_break_bp1 / total_map_reads_bp1, 2)
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 2)
    else:
        break2_ratio = 0

    max_break_ratio = max(break1_ratio, break2_ratio)
    breaks_dict[count_break_bp1] = total_map_reads_bp1
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_breaks = max(count_break_bp1, count_break_bp2)
    if break1_ratio >=0.4 and break2_ratio >= 0.4:
        genotype = "1/1"
    elif break1_ratio + break2_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio + break2_ratio < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
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
    genotype = "0/0"
    info_return = []
    breaks_dict = {}
    bp1_readsID, maq1 = sam2readsID(bp1_sam)
    bp2_readsID, maq2 = sam2readsID(bp2_sam)

    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    if max(bp1_tra, bp2_tra) > 0.90:
        genotype = "1/1"
    elif bp1_tra >= 0.8 and bp2_tra >= 0.8:
        genotype = "1/1"
    elif bp1_tra + bp2_tra < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    info_return.append(genotype)
    #print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return

def breaks2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
        return info_return
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio >= 0.85 and break2_ratio >= 0.85:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.10:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"***breakpoint to genotype TRA********* {sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2} ****************")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return


def supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate, maq_bp1 = 0,0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
        maq_bp2 = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    rate, bp1_map, bp2_map, maq = supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    info_return = []
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.05:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

def supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0

    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    info_return = []
    rate, bp1_map, bp2_map, maq = supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.625:
        genotype = '1/1'
    elif rate < 0.1:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

```

### ./PSV_Genotyper/back_2.Pop_lrSVGT_V1.py
```python
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
            genotype = little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift=shift)
        else:
            genotype = dupGT(sampleID, bp1_sam, bp2_sam, chrome, chrome, bp1, bp2, sv_size, min_maq, "DUP", shift=shift)
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
        chr2_s2, chr1_s1 = svid.split("_")[0].split(':'), svid.split("_")[1].split(":")
        chrome1, bp1 = chr1_s1[0], int(chr1_s1[1])
        chrome2, bp2 = chr2_s2[0], int(chr2_s2[1])
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

```

### ./PSV_Genotyper/sub_lr_SVGT.py
```python
import re
import subprocess
from collections import defaultdict
from math import ceil,floor
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types

def sam2readsID(region_sam):
    readsID = []
    maqs = 0
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        readID = row.query_name
        readsID.append(readID)
    if not readsID:
        return [], 0
    else:
        return readsID, ceil(maqs / len(readsID))

def sam_parser2Breaks(region_sam, min_maq):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq = row.mapping_quality
        maqs  += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (flag & 0x4) or maq < min_maq:
        if flag & 0x4:
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
            update_breakpoints(chrom, breakpoints_to_update)
    if total_map_reads >0:
        #print( breakpoints, total_map_reads, ceil(maqs / total_map_reads))
        return breakpoints, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, 0, 0


def sv_egde(current_svlen):
    if current_svlen < 100:
        window_size = 150 + current_svlen * 0.2
    elif 100 < current_svlen <= 500:
        window_size = 200 + current_svlen * 0.2
    elif 500 < current_svlen <= 1000:
        window_size = 300 + current_svlen * 0.2
    else:
        window_size = 500
    return window_size + 50



def sam_primary_parser2Breaks_Del(region_sam, min_maq, sv_size, sv_start, sv_end, minfilt):
    #sv_shift = sv_egde(sv_size)
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        del_size = 0
        readname = row.query_name
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):

            continue
        total_map_reads += 1
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start 
        align_end = row.reference_end
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS' and align_start > sv_end - 0.5*sv_size: ## right bp
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS' and align_end < sv_start + 0.5*sv_size: ## left breakpoint
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                if sv_start - 200 < current_start < sv_end + 200:
                    if sv_size*0.8 <= length <= sv_size*1.2:
                        deletions[readname] = length
                    if length >= minfilt:
                        del_size += length
                        if sv_size*0.8 <= del_size <= sv_size*1.2:
                            deletions[readname] = del_size
                current_start += length
    
        if 0.7*sv_size <= del_size <= 1.5*sv_size:
            deletions[readname] = del_size
    
    if total_map_reads >0:
        #print(f'{chrom}:{sv_start}\t{sv_size}:{deletions}')
        return breakpoints, deletions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_start, sv_end, minfilt):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    reads2ins_size = defaultdict(int) 
    total_map_reads = 0
    spans = 0
    maqs = 0
    effective_span = 0
    #sv_shift = sv_egde(sv_size)
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        ins_size = 0
        readname = row.query_name
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
            update_breakpoints(chrom, breakpoints_to_update)
        else:
            if (align_start + 1000 < sv_start)  and (align_end - 1000 > sv_end):
                effective_span += 1

        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I' and (sv_start - 300 < current_start <= sv_end + sv_size) and length >=minfilt: ## since window size in cluster is 500
                if 0.8*sv_size <= length <= 1.2*sv_size:
                    reads2ins_size[readname] = length
                else:
                    ins_size += length
                    if 0.8*sv_size <= length <= 1.2*sv_size:
                        reads2ins_size[readname] = ins_size
        if 0.7*sv_size <= ins_size < 1.5*sv_size:
            reads2ins_size[readname] = ins_size
        
    effective_span -= len(reads2ins_size)
    if total_map_reads >0:
        #print(f'{chrom}:{sv_start}-{sv_end}\t{sv_size}:{reads2ins_size}')
        return breakpoints, reads2ins_size, total_map_reads, ceil(maqs / total_map_reads), effective_span
    else:
        return {},{},0,0,effective_span
def sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size):
    ## some dup signal may hide in I cigar ##
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
            continue
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start
        align_end = row.reference_end
        total_map_reads += 1
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  
            elif ctype == 'I':  # Insertion
                if  0.7 * sv_size < length < 1.5 * sv_size  :  # Only count size equally 
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    else:
                        insertions[insert_key] += 1  # Count the insertion span
    if total_map_reads > 0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def determine_genotype(breaks, depth, homo_rate=0.75,ref_rate = 0.05):
    """
    depth and ratio base genotype
    """
    if depth == 0:
        return "./."
    if depth <= 5:
        if floor(homo_rate*depth) + 1 <= breaks:
            return "1/1"
        elif breaks / depth <= 0.2: ## 2 reads --> 5X
            return "0/0"
        else:
            return "0/1"

    if 5 < depth <= 10:
        if homo_rate * depth + 1 <=  breaks:
            return "1/1"
        elif 0.05 * depth + 1 > breaks:
            return "0/0"
        else:
            return "0/1"
    
    if depth > 10:
        if breaks / depth >= homo_rate: ## or 0.65   or 0.625
            return "1/1"
        elif breaks / depth < ref_rate:
            return "0/0"
        else:
            return "0/1"

def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift, minfilt):
    info_return = []
    genotype = "0/0"  # Default genotype
    shift = min(sv_size, 500)
    sv_start_shift = set(range(sv_s - shift, sv_s + shift ))
    sv_end_shift=set(range(sv_e-shift, sv_e+shift))
    dup_shift = set(range(sv_s - sv_size, sv_s + sv_size)) ## samll dup capture as ins
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, maq, effective_span = sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s, sv_e, minfilt)
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads},maq=0")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        count_break_and_Ins = len(inserts)
        if breakpoints:   #check breakpoints
            for breakpoint in breakpoints.get(chrome, {}).keys():
                if breakpoint in sv_start_shift or breakpoint in sv_end_shift: 
                    count_break_and_Ins += breakpoints[chrome][breakpoint]
        ins_ratio = round(count_break_and_Ins / total_map_reads, 3)
        genotype = determine_genotype(count_break_and_Ins, total_map_reads, homo_rate, ref_rate)
        if genotype == "1/1":
            if floor(0.1*total_map_reads) + 1 <=  effective_span:
                genotype = "0/1"
                #print(f"***************Correting SVINS {chrome}:{sv_s}-{sv_e} genotype to 0/1 since it has {effective_span} span reads*****************")
        #print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads},maq={maq}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return

def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift,minfilt):
    ############ SVDel Case ##############
    info_return = []
    breaks_dict = {}
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, maq_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq, sv_size, sv_s, sv_e, minfilt)
    breakpoints_r, deles_r, total_map_reads_r, maq_r = sam_primary_parser2Breaks_Del(right_sam, min_maq, sv_size, sv_s, sv_e, minfilt)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0,maq=0")
        info_return.append(f"deles_l_ratio=0,deles_r_ratio=0;DEL")
        return  info_return
    
    count_break_and_deles_l += len(deles_l)
    count_break_and_deles_r += len(deles_r)
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
    breaks_dict[count_break_and_deles_l ] = total_map_reads_l
    breaks_dict[count_break_and_deles_r ] = total_map_reads_r
    max_breaks = max(count_break_and_deles_l,count_break_and_deles_r)
    genotype = determine_genotype(max_breaks,breaks_dict[max_breaks], homo_rate, ref_rate)
    #print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l};total_map_reads_r={total_map_reads_r};maq={max(maq_l,maq_r)}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
        return info_return
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.8:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return

def little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, sv_type, shift=500):
    """
    Taking total local mapping  to genotyping the small(<7k) dup SV;
    1st Scan all insertion cigar; 2nd Capture all breakpoints;
    """
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome}:{bp1}", f"{chrome}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1  + shift))
    bp2_shift = set(range(bp2 - shift, bp2  + shift))
    ins_shift1 = set(range(bp1 - sv_size, bp1 + sv_size))
    ins_shift2 = set(range(bp2 - sv_size, bp2 + sv_size))


    breakpoints_bp, inserts, total_map_reads, maq = sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size)

    count_break_bp = 0
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts:
        for pos in inserts.keys():
            ins_s, _ = map(int, pos.split(":")[1].split("-"))
            if ins_s in ins_shift1 or ins_s in ins_shift2:
                count_break_bp += inserts[pos]
    if breakpoints_bp:
        for breakpoint in breakpoints_bp.get(chrome, {}).keys():
            if breakpoint in bp1_shift or breakpoint in bp2_shift:
                count_break_bp += breakpoints_bp[chrome][breakpoint]
    if total_map_reads:
        break_ratio =   round(count_break_bp / total_map_reads, 3)
    else:
        break_ratio = 0
    
    if break_ratio > 0.75:
        genotype = "1/1"
    elif break_ratio < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads={total_map_reads};\tbreakpoint_ratio={break_ratio}\tbp1={breakpoint1},bp2={breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads={total_map_reads};maq={maq}")
    info_return.append(f"bp1={breakpoint1},bp2={breakpoint2},bp_ratio={break_ratio};{sv_type}")
    return info_return


def dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints local mapping info to genotyping the big SV;
    
    """
    breaks_dict = {}
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift = set(range(bp2 - shift, bp2 + shift))
    ins_shift1 = set(range(bp1 - shift, bp1 + ceil(shift+ sv_size/2)))
    ins_shift2 = set(range(bp2 - shift - ceil(sv_size / 2), bp2 + shift))
    ins_shift = set(range(bp1 - shift, bp2 + shift))
    breakpoints_bp1, inserts_bp1, total_map_reads_bp1, maq1 = sam_primary_parser2Breaks_dup(bp1_sam, min_maq, sv_size)
    breakpoints_bp2, inserts_bp2, total_map_reads_bp2, maq2 = sam_primary_parser2Breaks_dup(bp2_sam, min_maq, sv_size)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts_bp1:  # Check if there are insertion entries
        for pos in inserts_bp1.keys():
            if pos in ins_shift:
                count_break_bp1 += inserts_bp1[pos]
    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]

    if inserts_bp2:  # Check if there are insertion entries
        for pos in inserts_bp2.keys():
            if pos in ins_shift:
                count_break_bp2 += inserts_bp2[pos]
    if breakpoints_bp2:
        for breakpoint in breakpoints_bp2.get(chrome2, {}).keys():
            if breakpoint in bp2_shift:
                count_break_bp2 += breakpoints_bp2[chrome2][breakpoint]

    if total_map_reads_bp1:
        break1_ratio =   round(count_break_bp1 / total_map_reads_bp1, 2)
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 2)
    else:
        break2_ratio = 0

    max_break_ratio = max(break1_ratio, break2_ratio)
    breaks_dict[count_break_bp1] = total_map_reads_bp1
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_breaks = max(count_break_bp1, count_break_bp2)
    if break1_ratio >=0.4 and break2_ratio >= 0.4:
        genotype = "1/1"
    elif break1_ratio + break2_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio + break2_ratio < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
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
    genotype = "0/0"
    info_return = []
    breaks_dict = {}
    bp1_readsID, maq1 = sam2readsID(bp1_sam)
    bp2_readsID, maq2 = sam2readsID(bp2_sam)

    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    if max(bp1_tra, bp2_tra) > 0.90:
        genotype = "1/1"
    elif bp1_tra >= 0.8 and bp2_tra >= 0.8:
        genotype = "1/1"
    elif bp1_tra + bp2_tra < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    info_return.append(genotype)
    #print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return

def breaks2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
        info_return.append(f"{breakpoint1}_ratio=0,{breakpoint2}_ratio=0;{sv_type}")
        return info_return
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio >= 0.85 and break2_ratio >= 0.85:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.10:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"***breakpoint to genotype TRA********* {sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2} ****************")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return


def supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate, maq_bp1 = 0,0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
        maq_bp2 = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    rate, bp1_map, bp2_map, maq = supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    info_return = []
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.05:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

def supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0

    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    info_return = []
    rate, bp1_map, bp2_map, maq = supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.625:
        genotype = '1/1'
    elif rate < 0.1:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

```

### ./PSV_Genotyper/SVGT_tab2vcf.py
```python
import argparse
from os.path import basename
import pandas as pd
def svindeltab2vcf(tab, outvcf, fai):
    header = f"""##fileformat=VCFv4.2
##source=PSVGT1.0\n"""
    with open(fai, 'r') as fai:
        lines = fai.readlines()
    for line in lines:
        genome_size = line.strip().split("\t")
        fline = f"##contig=<ID={genome_size[0]},length={genome_size[1]}>\n"
        header +=fline
    with open(tab, 'r') as f:
        lines = f.readlines()
    with open(outvcf, 'w') as outf:
        prehead = f"""{header}##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##INFO=<ID=CHR,Number=1,Type=String,Description="Chromosome">
##INFO=<ID=SV_START,Number=1,Type=Integer,Description="Start position of the structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=maq,Number=1,Type=Integer,Description="The mean mapping quality of SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length difference between the REF and ALT alleles">
##INFO=<ID=method,Number=1,Type=String,Description="Method used to call the variant">
##INFO=<ID=total_map_reads,Number=0,Type=Integer,Description="Total number of mapped reads supporting the variant">
##INFO=<ID=INS_rate,Number=1,Type=Float,Description="Insertion rate (e.g., read depth support)">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (e.g., INS means Insertion)">
##INFO=<ID=total_map_reads_l,Number=1,Type=Integer,Description="Total mapping reads at first breakpoints">
##INFO=<ID=total_map_reads_r,Number=1,Type=Integer,Description="Total mapping reads at second breakpoints">
##INFO=<ID=total_map_reads_bp1,Number=1,Type=Integer,Description="Total mapping reads at first breakpoints">
##INFO=<ID=total_map_reads_bp2,Number=1,Type=Integer,Description="Total mapping reads at second breakpoints">
##INFO=<ID=deles_l_ratio,Number=1,Type=Float,Description="The delestion ratio detected at left breakpoints">
##INFO=<ID=deles_r_ratio,Number=1,Type=Float,Description="The delestion ratio detected at right breakpoints">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{basename(tab)}"""
        print(prehead,file=outf)
        svlist = []
        for line in lines:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")
            chrom, refbase, sv_start, sv_end, svlen, gt, maps,rate,sv_type = info[0], "N", info[1], info[2], abs(int(info[3])), info[5], info[6], info[7].split(";")[0], info[7].split(";")[1]
            ID = f'{chrom}:{sv_start}-{sv_end}_{svlen}'
            outline = f"{chrom}\t{sv_start}\t{ID}\tN\t<{sv_type}>\t.\tPASS\tEND={sv_end};CHR={chrom};SV_START={sv_start};SVLEN={svlen};method=PSVGT;{maps};{rate};SVTYPE={sv_type}\tGT\t{gt}"
            svlist.append(outline.split("\t"))
        df = pd.DataFrame(svlist)
        df[1] = df[1].astype(int)
        df.sort_values(by=[0, 1, 2],inplace=True)
    df.to_csv(outvcf,header=None,sep="\t",index=None,mode='a')

def main():
    parser = argparse.ArgumentParser(description="Convert SVInDel tab-delimited file to VCF format")
    parser.add_argument('tab', help='Input tab-delimited file with SVInDel information')
    parser.add_argument('outvcf', help='Output VCF file path')
    parser.add_argument('fai', help='Output genome faidx file path')
    args = parser.parse_args()
    svindeltab2vcf(args.tab, args.outvcf,args.fai)
if __name__ == "__main__":
    main()


```

### ./PSV_Genotyper/2.Pop_srSVGT_V1.py
```python
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


def srSVGT(sv_table, mapf, name, outdir, min_maq, homo_rate, ref_rate, shift, span):
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

    with multiprocessing.Pool() as pool:
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
                    help="InDel record table can be generated by PSVGT step 1")
    IN.add_argument("-maq", dest="maq", default=1, type=int,
                    help="the mini mapping quality of the contigs, range from 0 - 60, default is 45")
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
    args = parser.parse_args()
    start_t = time()
    srSVGT(args.sv_info, args.mapf, args.ACC, args.dir, args.maq, args.homo_rate, args.ref_rate, args.shift, args.span)
    end_t = time()
    print(f"******************* Cost time {end_t - start_t}s *********************")

```

### ./PSV_Genotyper/2.Single_lrSVGT_V1.py
```python
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
    if entry_ratio > 0.65:
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
            out = [chrome, bp1, bp2, -sv_size, sv_size, dup_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'bp1={chrome}:{bp1},bp1_ratio={bp1_rate},bp2={chrome}:{bp2},bp2_ratio={bp2_rate};DUP']
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
            out = [chrome, bp1, bp2, -sv_size, sv_size, inv_gt, f'total_map_reads_bp1={left_cov};total_map_reads_bp2={right_cov}',f'bp1={chrome}:{bp1},bp1_ratio={bp1_rate},bp2={chrome}:{bp2},bp2_ratio={bp2_rate};INV']
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
        out = [chrome1, bp1, f"{chrome2}:{bp2}", sv_size, sv_size, tra_gt, f'total_map_reads_bp1={bp1_cov};total_map_reads_bp2={bp2_cov}',f'bp2={chrome2}:{bp2},bp2_ratio={tra_rate2},bp1={chrome1}:{bp2},bp1_ratio={tra_rate2};TRA']
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

```

### ./PSV_Genotyper/sub_lr_SVGT_force.py
```python
import re
import subprocess
from collections import defaultdict
from math import ceil,floor
def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigar):
    """Parses a CIGAR string into its component numbers and types."""
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    types = re.findall(r'[MIDNSHP=X]', cigar)
    return numbers, types

def sam2readsID(region_sam):
    readsID = []
    maqs = 0
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        readID = row.query_name
        readsID.append(readID)
    if not readsID:
        return [], 0
    else:
        return readsID, ceil(maqs / len(readsID))

def sam_parser2Breaks(region_sam, min_maq):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq = row.mapping_quality
        maqs  += maq
        align_start = int(row.reference_start)
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (flag & 0x4) or maq < min_maq:
        if flag & 0x4:
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
            update_breakpoints(chrom, breakpoints_to_update)
    if total_map_reads >0:
        #print( breakpoints, total_map_reads, ceil(maqs / total_map_reads))
        return breakpoints, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, 0, 0


def sam_primary_parser2Breaks_Del(region_sam, min_maq, sv_size):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    deletions   = defaultdict(int)    # To store deletions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = int(row.reference_start)
        chr = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
                if 0.5 * sv_size < length < 2 * sv_size:
                    deletion_start = current_start  # Position before deletion starts
                    deletion_end = current_start + length - 1  # Position before the next base
                    deletion_key = f"{chr}:{deletion_start}-{deletion_end}"
                    if deletion_key not in deletions:
                        deletions[deletion_key] = 1
                    else:
                        deletions[deletion_key] += 1  # Count the deletion span
                current_start += length  # Increment position past deletion
    if total_map_reads >0:
        return breakpoints, deletions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_start):
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    spans = 0
    maqs = 0
    effective_span = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chr = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
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
        else:
            if (align_start + 1000 < sv_start)  and (align_end - 1000 > sv_start):
                effective_span += 1

        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  # Increment position past deletion
            elif ctype == 'I' and (0.5 * sv_size < length < 2 * sv_size):
                insertion_start = current_start - 1
                insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                if insert_key not in insertions:
                    insertions[insert_key] = 1
                else:
                    insertions[insert_key] += 1
                effective_span -= 1

    if total_map_reads >0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads), effective_span
    else:
        return {},{},0,0,effective_span
def sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size):
    ## some dup signal may hide in I cigar ##
    breakpoints = defaultdict(lambda: defaultdict(int)) # To store breakpoints
    insertions  = defaultdict(int)   # To store insertions in span format
    total_map_reads = 0
    maqs = 0
    def update_breakpoints(chromosome, positions):
        for pos in positions:
            breakpoints[chromosome][pos] += 1
    for row in region_sam:
        flag = row.flag
        maq  = row.mapping_quality
        maqs += maq
        align_start = row.reference_start
        chrom = row.reference_name
        cigar = row.cigarstring
        #if (row.flag & 0x4) or (row.mapping_quality < min_maq):
        if (row.flag & 0x4):
            continue
        cigar_numbers, cigar_types = parse_cigar(cigar)
        current_start = align_start
        align_end = row.reference_end
        total_map_reads += 1
        breakpoints_to_update = []
        if 'H' in cigar_types or 'S' in cigar_types:
            if cigar_types[0] in 'HS':
                breakpoints_to_update.append(align_start)
            if cigar_types[-1] in 'HS':
                breakpoints_to_update.append(align_end)
            update_breakpoints(chrom, breakpoints_to_update)
        for i in range(len(cigar_numbers)):
            length = cigar_numbers[i]
            ctype = cigar_types[i]
            if ctype in ['M', '=', 'X']:  # Match or mismatch
                current_start += length  # Increment current position for these types
            elif ctype == 'D':
                current_start += length  
            elif ctype == 'I':  # Insertion
                if  0.7 * sv_size < length < 1.5 * sv_size  :  # Only count size equally 
                    insertion_start = current_start - 1
                    insert_key = f"{chr}:{insertion_start}-{insertion_start + 1}"
                    if insert_key not in insertions:
                        insertions[insert_key] = 1
                    else:
                        insertions[insert_key] += 1  # Count the insertion span
    if total_map_reads > 0:
        return breakpoints, insertions, total_map_reads, ceil(maqs / total_map_reads)
    else:
        return {}, {}, 0, 0

def determine_genotype(breaks, depth, homo_rate=0.75,ref_rate = 0.05):
    """
    depth and ratio base genotype
    """
    if depth == 0:
        return "./."
    if depth <= 5:
        if floor(homo_rate*depth) + 1 <= breaks:
            return "1/1"
        elif breaks / depth <= 0.2: ## 2 reads --> 5X
            return "0/0"
        else:
            return "0/1"

    if 5 < depth <= 10:
        if homo_rate * depth + 1 <=  breaks:
            return "1/1"
        elif 0.05 * depth + 1 > breaks:
            return "0/0"
        else:
            return "0/1"
    
    if depth > 10:
        if breaks / depth >= homo_rate: ## or 0.65   or 0.625
            return "1/1"
        elif breaks / depth < ref_rate:
            return "0/0"
        else:
            return "0/1"

def insGT(sampleID, region_sam, chrome, sv_s, sv_e,sv_size, min_maq, homo_rate, ref_rate, shift):
    info_return = []
    genotype = "0/0"  # Default genotype
    shift = min(sv_size, 500)
    sv_start_shift = set(range(sv_s - shift, sv_s + shift ))
    dup_shift = set(range(sv_s - 1000, sv_s - 1000)) ## samll dup capture as ins
    ############ SVIns Case #############
    breakpoints, inserts, total_map_reads, maq, effective_span = sam_primary_parser2Breaks_Ins(region_sam, min_maq, sv_size, sv_s)
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads={total_map_reads},maq=0")
        info_return.append(f"INS_rate=0;INS")
        return info_return
    else:
        count_break_and_Ins = 0
        if inserts:  
            for pos in inserts.keys():
                ins_s, ins_e = map(int, pos.split(":")[1].split("-"))
                if ins_s in sv_start_shift or ins_s in dup_shift:
                    count_break_and_Ins += inserts[pos]
        if breakpoints:   #check breakpoints
            for breakpoint in breakpoints.get(chrome, {}).keys():
                if breakpoint in sv_start_shift: 
                    count_break_and_Ins += breakpoints[chrome][breakpoint]
        ins_ratio = round(count_break_and_Ins / total_map_reads, 3)
        genotype = determine_genotype(count_break_and_Ins, total_map_reads, homo_rate, ref_rate)
        if genotype == "1/1":
            if floor(0.1*total_map_reads) + 1 <=  effective_span:
                genotype = "0/1"
                #print(f"***************Correting SVINS {chrome}:{sv_s}-{sv_e} genotype to 0/1 since it has {effective_span} span reads*****************")
        #print(f"INS\t{genotype}\t{sampleID}\ttotal_mapped_reads:{total_map_reads}\tIns_ratio:{ins_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
        info_return.append(genotype)
        info_return.append(f"total_map_reads={total_map_reads},maq={maq}")
        info_return.append(f"INS_rate={ins_ratio};INS")
    return info_return

def delGT(sampleID, left_sam, right_sam, chrome, sv_s, sv_e, sv_size, min_maq, homo_rate, ref_rate, shift):
    ############ SVDel Case ##############
    info_return = []
    breaks_dict = {}
    genotype = "0/0"  # Default genotype
    sv_start_shift = set(range(sv_s - shift, sv_s + shift))
    sv_end_shift   = set(range(sv_e - shift, sv_e + shift))
    breakpoints_l, deles_l, total_map_reads_l, maq_l = sam_primary_parser2Breaks_Del(left_sam,  min_maq, sv_size)
    breakpoints_r, deles_r, total_map_reads_r, maq_r = sam_primary_parser2Breaks_Del(right_sam, min_maq, sv_size)
    count_break_and_deles_l = 0
    count_break_and_deles_r = 0
    total_map_reads = total_map_reads_l + total_map_reads_r
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_l=0;total_map_reads_r=0,maq=0")
        info_return.append(f"deles_l_ratio=0,deles_r_ratio=0;DEL")
        return  info_return
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
    breaks_dict[count_break_and_deles_l ] = total_map_reads_l
    breaks_dict[count_break_and_deles_r ] = total_map_reads_r
    max_breaks = max(count_break_and_deles_l,count_break_and_deles_r)
    genotype = determine_genotype(max_breaks,breaks_dict[max_breaks], homo_rate, ref_rate)
    #print(f"DEL\t{genotype}\t{sampleID}\ttotal_mapped_reads_l={total_map_reads_l};total_mapped_reads_r={total_map_reads_r}\tdeles_l_ratio:{deles_l_ratio}\tdeles_r_ratio:{deles_r_ratio}\t{chrome}\t{sv_s}\t{sv_e}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_l={total_map_reads_l};total_map_reads_r={total_map_reads_r};maq={max(maq_l,maq_r)}")
    info_return.append(f"deles_l_ratio={deles_l_ratio},deles_r_ratio={deles_r_ratio};DEL")
    return info_return

def breaks2invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.8:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return

def little_dupGT(sampleID, region_sam, chrome, bp1, bp2, sv_size, min_maq, sv_type, shift=500):
    """
    Taking total local mapping  to genotyping the small(<7k) dup SV;
    1st Scan all insertion cigar; 2nd Capture all breakpoints;
    """
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome}:{bp1}", f"{chrome}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1  + shift))
    bp2_shift = set(range(bp2 - shift, bp2  + shift))
    ins_shift1 = set(range(bp1 - sv_size, bp1 + sv_size))
    ins_shift2 = set(range(bp2 - sv_size, bp2 + sv_size))

    breakpoints_bp, inserts, total_map_reads, maq = sam_primary_parser2Breaks_dup(region_sam, min_maq, sv_size)

    count_break_bp = 0
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
        return info_return
    if inserts:
        for pos in inserts.keys():
            ins_s, _ = map(int, pos.split(":")[1].split("-"))
            if ins_s in ins_shift1 or ins_s in ins_shift2:
                count_break_bp += inserts[pos]
    if breakpoints_bp:
        for breakpoint in breakpoints_bp.get(chrome, {}).keys():
            if breakpoint in bp1_shift or breakpoint in bp2_shift:
                count_break_bp += breakpoints_bp[chrome][breakpoint]
    if total_map_reads:
        break_ratio =   round(count_break_bp / total_map_reads, 3)
    else:
        break_ratio = 0
    
    if break_ratio > 0.75:
        genotype = "1/1"
    elif break_ratio < 0.05:
        genotype = "0/0"
    else:
        genotype = "0/1"
    
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads={total_map_reads};\tbreakpoint_ratio={break_ratio}\tbp1={breakpoint1},bp2={breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads={total_map_reads};maq={maq}")
    info_return.append(f"bp1={breakpoint1},bp2={breakpoint2},bp_ratio={break_ratio};{sv_type}")
    return info_return


def dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift):
    """
    Taking two breakpoints local mapping info to genotyping the big SV;
    
    """
    breaks_dict = {}
    info_return = []
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  # Default genotype
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift = set(range(bp2 - shift, bp2 + shift))
    ins_shift1 = set(range(bp1 - shift, bp1 + ceil(shift+ sv_size/2)))
    ins_shift2 = set(range(bp2 - shift - ceil(sv_size / 2), bp2 + shift))

    breakpoints_bp1, inserts_bp1, total_map_reads_bp1, maq1 = sam_primary_parser2Breaks_dup(bp1_sam, min_maq, sv_size)
    breakpoints_bp2, inserts_bp2, total_map_reads_bp2, maq2 = sam_primary_parser2Breaks_dup(bp2_sam, min_maq, sv_size)
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0,maq=0")
        info_return.append(f"bp1={breakpoint1},bp1_ratio=0,bp2={breakpoint2},bp2_ratio=0;{sv_type}")
    
    if inserts_bp1:  # Check if there are insertion entries
        for pos in inserts_bp1.keys():
            if pos in ins_shift1:
                count_break_bp1 += inserts_bp1[pos]
    if breakpoints_bp1:
        for breakpoint in breakpoints_bp1.get(chrome1, {}).keys():
            if breakpoint in bp1_shift:
                count_break_bp1 += breakpoints_bp1[chrome1][breakpoint]

    if inserts_bp2:  # Check if there are insertion entries
        for pos in inserts_bp2.keys():
            if pos in ins_shift2:
                count_break_bp2 += inserts_bp2[pos]
    if breakpoints_bp2:
        for breakpoint in breakpoints_bp2.get(chrome2, {}).keys():
            if breakpoint in bp2_shift:
                count_break_bp2 += breakpoints_bp2[chrome2][breakpoint]

    if total_map_reads_bp1:
        break1_ratio =   round(count_break_bp1 / total_map_reads_bp1, 2)
    else:
        break1_ratio = 0
    if total_map_reads_bp2:
        break2_ratio = round(count_break_bp2 / total_map_reads_bp2, 2)
    else:
        break2_ratio = 0

    max_break_ratio = max(break1_ratio, break2_ratio)
    breaks_dict[count_break_bp1] = total_map_reads_bp1
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_breaks = max(count_break_bp1, count_break_bp2)
    if break1_ratio >=0.4 and break2_ratio >= 0.4:
        genotype = "1/1"
    elif break1_ratio + break2_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio + break2_ratio < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"{sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2}")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"bp1={breakpoint1},bp1_ratio={break1_ratio},bp2={breakpoint2},bp2_ratio={break2_ratio};{sv_type}")
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
    genotype = "0/0"
    info_return = []
    breaks_dict = {}
    bp1_readsID, maq1 = sam2readsID(bp1_sam)
    bp2_readsID, maq2 = sam2readsID(bp2_sam)

    overlapID = list(set(bp1_readsID) & set(bp2_readsID))
    if bp1_readsID:
        bp1_tra = round(len(overlapID) /  len(bp1_readsID), 2)
    else:
        bp1_tra = 0
    if bp2_readsID:
        bp2_tra = round(len(overlapID) / len(bp2_readsID), 2)
    else:
        bp2_tra = 0
    if max(bp1_tra, bp2_tra) > 0.90:
        genotype = "1/1"
    elif bp1_tra >= 0.8 and bp2_tra >= 0.8:
        genotype = "1/1"
    elif bp1_tra + bp2_tra < 0.1:
        genotype = "0/0"
    else:
        genotype = "0/1"
    info_return.append(genotype)
    #print(f"************** TRA GT by reads name  ***************\nbp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    info_return.append(f'total_map_reads_bp1={len(bp1_readsID)};total_map_reads_bp2={len(bp2_readsID)};maq={max(maq1,maq2)}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_ratio={bp1_tra},bp2={chrome2}:{bp2},bp2_ratio={bp2_tra};TRA")
    return info_return

def breaks2traGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """
    Taking two breakpoints region mapping info to genotyping the SV;
    Here we design this func for genotype of INV, TRA, DUP.
    Although it will calculates the span of breakpoints,
    it isnt take these info into genotyping, future may discard or improve.
    """
    ############ two breakpoints Case #############
    info_return = []
    breaks_dict ={}
    breakpoint1, breakpoint2 = f"{chrome1}:{bp1}", f"{chrome2}:{bp2}"
    genotype = "0/0"  
    bp1_shift = set(range(bp1 - shift, bp1 + shift))
    bp2_shift   = set(range(bp2 - shift, bp2 + shift))
    breakpoints_bp1, total_map_reads_bp1, maq1 = sam_parser2Breaks(bp1_sam,  min_maq)
    breakpoints_bp2, total_map_reads_bp2, maq2 = sam_parser2Breaks(bp2_sam, min_maq)
    
    count_break_bp1 = 0
    count_break_bp2 = 0
    total_map_reads = total_map_reads_bp1 + total_map_reads_bp2
    if total_map_reads == 0:
        info_return.append("./.")
        info_return.append(f"total_map_reads_bp1=0;total_map_reads_bp2=0;maq=0")
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
    max_breaks = max(count_break_bp1, count_break_bp2)
    breaks_dict[count_break_bp1] = total_map_reads_bp1 
    breaks_dict[count_break_bp2] = total_map_reads_bp2
    max_break_ratio = max(break1_ratio, break2_ratio)
    if max_break_ratio >= 0.9:
        genotype = "1/1"
    elif break1_ratio >= 0.85 and break2_ratio >= 0.85:
        genotype = "1/1"
    elif break1_ratio+ break2_ratio< 0.10:
        genotype = "0/0"
    else:
        genotype = "0/1"
    #print(f"***breakpoint to genotype TRA********* {sv_type}\t{genotype}\t{sampleID}\ttotal_mapped_reads:bp1={total_map_reads_bp1};bp2={total_map_reads_bp2}\t{breakpoint1}_ratio={break1_ratio}\t{breakpoint2}_ratio={break2_ratio}\t{breakpoint1}\t{breakpoint2} ****************")
    info_return.append(genotype)
    info_return.append(f"total_map_reads_bp1={total_map_reads_bp1};total_map_reads_bp2={total_map_reads_bp2};maq={max(maq1,maq2)}")
    info_return.append(f"{breakpoint1}_ratio={break1_ratio},{breakpoint2}_ratio={break2_ratio};{sv_type}")
    return info_return


def supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0
    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate, maq_bp1 = 0,0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
        maq_bp2 = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def invGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    rate, bp1_map, bp2_map, maq = supp2INVGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    info_return = []
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.6:
        genotype = '1/1'
    elif rate < 0.05:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

def supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    bp1_map, bp2_map, genotype, maq1, maq2 = 0,0,0,0,0

    for line in bp1_sam:
        bp1_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq1 += maq0
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
        maq_bp1 = ceil(maq1 / bp1_map)
    else:
        bp1_rate = 0
    for line in bp2_sam:
        bp2_map += 1
        readname = line.query_name
        ref_chr = line.reference_name
        start = line.reference_start
        flag = line.flag
        maq0 = line.mapping_quality
        maq2 += maq0
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
        maq_bp2 = ceil(maq2 / bp2_map)
    else:
        bp2_rate = 0
    return bp1_rate + bp2_rate, bp1_map, bp2_map, max(maq_bp1, maq_bp2)

def supp_dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800):
    """"
    if a sv size if supper big, than supp aligns could be the signal round the breakpoints
    readsID truely mapping at chr1_bp1 and chr2_bp2 then count it as support not because it have breakpoints,
    a repeat or complex region may easy to cause breakpoints, but that kind of reads region is not a tra signal
    used: hifi and cr, big INV
    """
    info_return = []
    rate, bp1_map, bp2_map, maq = supp2dupGT(sampleID, bp1_sam, bp2_sam, chrome1, chrome2, bp1, bp2, sv_size, min_maq, sv_type, shift=800)
    if bp1_map + bp2_map == 0:
        genotype = "./."
    if rate >=0.625:
        genotype = '1/1'
    elif rate < 0.1:
        genotype = '0/0'
    else:
        genotype = '0/1'
    info_return.append(genotype)
    info_return.append(f'total_map_reads_bp1={bp1_map};total_map_reads_bp2={bp2_map};maq={maq}')
    info_return.append(f"bp1={chrome1}:{bp1},bp1_rate={round(rate,2)},bp2={chrome2}:{bp2},bp2_rate={round(rate,2)};{sv_type}")
    return info_return

```

## ./CapsPop
```shell
total 16
-rwxr-xr-x 1 lgb xinwang 1548 May 22 16:55 common_enzyme.list
-rwxr-xr-x 1 lgb xinwang 1686 May 22 16:55 mpileup_stdin4popcasp.py
-rwxr-xr-x 1 lgb xinwang 7470 May 22 16:55 pop_maf0.05_caps.py

```
### ./CapsPop/mpileup_stdin4popcasp.py
```python
import sys
import re
### samtools mpileup bam.list | python $1  > CapsMarker.input
def count_bases(sequence):
    # Initialize a dictionary to hold counts
    base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    # Loop through each character in the sequence
    for base in sequence.upper():  # Convert to uppercase to handle lowercase inputs
        if base in base_counts:
            base_counts[base] += 1
    return base_counts

for line in sys.stdin:
    line = re.sub(r'\^.',"", line)
    line = line.replace("$", "")
    info = line.strip().split("\t")
    num_base = count_bases(",".join(info[3:]).upper())
    num_samples = int(len(info[3:])/3)
    chroms, pos=info[0],info[1]
    most_freq, second_freq = sorted(num_base.values())[-1], sorted(num_base.values())[-2]
    marker_most_freq =  round(most_freq / num_samples,3)
    second_most_freq =  round(second_freq / num_samples,3)
    if most_freq != second_freq and second_freq >0 and marker_most_freq < 0.95 :
        most_freq_base =   [key for key, value in num_base.items() if value == most_freq][0]
        second_freq_base = [key for key, value in num_base.items() if value == second_freq][0]
        print(chroms, pos, most_freq_base, second_freq_base,most_freq, second_freq ,marker_most_freq, second_most_freq,sep="\t")
    elif most_freq == second_freq and most_freq >0 and marker_most_freq < 0.95:
        most_freq_base =   [key for key, value in num_base.items() if value == most_freq][0]
        second_freq_base = [key for key, value in num_base.items() if value == second_freq][1]
        print(chroms, pos, most_freq_base, second_freq_base,most_freq, second_freq,marker_most_freq ,second_most_freq,sep="\t") 

```

### ./CapsPop/pop_maf0.05_caps.py
```python
import argparse
import gzip
from collections import deque
from gzip import BadGzipFile
import re

hash_seq = {}
hash_site = {}
fa = {}

def parse_thermo_re(file):
    with open(file) as f:
        for line in f:
            enzyme, seq, site = line.strip().split("\t")
            hash_seq[enzyme] = seq
            hash_site[enzyme] = site
    return hash_seq,hash_site

def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

def find_primer_sequences(chrom, pos, con, mut, con_freq, mut_freq, lenp, fa):
    snp_seq_con = fa[chrom][pos - 6 : pos + 5]
    snp_seq_con = snp_seq_con[:5] + con + snp_seq_con[6:]

    snp_seq_mut = fa[chrom][pos - 6:pos + 5]
    snp_seq_mut = snp_seq_mut[:5] + mut + snp_seq_mut[6:]

    output_fa = []
    output_enzyme = []
    for enzyme, seq in hash_seq.items():
        if seq in snp_seq_con and seq not in snp_seq_mut:
            for i in range(lenp - 50, 49, -50):
                y = lenp - i
                snp_seq_con_primer = fa[chrom][pos - i:pos - i + lenp]
                snp_seq_con_primer = snp_seq_con_primer[:i-1] + con + snp_seq_con_primer[i:]
                enzyme_count = snp_seq_con_primer.count(seq)
                if enzyme_count <= 1:
                    output_enzyme.append(f"{chrom}\t{pos}\t{mut}\t{con}\tmut\tnot.{enzyme}\tcon\t{enzyme}.{enzyme_count}\t{hash_site[enzyme]}\t{snp_seq_mut}\t{snp_seq_con}\tmarker_freq:{con_freq}")
                    output_fa.append(f"{chrom}\t{pos}\t{enzyme}\t{i}.{y}\tcon\t{enzyme_count}\t{hash_site[enzyme]}\t{con_freq}")
                    output_fa.append(snp_seq_con_primer)
                elif enzyme_count == 2:
                    all_parts = snp_seq_con_primer.split(seq)
                    length_str = '.'.join(str(len(part)) for part in all_parts)
                    output_enzyme.append(f"{chrom}\t{pos}\t{mut}\t{con}\tmut\tnot.{enzyme}\tcon\t{enzyme}.{enzyme_count}.{i}.{y}-{length_str}\t{hash_site[enzyme]}\t{snp_seq_mut}\t{snp_seq_con}\tmarker_freq:{con_freq}")
                    output_fa.append(f"{chrom}\t{pos}\t{enzyme}\t{i}.{y}-{length_str}\tcon\t{enzyme_count}\t{hash_site[enzyme]}\t{con_freq}")
                    output_fa.append(snp_seq_con_primer)
        elif seq not in snp_seq_con and seq in snp_seq_mut:
            for i in range(lenp - 50, 49, -50):
                y = lenp - i
                snp_seq_mut_primer = fa[chrom][pos - i:pos - i + lenp]
                snp_seq_mut_primer = snp_seq_mut_primer[:i-1] + mut + snp_seq_mut_primer[i:]
                enzyme_count = snp_seq_mut_primer.count(seq)
                if enzyme_count <= 1:
                    output_enzyme.append(f"{chrom}\t{pos}\t{mut}\t{con}\tmut\t{enzyme}\tcon\tnot.{enzyme}\t{hash_site[enzyme]}\t{snp_seq_mut}\t{snp_seq_con}\tmarker_freq:{mut_freq}")
                    output_fa.append(f"{chrom}\t{pos}\t{enzyme}\t{i}.{y}\tmut\t{enzyme_count}\t{hash_site[enzyme]}\t{mut_freq}")
                    output_fa.append(snp_seq_mut_primer)
                elif enzyme_count == 2:
                    all_parts = snp_seq_mut_primer.split(seq)
                    length_str = '.'.join(str(len(part)) for part in all_parts)
                    output_enzyme.append(f"{chrom}\t{pos}\t{mut}\t{con}\tmut\t{enzyme}\tcon\tnot.{enzyme}\t{enzyme}.{enzyme_count}.{i}.{y}-{length_str}\t{hash_site[enzyme]}\t{snp_seq_mut}\t{snp_seq_con}\tmarker_freq:{mut_freq}")
                    output_fa.append(f"{chrom}\t{pos}\t{enzyme}\t{i}.{y}-{length_str}\tmut\t{enzyme_count}\t{hash_site[enzyme]}\t{mut_freq}")
                    output_fa.append(snp_seq_mut_primer)
    return output_fa, output_enzyme

def main():
    parser = argparse.ArgumentParser(description='SNP Primer Finder')
    parser.add_argument('thermo_re', help='Thermo_re.txt file')
    parser.add_argument('ref_fasta', help='Reference FASTA file')
    parser.add_argument('snp_file', help='SNP input file')
    parser.add_argument('out_file', help='Output file')
    parser.add_argument('primer_length', type=int, help='Primer length')

    args = parser.parse_args()
    hash_seq, hash_site = parse_thermo_re(args.thermo_re)
    print(hash_seq)
    fa = readfa2Dict(args.ref_fasta)
    print(fa.keys())
    #print(fa['Lsat_1_v11_chr1'])
    with open(args.snp_file,"r") as snp:
        snplines = snp.readlines()
    out = open(args.out_file, 'w')
    outseq = open(args.out_file + ".seq.txt", "w")
    print("chr\tpos\tmut\tcon\tmut\tmut_enzyme\tcon\tcon_enzyme\tenzyme_sequence\tmut_snp_substr\tcon_snp_substr\tmarker_pop_freq",file=out)
    from tqdm import tqdm
    for line in tqdm(snplines):
        if len(line.strip().split("\t")) != 8:  ##### formatiing ####
            continue
        chrom, pos = line.strip().split("\t")[0], int(line.strip().split("\t")[1])
        if pos -20 < args.primer_length or len(fa[chrom]) - pos < args.primer_length + 20: ## skip row that snp start too start and end too end
            continue
        con, mut = line.strip().split("\t")[2], line.strip().split("\t")[3]
        con_freq, mut_freq = line.strip().split("\t")[6], line.strip().split("\t")[7]
        if con in ['A', 'T', 'C', 'G'] and mut in ['A', 'T', 'C', 'G']:
            output_fa,output_enzyme = find_primer_sequences(chrom, pos, con, mut, con_freq, mut_freq, args.primer_length,fa)
            for output_line in output_fa:
                print(output_line,file =outseq)
            for outline in output_enzyme:
                print(outline,file=out)
    out.close()
    outseq.close()
    
    with open(args.out_file + ".seq.txt", 'r') as file:
        lines = file.readlines()
    combined_lines = []
    for i in range(0, len(lines), 2):
        if i+1 < len(lines):  # Check if there's a next line to combine
            combined_lines.append(lines[i].strip() + '\t' + lines[i+1].strip())
    with open(args.out_file + ".oneline_seq.txt", 'w') as file:
        file.write('\n'.join(combined_lines))

if __name__ == '__main__':
    main()

```

## ./SVInDel_Anno
```shell
total 72392
-rw-r--r-- 1 lgb xinwang    44685 May 22 16:55 Out_SVInDel_each_gene.txt
-rw-r--r-- 1 lgb xinwang   342678 May 22 16:55 Out_SVInDel_pos.txt
-rwxr-xr-x 1 lgb xinwang      181 May 22 16:55 run.sh
-rwxr-xr-x 1 lgb xinwang    14652 Jun 11 09:34 SV_Features_Annotation.py
-rwxr-xr-x 1 lgb xinwang     4551 May 22 16:55 SV_Features_Position.py
-rwxr-xr-x 1 lgb xinwang    13390 May 22 16:55 SVInDels_Feature_Annotation.py
-rwxr-xr-x 1 lgb xinwang     4551 May 22 16:55 SVInDels_Features_Position.py
-rw-r--r-- 1 lgb xinwang   145507 May 22 16:55 SVInDels_Lead_Gene_Variants.txt
-rw-r--r-- 1 lgb xinwang 26875163 May 22 16:55 TAIR10_gene.gff
-rw-r--r-- 1 lgb xinwang 46656028 May 22 16:55 TAIR10_gene.gtf
-rw-r--r-- 1 lgb xinwang       50 May 22 16:55 test

```
### ./SVInDel_Anno/run.sh
```shell
################# example ##################
python SVInDels_Feature_Annotation.py -g TAIR10_gene.gff -s ../final.gt.txt -m ID -c Parent -o test
python SVInDels_Features_Position.py TAIR10_gene.gff ../final.gt.txt Out_SVInDel

```

### ./SVInDel_Anno/SV_Features_Position.py
```python
import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

def parse_args():
    parser = argparse.ArgumentParser(description='Process GFF and indel files')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file or gtf')
    parser.add_argument('pos_file', type=str, help='Path to the pos file, which must have the 3 columns (#Target_name	Target_start	Target_end)')
    parser.add_argument('output_file', type=str, help='the output prefix name')
    return parser.parse_args()

def load_gff(gff_file):
    gff = pd.read_csv(gff_file, header=None,comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = list(set(gff[0].tolist()))
    return gff, chroms


def extract_gene_cds(gff, gff_file):
    if "gff" in str(gff_file):
        gene = gff[gff[2]=="gene"]
        cds = gff[gff[2]=="CDS"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'ID=([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'Parent=([^;]+)')
    elif "gtf" in str(gff_file):
        gene = gff[gff[2]=="transcript"]
        cds = gff[gff[2]=="exon"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'gene_name ([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'gene_name ([^;]+)')

    return gene, cds

def get_3k_promoter(gene):
    pro3k = gene.copy()
    for i in range(pro3k.shape[0]):
        if pro3k.iloc[i,6] == "+":
            pro3k.iloc[i, 4] = pro3k.iloc[i, 3] - 1
            pro3k.iloc[i, 3] = pro3k.iloc[i, 3] - 3000
        else:
            pro3k.iloc[i,3] = pro3k.iloc[i,4] + 1
            pro3k.iloc[i,4] = pro3k.iloc[i,4] + 3000
        pro3k.iloc[i,2] = "promoter"
    return pro3k

def load_indel(indel_file):
    ori_indel = pd.read_csv(indel_file, header=0, sep="\t", index_col=None)
    ori_indel["gene"] = ''
    ori_indel["CDS"] = ''
    ori_indel["3k_promoter"] = ''
    return ori_indel[["#Target_name", "Target_start", "Target_end"]], ori_indel
def overlapper(indel_pos, region, struc, chrom, ori_indel, outeach):
    indel_pos = indel_pos[indel_pos["#Target_name"] == chrom]
    region.columns = ["chr","struc","start","end","ID"]
    region = region[region["chr"] == chrom]
    rtree = IntervalTree()
    
    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        rtree[rStart:rEnd] = region.loc[r, "ID"].replace('"', "")

    for i in indel_pos.index:
        pos1 = indel_pos.loc[i, "Target_start"]
        pos2 = indel_pos.loc[i, "Target_end"]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    print(f'GenePAV\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")

                else:
                    print(f'{struc}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")
def main():
    args = parse_args()
    gff_file = args.gff_file
    indel_file = args.pos_file
    output_file = args.output_file
    gff,chroms = load_gff(gff_file)
    gene, cds = extract_gene_cds(gff, gff_file)
    pro3k = get_3k_promoter(gene)
    
    gene = gene[[0, 2,3,4,"ID"]]
    gene.index = range(gene.shape[0])
    cds = cds[[0, 2,3,4,"ID"]]
    cds.index = range(cds.shape[0])
    pro3k = pro3k[[0,2,3,4,"ID"]]
    pro3k.index = range(pro3k.shape[0])

    indel, ori_indel = load_indel(indel_file)
    outeach = open(f'{output_file}_each_gene.txt', "w")
    for chrom in tqdm(chroms):
        overlapper(indel, gene, "gene", chrom, ori_indel, outeach)
        overlapper(indel, cds, "CDS", chrom, ori_indel, outeach)
        overlapper(indel, pro3k, "3k_promoter", chrom, ori_indel, outeach)

    ori_indel.to_csv(f"{output_file}_pos.txt", header=1, index=None, sep="\t", quoting=None)
    outeach.close()
if __name__ == '__main__':
    main()

```

### ./SVInDel_Anno/SV_Features_Annotation.py
```python
import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from copy import deepcopy

def del_overlapper(del_pos, region, chrom, ori_indel, outeach):
    chr_del = del_pos[del_pos["#Target_name"] == chrom]
    dele_pos = chr_del[chr_del["Target_end"] - chr_del["Target_start"] >40]
    print(dele_pos.head())
    print(chr_del.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in dele_pos.index:
        pos1 = dele_pos.loc[i, "Target_start"]
        pos2 = dele_pos.loc[i, "Target_end"]
        svsize = dele_pos.loc[i,"Target_size" ]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    if "cds" in intervalX.data:
                        impact = "AS"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_lost"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                elif pos1 < intervalX.begin and pos2 <= intervalX.end and pos2 >intervalX.begin:
                    delLeft = pos2 - intervalX.begin + 1
                    if "cds" in intervalX.data and delLeft % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delLeft % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    IDlst.append(f"{intervalX.data}:{impact}")
                elif pos1 > intervalX.begin and pos1 <  intervalX.end and pos2 >intervalX.end:
                    delright = intervalX.end - pos1 + 1
                    if "cds" in intervalX.data and delright % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delright % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smalller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRignt_Lost:{delright}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRight_Lost:{delright}bp\t{impact}')
                elif pos1 >= intervalX.begin and pos2 <= intervalX.end:
                    delcenter = pos2 - pos1 + 1
                    if "cds" in intervalX.data and delcenter % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delcenter % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}')


def ins_overlapper(indel_pos, region, chrom, ori_indel, outeach):
    chr_indel = indel_pos[indel_pos["#Target_name"] == chrom]
    ins_pos = chr_indel[chr_indel["Target_end"] - chr_indel["Target_start"] == 1]
    print(chr_indel.head())
    print(ins_pos.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in ins_pos.index:
        pos1 = ins_pos.loc[i, "Target_start"]
        pos2 = ins_pos.loc[i, "Target_end"]
        svsize = abs(ins_pos.loc[i, "Target_size"])
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if "cds" in intervalX.data and svsize % 3 ==0:
                    impact = "AS"
                elif "cds" in intervalX.data and svsize % 3 !=0:
                    impact = "frameshift"
                elif "utr" in intervalX.data:
                    impact = "Exp"
                elif "intron" in intervalX.data:
                    impact = "intron_exten"
                IDlst.append(f"{intervalX.data}:{impact}")
                print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}', file=outeach)
                #print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}')
            ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")

import pandas as pd

def extract_gene_cds(gff_file, geneid, cdsid):
    gff = pd.read_csv(gff_file, header=None, comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = gff[0].unique()
    if "gff" in str(gff_file):
        gene = gff[gff[2] == "mRNA"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid}=([^;]+)')  # Using f-string for regex
        gene["ID"] = gene["ID"].str.replace('"', '')
        cds["ID"] = cds[8].str.extract(fr'{cdsid}=([^;]+)') 
        cds["ID"] = cds["ID"].str.replace('"', '')
        
    elif "gtf" in str(gff_file):
        gene = gff[gff[2] == "transcript"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid} ([^;]+)')  
        cds["ID"] = cds[8].str.extract(fr'{cdsid} ([^;]+)')  
    return gene, cds, chroms

def gene_feature(geneid, genedf, cdsdf):
    genes_dict = {}
    gene_intron_regions = {}
    gene_cds_regions = {}
    utr5_region = {}
    utr3_region = {}
    geneIDdf = genedf[genedf["ID"] == geneid]
    geneIDcds = cdsdf[cdsdf["ID"] == geneid]
    cds_regions = []
    ### to build interval tree we dont take region strandness ##
    gene_chr = geneIDdf.iloc[0,0]
    gene_start = geneIDdf.iloc[0, 3]
    gene_end = geneIDdf.iloc[0, 4]
    # Collect CDS regions
    for index in geneIDcds.index:
        cds_region = [geneIDcds.loc[index, 3], geneIDcds.loc[index, 4]]
        cds_regions.append(cds_region)
    gene_cds_regions[geneid] = cds_regions
    # Calculate UTR regions
    first_cds_start = geneIDcds.iloc[0, 3]
    last_cds_end = geneIDcds.iloc[-1, 4]
    gene_strand = geneIDdf.iloc[0, 6]
    utr5 = ""
    utr3 = ""
    if gene_strand == "+":
        if int(gene_start) < int(first_cds_start):
            utr5 = [gene_start, first_cds_start]
        else:
            utr5 = "UTR5 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr3 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5
    elif gene_strand == "-":
        # Calculate UTR regions
        if int(gene_start) < int(first_cds_start):
            utr3 = [gene_start, first_cds_start]
        else:
            utr3 = "UTR3 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr5 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
    # Store UTR regions
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5

    # Calculate intron regions
    intron_regions = []
    if len(cds_regions) > 1:
        for i in range(len(cds_regions) - 1):
            if gene_strand == "+":
                # Intron is between the end of the current CDS and the start of the next CDS
                intron_start = cds_regions[i][1] + 1
                intron_end = cds_regions[i + 1][0] - 1
            else:  # "-" strand
                # Intron is between the end of the next CDS and the start of the current CDS
                intron_start = cds_regions[i + 1][1] + 1
                intron_end = cds_regions[i][0] - 1
            if intron_start <= intron_end:  # Only add valid intron regions
                intron_regions.append([intron_start, intron_end])
    gene_intron_regions[geneid] = intron_regions
    return gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr

############## vcf format to below df ########################
# head ../final.gt.txt 
# POS	#Target_name	Target_start	Target_end	Target_size	Query_sizeERR11436000_1.fastq_dedup	ERR11436001_1.fastq_dedup
# Chr1:17145-17369	Chr1	17145	17369	225	-225	0/0	0/0
# Chr1:132223-132224	Chr1	132223	132224	-127	127	1/1	1/1

def main(gff_file, sv_file, output_file, geneid, cdsid):
    vcf = pd.read_csv(sv_file,sep="\t",header=0,index_col=None)
    use = ['ID', '#CHROM','POS','ALT']
    for acc in vcf.columns[9:]:
        use.append(acc)
    
    sv = vcf[use]
    ## insert target_end
    sv.insert(3,'Target_end',sv['ID'].str.split(":",expand=True)[1])
    sv['Target_end'] = sv['Target_end'].str.split("-|_",expand=True)[1].astype(int)
    ## Target size and query_size
    sv.insert(4,'Target_size',sv['ID'].str.split(":",expand=True)[1])
    sv["Target_size"] = sv["Target_size"].str.split("-|_",expand=True)[2].astype(int)
    sv.insert(5,'Query_size',sv['ID'].str.split(":",expand=True)[1])
    sv["Query_size"] = sv["Query_size"].str.split("-|_",expand=True)[2].astype(int)
    for index in sv.index:
        if sv.loc[index,'ALT'] == "<DEL>":
            sv.loc[index,'Query_size'] = -sv.loc[index,'Query_size']
        elif sv.loc[index,'ALT'] == "<INS>":
            sv.loc[index,'Target_size'] = -sv.loc[index,'Target_size']
    print(sv.head())
    cols = ['POS','#Target_name','Target_start','Target_end','Target_size','Query_size','SV'] 
    for acc in vcf.columns[9:]:
        cols.append(acc)
    sv.columns = cols
    print(sv.head())
    sv.to_csv(f'{sv_file}_tmp.tab',header=True,sep="\t",index=None)

    #print(sv.iloc[0,[3]])
    genedf, cdsdf, chroms = extract_gene_cds(gff_file, geneid, cdsid)
    genedf["ID"] = genedf["ID"].str.replace('"', "")
    cdsdf["ID"] = cdsdf["ID"].str.replace('"', "")
    print(cdsdf.head())
    print(genedf.head())
    # Prepare output DataFrame for SV
    outsv = deepcopy(sv)  # Assuming sv is defined globally or passed as an argument
    outsv["struc"] = ""

    with open(output_file, "w") as outeach:
        print('SV\tSVsize\tChr:Start-End\tGeneFeature\tGeneID\tImpact', file=outeach)
        for chrom in chroms:
            bdict = {}
            gdict = {}
            genedf_chr = genedf[genedf[0] == chrom]
            cdsdf_chr = cdsdf[cdsdf[0] == chrom]
            for gene in list(set(cdsdf_chr["ID"].tolist())):
                gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr = gene_feature(gene, genedf_chr, cdsdf_chr)
                for i, cds_region in enumerate(gene_cds_regions[gene]):
                    bdict[f"{gene}__cds{i + 1}"] = cds_region
                for i, intron_region in enumerate(gene_intron_regions[gene]):
                    bdict[f'{gene}__intron{i + 1}'] = intron_region
                    if gene in utr3_region.keys():
                        if utr3_region[gene] != "UTR3 Info Lost":
                            bdict[f'{gene}__utr3'] = utr3_region[gene]
                    if gene in  utr5_region.keys():
                        if utr5_region[gene] != "UTR5 Info Lost":
                            bdict[f'{gene}__utr5'] = utr5_region[gene]
            if bdict:
                gene_body = pd.DataFrame(bdict).T
                gene_body["gene"] = gene_body.index
            
                gene_body["gene"] = gene_body["gene"].str.split("__", expand=True)[0]
                gene_body.dropna(inplace=True)
                print(f"********************* the gene_body format df **********************\n{gene_body.head()}")
                del_overlapper(del_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
                ins_overlapper(indel_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
    outsv.to_csv(f"{args.sv}_anno.txt",header=True,index=None,sep="\t")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process gene and structural variant data.")
    parser.add_argument('-g', '--gff', required=True, help='Path to the GFF file (e.g., TAIR10_gene.gff)')
    parser.add_argument('-s', '--sv', required=True, help='Path to the sv vcf file')
    parser.add_argument('-m', '--mRNA', default="ID", help='to extract GeneID, some gff is not in ID=Geneid, here default is ID')
    parser.add_argument('-c', '--cds', default="Parent", help='to extract cds feature coordinated GeneID, some gff is not in Parent=Geneid, here default is Parent')
    parser.add_argument('-o', '--output', default='SVInDels_Lead_Gene_Variants.txt', help='Output file name (default: SVInDels_Lead_Gene_Variants.txt)')
    args = parser.parse_args()
    main(args.gff, args.sv, args.output, args.mRNA, args.cds)


```

### ./SVInDel_Anno/SVInDels_Feature_Annotation.py
```python
import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from copy import deepcopy

def del_overlapper(del_pos, region, chrom, ori_indel, outeach):
    chr_del = del_pos[del_pos["#Target_name"] == chrom]
    dele_pos = chr_del[chr_del["Target_end"] - chr_del["Target_start"] >2]
    print(dele_pos.head())
    print(chr_del.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in dele_pos.index:
        pos1 = dele_pos.loc[i, "Target_start"]
        pos2 = dele_pos.loc[i, "Target_end"]
        svsize = dele_pos.loc[i,"Target_size" ]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    if "cds" in intervalX.data:
                        impact = "AS"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_lost"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                elif pos1 < intervalX.begin and pos2 <= intervalX.end and pos2 >intervalX.begin:
                    delLeft = pos2 - intervalX.begin + 1
                    if "cds" in intervalX.data and delLeft % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delLeft % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    IDlst.append(f"{intervalX.data}:{impact}")
                elif pos1 > intervalX.begin and pos1 <  intervalX.end and pos2 >intervalX.end:
                    delright = intervalX.end - pos1 + 1
                    if "cds" in intervalX.data and delright % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delright % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smalller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRignt_Lost:{delright}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRight_Lost:{delright}bp\t{impact}')
                elif pos1 >= intervalX.begin and pos2 <= intervalX.end:
                    delcenter = pos2 - pos1 + 1
                    if "cds" in intervalX.data and delcenter % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delcenter % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}')


def ins_overlapper(indel_pos, region, chrom, ori_indel, outeach):
    chr_indel = indel_pos[indel_pos["#Target_name"] == chrom]
    ins_pos = chr_indel[chr_indel["Target_end"] - chr_indel["Target_start"] == 1]
    print(chr_indel.head())
    print(ins_pos.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in ins_pos.index:
        pos1 = ins_pos.loc[i, "Target_start"]
        pos2 = ins_pos.loc[i, "Target_end"]
        svsize = abs(ins_pos.loc[i, "Target_size"])
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if "cds" in intervalX.data and svsize % 3 ==0:
                    impact = "AS"
                elif "cds" in intervalX.data and svsize % 3 !=0:
                    impact = "frameshift"
                elif "utr" in intervalX.data:
                    impact = "Exp"
                elif "intron" in intervalX.data:
                    impact = "intron_exten"
                IDlst.append(f"{intervalX.data}:{impact}")
                print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}', file=outeach)
                #print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}')
            ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")

import pandas as pd

def extract_gene_cds(gff_file, geneid, cdsid):
    gff = pd.read_csv(gff_file, header=None, comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = gff[0].unique()
    if "gff" in str(gff_file):
        gene = gff[gff[2] == "mRNA"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid}=([^;]+)')  # Using f-string for regex
        gene["ID"] = gene["ID"].str.replace('"', '')
        cds["ID"] = cds[8].str.extract(fr'{cdsid}=([^;]+)') 
        cds["ID"] = cds["ID"].str.replace('"', '')
        
    elif "gtf" in str(gff_file):
        gene = gff[gff[2] == "transcript"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid} ([^;]+)')  
        cds["ID"] = cds[8].str.extract(fr'{cdsid} ([^;]+)')  
    return gene, cds, chroms

def gene_feature(geneid, genedf, cdsdf):
    genes_dict = {}
    gene_intron_regions = {}
    gene_cds_regions = {}
    utr5_region = {}
    utr3_region = {}
    geneIDdf = genedf[genedf["ID"] == geneid]
    geneIDcds = cdsdf[cdsdf["ID"] == geneid]
    cds_regions = []
    ### to build interval tree we dont take region strandness ##
    gene_chr = geneIDdf.iloc[0,0]
    gene_start = geneIDdf.iloc[0, 3]
    gene_end = geneIDdf.iloc[0, 4]
    # Collect CDS regions
    for index in geneIDcds.index:
        cds_region = [geneIDcds.loc[index, 3], geneIDcds.loc[index, 4]]
        cds_regions.append(cds_region)
    gene_cds_regions[geneid] = cds_regions
    # Calculate UTR regions
    first_cds_start = geneIDcds.iloc[0, 3]
    last_cds_end = geneIDcds.iloc[-1, 4]
    gene_strand = geneIDdf.iloc[0, 6]
    utr5 = ""
    utr3 = ""
    if gene_strand == "+":
        if int(gene_start) < int(first_cds_start):
            utr5 = [gene_start, first_cds_start]
        else:
            utr5 = "UTR5 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr3 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5
    elif gene_strand == "-":
        # Calculate UTR regions
        if int(gene_start) < int(first_cds_start):
            utr3 = [gene_start, first_cds_start]
        else:
            utr3 = "UTR3 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr5 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
    # Store UTR regions
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5

    # Calculate intron regions
    intron_regions = []
    if len(cds_regions) > 1:
        for i in range(len(cds_regions) - 1):
            if gene_strand == "+":
                # Intron is between the end of the current CDS and the start of the next CDS
                intron_start = cds_regions[i][1] + 1
                intron_end = cds_regions[i + 1][0] - 1
            else:  # "-" strand
                # Intron is between the end of the next CDS and the start of the current CDS
                intron_start = cds_regions[i + 1][1] + 1
                intron_end = cds_regions[i][0] - 1
            if intron_start <= intron_end:  # Only add valid intron regions
                intron_regions.append([intron_start, intron_end])
    gene_intron_regions[geneid] = intron_regions
    return gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr

############## final.gt.txt format ########################
# head ../final.gt.txt 
# POS	#Target_name	Target_start	Target_end	Target_size	Query_sizeERR11436000_1.fastq_dedup	ERR11436001_1.fastq_dedup
# Chr1:17145-17369	Chr1	17145	17369	225	-225	0/0	0/0
# Chr1:132223-132224	Chr1	132223	132224	-127	127	1/1	1/1

def main(gff_file, sv_file, output_file, geneid, cdsid):
    sv = pd.read_csv(sv_file,sep="\t",header=0,index_col=None)
    print(sv.iloc[0,[3]])
    genedf, cdsdf, chroms = extract_gene_cds(gff_file, geneid, cdsid)
    genedf["ID"] = genedf["ID"].str.replace('"', "")
    cdsdf["ID"] = cdsdf["ID"].str.replace('"', "")
    print(cdsdf.head())
    print(genedf.head())
    # Prepare output DataFrame for SV
    outsv = deepcopy(sv)  # Assuming sv is defined globally or passed as an argument
    outsv["struc"] = ""

    with open(output_file, "w") as outeach:
        print('SV\tSVsize\tChr:Start-End\tGeneFeature\tGeneID\tImpact', file=outeach)
        for chrom in chroms:
            bdict = {}
            gdict = {}
            genedf_chr = genedf[genedf[0] == chrom]
            cdsdf_chr = cdsdf[cdsdf[0] == chrom]
            for gene in list(set(cdsdf_chr["ID"].tolist())):
                gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr = gene_feature(gene, genedf_chr, cdsdf_chr)
                for i, cds_region in enumerate(gene_cds_regions[gene]):
                    bdict[f"{gene}__cds{i + 1}"] = cds_region
                for i, intron_region in enumerate(gene_intron_regions[gene]):
                    bdict[f'{gene}__intron{i + 1}'] = intron_region
                    if utr3_region[gene] != "UTR3 Info Lost":
                        bdict[f'{gene}__utr3'] = utr3_region[gene]
                    if utr5_region[gene] != "UTR5 Info Lost":
                        bdict[f'{gene}__utr5'] = utr5_region[gene]
            gene_body = pd.DataFrame(bdict).T
            gene_body["gene"] = gene_body.index
            gene_body["gene"] = gene_body["gene"].str.split("__", expand=True)[0]
            gene_body.dropna(inplace=True)
            print(f"********************* the gene_body format df **********************\n{gene_body.head()}")
            del_overlapper(del_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
            ins_overlapper(indel_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
    outsv.to_csv(f"{args.sv}_anno.txt",header=True,index=None,sep="\t")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process gene and structural variant data.")
    parser.add_argument('-g', '--gff', required=True, help='Path to the GFF file (e.g., TAIR10_gene.gff)')
    parser.add_argument('-s', '--sv', required=True, help='Path to the SV file (e.g., SVInDel.gt.txt) <POS   #Target_name    Target_start    Target_end      Target_size     Query_size GT>')
    parser.add_argument('-m', '--mRNA', default="ID", help='to extract GeneID, some gff is not in ID=Geneid, here default is ID')
    parser.add_argument('-c', '--cds', default="Parent", help='to extract cds feature coordinated GeneID, some gff is not in Parent=Geneid, here default is Parent')
    parser.add_argument('-o', '--output', default='SVInDels_Lead_Gene_Variants.txt', help='Output file name (default: SVInDels_Lead_Gene_Variants.txt)')
    args = parser.parse_args()
    main(args.gff, args.sv, args.output, args.mRNA, args.cds)


```

### ./SVInDel_Anno/SVInDels_Features_Position.py
```python
import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

def parse_args():
    parser = argparse.ArgumentParser(description='Process GFF and indel files')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file or gtf')
    parser.add_argument('pos_file', type=str, help='Path to the pos file, which must have the 3 columns (#Target_name	Target_start	Target_end)')
    parser.add_argument('output_file', type=str, help='the output prefix name')
    return parser.parse_args()

def load_gff(gff_file):
    gff = pd.read_csv(gff_file, header=None,comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = list(set(gff[0].tolist()))
    return gff, chroms


def extract_gene_cds(gff, gff_file):
    if "gff" in str(gff_file):
        gene = gff[gff[2]=="gene"]
        cds = gff[gff[2]=="CDS"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'ID=([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'Parent=([^;]+)')
    elif "gtf" in str(gff_file):
        gene = gff[gff[2]=="transcript"]
        cds = gff[gff[2]=="exon"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'gene_name ([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'gene_name ([^;]+)')

    return gene, cds

def get_3k_promoter(gene):
    pro3k = gene.copy()
    for i in range(pro3k.shape[0]):
        if pro3k.iloc[i,6] == "+":
            pro3k.iloc[i, 4] = pro3k.iloc[i, 3] - 1
            pro3k.iloc[i, 3] = pro3k.iloc[i, 3] - 3000
        else:
            pro3k.iloc[i,3] = pro3k.iloc[i,4] + 1
            pro3k.iloc[i,4] = pro3k.iloc[i,4] + 3000
        pro3k.iloc[i,2] = "promoter"
    return pro3k

def load_indel(indel_file):
    ori_indel = pd.read_csv(indel_file, header=0, sep="\t", index_col=None)
    ori_indel["gene"] = ''
    ori_indel["CDS"] = ''
    ori_indel["3k_promoter"] = ''
    return ori_indel[["#Target_name", "Target_start", "Target_end"]], ori_indel
def overlapper(indel_pos, region, struc, chrom, ori_indel, outeach):
    indel_pos = indel_pos[indel_pos["#Target_name"] == chrom]
    region.columns = ["chr","struc","start","end","ID"]
    region = region[region["chr"] == chrom]
    rtree = IntervalTree()
    
    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        rtree[rStart:rEnd] = region.loc[r, "ID"].replace('"', "")

    for i in indel_pos.index:
        pos1 = indel_pos.loc[i, "Target_start"]
        pos2 = indel_pos.loc[i, "Target_end"]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    print(f'GenePAV\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")

                else:
                    print(f'{struc}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")
def main():
    args = parse_args()
    gff_file = args.gff_file
    indel_file = args.pos_file
    output_file = args.output_file
    gff,chroms = load_gff(gff_file)
    gene, cds = extract_gene_cds(gff, gff_file)
    pro3k = get_3k_promoter(gene)
    
    gene = gene[[0, 2,3,4,"ID"]]
    gene.index = range(gene.shape[0])
    cds = cds[[0, 2,3,4,"ID"]]
    cds.index = range(cds.shape[0])
    pro3k = pro3k[[0,2,3,4,"ID"]]
    pro3k.index = range(pro3k.shape[0])

    indel, ori_indel = load_indel(indel_file)
    outeach = open(f'{output_file}_each_gene.txt', "w")
    for chrom in tqdm(chroms):
        overlapper(indel, gene, "gene", chrom, ori_indel, outeach)
        overlapper(indel, cds, "CDS", chrom, ori_indel, outeach)
        overlapper(indel, pro3k, "3k_promoter", chrom, ori_indel, outeach)

    ori_indel.to_csv(f"{output_file}_pos.txt", header=1, index=None, sep="\t", quoting=None)
    outeach.close()
if __name__ == '__main__':
    main()

```

## ./SVInDel_Imputation
```shell
total 164
-rwxr-xr-x 1 lgb xinwang   6067 May 22 16:55 Imputing_SVInDel.py
-rw-r--r-- 1 lgb xinwang 157866 May 22 16:55 test.SVInDel.tab

```
### ./SVInDel_Imputation/Imputing_SVInDel.py
```python
import pandas as pd
import numpy as np
from copy import deepcopy
import argparse
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
def calculate_missing_rate(row):
    # Count the number of missing genotypes
    missing_count = sum(genotype in {'.', './.'} for genotype in row)
    total_count = len(row)
    return missing_count / total_count if total_count > 0 else 0

def gtdf_to_haplotype_array_one_hot(df):
    haplotype_dict = {}
    sample_ids = df.columns[1:]  # Assuming the first column contains something like IDs or a placeholder
    for index, row in df.iterrows():
        genotypes = row[1:]  # Skip the first column if it doesn't contain genotype data
        for sample_index, genotype in enumerate(genotypes):
            sample_id = sample_ids[sample_index]
            if sample_id not in haplotype_dict:
                haplotype_dict[sample_id] = []
            if genotype == '0/0':
                haplotype_dict[sample_id].append([1, 0, 0])  # Homozygous reference
            elif genotype == '1/1':
                haplotype_dict[sample_id].append([0, 1, 0])  # Homozygous alternate
            elif genotype in ['0/1', '1/0']:
                haplotype_dict[sample_id].append([0, 0, 1])  # Heterozygous
            elif genotype in ['.', './.']:  # Handle missing genotypes
                haplotype_dict[sample_id].append([0, 0, 0])  # Missing genotype representation
            else:
                haplotype_dict[sample_id].append([0, 0, 0])  # Unknown genotype
    haplotypes_array = np.array([np.array(haplotype_dict[sample_id]).flatten() for sample_id in haplotype_dict])
    return haplotypes_array, sample_ids.tolist()

def blood_hood(gtdf):
    haplotypes_array, sample_ids = gtdf_to_haplotype_array_one_hot(gtdf)
    unique_haplotypes, unique_indices = np.unique(haplotypes_array, axis=0, return_inverse=True)
    n_clusters = len(unique_haplotypes)   ## this values matter a lots
    distance_matrix = pairwise_distances(unique_haplotypes, metric='hamming')
    # Clustering using Agglomerative Clustering
    agg_cluster = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
    clusters = agg_cluster.fit_predict(distance_matrix)
    # Create a DataFrame to map unique haplotypes to their cluster assignments
    cluster_mapping = pd.DataFrame({
        'Unique_Haplotype': [list(hap) for hap in unique_haplotypes],
        'Cluster': clusters
    })
    # Map each sample to its respective cluster based on the unique haplotypes
    sample_clusters = pd.DataFrame({
        'Sample_ID': sample_ids,
        'Cluster': clusters[unique_indices]
    })
    #print("Number of unique clusters:", n_clusters)
    #print("Cluster assignments for each sample:")
    #print(sample_clusters)
    return sample_clusters

def hap_cluster2impute(rowGT: pd.Series, hap_clusters: pd.DataFrame) -> pd.Series:
    """
    Process genotype data and aggregate by haplotype clusters, maintaining original rowGT format.
    Parameters:
    - rowGT: Series containing genotype data for a specific genomic region.
    - hap_clusters: DataFrame containing sample IDs and their corresponding clusters.
    Returns:
    - A Series with aggregated genotype information for each sample, structured like rowGT.
    """
    result_series = pd.Series(index=rowGT.index , dtype=str)
    for sample in rowGT.index:
        result_series[sample] = rowGT[sample]  # Start with the original value
    for cluster_id in hap_clusters['Cluster'].unique():
        cluster_samples = hap_clusters[hap_clusters['Cluster'] == cluster_id]['Sample_ID'].tolist()
        cluster_genotypes = rowGT[cluster_samples]
        non_missing_gt = cluster_genotypes[cluster_genotypes != './.'].unique()
        if len(non_missing_gt) > 0:
            # If there's a valid genotype, use it to fill in missing genotypes
            valid_gt = non_missing_gt[0]
            for sample in cluster_samples:
                if result_series[sample] == './.':
                    result_series[sample] = valid_gt
    return result_series

def imputing(args):
    gt = pd.read_csv(args.SVInDelGT,header=0,sep="\t", index_col=0)
    gt.sort_values(by=["#Target_name", "Target_start"], inplace=True)
    pos = gt.index
    samples = gt.columns[5:]
    gt.index = range(gt.shape[0]) ## avoid same SVInDel break record of Insertion
    out_imputed_gt = deepcopy(gt)
    for chrom in gt["#Target_name"].unique():
        chromgt = gt[gt["#Target_name"] == chrom]
        for i, index in enumerate(chromgt.index):
            rowGT = chromgt.loc[index, :]
            missing = calculate_missing_rate(rowGT)
            if missing <= args.miss and missing > 0 :
                lower_bound = max(0, i - args.Kneigbor)  # Ensure we don't go below 0
                upper_bound = min(len(chromgt) - 1, i + args.Kneigbor)  # Ensure we don't exceed the bounds
                Local_df = chromgt.iloc[lower_bound:upper_bound + 1]  # Include upper bound
                Local_gt = Local_df[samples]
                haplotype_clusters = blood_hood(Local_gt)
                target_gt = Local_gt[samples].loc[index, :]
                imputed_gt = hap_cluster2impute(target_gt, haplotype_clusters)
                out_imputed_gt.loc[index,samples] = imputed_gt
    out_imputed_gt.index = pos
    out_imputed_gt.to_csv(f"{args.SVInDelGT}_imputed_miss{args.miss}_Kneigbor{args.Kneigbor}.txt",sep="\t",header=True,index=True)
    return out_imputed_gt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Haplotypes Cluster to Impute the SVInDel Genotype file")
    parser.add_argument('SVInDelGT', type=str, help='the Genotype File Generated by SVInDel Program.')
    parser.add_argument('--Kneigbor', type=int, help='to Impute the Missing, we use the up down Neigbor SVs to Cluster Haplotype default is 5', default=5)
    parser.add_argument('--miss', type=float,default=0.8, help='missing rate to impute,  0.8 recommended, default is 0.8')
    args = parser.parse_args()
    imputing(args)

```

## ./utils
```shell
total 4
-rw-r--r-- 1 lgb xinwang 680 Aug 14 17:26 longreads_265M_cigar_out.py

```
### ./utils/longreads_265M_cigar_out.py
```python
import re
import sys
insam = sys.argv[1]
outsam = sys.argv[2]
MAX_HS_LEN = 265000000
cigar_pattern = re.compile(r'(\d+)([MIDNSHPX=])')
with open(insam, 'r') as fin, open(outsam, 'w') as fout:
    for line in fin:
        if line.startswith('@'):
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        cigar = parts[5]
        valid = True
        for length_str, op in cigar_pattern.findall(cigar):
            if op in ['H', 'S']:
                if int(length_str) > MAX_HS_LEN:
                    valid = False
                    break
        if valid:
            fout.write(line)

```

