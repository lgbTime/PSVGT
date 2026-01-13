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
