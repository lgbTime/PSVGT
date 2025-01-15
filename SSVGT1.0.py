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
