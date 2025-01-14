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
if __name__ == "__main__":
    import pyfiglet
    def print_large_pro(pro):
        ascii_art = pyfiglet.figlet_format(f'{pro}', font="slant")
        print(ascii_art)
    print_large_pro()
    import argparse
    from time import time
    parser = argparse.ArgumentParser(description="A versitile tools from PSVGT", formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outdir", default="PopSVGT_out",help="output dir of population mode")
    parser.add_argument("-sr", "--srdir", help="a directory contain pair end short reads data only")
    parser.add_argument("-hifi", "--hifidir", help="a directory contain hifi reads file")
    parser.add_argument("-ont", "--ontdir", help="a directory contain nanopore ont reads file")
    parser.add_argument("-pb", "--pbdir", help="a directory contain PacBio CLR genomic reads file")
    parser.add_argument("-cr", "--crdir", help="a directory contain the fasta file of assembly genome in contig level or chromosome level")
    parser.add_argument("-w", "--max_workers", default = 4,type=int, help="the max workers thread pool excutor, 4 means run for samples at a time")
    parser.add_argument("-t", "--threads", default = 10, help="the cpu use to assembly contig or bwa mapping")
    parser.add_argument("-minimapCPU", "--minimapCPU", default = 10, help="the cpu in minimap mapping")
    parser.add_argument("-r", "--refGenome", required=True, help="the reference genome use in mapping")
    parser.add_argument("-g", "--gff", help="gff file to annotate the SV genotyping")
    parser.add_argument("-m",  "--min", default=50, help= "The min length of an SVIndel ")
    parser.add_argument("-e",  "--popcaps",default="no", help= "population caps analysis, the caps marker has a maf >= 0.05 will be output, input yes PopCaps will perform the analysis")
    parser.add_argument("-p",  "--popInDel",default="yes", help= "using the primer3 to design the primer for each SVInDel")
    parser.add_argument("-b",  "--breaker",default="no", help= "using the break points info to support the SVInDel Genotyping, this will perform bwa mapping process and breakpoints genotype")
    parser.add_argument("-maq",  "--maq",default=45,type=int, help= "the mapping quality to caculate break points and mapping coverge range from 30-60")
    parser.add_argument("-support",  "--support",default=0.2, type=float, help= "the percent of reads that support a candidate SVInDel(0.15 means 20X data should have at least 3 reads support SVINDel), this parameter is for the variaty depth of HIFI/ONT/PB samples")
    parser.add_argument("-msv","--msv_mode",default="no", help= "msv mode is to use the supplementary alignment from ONT/HIFI/CLR reads to detect INS,DEL,INV,DUP,TRA, default is SVInDel mode Only detect SVInDel by unique mapping")

def main():
    parser = argparse.ArgumentParser(description='versite tools from PSVGT')
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', dest='subcommand')

    task1_parser = subparsers.add_parser('vcf2InDelprimer', help='taken the PSVGT SV format file to design the maf setted primers of InDel')
    task1_parser.add_argument('-vcf',dest='vcf', type=str, help='vcf in PSVGT format', required=True)
    task1_parser.add_argument('-r',dest='ref', type=str, help='the reference genome', required=True)
    task1_parser.add_argument('-min',dest='min', type=str, help='the min size of InDel', default=50, type=int)
    task1_parser.add_argument("-max", dest= "max", help="the max size of InDel", default=500, type=int)
    task1_parser.add_argument("-maf", dest= "maf", help="the minor allele frequency of the InDel", default=0.01, type = float)

    task2_parser = subparsers.add_parser('fine_gwas', help='from qqnorm pheno file to emmax gwas result,kinship;genotype already done before')
    task2_parser.add_argument('-p',dest='pheno',type=str, help='qqnorm of phenotype file', required=True)
    task2_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task2_parser.add_argument("-c", dest="pca", help="pca file format for emmax")
    task2_parser.add_argument("-sig", dest= "sig", help="significance line height which equal to -log10(1 or 0.05 /(effective snps))")
    task2_parser.add_argument("-kin", dest= "kin", help="BN or IBS kin-ship in emmax association analysis", default="BN")
    task2_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")
    task2_parser.add_argument("-t", dest="binary",help="if binary trait 0/1 fill 'yes',else,default is no", default="no")

    task3_parser = subparsers.add_parser('emmax-kin', help='generate emmax-kin ship matrix')
    task3_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task3_parser.add_argument("-kin", dest= "kin", help="BN or IBS", default="BN")
    task3_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task6_parser = subparsers.add_parser('emmax_genotype', help='generate emmax genotyoe file')
    task6_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task6_parser.add_argument('-keep', dest='keeplst',help='a fam ID list has two columns, first col is the same as second col ')
    task6_parser.add_argument('-gt', dest='outgt',help='the prefix of the output genotype file')
    task6_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task7_parser = subparsers.add_parser('pca', help='flashpca2 to perform pca analysis of genotype file')
    task7_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task7_parser.add_argument('-mp', dest='pca',help='at least pca number to select for emmax analysis',default=2,type=int)
    task7_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task4_parser = subparsers.add_parser('ps2plotdata', help='to only perform qq plot and manhatan plot task we should convert emmax_result.ps file to snpID chr position pvalue format')
    task4_parser.add_argument('-ps',dest= "ps", type=str, help='emmax out result')
    task4_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task5_parser = subparsers.add_parser('plot', help='perform qq plot and manhatan plot use my.qqman file')
    task5_parser.add_argument('-qqman',dest= "qqman", type=str, help='emmax out result in my.qqman format')
    task5_parser.add_argument("-sig", dest= "sig", help="significance line height which equal to -log10(1 or 0.05 /(effective snps))",required=True)
    task5_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    args = parser.parse_args()
    start_t = time()
    all_log = open("log4SVGT.txt", "w")
    check_dir(args.outdir)
    fa = readfa2Dict(args.refGenome)
    if not os.path.isfile(f'{args.refGenome}.fai'):
        run_command(f"samtools faidx {args.refGenome}")
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
    Pop_SVInDel_Analysor_cmds = []
    done_analysor = file_capture(args.outdir, ".record.txt")
    def add_commands(files, dtype):
        for contig in files:
            done_name = f"{args.outdir}/0_tmp_{basename(contig)}.record.txt"
            if done_name not in done_analysor:
                maq = min(args.maq + 5, 60) ## biger than illumina breakpoints quality
                cmd = f'python {PSVGT_folder}/SV_Genotyper/0.sospop_signaling.py -i {contig} -dtype {dtype} -r {args.refGenome} -m {args.min} -maq {maq} -o {args.outdir} -minimapCPU {args.minimapCPU} -msv {args.msv_mode}'
                if dtype not in ['sr','cr']:  # Only add filter for non-cr types
                    cmd += f' && python {PSVGT_folder}/SV_Genotyper/0.sv_signal_filter.py -f {done_name} -c {args.support}'
                print(cmd)
                Pop_SVInDel_Analysor_cmds.append(cmd)
    # Capture files for different types
    if args.srdir:
        add_commands(contigs, "sr") #### the short reads assembly reads use ont mode to mapping ####
    if args.hifidir:
        add_commands(file_capture(args.hifidir, ".gz"), "hifi")
    if args.ontdir:
        add_commands(file_capture(args.ontdir, ".gz"), "ont")
    if args.pbdir:
        add_commands(file_capture(args.pbdir, ".gz"), "pb")
    if args.crdir:
        add_commands(file_capture(args.crdir, "a"), "cr")  ## to get fasta or fa
    # Execute commands using ThreadPool
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures1 = [executor.submit(execute_commands, cmd) for cmd in Pop_SVInDel_Analysor_cmds]
        results1 =[]
        for future in as_completed(futures1):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=all_log)
            else:
                print(stdout,file=all_log)
            results1.append((stdout,stderr,returncode,cmd))
    
    ## step1 to get uniq population SV records and clustering the signal by breakpoints shift ##
    run_command(f"python {PSVGT_folder}/SV_Genotyper/1.popSV_signal_cluster.py -d {args.outdir} -s 50")
    

    ## step2 genotypiing by long seq mapping map ##
    mapinfo_files = file_capture(f"./{args.outdir}", ".bam")
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

            cmd = f"python {PSVGT_folder}/SV_Genotyper/2.Pop_lrSVGT_V1.py -i {args.outdir}/PopSV_clustered_Record.txt -mapf {mapinfo_file}  -n {acc_name} -o {args.outdir} && python {PSVGT_folder}/SV_Genotyper/sospop_tab2vcf.py {if_done_name} {if_done_name.replace('.txt', '')}.vcf"
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
                bpgt_cmd =  f"python {PSVGT_folder}/SV_Genotyper/2.Pop_srSVGT_V1.py -i {args.outdir}/PopSV_clustered_Record.txt -mapf {bam} -s 25 -n {sampleID} -o {args.outdir} && python {PSVGT_folder}/SV_Genotyper/sospop_tab2vcf.py {args.outdir}/2_tmp_{sampleID}_bpgenotype.txt {args.outdir}/2_tmp_{sampleID}_bpgenotype.vcf"
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
        print("python {PSVGT_folder}/SVInDel_Primer/vcf2primer.py {args.outdir}/PSVGT_all.vcf2.SVInDel {args.refGenome} 50 500 500 > {args.outdir}/PSVInDel_Primer4Pop.txt ")
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
