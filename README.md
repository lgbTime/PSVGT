![PSVGT](https://github.com/lgbTime/PSVGT/logo.png)                         
## PSVGT is a versatile program for SV detection and genotyping at population scale, supporting all sequence length (>150+ pbs) and depth free, providing CAPS & InDel development Modules and SV annotation function
The short reads samples will be de novo into contigs, while assemble samples and long reads sequencing data could be directory used in the SV detection and Genotyping.
## 💻 Installation
PSVGT required **megahit**, **minimap2**, **samtools**, and has been tested at python3.11
the  folow python dependency is required.
```sh
git clone https://github.com/lgbTime/PSVGT.git
pip install primer3-py
pip install pyfiglet
pip install pandas 
pip install tqdm
pip install pysam
pip install intervaltree
```

## Performance of PSVGT at benchmark 
```
Detect & GT	PSVGT			DeBreak			cuteSV			Sniffles2
	Rec	Pre	F1	Rec	Pre	F1	Rec	Pre	F1	Rec	Pre	F1
HiFi												
DEL	99.86	99.91	99.88	96.97	97.03	97.00	99.51	99.85	99.68	99.89	99.95	99.92
INS	98.46	99.62	99.03	83.82	84.89	84.35	97.90	99.69	98.79	94.49	94.20	94.34
DUP	92.60	94.39	93.49	90.20	90.02	90.11	65.70	69.82	67.70	78.00	66.38	71.72
INV	96.30	98.27	97.27	84.60	85.63	85.11	97.20	96.14	96.67	71.40	89.92	79.60
TRA	98.63	99.59	99.11	64.93	64.93	64.93	98.49	96.64	97.56	86.30	86.07	86.18
ONT												
DEL	99.89	99.89	99.89	97.01	96.99	97.00	94.97	99.78	97.32	99.87	99.39	99.63
INS	98.47	99.05	98.76	84.17	84.98	84.57	93.19	99.46	96.22	94.44	93.70	94.07
DUP	95.40	96.17	95.78	86.30	89.80	88.02	63.20	74.35	68.32	78.50	71.43	74.80
INV	95.40	98.45	96.90	86.80	88.93	87.85	90.90	95.58	93.18	70.50	88.01	78.29
TRA	98.49	99.31	98.90	66.30	66.12	66.21	98.77	96.01	97.37	86.99	86.87	86.93
PacBio CLR												
DEL	99.79	99.92	99.85	97.87	98.04	97.95	99.16	99.61	99.38	99.58	99.84	99.71
INS	97.62	99.24	98.42	86.01	87.93	86.96	96.20	99.12	97.64	93.89	94.20	94.05
DUP	96.80	97.98	97.38	89.80	93.15	91.45	61.80	66.81	64.21	77.90	71.47	74.55
INV	97.20	99.18	98.18	88.20	88.73	88.47	97.10	95.29	96.19	72.20	83.86	77.59
TRA	98.36	99.58	98.97	66.30	65.41	65.85	99.04	96.14	97.57	98.90	89.25	93.83
Short Reads(150bp)												
DEL	34.33	96.14	50.59	*	*	*	*	*	*	*	*	*
INS	91.32	99.84	95.39	*	*	*	*	*	*	*	*	*
DUP	*	*	*	*	*	*	*	*	*	*	*	*
INV	*	*	*	*	*	*	*	*	*	*	*	*
TRA	*	*	*	*	*	*	*	*	*	*	*	*
Assembly(two haplotype genome)												
DEL	98.61	99.97	99.29	*	*	*	*	*	*	*	*	*
INS	98.64	99.99	99.31	*	*	*	*	*	*	*	*	*
DUP	96.50	100.00	98.22	*	*	*	*	*	*	*	*	*
INV	*	*	*	*	*	*	*	*	*	*	*	*
TRA	*	*	*	*	*	*	*	*	*	*	*	*
```									
																								
## force GT using hifi signal
```
Force GT	PSVGT			DeBreak			cuteSV			Sniffles2		
	Rec	Pre	F1	Rec	Pre	F1	Rec	Pre	F1	Rec	Pre	F1
Short Reads(150bp)												
DEL GT	99.55	99.68	99.61	*	*	*	*	*	*	*	*	*
INS GT	98.14	99.29	98.71	*	*	*	*	*	*	*	*	*
DUP GT	95.34	100.00	97.61	*	*	*	*	*	*	*	*	*
INV GT	97.50	99.49	98.48	*	*	*	*	*	*	*	*	*
TRA GT	98.49	99.86	99.17	*	*	*	*	*	*	*	*	*
Assembly Genome(two haplotype genome)											
DEL GT	99.94	99.99	99.96	*	*	*	*	*	*	*	*	*
INS GT	98.64	99.99	99.31	*	*	*	*	*	*	*	*	*
INV GT	98.00	100.00	98.99	*	*	*	*	*	*	*	*	*
TRA GT	98.63	100.00	99.31	*	*	*	*	*	*	*	*	*
``` 

## ⏩ Quick Start
### One-Step PSVGT

```shell
python PSVGT1.0/PSVGT1.0.py -hifi test_hifi \
			    -ont test_ont  \
			    -pb test_pb  \
			    -sr test_sr \
			    -cr test_cr \
			    -r Db-1_ref/Db-1_genome.fa -b yes -o out4PSVGT -msv yes
```

For specific sequence types, PSVGT uses a coordinated mapping strategy and program. In the folders `test_hifi`, `test_ont`, and `test_pb`, the sequencing data should be in gz-compressed FASTQ or mapped BAM format. The samples in the `test_cr` folder should be FASTA files with `.fasta` or `.fa` file extensions. The `test_sr` folder allows paired-end short reads data in gz-compressed FASTQ format.

- `-r`: Specifies the reference genome.
- `-o`: Specifies the output folder.
- `-msv`: Detects all SVs instead of just SVInDels.
- `-cr`: the samples in test_cr folder should be assembly contigs or genome.
- `-sr`: the samples in test_sr folder should be short reads paired data.


## Step by Step PSVGT from pre-align BAM file
### demo data
```sh
(base) ➜  demo ls
0_tmp_hifi_5.gz.bam  0_tmp_hifi_5.gz.bam.bai  PSVGT1.0 Db-1_genome.fa.fai
 
```

### 1. SV sigaling from bam file
```sh
(base) ➜  demo python PSVGT1.0/PSV_Signal/0.Signal4bam_PSVGT.py -h
usage: SV signal extract from sam file [-h] -b BAM -o OUT [-m MIN] [-maq MAQ] -dtype DTYPE [-M MAX]
                                       [-fai FAI] [-msv MSV]
options:
  -h, --help    show this help message and exit

Signal Capture:
  -b BAM        sorted bam and index bam only (default: None)
  -o OUT        Output SV signal info file chromsomely (default: None)
  -m MIN        SV min length (default: 45)
  -maq MAQ      The min mapping quality (default: 50)
  -dtype DTYPE  The sequence type (ont, hifi, pb, cr, sr) (default: None)
  -M MAX        SV max length (default: 10000000)
  -fai FAI      Chromosome fai index file (default: None)
  -msv MSV      Detecting complex SVs (INV, DP, TRA, INS, DEL) from supplementary alignment (default:
                no)

python PSVGT1.0/PSV_Signal/0.Signal4bam_PSVGT.py -b 0_tmp_hifi_5.gz.bam -o 5X -m 40 -maq 30 -dtype hifi -fai Db-1_genome.fa.fai -msv yes

## out result ##
(base) ➜  demo python PSVGT1.0/PSV_Signal/0.Signal4bam_PSVGT.py -b 0_tmp_hifi_5.gz.bam -o 5X -m 40 -maq 30 -dtype hifi -fai ../../../../Db-1_ref/Db-1_genome.fa.fai -msv yes
**************************************** done SV searching ****************************************
cost time: 3.157372236251831

(base) ➜  demo ll
lrwxrwxrwx 1 lgb xinwang   59 Jan 28 20:16 0_tmp_hifi_5.gz.bam 
lrwxrwxrwx 1 lgb xinwang   63 Jan 28 20:16 0_tmp_hifi_5.gz.bam.bai
-rw-r--r-- 1 lgb xinwang 2.2M Jan 28 20:24 5X_Db-Chr1.record.txt
-rw-r--r-- 1 lgb xinwang 2.5M Jan 28 20:24 5X_Db-Chr1.record.txt.cov
-rw-r--r-- 1 lgb xinwang 228K Jan 28 20:24 5X_Db-Chr1.record.txt.suppAlign
-rw-r--r-- 1 lgb xinwang 1.5M Jan 28 20:24 5X_Db-Chr2.record.txt
-rw-r--r-- 1 lgb xinwang 1.7M Jan 28 20:24 5X_Db-Chr2.record.txt.cov
-rw-r--r-- 1 lgb xinwang 167K Jan 28 20:24 5X_Db-Chr2.record.txt.suppAlign
-rw-r--r-- 1 lgb xinwang 1.6M Jan 28 20:24 5X_Db-Chr3.record.txt
-rw-r--r-- 1 lgb xinwang 2.0M Jan 28 20:24 5X_Db-Chr3.record.txt.cov
-rw-r--r-- 1 lgb xinwang 267K Jan 28 20:24 5X_Db-Chr3.record.txt.suppAlign
-rw-r--r-- 1 lgb xinwang 1.6M Jan 28 20:24 5X_Db-Chr4.record.txt
-rw-r--r-- 1 lgb xinwang 1.6M Jan 28 20:24 5X_Db-Chr4.record.txt.cov
-rw-r--r-- 1 lgb xinwang 183K Jan 28 20:24 5X_Db-Chr4.record.txt.suppAlign
-rw-r--r-- 1 lgb xinwang 1.8M Jan 28 20:24 5X_Db-Chr5.record.txt
-rw-r--r-- 1 lgb xinwang 2.3M Jan 28 20:24 5X_Db-Chr5.record.txt.cov
-rw-r--r-- 1 lgb xinwang 230K Jan 28 20:24 5X_Db-Chr5.record.txt.suppAlign
```


### 2. SV Signal Cluster and Local Depth Based Filter:
```sh
python PSV_Signal/0.KLookCluster_LocalDepthPASS.py   
usage: signal filtering through support reads ratio [-h] -f RAW_SIGNAL [-s SHIFT] [-M MAX] -dtype DTYPE [--cov COVFILE] [--b BAM]
signal filtering through support reads ratio: error: the following arguments are required: -f, -dtype

options:
  -h, --help     show this help message and exit

Input File :
  -f RAW_SIGNAL  the raw sv signal record file from 0 step signalling (default: None)
  -s SHIFT       the distance shift of breakpoint to cluster the TRA/INV/DUP signal (default: 1000)
  -M MAX         the max SV length (default: 6868686)
  -dtype DTYPE   the sequencing type of samples (default: None)
  --cov COVFILE  Coverage File (default: None)
  --b BAM        the bam file of Individual (default: None)

python PSVGT1.0/PSV_Signal/0.KLookCluster_LocalDepthPASS.py  -f 5X_Db-Chr1.record.txt -s 800 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/PSV_Signal/0.KLookCluster_LocalDepthPASS.py  -f 5X_Db-Chr2.record.txt -s 800 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/PSV_Signal/0.KLookCluster_LocalDepthPASS.py  -f 5X_Db-Chr3.record.txt -s 800 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/PSV_Signal/0.KLookCluster_LocalDepthPASS.py  -f 5X_Db-Chr4.record.txt -s 800 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/PSV_Signal/0.KLookCluster_LocalDepthPASS.py  -f 5X_Db-Chr5.record.txt -s 800 -dtype hifi --b 0_tmp_hifi_5.gz.bam

### Signal files will be generated ###
5X_Db-Chr1.record.txt_DEL.signal  5X_Db-Chr2.record.txt_DUP.signal  5X_Db-Chr3.record.txt_INS.signal  5X_Db-Chr4.record.txt_INV.signal
5X_Db-Chr1.record.txt_DUP.signal  5X_Db-Chr2.record.txt_INS.signal  5X_Db-Chr3.record.txt_INV.signal  5X_Db-Chr4.record.txt_TRA.signal
5X_Db-Chr1.record.txt_INS.signal  5X_Db-Chr2.record.txt_INV.signal  5X_Db-Chr3.record.txt_TRA.signal  5X_Db-Chr5.record.txt_DEL.signal
5X_Db-Chr1.record.txt_INV.signal  5X_Db-Chr2.record.txt_TRA.signal  5X_Db-Chr4.record.txt_DEL.signal  5X_Db-Chr5.record.txt_DUP.signal
5X_Db-Chr1.record.txt_TRA.signal  5X_Db-Chr3.record.txt_DEL.signal  5X_Db-Chr4.record.txt_DUP.signal  5X_Db-Chr5.record.txt_INS.signal
5X_Db-Chr2.record.txt_DEL.signal  5X_Db-Chr3.record.txt_DUP.signal  5X_Db-Chr4.record.txt_INS.signal  5X_Db-Chr5.record.txt_INV.signal

```

### 3. Merge Individual Chromosome SV Signal
```sh
python PSV_Signal/1.ACCSV_Signal_Cluster.py -h
usage: signal filtering through support reads ratio [-h] -preffix PREFFIX -fai FAI [-s SHIFT] [-M MAX]
python PSV_Signal/1.ACCSV_Signal_Cluster.py -preffix 5X -fai Db-1_genome.fa.fai
``` 
Output File: 5X_Cluster_Record.txt

### 4. Merge all Individual SV Signals
```sh
python  PSVGT1.0/PSV_Signal/1.PSV_signal_cluster.py -h
usage: signal filtering through support reads ratio [-h] -d SV_DIR [-s SHIFT] [-M MAX]

options:
  -h, --help  show this help message and exit

Input file :
  -d SV_DIR   the PSVGT output directory (default: None)
  -s SHIFT    the distance of shifting the breakpoints (default: 30)
  -M MAX      the max SV length (default: 6868886)

(base) ➜  demo python  PSVGT1.0/PSV_signal/1.PSV_signal_cluster.py -d ./ -s 40
************************** chromsomes numer is 5 ***************************
The Db-Chr3 DEL data (805, 11) after clustering by shift:30 is 805
The Db-Chr4 DEL data (705, 11) after clustering by shift:30 is 705
The Db-Chr2 DEL data (728, 11) after clustering by shift:30 is 728
The Db-Chr5 DEL data (801, 11) after clustering by shift:30 is 801
The Db-Chr1 DEL data (1047, 11) after clustering by shift:30 is 1047
The Db-Chr4 INS data (670, 11) after clustering by shift:30 is 670

Here is only one sample, so the cluster results will not be changed
## Final Candidate SV signal file is PopSV_clustered_Record.txt
```
## 4. SV Genotyping
```sh
python PSVGT1.0/PSV_Genotyper/2.Pop_lrSVGT_V1.py -h
usage: the SVGT for long reads samples [-h] -i SV_INFO [-m MAQ] [-s SHIFT] -mapf MAPF -n ACC -o DIR

options:
  -h, --help  show this help message and exit

Input file :
  -i SV_INFO  sv call set (default: None)
  -m MAQ      the mini mapping quality of the contigs, range from 45-60, default is 50 (default: 50)
  -s SHIFT    breakpoints may shifting, here we shift 100bp as default (default: 100)
  -mapf MAPF  the mapping file of sample in sam or bam format (default: None)
  -n ACC      Accession name of the Individual (default: None)
  -o DIR      the output dir (default: None)

python PSVGT1.0/PSV_Genotyper/2.Pop_lrSVGT_V1.py -i PopSV_clustered_Record.txt -m 30 -s 100 -mapf 0_tmp_hifi_5.gz.bam -n 5X_hifi -o

## Final Genotype Table ##
cat 2_tmp_5X_hifi_genotype.txt |head                                                                          
#Target_name	Target_start	Target_end	Target_size	Query_size	5X_hifi	Total_Map_Reads	SV_support
Db-Chr1	10191632	10191633	-172	172	1/1	total_map_reads=8	INS_rate=1.0;INS
Db-Chr1	10075481	10075482	-4035	4035	1/1	total_map_reads=3	INS_rate=1.0;INS
Db-Chr1	10116536	10116537	-79	79	1/1	total_map_reads=7	INS_rate=1.0;INS
```

## PSVGT on genome assemblies datasets ##
### haploid assembly
```
## the pre-aligned and indexed file of genome should be in folder genome_bam 
python PSVGT1.0/PSVGT1.0.py -cr genome_bam -r ref.fa -msv yes -m 50 -o outfolder 
```
### diploid assemblies
```
cat diploid.info 
HG00171_clr_hap1.fasta.bam	HG00171_clr_hap2.fasta.bam	HG00171_clr
# two haplotype mapped bams or fastas file should be in genome_bam folder
python PSVGT1.0/PSVGT1.0.py -cr genome_bam -r ref.fa -msv yes -m 50 -o outfolder --diploid diploid.info
```
### polyploid assemblies ###
```
cat polyploid.info
C_hap1_genome.fasta.sorted.bam	C_hap2_genome.fasta.sorted.bam	C_hap3_genome.fasta.sorted.bam	C_hap4_genome.fasta.sorted.bam	Eig
python PSVGT1.0.py -cr genome_bam -r ref.fa -msv no -m 50 -o outfolder --polyploid polyploid.info
```

## 🧬 PSVGT ToolKits Commands
### A more detail usage will be updated sooner or later
 - [CAPSPop] Minor Allele Frequency CAPS Marker Developement 
```shell
## mpileup to get maf 0.05 SNPs site
samtools mpileup -b bam_lst.txt -q 55 -Q 30  | python PSVGT1.0/CapsPop/mpileup_stdin4popcasp.py > PopCaps_input.txt
## CAPS development
python PSVGT1.0/CapsPop/pop_maf0.05_caps.py PSVGT1.0/CapsPop/common_enzyme.list reference.fa PopCaps_input.txt Out_PopCaps_maf0.05.txt 300 
rm PopCaps_input.txt
```

 - [SVInDel_Anno]   - Based on gff to annotate the SV impact to GenePAV, Expression, IntronLost, GeneAS, Frameshift
```shell
## The PSVGT_all.vc2.SVInDel format 
head -n 5 out0129/PSVGT_all.vcf2.SVInDel        
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2_tmp_hifi_5.gz_genotype.txt
Db-Chr1	10038197	Db-Chr1:10038197-10038740_544	N	<DEL>	.	PASS	.	GT	1/1
Db-Chr1	10038841	Db-Chr1:10038841-10038842_58	N	<INS>	.	PASS	.	GT	1/1
Db-Chr1	10116536	Db-Chr1:10116536-10116537_79	N	<INS>	.	PASS	.	GT	1/1
Db-Chr1	10118152	Db-Chr1:10118152-10118663_512	N	<DEL>	.	PASS	.	GT	1/1

## The gff format mRNA or Gene Line has ID= and CDS Line has Parent=
Db-Chr1	TAIR10	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;
Db-Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1;
Db-Chr1	TAIR10	CDS	3760	3913	.	+	0	Parent=AT1G01010.1;
Db-Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1;
Db-Chr1	TAIR10	CDS	3996	4276	.	+	2	Parent=AT1G01010.1;
Db-Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1;
Db-Chr1	TAIR10	CDS	4486	4605	.	+	0	Parent=AT1G01010.1;
Db-Chr1	TAIR10	exon	5174	5326	.	+	.	Parent=AT1G01010.1;
Db-Chr1	TAIR10	CDS	5174	5326	.	+	0	Parent=AT1G01010.1;
Db-Chr1	TAIR10	exon	5439	5899	.	+	.	Parent=AT1G01010.1;
Db-Chr1	TAIR10	CDS	5439	5630	.	+	0	Parent=AT1G01010.1;
Db-Chr1	TAIR10	mRNA	5928	8737	.	-	.	ID=AT1G01020.1;
Db-Chr1	TAIR10	exon	5928	6263	.	-	.	Parent=AT1G01020.1;
Db-Chr1	TAIR10	exon	6437	7069	.	-	.	Parent=AT1G01020.1;

python SVInDel_Anno/SV_Features_Annotation.py -g test.gff3 -s  PSVGT_all.vcf2.SVInDel -m ID -c Parent -o SVInDels_Lead_Gene_Variant.txt
python SVInDel_Anno/SV_Features_Position.py test.gff3 PSVGT_all.vcf2.SVInDel_tmp.tab PSVInDel
```

 - [SVInDel_Primer] - High quality primers development for SVInDel.

```sh
head -n 6 PSVGT_all.vcf2.SVInDel 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2_tmp_hifi_5.gz_genotype.txt	2_tmp_sr_40X.contigs.fa.gz_genotype.txt	2_tmp_sr_40X.contigs.fa_genotype.txt
Db-Chr1	137034	Db-Chr1:137034-137160_127	N	<DEL>	.	PASS	.	GT	1/1	1/1	1/1
Db-Chr1	256076	Db-Chr1:256076-256296_221	N	<DEL>	.	PASS	.	GT	1/1	1/1	1/1
Db-Chr1	356036	Db-Chr1:356036-357213_1178	N	<DEL>	.	PASS	.	GT	1/1	1/1	1/1
Db-Chr1	359380	Db-Chr1:359380-359381_88	N	<INS>	.	PASS	.	GT	1/1	1/1	1/1
Db-Chr1	381761	Db-Chr1:381761-381762_256	N	<INS>	.	PASS	.	GT	1/1	1/1	1/1
Db-Chr1	399719	Db-Chr1:399719-399720_91	N	<INS>	.	PASS	.	GT	1/1	1/1	1/1

python SVInDel_Primer/vcf2primer.py PSVGT_all.vcf2.SVInDel Db-1_genome.fa 50 500 500 > PSVInDel_Primer4Pop.txt
```

## 🔎 Note !!!!!!! Information
The reference genome should not has thousands of chromosome, please remove the small fragment contigs from the genome  
Please emails me 13414960404@163.com
