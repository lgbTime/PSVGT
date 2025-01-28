![PSVGT](https://github.com/lgbTime/PSVGT/logo.png)  
    ____                 ______    ______________
   / __ \____  ____     / ___/ |  / / ____/_  __/
  / /_/ / __ \/ __ \    \__ \| | / / / __  / /   
 / ____/ /_/ / /_/ /   ___/ /| |/ / /_/ / / /    
/_/    \____/ .___/   /____/ |___/\____/ /_/     
           /_/                              
## PSVGT is a versatile program for SV detection and genotyping at population scale, supporting all sequence length (>150+ pbs) and depth free, providing CAPS & InDel development Modules and SV annotation function
The short reads samples will be de novo into contigs, while assemble samples and long reads sequencing data could be directory used in the SV detection and Genotyping.
## 💻 Installation
PSVGT required **megahit**, **minimap2**, **samtools**, and has been tested at python3.11
the  folow python dependency is required.
```sh
pip install primer3-py
pip install pyfiglet
pip install pandas 
pip install tqdm
pip install pysam
pip install intervaltree
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


## Step by Step PSVGT from aready mapped BAM file
### demo data
```sh
(base) ➜  demo ls
0_tmp_hifi_5.gz.bam  0_tmp_hifi_5.gz.bam.bai  PSVGT1.0 Db-1_genome.fa.fai
 
```

### 1. SV sigaling from bam file
```sh
(base) ➜  demo python PSVGT1.0/SV_Genotyper/0.Signal4bam_PSVGT.py -h
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

python PSVGT1.0/SV_Genotyper/0.Signal4bam_PSVGT.py -b 0_tmp_hifi_5.gz.bam -o 5X -m 40 -maq 30 -dtype hifi -fai Db-1_genome.fa.fai -msv yes

## out result ##
(base) ➜  demo python PSVGT1.0/SV_Genotyper/0.Signal4bam_PSVGT.py -b 0_tmp_hifi_5.gz.bam -o 5X -m 40 -maq 30 -dtype hifi -fai ../../../../Db-1_ref/Db-1_genome.fa.fai -msv yes
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
python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -h
usage: signal filtering through support reads ratio [-h] -f RAW_SIGNAL [-s SHIFT] [-M MAX] [-dtype DTYPE] [-csv CSV] [--cov COVFILE] [--b BAM]

options:
  -h, --help     show this help message and exit

Input File :
  -f RAW_SIGNAL  the raw sv signal record file from 0 step signalling (default: None)
  -s SHIFT       the distance shift of breakpoint to cluster the TRA/INV/DUP signal (default: 1000)
  -M MAX         the max SV length (default: 6868686)
  -dtype DTYPE   the sequencing type of samples (default: None)
  -csv CSV       the paramter to filter the sv signal / local_total < 0.2 (default: 0.2)
  --cov COVFILE  Coverage File (default: None)
  --b BAM        the bam file of Individual (default: None)

python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f 5X_Db-Chr1.record.txt -s 1000 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f 5X_Db-Chr2.record.txt -s 1000 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f 5X_Db-Chr3.record.txt -s 1000 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f 5X_Db-Chr4.record.txt -s 1000 -dtype hifi --b 0_tmp_hifi_5.gz.bam
python PSVGT1.0/SV_Genotyper/0.SignalCluster_LocalDepthFil.py -f 5X_Db-Chr5.record.txt -s 1000 -dtype hifi --b 0_tmp_hifi_5.gz.bam

### Signal files will be generated ###
5X_Db-Chr1.record.txt_DEL.signal  5X_Db-Chr2.record.txt_DUP.signal  5X_Db-Chr3.record.txt_INS.signal  5X_Db-Chr4.record.txt_INV.signal
5X_Db-Chr1.record.txt_DUP.signal  5X_Db-Chr2.record.txt_INS.signal  5X_Db-Chr3.record.txt_INV.signal  5X_Db-Chr4.record.txt_TRA.signal
5X_Db-Chr1.record.txt_INS.signal  5X_Db-Chr2.record.txt_INV.signal  5X_Db-Chr3.record.txt_TRA.signal  5X_Db-Chr5.record.txt_DEL.signal
5X_Db-Chr1.record.txt_INV.signal  5X_Db-Chr2.record.txt_TRA.signal  5X_Db-Chr4.record.txt_DEL.signal  5X_Db-Chr5.record.txt_DUP.signal
5X_Db-Chr1.record.txt_TRA.signal  5X_Db-Chr3.record.txt_DEL.signal  5X_Db-Chr4.record.txt_DUP.signal  5X_Db-Chr5.record.txt_INS.signal
5X_Db-Chr2.record.txt_DEL.signal  5X_Db-Chr3.record.txt_DUP.signal  5X_Db-Chr4.record.txt_INS.signal  5X_Db-Chr5.record.txt_INV.signal

```
### 3. Merge All Signal at Population Scale 
```sh
python  PSVGT1.0/SV_Genotyper/1.PSV_signal_cluster.py -h
usage: signal filtering through support reads ratio [-h] -d SV_DIR [-s SHIFT] [-M MAX]

options:
  -h, --help  show this help message and exit

Input file :
  -d SV_DIR   the PSVGT output directory (default: None)
  -s SHIFT    the distance of shifting the breakpoints (default: 30)
  -M MAX      the max SV length (default: 6868886)

(base) ➜  demo python  PSVGT1.0/SV_Genotyper/1.PSV_signal_cluster.py -d ./ -s 30
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
python PSVGT1.0/SV_Genotyper/2.Pop_lrSVGT_V1.py -h
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

python PSVGT1.0/SV_Genotyper/2.Pop_lrSVGT_V1.py -i PopSV_clustered_Record.txt -m 30 -s 100 -mapf 0_tmp_hifi_5.gz.bam -n 5X_hifi -o

## Final Genotype File ##
cat 2_tmp_5X_hifi_genotype.txt |head                                                                          
#Target_name	Target_start	Target_end	Target_size	Query_size	5X_hifi	Total_Map_Reads	SV_support
Db-Chr1	10191632	10191633	-172	172	1/1	total_map_reads=8	INS_rate=1.0;INS
Db-Chr1	10075481	10075482	-4035	4035	1/1	total_map_reads=3	INS_rate=1.0;INS
Db-Chr1	10116536	10116537	-79	79	1/1	total_map_reads=7	INS_rate=1.0;INS
```


## 🧬 PSVGT Commands
## The usage will be update sooner or later
 - [CAPSPop] Minor Allele Frequency CAPS Marker developement 
 - [SVInDel_Anno]   - Based on gff to annotate the SV impact to GenePAV, Expression, IntronLost, GeneAS, Frameshift
 - [SVInDel_Primer] - High quality primers development for SVInDel.
## 🔎 More Information
Please emails me 13414960404@163.com
