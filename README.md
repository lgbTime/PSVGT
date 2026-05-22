## 🧬 PSVGT: Population-Scale Structural Variant Calling and Genotyping Computation Framework

**PSVGT** is a versatile toolkit for **structural variant (SV) detection and genotyping at population scale**.It supports sequencing reads of both long/short platform, and assemblies. the performance has been tested on:

- Haplotype assemblies
- Diploid and tetraploid assemblies
- Allodiploid assemblies
- Common genome assemblies
- Long-read data (HiFi, ONT, CLR)
- Illumina short reads

PSVGT also provides integrated modules for:

- **CAPS marker development**
- **InDel primer design**
- **SV functional annotation**

Short-read samples are assembled into contigs before SV detection, while haplotype assemblies and long-read data can be directly used for SV detection and genotyping.

---

## 💻 Installation

### Requirements

- megahit
- minimap2
- samtools
- Python ≥ 3.11

### Install

```bash
git clone https://github.com/lgbTime/PSVGT.git
pip install primer3-py pyfiglet pandas tqdm pysam intervaltree
```

---

## 🔨 General Usage

```bash
usage: PSVGT1.0.py [-h] [-o OUTDIR] [-sr SRDIR] [-hifi HIFIDIR] [-ont ONTDIR]
                   [-pb PBDIR] [-cr CRDIR] [-diploid DIPLOID] [-polyploid POLYPLOID]
                   [-r REFGENOME] [--num_hap NUM_HAP] 
                   [-msv MSV_MODE]
```



> **Note:** `-r / --refGenome` is required.

### Polyploid SV Detection

Set `--num_hap` according to ploidy:

* Diploid: `--num_hap 2`
* Tetraploid: `--num_hap 4`

### SV Genotyping

PSVGT performs **forced genotyping**, using candidate SVs, alignments, CIGAR strings, breakpoints, and spanning reads to infer genotypes (0/0, 0/1, 1/1).

Recommended long-read parameters:

* Tetraploid long reads: `lr_homo_rate` in range **0.75–0.95**
* Assembly long reads: `lr_homo_rate` in range **0.75–0.95**

---

## ⏩ One-Step PSVGT (Quick Start)

This is the simplest way to run PSVGT on mixed sequencing data.

```bash
python PSVGT1.0/PSVGT1.0.py \
    -hifi test_hifi \
    -ont test_ont \
    -pb test_pb \
    -sr test_sr \
    -cr test_cr \
    -r Db-1_ref/Db-1_genome.fa \
    -b yes \
    -msv yes \
    -o out4PSVGT
```

### Input Folder Requirements

| Folder                         | Content Type                  |
| ------------------------------ | ----------------------------- |
| test_hifi / test_ont / test_pb | gzipped FASTQ or mapped BAM   |
| test_cr                        | Assembly FASTA (.fa / .fasta) |
| test_sr                        | Paired-end gzipped FASTQ      |

### Key Parameters

* `-r` Reference genome
* `-o` Output directory
* `-msv yes` Detect all SV types (not only InDels)
* `-cr` Assembly input
* `-sr` Short-read input

---

## 📊 Benchmark Performance

### Highest sensitivity and genotype accuracy on long-read benchmark

![downsample](./img/downsample.png)

### Performance on simulated long-read benchmark

![simulated](./img/simulated_benchmark.png)

---

## 🧩 Step-by-Step: Starting from Pre-aligned BAM

### 1. SV Signal Extraction

```bash
python PSVGT1.0/PSV_Signal/0.Signal4bam_PSVGT.py \
    -b 0_tmp_hifi_5.gz.bam \
    -o 5X \
    -m 40 \
    -maq 30 \
    -dtype hifi \
    -fai Db-1_genome.fa.fai \
    -msv yes
```

---

### 2. Signal Clustering and Depth Filtering

```bash
python PSVGT1.0/PSV_Signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive.py \
    -f 5X_Db-Chr1.record.txt \
    -s 800 \
    -dtype hifi \
    --b 0_tmp_hifi_5.gz.bam
```

---

### 3. Merge Chromosome-Level Signals

```bash
python PSVGT1.0/PSV_Signal/1.ACCSV_Signal_Cluster.py \
    -preffix 5X \
    -fai Db-1_genome.fa.fai
```

Output:

```
5X_Cluster_Record.txt
```

---

### 4. Merge Population SV Signals

```bash
python PSVGT1.0/PSV_signal/1.PSV_signal_cluster.py -d ./ -s 40
```

Final candidate SV file:

```
PopSV_clustered_Record.txt
```

---

## 🧬 SV Genotyping

```bash
python PSVGT1.0/PSV_Genotyper/2.Pop_lrSVGT_V1.py \
    -i PopSV_clustered_Record.txt \
    -m 30 \
    -s 100 \
    -mapf 0_tmp_hifi_5.gz.bam \
    -n 5X_hifi \
    -o output_dir
```

---

## 🧬 For User Interest SV Genotyping

To genotyping your interest SV sets, you have to have the mapping bam files and provie the SV table as below:
**1st columns is chr, 2nd is position, 3rd is SV size or translocation chrx, 4th is SV type or position of chrx**

```sh
cat sv_table.txt
Chr1    10025486    2129    DEL
Chr2    12548600    1290    DEL
Chr1    10041503    1140    INS
Chr1    10141999    5440    INS
Chr2    61159250    65543   DUP
Chr2    92422310    1223    DUP
Chr1    14132000    1656    INV
Chr5    19413200    2676    INV
Chr2    5100000 Chr1    155500
Chr3    2000    Chr4    3030320
```

**Convert SV Table to PSVGT Input Sigals**

```sh
python PSVGT1.0/PSV_Signal/SVinfo2PSVGT_Candidate.py -i sv_table.txt -o PSV.signal.txt

# Genotyping your SV in the  samples of long reads or genome mapping
python PSVGT1.0/PSV_Genotyper/2.Pop_lrSVGT_V1.py \
    -i PSV.signal.txt \
    -m 30 -s 100 \
    -mapf 0_tmp_hifi_5.gz.bam \
    -n 5X_hifi -o output_dir

# Genotyping your SV in the  samples of short reads mapping
python PSVGT1.0/PSV_Genotyper/2.Pop_srSVGT_V1.py  \
    -i PSV.signal.txt \
    -m 30 -s 30 \
    -mapf samplexx_illumina_short_reads.bam \
    -n samplexx -o output_dir
```

---

## 🧬 PSVGT on Genome Assemblies

#### 🧬 Haploid Assembly

```bash
python PSVGT1.0/PSVGT1.0.py -cr genome_bam -r ref.fa -msv yes -m 50 -o outfolder
```

#### 🧬🧬 Diploid Assemblies

```bash
cat diploid.info
HG00171_clr_hap1.fasta.bam HG00171_clr_hap2.fasta.bam HG00171_clr

python PSVGT1.0/PSVGT1.0.py -cr genome_bam -r ref.fa -msv yes -m 50 -o outfolder --diploid diploid.info
```

#### 🧬🧬🧬🧬 Tetraploid Assemblies

```bash
cat polyploid.info
C_hap1.bam C_hap2.bam C_hap3.bam C_hap4.bam SampleName

python PSVGT1.0.py -cr genome_bam -r ref.fa -msv no -m 50 -o outfolder --polyploid polyploid.info --num_hap 4
```

---

## 🧬 PSVGT on Complex Rearrangement Genomics

#### For Geomme Assemblies (Inversion host nested SV)

###### PSVGT using multi-round recalling on a modified genome based on previous step found inversion.

```
usage: genome_nsv_calling.py [-h] -r REF -q QUERY [--psvgt_dir PSVGT_DIR]

Recursive PSVGT Inversion Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Initial reference FASTA file
  -q QUERY, --query QUERY
                        Input query folder (containing genome/reads)
  --psvgt_dir PSVGT_DIR
                        Base directory for PSVGT scripts

## example run in yeast genome
python ~/bin/PSVGT1.0/NSV_Genome/genome_nsv_calling.py -r CBS432.genome.fa  -q yeast_asm/N44 --psvgt_dir ~/bin/PSVGT1.0
```

#### For long-read datasets calling nested SV such as INV host INS/DEL/DUP, DUP host DEL/INV.

###### Callbreak of PSVGT can be used to capture low density CIGAR nested SV based on breakpoint counts

```
## a human genomics Nested SV calling case
for i in {1..22}; do python PSVGT1.0/PSV_signal/0.KLOOK_Cluster_Flexible_Window_Break_Depth_Adaptive_callbreak.py -dtype hifi -f 0_tmp_Nested_hifi.fastq.gz_${i}.record.txt --callbreak yes --b 0_tmp_Nested_hifi.fastq.gz.bam  --nreads 4 >log${i} &
done
```

---

## 🧰 PSVGT Toolkits

### CAPS Marker Development

```bash
samtools mpileup -b bam_lst.txt -q 55 -Q 30 \
 | python PSVGT1.0/CapsPop/mpileup_stdin4popcasp.py > PopCaps_input.txt

python PSVGT1.0/CapsPop/pop_maf0.05_caps.py \
    PSVGT1.0/CapsPop/common_enzyme.list \
    reference.fa \
    PopCaps_input.txt \
    Out_PopCaps_maf0.05.txt \
    300
```

---

### SV Functional Annotation

```bash
python SVInDel_Anno/SV_Features_Annotation.py \
    -g test.gff3 \
    -s PSVGT_all.vcf2.SVInDel \
    -m ID \
    -c Parent \
    -o SVInDels_Lead_Gene_Variant.txt
```

---

### SVInDel Primer Design

#### This module is to develop high quality markers for agarose gel electrophoresis experiment and marker association selection.

###### it will report the best pair primers, marker MAF, GC, TM, and scores.

##### Version 1

```bash
python vcf2primer_v1.py -h

usage: MAF InDel Marker Analysis For PSVGT InDel VCF [-h] [--min MIN] [--max MAX] [--frank FRANK] [--maf MAF] vcf_file ref
positional arguments:
  vcf_file       the vcf file from PSVGT
  ref            the reference genome to extract sequence design primers

options:
  -h, --help     show this help message and exit
  --min MIN      the min size of InDel (default: 50)
  --max MAX      the max size of InDel (default: 600)
  --frank FRANK  the exten size from the target region (default: 300)
  --maf MAF      the minor allele frequency setting for the InDel Primer (default: 0.01)

## example run v1
python vcf2primer_v1.py PSVGT_all.vcf2.SVInDel ref_Col-0-CEN/Col-CEN_v1.2_Chr.fasta | head -n 3         
#CHROM	POS	ID	REF	ALT	2_tmp_Db-1.fa_genotype.txt	MAF	Seqid	Seq.	Forward primers	Reverse primers	MiniPCR Product Primer	Ref PCR size	Query PCR size
Chr1	10038075	Chr1:10038075-10038076_106	N	<INS>	1/1	1.0	Chr1:10037775-10038376	AGCT**	Forward primers: {'TACCATTCTTCTGTGCACCGG_start_at_9', 'TCTTCTGTGCACCGGAAATCT_start_at_15'}	Reverse primers: {'TAATCGAAAGGTGACGCGACG_star_at_544', 'TTAATCGAAAGGTGACGCGAC_star_at_545', 'CGGAGAACAAACGACGGTGAA_star_at_525'}	MiniPCR_Pair_Primer:('TCTTCTGTGCACCGGAAATCT_start_at_15', 'CGGAGAACAAACGACGGTGAA_star_at_525')  510 616
Chr1	10042709	Chr1:10042709-10043267_559	N	<DEL>	1/1	1.0	Chr1:10042409-10043567	AGCT***	Forward primers: {'GTCCCCACGTCTTGTTGAGAT_start_at_60', 'ATTGGTGAGCATGTTGAACGC_start_at_79'}	Reverse primers: {'CTCCGCTCAAACCGTCAATTG_star_at_1006', 'CCGCTCAAACCGTCAATTGTT_star_at_1004', 'GAAAGATCTTCCCACCTCCGC_star_at_1021'}	MiniPCR_Pair_Primer:('ATTGGTGAGCATGTTGAACGC_start_at_79', 'CCGCTCAAACCGTCAATTGTT_star_at_1004') 925 366
```

##### Version 2

```
python vcf2primer_v2.py -h
usage: vcf2primer_v2.py [-h] [--min MIN] [--max MAX] [--frank FRANK] [--maf MAF] [--min_product MIN_PRODUCT] [--max_product MAX_PRODUCT] [--min_size_diff MIN_SIZE_DIFF]
                        [--single_sample] [--best_only]
                        vcf_file ref

Agarose-gel InDel marker primer designer
positional arguments:
  vcf_file              Input VCF with InDels
  ref                   Reference genome FASTA (can be gzipped)

options:
  -h, --help            show this help message and exit
  --min MIN             Minimum InDel size (bp) (default: 50)
  --max MAX             Maximum InDel size (bp) (default: 600)
  --frank FRANK         Flanking region size (bp) (default: 300)
  --maf MAF             Minor allele frequency threshold (multi-sample) (default: 0.01)
  --min_product MIN_PRODUCT
                        Minimum total Ref product size (bp) (default: 100)
  --max_product MAX_PRODUCT
                        Maximum total Ref product size (bp) (default: 1000)
  --min_size_diff MIN_SIZE_DIFF
                        Minimum size difference between alleles (bp) (default: 30)
  --single_sample       Ignore MAF filter (single sample mode) (default: False)
  --best_only           Output only the best primer pair (no verbose list) (default: False)

## example run 
python vcf2primer_v2.py out_cr/PSVGT_all.vcf2.SVInDel Col-CEN_v1.2_Chr.fasta  --min 50 --max 600 --min_product 100 --max_product 1000 --min_size_diff 30 --frank 300 --best_only 

############################ out example #######################
----------------------------------------------------------
CHROM	POS	ID	SV_TYPE	SV_SIZE	MAF	REF_SIZE	ALT_SIZE	F_PRIMER	R_PRIMER	F_TM	R_TM	F_GC	R_GC	HETERO_TM	TM_DIFF	SCORE
Chr1	10038075	Chr1:10038075-10038076_106	INS	106	0.75	537	643	TACCATTCTTCTGTGCACCGG	TTAATCGAAAGGTGACGCGAC	60.3	59.0	52.4	47.-9.5	1.3	100.0
Chr1	10042709	Chr1:10042709-10043268_559	DEL	559	0.38	945	386	GTCCCCACGTCTTGTTGAGAT	CCGCTCAAACCGTCAATTGTT	60.0	60.0	52.4	47.-9.0	0.0	100.0
```

---

## ⚠️ Notes

* The reference genome should **not contain thousands of small contigs**.
  Remove very small scaffolds before running PSVGT.

---

## 📚 Preprint

This project is described in the following preprint:

**Guangbao Luo**, Li Xiao, Zhangjun Fei, Xin Wang
*A sensitive and accurate framework for population-scale structural variant discovery and genotyping across sequence types.*
**bioRxiv 2026.01.10.698766v2**
[https://doi.org/10.1101/2026.01.10.698766v2](https://doi.org/10.1101/2026.01.10.698766v2)

https://www.biorxiv.org/content/10.64898/2026.01.10.698766v2
------------------------------------------------------------

## 📧 Contact

For questions or support:
**Email:** [13414960404@163.com]

