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
parser.add_argument("--frank", help="the exten size from the target region", default=600, type=int)
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
