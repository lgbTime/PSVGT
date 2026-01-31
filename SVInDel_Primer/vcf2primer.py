#!/usr/bin/env python3

import argparse
import primer3
import gzip
import re
from collections import deque
from gzip import BadGzipFile

################################
# FASTA reader
################################

def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()

    with gzip.open(fa, 'rb') as f:
        try:
            f.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False

    fin = gzip.open(fa, 'rt') if isgzip else open(fa, 'r')

    for line in fin:
        if line.startswith(">"):
            if geneID:
                bigFa[geneID] = ''.join(geneSeq)
            geneID = re.split(r'\s+', line[1:].strip())[0]
            geneSeq = deque()
        else:
            geneSeq.append(line.strip())

    if geneID:
        bigFa[geneID] = ''.join(geneSeq)

    fin.close()
    return bigFa


def is_valid_dna_sequence(seq):
    return all(c in "acgtACGTN" for c in seq)


################################
# Primer + PCR logic
################################

def my_primer(vcf_fields, seq, seqid,
              seq_args_left, seq_args_right,
              global_args, sv_size):

    resL = primer3.bindings.design_primers(seq_args_left, global_args)
    resR = primer3.bindings.design_primers(seq_args_right, global_args)

    if resL["PRIMER_LEFT_NUM_RETURNED"] == 0 or resR["PRIMER_RIGHT_NUM_RETURNED"] == 0:
        return

    # Best-ranked primers
    f_seq = resL["PRIMER_LEFT_0_SEQUENCE"]
    f_start = resL["PRIMER_LEFT_0"][0] + 1

    r_seq = resR["PRIMER_RIGHT_0_SEQUENCE"]
    r_start_local = resR["PRIMER_RIGHT_0"][0]
    r_len = resR["PRIMER_RIGHT_0"][1]

    # Map right primer to full-sequence coordinate
    right_offset = len(seq) - seq_args_right["SEQUENCE_INCLUDED_REGION"][1]
    #r_end = right_offset + r_start_local + r_len - 1
    r_end = right_offset + r_start_local + 1
    if r_end <= f_start:
        return

    RefPCRsize = r_end - f_start + 1
    QueryPCRsize = RefPCRsize - sv_size

    print(
        "\t".join(vcf_fields) + "\t" +
        "\t".join([
            seq,
            seqid,
            f"{f_seq}@{f_start}",
            f"{r_seq}@{r_end}",
            str(RefPCRsize),
            str(QueryPCRsize)
        ])
    )


################################
# Arguments
################################

parser = argparse.ArgumentParser(
    "MAF InDel Marker Analysis For PSVGT InDel VCF",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("vcf_file")
parser.add_argument("ref")
parser.add_argument("--min", default=50, type=int)
parser.add_argument("--max", default=600, type=int)
parser.add_argument("--frank", default=300, type=int)
parser.add_argument("--maf", default=0.01, type=float)

args = parser.parse_args()

global_args = {
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_OPT_SIZE': 21,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 53.0,
    'PRIMER_MAX_TM': 65.0,
    'PRIMER_MIN_GC': 40.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
}

################################
# Main
################################

fa = readfa2Dict(args.ref)

with open(args.vcf_file) as f:
    for line in f:
        line = line.rstrip()

        # Meta-header
        if line.startswith("##"):
            print(line)
            continue

        # Main header
        if line.startswith("#CHROM"):
            print(
                line +
                "\tRefSeq\tSeqID\tForwardPrimer\tReversePrimer\tRefPCRsize\tQueryPCRsize"
            )
            continue

        # Variant line
        fields = line.split("\t")

        SVID = fields[2]
        chrom = SVID.split(":")[0]
        startI = int(SVID.split(":")[1].split("-")[0])
        endI = int(SVID.split(":")[1].split("-")[1].split("_")[0])
        sv_size = int(SVID.split("_")[1])

        genotypes = fields[9:]
        alt_maf = genotypes.count("1/1") / len(genotypes)

        if not (args.maf < alt_maf < 1 - args.maf):
            continue
        if not (args.min < sv_size < args.max):
            continue
        if startI - args.frank < 0:
            continue

        ## convert INS sive ##
        if "INS" in line:
            sv_size = -sv_size
        seq = fa[chrom][startI - args.frank:endI + args.frank].upper()
        if not is_valid_dna_sequence(seq):
            continue

        pos = args.frank
        seqid = f"{chrom}:{startI - args.frank}-{endI + args.frank}"

        seq_args_left = {
            "SEQUENCE_ID": seqid + "_L",
            "SEQUENCE_TEMPLATE": seq[:pos],
            "SEQUENCE_INCLUDED_REGION": [0, pos]
        }

        seq_args_right = {
            "SEQUENCE_ID": seqid + "_R",
            "SEQUENCE_TEMPLATE": seq[-pos:],
            "SEQUENCE_INCLUDED_REGION": [0, pos]
        }

        my_primer(
            fields,
            seq,
            seqid,
            seq_args_left,
            seq_args_right,
            global_args,
            sv_size
        )

