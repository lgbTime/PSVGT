#!/usr/bin/env python3
"""
Design PCR primers flanking an InDel for agarose-gel genotyping.
Input: VCF with InDels (SVLEN/END in INFO, type in ALT), FASTA reference.
Output: best primer pair(s) per variant, with expected Ref/Alt product sizes.
"""

import argparse
import primer3
import gzip
import sys
from collections import deque
from gzip import BadGzipFile
import re

# ---------- reverse complement ----------
def reverse_complement(primer):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(primer))

# ---------- find primer positions ----------
def find_primer_positions(sequence, primer, strand):
    """Return start (forward) or end (reverse) coordinate on the full sequence."""
    if strand == "+":
        start_f = sequence.find(primer)
        if start_f == -1:
            return None
        return start_f
    else:
        rc = reverse_complement(primer)
        start_r = sequence.find(rc)
        if start_r == -1:
            return None
        end_r = start_r + len(primer) - 1
        return end_r

# ---------- FASTA reader ----------
def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa, 'rb') as f_in:
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
                            geneSeq = deque()
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
                            geneID = re.split(r'\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split(r'\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
        
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa

# ---------- core primer design & scoring ----------
def design_best_primer_pair(info, seq, seqid,
                            seq_args_left, seq_args_right,
                            global_args, sv_size, args):
    """Return dict with best primer pair or None."""
    left_res = primer3.bindings.design_primers(seq_args=seq_args_left, global_args=global_args)
    right_res = primer3.bindings.design_primers(seq_args=seq_args_right, global_args=global_args)

    # Collect left (forward) primers
    left_primers = []
    for i in range(left_res.get('PRIMER_LEFT_NUM_RETURNED', 0)):
        f_seq = left_res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        if not f_seq:
            continue
        f_tm = left_res.get(f"PRIMER_LEFT_{i}_TM", None)
        f_gc = left_res.get(f"PRIMER_LEFT_{i}_GC_PERCENT", None)
        pos = find_primer_positions(seq, f_seq, "+")
        if pos is None:
            continue
        left_primers.append({
            'seq': f_seq,
            'tm': f_tm,
            'gc': f_gc,
            'pos': pos
        })

    # Collect right (reverse) primers
    right_primers = []
    for i in range(right_res.get('PRIMER_RIGHT_NUM_RETURNED', 0)):
        r_seq = right_res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        if not r_seq:
            continue
        r_tm = right_res.get(f"PRIMER_RIGHT_{i}_TM", None)
        r_gc = right_res.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", None)
        pos = find_primer_positions(seq, r_seq, "-")
        if pos is None:
            continue
        right_primers.append({
            'seq': r_seq,
            'tm': r_tm,
            'gc': r_gc,
            'pos': pos
        })

    if len(left_primers) == 0 or len(right_primers) == 0:
        print(f"# Warning: No primers for {info[2]} (left {len(left_primers)}, right {len(right_primers)})",
              file=sys.stderr)
        return None

    # Score all left × right pairs
    best_pair = None
    best_score = -999999
    best_details = {}

    for left in left_primers:
        for right in right_primers:
            start_f = left['pos']
            end_r = right['pos']
            if end_r <= start_f:
                continue
            ref_size = end_r - start_f + 1
            alt_size = ref_size - sv_size  # sv_size pos for DEL, neg for INS

            # agarose gel filters
            if not (args.min_product <= ref_size <= args.max_product):
                continue
            if not (args.min_product <= alt_size <= args.max_product):
                continue
            if abs(sv_size) < args.min_size_diff:
                continue

            # compatibility scoring
            hetero = primer3.bindings.calc_heterodimer(left['seq'], right['seq'])
            hetero_penalty = max(0, hetero.tm - 45) * 2 if hetero.tm else 0
            tm_diff = abs(left['tm'] - right['tm']) if left['tm'] and right['tm'] else 10
            tm_penalty = max(0, tm_diff - 3) * 5
            gc_penalty = 0
            for gc in (left['gc'], right['gc']):
                if gc is not None:
                    if gc < 40:
                        gc_penalty += (40 - gc) * 0.5
                    elif gc > 60:
                        gc_penalty += (gc - 60) * 0.5

            score = 100 - hetero_penalty - tm_penalty - gc_penalty

            if score > best_score:
                best_score = score
                best_pair = (left, right)
                best_details = {
                    'ref_size': ref_size,
                    'alt_size': alt_size,
                    'hetero_tm': hetero.tm,
                    'hetero_dg': hetero.dg,
                    'tm_diff': tm_diff
                }

    if best_pair is None:
        return None

    return {
        'forward_seq': best_pair[0]['seq'],
        'reverse_seq': best_pair[1]['seq'],
        'forward_tm': best_pair[0]['tm'],
        'reverse_tm': best_pair[1]['tm'],
        'forward_gc': best_pair[0]['gc'],
        'reverse_gc': best_pair[1]['gc'],
        'ref_size': best_details['ref_size'],
        'alt_size': best_details['alt_size'],
        'hetero_tm': best_details['hetero_tm'],
        'tm_diff': best_details['tm_diff'],
        'score': best_score
    }

# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(
        description="Agarose-gel InDel marker primer designer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf_file", help="Input VCF with InDels")
    parser.add_argument("ref", help="Reference genome FASTA (can be gzipped)")
    parser.add_argument("-o", "--output", type=str, default="designed_markers_output.txt", 
                        help="Path to save the clean tab-delimited output results file")
    
    parser.add_argument("--min", type=int, default=50, help="Minimum InDel size (bp)")
    parser.add_argument("--max", type=int, default=600, help="Maximum InDel size (bp)")
    parser.add_argument("--frank", type=int, default=300, help="Flanking region size (bp)")
    parser.add_argument("--maf", type=float, default=0.01, help="Minor allele frequency threshold (multi-sample)")
    parser.add_argument("--min_product", type=int, default=100, help="Minimum total Ref product size (bp)")
    parser.add_argument("--max_product", type=int, default=1000, help="Maximum total Ref product size (bp)")
    parser.add_argument("--min_size_diff", type=int, default=30, help="Minimum size difference between alleles (bp)")
    parser.add_argument("--single_sample", action="store_true", help="Ignore MAF filter (single sample mode)")
    parser.add_argument("--best_only", action="store_true", help="Output only the best primer pair (no verbose list)")

    args = parser.parse_args()

    global_args = {
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

    fa = readfa2Dict(args.ref)

    # Open the dedicated output file to write clean results
    print(f"Initializing run. Output will be saved cleanly to: {args.output}")
    f_out = open(args.output, "w")

    # Write header directly to file
    out_fields = [
        "CHROM", "POS", "ID", "SV_TYPE", "SV_SIZE",
        "MAF", "REF_SIZE", "ALT_SIZE",
        "F_PRIMER", "R_PRIMER",
        "F_TM", "R_TM", "F_GC", "R_GC",
        "HETERO_TM", "TM_DIFF", "SCORE"
    ]
    f_out.write("\t".join(out_fields) + "\n")

    success_count = 0
    total_variants = 0

    with open(args.vcf_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("##") or line.startswith("#CHROM"):
                continue

            fields = line.split("\t")
            chrom = fields[0]
            pos = fields[1]              
            startI = int(pos)            
            svid = fields[2]
            ref_seq_vcf = fields[3]
            alt_seq_vcf = fields[4]
            info = fields[7]
            genotypes = fields[9:]

            total_variants += 1

            if "<DEL>" in alt_seq_vcf or "DEL" in alt_seq_vcf:
                sv_type = "DEL"
            elif "<INS>" in alt_seq_vcf or "INS" in alt_seq_vcf:
                sv_type = "INS"
            else:
                continue

            endI = startI
            if "END=" in info:
                try:
                    endI = int(info.split("END=")[1].split(";")[0])
                except (IndexError, ValueError):
                    pass

            sv_size = None
            if "SVLEN=" in info:
                try:
                    sv_size = abs(int(info.split("SVLEN=")[1].split(";")[0]))
                except (IndexError, ValueError):
                    pass

            if sv_size is None:
                sv_size = abs(len(ref_seq_vcf) - len(alt_seq_vcf))
                if sv_size == 0 and endI != startI:
                    sv_size = abs(endI - startI)

            if sv_size is None:
                continue

            if sv_type == "INS":
                sv_size = -sv_size   

            if not args.single_sample and len(genotypes) > 1:
                alt_maf = genotypes.count("1/1") / len(genotypes)
                if alt_maf <= args.maf or alt_maf >= 1 - args.maf:
                    continue

            if abs(sv_size) < args.min or abs(sv_size) > args.max:
                continue

            if startI - args.frank < 0 or chrom not in fa:
                continue

            seq = fa[chrom][startI - args.frank : endI + args.frank].upper()
            if 'N' in seq or not seq:
                continue

            flank = args.frank   
            seqid = f"{chrom}:{startI - flank}-{endI + flank}"

            left_seq = seq[:flank]
            right_seq = seq[len(seq)-flank:]

            seq_args_left = {
                'SEQUENCE_ID': seqid + '_left',
                'SEQUENCE_TEMPLATE': left_seq,
                'SEQUENCE_INCLUDED_REGION': [0, len(left_seq)]
            }
            seq_args_right = {
                'SEQUENCE_ID': seqid + '_right',
                'SEQUENCE_TEMPLATE': right_seq,
                'SEQUENCE_INCLUDED_REGION': [0, len(right_seq)]
            }

            best = design_best_primer_pair(
                info=fields,
                seq=seq,
                seqid=seqid,
                seq_args_left=seq_args_left,
                seq_args_right=seq_args_right,
                global_args=global_args,
                sv_size=sv_size,
                args=args
            )

            if best is None:
                continue

            maf_str = "N/A"
            if len(genotypes) > 1:
                maf_str = f"{genotypes.count('1/1') / len(genotypes):.2f}"

            out = [
                chrom, pos, svid, sv_type, str(abs(sv_size)),
                maf_str,
                str(best['ref_size']), str(best['alt_size']),
                best['forward_seq'], best['reverse_seq'],
                f"{best['forward_tm']:.1f}" if best['forward_tm'] else "NA",
                f"{best['reverse_tm']:.1f}" if best['reverse_tm'] else "NA",
                f"{best['forward_gc']:.1f}" if best['forward_gc'] else "NA",
                f"{best['reverse_gc']:.1f}" if best['reverse_gc'] else "NA",
                f"{best['hetero_tm']:.1f}" if best['hetero_tm'] else "NA",
                f"{best['tm_diff']:.1f}",
                f"{best['score']:.1f}",
                f"{seq}"
            ]
            
            # Write row directly to output file
            f_out.write("\t".join(out) + "\n")
            success_count += 1

    f_out.close()
    print(f"Finished processing. Successfully designed assays for {success_count}/{total_variants} checked variants.")

if __name__ == "__main__":
    main()
