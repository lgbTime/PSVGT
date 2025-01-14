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
