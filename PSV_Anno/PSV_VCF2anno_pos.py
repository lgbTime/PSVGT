#!/usr/bin/env python3
"""
Process GFF/GTF structural layers and correlate them against VCF Structural Variants.
Pools features dynamically by chromosome to optimize multi-interval lookups.
"""

import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
import gzip
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Process GFF/GTF annotations against VCF SV files')
    parser.add_argument('gff_file', type=str, help='Path to the GFF or GTF file')
    parser.add_argument('vcf_file', type=str, help='Path to the input VCF/VCF.GZ file containing SVs')
    parser.add_argument('output_file', type=str, help='Output base prefix name')
    return parser.parse_args()

def load_gff(gff_file):
    print("Loading genomic annotations...")
    gff = pd.read_csv(gff_file, header=None, comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = sorted(list(set(gff[0].tolist())))
    return gff, chroms

def extract_gene_cds(gff, gff_file):
    # Parse GFF structures
    if "gff" in str(gff_file).lower():
        gene = gff[gff[2] == "gene"].copy()
        cds = gff[gff[2] == "CDS"].copy()
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'ID=([^;]+)')
        cds["ID"] = cds[8].str.extract(r'Parent=([^;]+)')
    # Parse GTF structures
    elif "gtf" in str(gff_file).lower():
        gene = gff[gff[2] == "transcript"].copy()
        cds = gff[gff[2] == "exon"].copy()
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'gene_name "([^"]+)"')
        if gene["ID"].isna().all():
            gene["ID"] = gene[8].str.extract(r'gene_name ([^;]+)')
        cds["ID"] = cds[8].str.extract(r'gene_name "([^"]+)"')
        if cds["ID"].isna().all():
            cds["ID"] = cds[8].str.extract(r'gene_name ([^;]+)')
    return gene, cds

def get_3k_promoter(gene):
    pro3k = gene.copy()
    # Optimized vectorized assignment logic to prevent slow iloc iterations
    plus_strand = pro3k[6] == "+"
    minus_strand = pro3k[6] != "+"

    # Forward strand coordinates
    pro3k.loc[plus_strand, 4] = pro3k.loc[plus_strand, 3] - 1
    pro3k.loc[plus_strand, 3] = pro3k.loc[plus_strand, 3] - 3000

    # Reverse strand coordinates
    pro3k.loc[minus_strand, 3] = pro3k.loc[minus_strand, 4] + 1
    pro3k.loc[minus_strand, 4] = pro3k.loc[minus_strand, 4] + 3000

    pro3k[2] = "promoter"
    # Bound check coordinates so they don't drop below 1
    pro3k.loc[pro3k[3] < 1, 3] = 1
    return pro3k

def load_vcf_variants(vcf_file):
    """Parses structural variant rows from standardized VCF and VCF.GZ formats."""
    print("Parsing VCF structural variant indices...")
    variants = []
    
    # Handle compressed or raw file streams seamlessly
    is_gzip = vcf_file.endswith(".gz")
    f = gzip.open(vcf_file, "rt") if is_gzip else open(vcf_file, "r")
    
    try:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            sv_id = fields[2]
            ref_allele = fields[3]
            alt = fields[4]
            if "TRA" in alt:
                continue
            info = fields[7]
            
            # Extract standard structural Variant End coordinate
            end_pos = None
            if "END=" in info:
                try:
                    end_pos = int(info.split("END=")[1].split(";")[0])
                except (IndexError, ValueError):
                    pass
            
            # Fallback handling if no structural END tags exist
            if end_pos is None:
                end_pos = pos + len(ref_allele)
            
            variants.append({
                "CHROM": chrom,
                "START": pos,
                "END": end_pos,
                "ID": sv_id,
                "VCF_LINE": fields,
                "gene": [],
                "CDS": [],
                "3k_promoter": []
            })
    finally:
        f.close()
        
    print(f"Loaded {len(variants)} variants from VCF.")
    return variants

def main():
    args = parse_args()
    
    # Load and map genomic structural fields
    gff, chroms = load_gff(args.gff_file)
    gene, cds = extract_gene_cds(gff, args.gff_file)
    pro3k = get_3k_promoter(gene)
    
    # Filter working data structures down to essential indexed columns
    gene = gene[[0, 3, 4, "ID"]].dropna()
    cds = cds[[0, 3, 4, "ID"]].dropna()
    pro3k = pro3k[[0, 3, 4, "ID"]].dropna()
    
    # Load target structural variants from VCF
    variants = load_vcf_variants(args.vcf_file)
    
    # Pool structural variants by chromosome for high performance
    vcf_pools = {}
    for var in variants:
        vcf_pools.setdefault(var["CHROM"], []).append(var)
        
    outeach = open(f'{args.output_file}_each_gene.txt', "w")
    
    print("Beginning pooled intersection calculations across chromosomes...")
    # Iterate across structural features chromosome-by-chromosome
    for chrom in tqdm(chroms):
        if chrom not in vcf_pools:
            continue
            
        current_vcf_batch = vcf_pools[chrom]
        
        # 1. Build Interval Trees for the active chromosome
        gene_tree = IntervalTree()
        for _, row in gene[gene[0] == chrom].iterrows():
            if row[3] < row[4]:
                gene_tree[row[3]:row[4]] = str(row["ID"]).replace('"', '')
                
        cds_tree = IntervalTree()
        for _, row in cds[cds[0] == chrom].iterrows():
            if row[3] < row[4]:
                cds_tree[row[3]:row[4]] = str(row["ID"]).replace('"', '')
                
        pro3k_tree = IntervalTree()
        for _, row in pro3k[pro3k[0] == chrom].iterrows():
            if row[3] < row[4]:
                pro3k_tree[row[3]:row[4]] = str(row["ID"]).replace('"', '')
                
        # 2. Map coordinates simultaneously against the active trees
        for var in current_vcf_batch:
            p1, p2 = var["START"], var["END"]
            if p1 >= p2:
                p2 = p1 + 1 # Validation fix for empty interval protections
                
            # Intersect with Genes
            gene_hits = gene_tree.overlap(p1, p2)
            for hit in gene_hits:
                var["gene"].append(hit.data)
                # Output tracking logic
                if p1 <= hit.begin and p2 >= hit.end:
                    print(f'GenePAV\t{chrom}:{p1}-{p2}\t{hit.data}', file=outeach)
                else:
                    print(f'gene\t{chrom}:{p1}-{p2}\t{hit.data}', file=outeach)
                    
            # Intersect with CDS regions
            cds_hits = cds_tree.overlap(p1, p2)
            for hit in cds_hits:
                var["CDS"].append(hit.data)
                if p1 <= hit.begin and p2 >= hit.end:
                    print(f'GenePAV_CDS\t{chrom}:{p1}-{p2}\t{hit.data}', file=outeach)
                else:
                    print(f'CDS\t{chrom}:{p1}-{p2}\t{hit.data}', file=outeach)
                    
            # Intersect with 3k Promoters
            pro_hits = pro3k_tree.overlap(p1, p2)
            for hit in pro_hits:
                var["3k_promoter"].append(hit.data)
                print(f'3k_promoter\t{chrom}:{p1}-{p2}\t{hit.data}', file=outeach)

    outeach.close()

    # 3. Format pooled outputs directly to flat records file
    print(f"Writing complete results structure to prefix: {args.output_file}_pos.txt")
    with open(f"{args.output_file}_pos.txt", "w") as fout:
        # Structure headers
        headers = ["#CHROM", "START", "END", "ID", "Gene_Overlap", "CDS_Overlap", "Promoter_Overlap"]
        fout.write("\t".join(headers) + "\n")
        
        for var in variants:
            # Drop duplicates and join multiple values cleanly into single comma strings
            g_str = ",".join(sorted(list(set(var["gene"])))) if var["gene"] else "."
            c_str = ",".join(sorted(list(set(var["CDS"])))) if var["CDS"] else "."
            p_str = ",".join(sorted(list(set(var["3k_promoter"])))) if var["3k_promoter"] else "."
            
            row = [
                str(var["CHROM"]), str(var["START"]), str(var["END"]), str(var["ID"]),
                g_str, c_str, p_str
            ]
            fout.write("\t".join(row) + "\n")

    print("Success. All operations complete.")

if __name__ == '__main__':
    main()
