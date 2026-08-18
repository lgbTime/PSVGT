#!/usr/bin/env python3
"""
Annotate VCF files directly with precise functional SV mutations and overlapping Gene IDs.
Generates group-level and population-level summary matrices and dual publication-ready plots.
Works with or without an explicit population group structure file.
"""

import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
import gzip
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# Adobe-editable Vector Font Configurations
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

def parse_args():
    parser = argparse.ArgumentParser(description='Population-level and group-level structural variant dual pie annotator with Gene ID integration.')
    parser.add_argument('-a', '--annotation', required=True, help='Path to GFF or GTF file')
    parser.add_argument('-v', '--vcf', required=True, help='Path to input VCF or VCF.GZ file')
    parser.add_argument('-g', '--groups', required=False, default=None, 
                        help='Optional tab-delimited group configuration file (Columns: SampleID\\tGroup)')
    parser.add_argument('-o', '--output_prefix', default='lettuce_sv_pop', help='Output base file prefix name')
    parser.add_argument('-f', '--min_freq', type=float, default=0.05, help='Frequency threshold (Default: 0.05)')
    parser.add_argument('-p', '--promoter_length', type=int, default=3000,
                        help='Upstream promoter region length in bp (Default: 3000)')
    return parser.parse_args()

def load_group_map(group_file):
    if not group_file:
        return None
    group_map = {}
    with open(group_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fields = line.rstrip().split("\t")
            if len(fields) < 2: continue
            group_map[fields[0]] = fields[1]
    return group_map

def load_gff(gff_file):
    print("Loading genomic annotations...")
    gff = pd.read_csv(gff_file, header=None, comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    return gff, sorted(list(set(gff[0].tolist())))

def extract_gene_cds(gff, gff_file):
    if "gff" in str(gff_file).lower():
        gene = gff[gff[2] == "gene"].copy()
        cds = gff[gff[2] == "CDS"].copy()
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'ID=([^;]+)')
        cds["ID"] = cds[8].str.extract(r'Parent=([^;]+)')
    elif "gtf" in str(gff_file).lower():
        gene = gff[gff[2] == "transcript"].copy()
        cds = gff[gff[2] == "exon"].copy()
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'gene_name "([^"]+)"')
        if gene["ID"].isna().all(): gene["ID"] = gene[8].str.extract(r'gene_name ([^;]+)')
        cds["ID"] = cds[8].str.extract(r'gene_name "([^"]+)"')
        if cds["ID"].isna().all(): cds["ID"] = cds[8].str.extract(r'gene_name ([^;]+)')
    return gene, cds

def get_promoter(gene, promoter_length=3000):
    pro = gene.copy()
    plus_strand = pro[6] == "+"
    minus_strand = pro[6] != "+"
    pro.loc[plus_strand, 4] = pro.loc[plus_strand, 3] - 1
    pro.loc[plus_strand, 3] = pro.loc[plus_strand, 3] - promoter_length
    pro.loc[minus_strand, 3] = pro.loc[minus_strand, 4] + 1
    pro.loc[minus_strand, 4] = pro.loc[minus_strand, 4] + promoter_length
    pro[2] = "promoter"
    pro.loc[pro[3] < 1, 3] = 1
    return pro

def compute_introns(cds_df):
    """Compute intron intervals as gaps between consecutive exons, grouped by gene ID."""
    intron_rows = []
    if cds_df.empty:
        return pd.DataFrame(columns=[0, 3, 4, "ID"])
    for (chrom, gid), group in cds_df.groupby([0, "ID"]):
        sorted_exons = group.sort_values(3)
        starts = sorted_exons[3].tolist()
        ends = sorted_exons[4].tolist()
        for i in range(len(starts) - 1):
            gap_start = int(ends[i]) + 1
            gap_end = int(starts[i + 1]) - 1
            if gap_start <= gap_end:
                intron_rows.append([chrom, gap_start, gap_end, gid])
    if intron_rows:
        return pd.DataFrame(intron_rows, columns=[0, 3, 4, "ID"])
    return pd.DataFrame(columns=[0, 3, 4, "ID"])

def process_vcf_lines_and_groups(vcf_file, group_map, min_freq):
    print("Parsing structural VCF streaming matrices...")
    variants = []
    sample_to_group = {}
    
    is_gzip = vcf_file.endswith(".gz")
    f = gzip.open(vcf_file, "rt") if is_gzip else open(vcf_file, "r")
    
    header_lines = []
    all_groups = []
    
    for line in f:
        if line.startswith("##"):
            header_lines.append(line.rstrip())
            continue
        if line.startswith("#CHROM"):
            header_lines.append('##INFO=<ID=ANNOTATION,Number=1,Type=String,Description="Functional SV Annotation Class determined via interval hierarchy cascades">')
            header_lines.append('##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Overlapping Feature Gene Identifier">')
            header_lines.append(line.rstrip())
            
            header = line.rstrip().split("\t")
            vcf_samples = header[9:]
            
            if group_map is None:
                sample_to_group = {idx: "All_Samples" for idx, s in enumerate(vcf_samples)}
                all_groups = ["All_Samples"]
            else:
                for idx, sample in enumerate(vcf_samples):
                    if sample in group_map:
                        sample_to_group[idx] = group_map[sample]
                all_groups = sorted(list(set(group_map.values())))
            continue
            
        fields = line.rstrip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref_allele = fields[3]
        alt = fields[4]
        if "TRA" in alt: continue
        
        info = fields[7]
        sv_type = "UNK"
        if "SVTYPE=" in info:
            sv_type = info.split("SVTYPE=")[1].split(";")[0]
        else:
            if "DEL" in alt: sv_type = "DEL"
            elif "DUP" in alt: sv_type = "DUP"
            elif "INV" in alt: sv_type = "INV"
            elif "INS" in alt: sv_type = "INS"
        
        end_pos = None
        if "END=" in info:
            try: end_pos = int(info.split("END=")[1].split(";")[0])
            except (IndexError, ValueError): pass
        if end_pos is None: end_pos = pos + len(ref_allele)

        svlen = None
        if "SVLEN=" in info:
            try: svlen = abs(int(info.split("SVLEN=")[1].split(";")[0]))
            except (IndexError, ValueError): pass
            
        group_counts = {g: [0, 0] for g in all_groups}
        genotypes = fields[9:]
        global_alt, global_tot = 0, 0
        
        for idx, gt_str in enumerate(genotypes):
            gt = gt_str.split(":")[0]
            if "./." in gt or "." in gt: continue
            alleles = gt.replace("|", "/").split("/")
            for a in alleles:
                global_tot += 1
                if a == "1": global_alt += 1
                
                if idx in sample_to_group:
                    group = sample_to_group[idx]
                    group_counts[group][1] += 1
                    if a == "1": group_counts[group][0] += 1
        
        active_groups = []
        for group, counts in group_counts.items():
            alt_c, tot_c = counts[0], counts[1]
            freq = alt_c / tot_c if tot_c > 0 else 0.0
            if freq >= min_freq: active_groups.append(group)
            
        global_freq = global_alt / global_tot if global_tot > 0 else 0.0
        if global_freq >= min_freq:
            active_groups.append("Global_Population")
        
        variants.append({
            "CHROM": chrom, "START": pos, "END": end_pos, "FIELDS": fields,
            "SV_TYPE": sv_type, "SVLEN": svlen, "ACTIVE_GROUPS": active_groups,
            "ANNOTATION": "Intergenic", "GENE_ID": "None"
        })
        
    f.close()
    return header_lines, variants, all_groups

def generate_plots(detailed_df, broad_df, plot_targets, prefix):
    print("Plotting Category 1: High-Resolution Detailed Donut Charts (Excluding Intergenic)...")
    detailed_plot_data = detailed_df.drop(index="Intergenic", errors="ignore")
    detailed_plot_data = detailed_plot_data[detailed_plot_data.sum(axis=1) > 0]
    
    color_palette_1 = sns.color_palette("Set3", n_colors=len(detailed_plot_data))
    
    for target in plot_targets:
        if target not in detailed_plot_data.columns: continue
        counts = detailed_plot_data[target]
        if counts.sum() == 0: continue
            
        fig, ax = plt.subplots(figsize=(10.5, 7.5))
        wedges, texts, autotexts = ax.pie(
            counts, labels=None, autopct='%1.1f%%', startangle=140,
            colors=color_palette_1, pctdistance=0.85,
            textprops=dict(color="black", weight="bold", size=9)
        )
        
        centre_circle = plt.Circle((0,0), 0.70, fc='white')
        ax.add_artist(centre_circle)
        
        legend_labels = [f"{label} ({int(counts[label]):,})" for label in counts.index]
        ax.legend(wedges, legend_labels, title="Detailed Structural Impact", loc="center left", bbox_to_anchor=(1, 0.5), frameon=True)
        
        title_label = "Global Population Pool" if target == "Global_Population" else f"Cohort Variety Group: {target}"
        ax.set_title(f"Detailed Structural Mutation Profile\\n{title_label}\\n(Excluding Intergenic Deserts)", fontsize=12, weight="bold", pad=15)
        plt.tight_layout()
        
        plt.savefig(f"{prefix}_detailed_donut_{target}.pdf", dpi=300, bbox_inches='tight')
        plt.close()

    print("Plotting Category 2: Traditional Landscape Pie Charts (Including Intergenic)...")
    color_palette_2 = ["#729ece", "#ff9e4a", "#67bf5c", "#ed665d"]
    
    for target in plot_targets:
        if target not in broad_df.columns: continue
        counts = broad_df[target]
        if counts.sum() == 0: continue
            
        fig, ax = plt.subplots(figsize=(8.5, 7))
        wedges, texts, autotexts = ax.pie(
            counts, labels=None, autopct='%1.1f%%', startangle=90,
            colors=color_palette_2, pctdistance=0.75,
            wedgeprops={"edgecolor": "black", "linewidth": 1.2, "antialiased": True},
            textprops=dict(color="black", weight="bold", size=10)
        )
        
        legend_labels = [f"{label} ({int(counts[label]):,})" for label in counts.index]
        ax.legend(wedges, legend_labels, title="Genomic Region", loc="center left", bbox_to_anchor=(1, 0.5), frameon=True)
        
        title_label = "Global Population Pool" if target == "Global_Population" else f"Cohort Variety Group: {target}"
        ax.set_title(f"Genomic Landscape of SV Mutations\\n{title_label}\\n(Complete Sequence Proportions)", fontsize=13, weight="bold", pad=15)
        plt.tight_layout()
        
        plt.savefig(f"{prefix}_traditional_landscape_pie_{target}.pdf", dpi=300, bbox_inches='tight')
        plt.close()

def print_glossary(promoter_length):
    """Print a glossary of annotation terms."""
    plen = promoter_length
    plabel = f"{plen // 1000}kb" if plen >= 1000 and plen % 1000 == 0 else f"{plen}bp"
    glossary = [
        ("GenePAV",              "DEL covers >=60% of a gene body (transcript)"),
        ("GenePAV_CDS",          "DEL fully encompasses a CDS/exon, but <60% of the gene body"),
        ("Gene_Copy",            "DUP covers >=80% of a gene body"),
        ("Gene_Copy_CDS",        "DUP fully encompasses a CDS/exon, but <80% of the gene body"),
        ("Exonic_Frameshift",    "SV in CDS/exon where the affected CDS bases are not a multiple of 3"),
        ("Exonic_Disruption_Partial_DEL", "DEL partially overlaps CDS/exon, in-frame (CDS overlap % 3 == 0)"),
        ("Exonic_Duplication_Partial_DUP","DUP partially overlaps CDS/exon, in-frame"),
        ("Exonic_Insertion_Disruption",   "INS in CDS/exon, in-frame (SVLEN % 3 == 0)"),
        ("Exonic_Inversion_Breakpoint",   "INV partially overlaps CDS/exon (no whole-gene cover)"),
        ("Whole_Exon_Inversion",          "INV fully encompasses a CDS/exon, but not the entire gene body"),
        ("CDS_Overlap",                   "Other SV type overlapping CDS/exon"),
        ("Intron_Lost",           "DEL fully encompasses an intron (both intron boundaries within DEL)"),
        ("Intronic_Disruption_Partial_DEL","DEL overlaps intron but does not fully encompass it"),
        ("Intronic_Duplication_Partial_DUP","DUP overlaps intron"),
        ("Intronic_Inversion_Breakpoint",  "INV partially overlaps gene body"),
        ("Whole_Gene_Inversion",           "INV fully encompasses a gene body (takes priority over exon-level)"),
        ("Intronic_Insertion",             "INS in intron"),
        ("Gene_Body_Overlap",              "Other SV type overlapping gene body"),
        ("Promoter_Deletion",     f"DEL in {plabel} upstream promoter region"),
        ("Promoter_Duplication",  f"DUP in {plabel} upstream promoter region"),
        ("Promoter_Insertion",    f"INS in {plabel} upstream promoter region"),
        ("Promoter_Inversion",    f"INV in {plabel} upstream promoter region"),
        ("Promoter_Region_3kb_Upstream", f"Other SV type in {plabel} promoter region"),
        ("Intergenic",           "SV does not overlap any gene, CDS, or promoter"),
    ]
    print()
    print("=" * 72)
    print("ANNOTATION GLOSSARY")
    print("=" * 72)
    for term, definition in glossary:
        print(f"  {term:<36} {definition}")
    print("=" * 72)

def main():
    args = parse_args()
    group_map = load_group_map(args.groups)
    gff, chroms = load_gff(args.annotation)
    gene, cds = extract_gene_cds(gff, args.annotation)
    pro3k = get_promoter(gene, args.promoter_length)
    
    gene = gene[[0, 3, 4, "ID"]].dropna()
    cds = cds[[0, 3, 4, "ID"]].dropna()
    pro3k = pro3k[[0, 3, 4, "ID"]].dropna()
    introns = compute_introns(cds)
    
    headers, variants, basic_groups = process_vcf_lines_and_groups(args.vcf, group_map, args.min_freq)
    all_tracking_targets = basic_groups + ["Global_Population"]
    
    vcf_pools = {}
    for var in variants:
        vcf_pools.setdefault(var["CHROM"], []).append(var)
        
    print("Running interval tree overlaps and checking PAV/Copy rules...")
    for chrom in tqdm(chroms):
        if chrom not in vcf_pools: continue
        current_vcf_batch = vcf_pools[chrom]
        
        # Store metadata inside structural tree objects: tuple format (start, end, data_payload)
        gene_intervals = [(int(r[3]), int(r[4]), str(r["ID"])) for _, r in gene[gene[0] == chrom].iterrows() if r[3] < r[4]]
        cds_intervals = [(int(r[3]), int(r[4]), str(r["ID"])) for _, r in cds[cds[0] == chrom].iterrows() if r[3] < r[4]]
        pro_intervals = [(int(r[3]), int(r[4]), str(r["ID"])) for _, r in pro3k[pro3k[0] == chrom].iterrows() if r[3] < r[4]]
        
        intron_intervals = [(int(r[3]), int(r[4]), str(r["ID"])) for _, r in introns[introns[0] == chrom].iterrows() if r[3] < r[4]]

        gene_tree = IntervalTree.from_tuples(gene_intervals) if gene_intervals else IntervalTree()
        cds_tree = IntervalTree.from_tuples(cds_intervals) if cds_intervals else IntervalTree()
        pro3k_tree = IntervalTree.from_tuples(pro_intervals) if pro_intervals else IntervalTree()
        intron_tree = IntervalTree.from_tuples(intron_intervals) if intron_intervals else IntervalTree()
        
        for var in current_vcf_batch:
            p1, p2 = var["START"], var["END"]
            if p1 >= p2: p2 = p1 + 1
            sv_type = var["SV_TYPE"]
            
            cds_hits = cds_tree.overlap(p1, p2)
            gene_hits = gene_tree.overlap(p1, p2)
            pro_hits = pro3k_tree.overlap(p1, p2)

            # ---------------------------------------------------------
            # Priority 1: GenePAV (DEL >= 60% of gene body) / Gene_Copy (DUP >= 80% of gene body)
            # ---------------------------------------------------------
            gene_pav = False
            gene_pav_ids = []
            gene_copy = False
            gene_copy_ids = []
            for hit in gene_hits:
                overlap_len = min(p2, hit.end) - max(p1, hit.begin) + 1
                gene_len = hit.end - hit.begin + 1
                if gene_len > 0:
                    if sv_type == "DEL" and overlap_len / gene_len >= 0.6:
                        gene_pav = True
                        gene_pav_ids.append(hit.data)
                    elif sv_type == "DUP" and overlap_len / gene_len >= 0.8:
                        gene_copy = True
                        gene_copy_ids.append(hit.data)

            if gene_pav and sv_type == "DEL":
                var["GENE_ID"] = ",".join(list(set(gene_pav_ids)))
                var["ANNOTATION"] = "GenePAV"
            elif gene_copy and sv_type == "DUP":
                var["GENE_ID"] = ",".join(list(set(gene_copy_ids)))
                var["ANNOTATION"] = "Gene_Copy"
            # ---------------------------------------------------------
            # Priority 2: CDS/exon-level annotations
            # ---------------------------------------------------------
            elif cds_hits:
                var["GENE_ID"] = ",".join(list(set([hit.data for hit in cds_hits])))
                svlen = var.get("SVLEN")
                if sv_type == "DEL":
                    is_whole_exon = any(p1 <= hit.begin and p2 >= hit.end for hit in cds_hits)
                    if is_whole_exon:
                        # Frameshift if any fully-encompassed exon length is not a multiple of 3
                        whole_exon_fs = any((hit.end - hit.begin + 1) % 3 != 0 for hit in cds_hits if p1 <= hit.begin and p2 >= hit.end)
                        var["ANNOTATION"] = "Exonic_Frameshift" if whole_exon_fs else "GenePAV_CDS"
                    else:
                        # Total CDS overlap across all affected exons determines frameshift
                        total_ov = sum(min(p2, hit.end) - max(p1, hit.begin) + 1 for hit in cds_hits)
                        var["ANNOTATION"] = "Exonic_Frameshift" if total_ov % 3 != 0 else "Exonic_Disruption_Partial_DEL"
                elif sv_type == "DUP":
                    is_whole_exon = any(p1 <= hit.begin and p2 >= hit.end for hit in cds_hits)
                    if is_whole_exon:
                        # Frameshift if any fully-encompassed exon length is not a multiple of 3
                        whole_exon_fs = any((hit.end - hit.begin + 1) % 3 != 0 for hit in cds_hits if p1 <= hit.begin and p2 >= hit.end)
                        var["ANNOTATION"] = "Exonic_Frameshift" if whole_exon_fs else "Gene_Copy_CDS"
                    else:
                        # Total CDS overlap across all affected exons determines frameshift
                        total_ov = sum(min(p2, hit.end) - max(p1, hit.begin) + 1 for hit in cds_hits)
                        var["ANNOTATION"] = "Exonic_Frameshift" if total_ov % 3 != 0 else "Exonic_Duplication_Partial_DUP"
                elif sv_type == "INS":
                    if svlen is not None and svlen % 3 != 0:
                        var["ANNOTATION"] = "Exonic_Frameshift"
                    else:
                        var["ANNOTATION"] = "Exonic_Insertion_Disruption"
                elif sv_type == "INV":
                    # Priority: whole-gene inversion > whole-exon inversion
                    # Check if INV fully covers any complete gene body
                    is_whole_gene = any(p1 <= hit.begin and p2 >= hit.end for hit in gene_hits)
                    if is_whole_gene:
                        var["GENE_ID"] = ",".join(list(set([hit.data for hit in gene_hits if p1 <= hit.begin and p2 >= hit.end])))
                        var["ANNOTATION"] = "Whole_Gene_Inversion"
                    else:
                        is_whole_exon = any(p1 <= hit.begin and p2 >= hit.end for hit in cds_hits)
                        var["ANNOTATION"] = "Whole_Exon_Inversion" if is_whole_exon else "Exonic_Inversion_Breakpoint"
                else: var["ANNOTATION"] = "CDS_Overlap"
            elif gene_hits:
                var["GENE_ID"] = ",".join(list(set([hit.data for hit in gene_hits])))
                if sv_type == "DEL":
                    intron_hits = intron_tree.overlap(p1, p2)
                    is_intron_loss = any(p1 <= hit.begin and p2 >= hit.end for hit in intron_hits)
                    var["ANNOTATION"] = "Intron_Lost" if is_intron_loss else "Intronic_Disruption_Partial_DEL"
                elif sv_type == "DUP":
                    var["ANNOTATION"] = "Intronic_Duplication_Partial_DUP"
                elif sv_type == "INV":
                    is_whole = any(p1 <= hit.begin and p2 >= hit.end for hit in gene_hits)
                    var["ANNOTATION"] = "Whole_Gene_Inversion" if is_whole else "Intronic_Inversion_Breakpoint"
                elif sv_type == "INS":
                    var["ANNOTATION"] = "Intronic_Insertion"
                else: var["ANNOTATION"] = "Gene_Body_Overlap"
            elif pro_hits:
                var["GENE_ID"] = ",".join(list(set([hit.data for hit in pro_hits])))
                if sv_type == "INV": var["ANNOTATION"] = "Promoter_Inversion"
                elif sv_type == "INS": var["ANNOTATION"] = "Promoter_Insertion"
                elif sv_type == "DEL": var["ANNOTATION"] = "Promoter_Deletion"
                elif sv_type == "DUP": var["ANNOTATION"] = "Promoter_Duplication"
                else: var["ANNOTATION"] = "Promoter_Region_3kb_Upstream"
            else:
                var["ANNOTATION"] = "Intergenic"
                var["GENE_ID"] = "None"

    out_vcf_name = f"{args.output_prefix}_annotated.vcf"
    print(f"Writing updated variants into customized output VCF: {out_vcf_name}")
    with open(out_vcf_name, 'w') as out_vcf:
        for line in headers: out_vcf.write(line + "\n")
        for var in variants:
            fields = var["FIELDS"]
            # Prepend both keys sequentially into info fields column
            fields[7] = f"ANNOTATION={var['ANNOTATION']};GENE_ID={var['GENE_ID']};{fields[7]}"
            out_vcf.write("\t".join(fields) + "\n")

    # -------------------------------------------------------------
    # GENERATE MATRIX 1: FULL DETAIL COUNTS
    # -------------------------------------------------------------
    categories = [
        "GenePAV_CDS", "Exonic_Disruption_Partial_DEL", "Exonic_Frameshift",
        "Gene_Copy_CDS", "Exonic_Duplication_Partial_DUP",
        "Exonic_Inversion_Breakpoint", "Whole_Exon_Inversion", "Exonic_Insertion_Disruption", "CDS_Overlap",
        "GenePAV", "Intron_Lost", "Intronic_Disruption_Partial_DEL",
        "Gene_Copy", "Intronic_Duplication_Partial_DUP",
        "Intronic_Inversion_Breakpoint", "Whole_Gene_Inversion", "Intronic_Insertion", "Gene_Body_Overlap",
        "Promoter_Inversion", "Promoter_Insertion", "Promoter_Deletion", "Promoter_Duplication",
        "Promoter_Region_3kb_Upstream", "Intergenic"
    ]
    
    stats_dict = {g: {cat: 0 for cat in categories} for g in all_tracking_targets}
    for var in variants:
        cat = var["ANNOTATION"]
        if cat in stats_dict["Global_Population"]:
            for group in var["ACTIVE_GROUPS"]:
                stats_dict[group][cat] += 1
                
    detailed_count_df = pd.DataFrame(stats_dict)
    detailed_counts_file = f"{args.output_prefix}_detailed_counts_matrix.txt"
    detailed_count_df.to_csv(detailed_counts_file, sep="\t")
    print(f"Detailed counts matrix generated at: {detailed_counts_file}")
    
    # -------------------------------------------------------------
    # GENERATE MATRIX 2: COLLAPSED GENOMIC REGIONS
    # -------------------------------------------------------------
    category_mapping = {
        "GenePAV_CDS": "Exon", "Exonic_Disruption_Partial_DEL": "Exon",
        "Exonic_Frameshift": "Exon",
        "Gene_Copy_CDS": "Exon", "Exonic_Duplication_Partial_DUP": "Exon",
        "Exonic_Inversion_Breakpoint": "Exon", "Whole_Exon_Inversion": "Exon",
        "Exonic_Insertion_Disruption": "Exon", "CDS_Overlap": "Exon",

        "GenePAV": "Intron", "Intron_Lost": "Intron",
        "Intronic_Disruption_Partial_DEL": "Intron",
        "Gene_Copy": "Intron", "Intronic_Duplication_Partial_DUP": "Intron",
        "Intronic_Inversion_Breakpoint": "Intron", "Whole_Gene_Inversion": "Intron",
        "Intronic_Insertion": "Intron", "Gene_Body_Overlap": "Intron",

        "Promoter_Inversion": "Promoter", "Promoter_Insertion": "Promoter",
        "Promoter_Deletion": "Promoter", "Promoter_Duplication": "Promoter",
        "Promoter_Region_3kb_Upstream": "Promoter",

        "Intergenic": "Intergenic"
    }
    
    collapsed_df = detailed_count_df.copy()
    collapsed_df['Major_Category'] = collapsed_df.index.map(category_mapping)
    broad_count_df = collapsed_df.groupby('Major_Category').sum()
    broad_count_df = broad_count_df.reindex(["Intergenic", "Exon", "Intron", "Promoter"]).fillna(0).astype(int)
    
    broad_counts_file = f"{args.output_prefix}_major_landscape_matrix.txt"
    broad_count_df.to_csv(broad_counts_file, sep="\t")
    print(f"Broad landscape matrix generated at: {broad_counts_file}")
    
    generate_plots(detailed_count_df, broad_count_df, all_tracking_targets, args.output_prefix)
    print("Execution complete. Output figures saved successfully.")
    print_glossary(args.promoter_length)

if __name__ == '__main__':
    main()

