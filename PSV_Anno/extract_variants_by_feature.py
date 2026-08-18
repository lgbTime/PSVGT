#!/usr/bin/env python3
"""
extract_variants_by_feature.py  (v2)
====================================
Extract/annotate variants from a PAV matrix or VCF file based on their genomic
overlap with annotation features from a GTF file.

Features supported:
  exon      - overlapping any annotated exon
  cds       - overlapping any annotated CDS
  intron    - overlapping computed intronic gaps (between exons of same gene)
  gene_body - overlapping any annotated transcript region
  promoter  - overlapping upstream promoter region (± --promoter-length bp)

Input auto-detection
--------------------
  *.vcf / *.vcf.gz   → treated as VCF  (produces filtered + annotated VCFs)
  *.txt / *.tsv       → treated as PAV matrix (produces per-feature tab files)
  Use --input-format vcf|pav to override.

VCF output (for VCF input)
--------------------------
  1. {prefix}.{feature}.filtered.vcf
     Only variants that overlap the requested feature(s).
     One filtered VCF per feature.

  2. {prefix}.annotated.vcf
     ALL variants from the input VCF, with ANNOTATION=<feature> and
     GENE_ID=<ids> prepended to the INFO column.  Variants matching
     multiple requested features have a comma-separated ANNOTATION.
     Intergenic variants get ANNOTATION=Intergenic;GENE_ID=None.

PAV output (for PAV input)
--------------------------
  1. {prefix}.{feature}.txt
     PAV rows whose interval overlaps the requested feature.
     Two extra columns appended: FeatureType, OverlappingGeneIDs.

Usage
-----
  # VCF input — produces filtered + annotated VCFs
  python extract_variants_by_feature.py \\
      -g annotation.gtf -i variants.vcf -f exon intron cds -o out

  # PAV input — produces per-feature tab files
  python extract_variants_by_feature.py \\
      -g annotation.gtf -i 5.PSVIndel.PAV.txt -f exon intron gene_body -o out

  # Only the annotated VCF (skip per-feature filtered VCFs)
  python extract_variants_by_feature.py \\
      -g annotation.gtf -i variants.vcf -f exon intron \\
      --no-filtered -o out

Author : auto-generated
Date   : 2026-06-06
"""

import argparse
import sys
import gzip
from collections import defaultdict
from intervaltree import IntervalTree

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

VALID_FEATURES = {"exon", "cds", "intron", "gene_body", "promoter"}

# Mapping from requested feature to the GTF / computed source
FEATURE_SOURCES = {
    "gene_body": "transcript",
    "exon":      "exon",
    "cds":       "CDS",
}


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract/annotate variants overlapping genomic features "
                    "from a PAV matrix or VCF file."
    )
    p.add_argument("-g", "--gtf", required=True,
                   help="GTF annotation file.")
    p.add_argument("-i", "--input", required=True,
                   help="Input variant file (PAV tab-delimited or VCF).")
    p.add_argument("-f", "--features", nargs="+", required=True,
                   choices=sorted(VALID_FEATURES),
                   help="Feature type(s) to extract/annotate.")
    p.add_argument("-o", "--output-prefix", default="extracted",
                   help="Output file prefix (default: extracted).")
    p.add_argument("-p", "--promoter-length", type=int, default=3000,
                   help="Upstream promoter region length in bp (default: 3000).")
    p.add_argument("--input-format", choices=("vcf", "pav"), default=None,
                   help="Force input format; by default auto-detected from "
                        "file extension.")
    p.add_argument("--no-filtered", action="store_true",
                   help="For VCF input: skip per-feature filtered VCF files; "
                        "only write the full annotated VCF.")
    p.add_argument("--no-annotated", action="store_true",
                   help="For VCF input: skip the full annotated VCF; "
                        "only write per-feature filtered VCFs.")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="Suppress progress messages.")
    return p.parse_args()


# ---------------------------------------------------------------------------
#  GTF parsing
# ---------------------------------------------------------------------------

def _extract_gene_id(attr_string: str) -> str:
    """Extract gene_name from GTF attributes; fall back to gene_id."""
    for field in attr_string.rstrip(";").split("; "):
        if field.startswith('gene_name "'):
            return field.split('gene_name "')[1].rstrip('"')
    for field in attr_string.rstrip(";").split("; "):
        if field.startswith("gene_name "):
            return field.split("gene_name ")[1]
    for field in attr_string.rstrip(";").split("; "):
        if field.startswith('gene_id "'):
            return field.split('gene_id "')[1].rstrip('"')
    for field in attr_string.rstrip(";").split("; "):
        if field.startswith('transcript_id "'):
            return field.split('transcript_id "')[1].rstrip('"')
    return "unknown"


def load_gtf(gtf_path: str) -> dict:
    """
    Parse a GTF file.

    Returns a dict keyed by chromosome.  Each value is another dict:
        {"transcript": [(start, end, gene_id, strand), ...],
         "exon":       [...],
         "CDS":        [...]}
    Special key '__chroms' holds the sorted list of all chromosome names.
    """
    log(f"[1/5] Parsing GTF: {gtf_path}")
    raw = defaultdict(list)
    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attr = parts
            if ftype not in ("transcript", "exon", "CDS"):
                continue
            gid = _extract_gene_id(attr)
            raw[chrom].append((ftype, int(start), int(end), gid, strand))

    chrom_data: dict = {}
    for chrom, entries in raw.items():
        chrom_data[chrom] = {"transcript": [], "exon": [], "CDS": []}
        for ftype, start, end, gid, strand in entries:
            chrom_data[chrom][ftype].append((start, end, gid, strand))

    chroms = sorted(chrom_data.keys())
    chrom_data["__chroms"] = chroms

    n_t = sum(len(v["transcript"]) for k, v in chrom_data.items() if k != "__chroms")
    n_e = sum(len(v["exon"])       for k, v in chrom_data.items() if k != "__chroms")
    n_c = sum(len(v["CDS"])        for k, v in chrom_data.items() if k != "__chroms")
    log(f"       → {n_t:,} transcripts, {n_e:,} exons, {n_c:,} CDS "
        f"across {len(chroms)} chromosomes")
    return chrom_data


# ---------------------------------------------------------------------------
#  Derived features — introns & promoters
# ---------------------------------------------------------------------------

def compute_introns(chrom_data: dict) -> dict:
    """
    Compute intron intervals as gaps between consecutive exons of the same
    gene (gene_id).  Returns {chrom: [(start, end, gene_id), ...]}
    """
    log("[2/5] Computing intron intervals …")
    introns: dict = defaultdict(list)
    for chrom, features in chrom_data.items():
        if chrom == "__chroms":
            continue
        exons_by_gene = defaultdict(list)
        for s, e, gid, _ in features.get("exon", []):
            exons_by_gene[gid].append((s, e))
        for gid, exons in exons_by_gene.items():
            exons.sort(key=lambda x: x[0])
            for i in range(len(exons) - 1):
                gap_start = exons[i][1] + 1
                gap_end   = exons[i + 1][0] - 1
                if gap_start <= gap_end:
                    introns[chrom].append((gap_start, gap_end, gid))
    total = sum(len(v) for v in introns.values())
    log(f"       → {total:,} intron intervals")
    return dict(introns)


def compute_promoters(chrom_data: dict, promoter_len: int = 3000) -> dict:
    """
    Compute upstream promoter regions.
    + strand: [start - promoter_len, start - 1]
    - strand: [end + 1, end + promoter_len]
    Returns {chrom: [(start, end, gene_id), ...]}
    """
    log(f"[2/5] Computing promoter regions (±{promoter_len} bp) …")
    promoters: dict = defaultdict(list)
    for chrom, features in chrom_data.items():
        if chrom == "__chroms":
            continue
        for s, e, gid, strand in features.get("transcript", []):
            if strand == "+":
                p_s, p_e = max(1, s - promoter_len), s - 1
            else:
                p_s, p_e = e + 1, e + promoter_len
            if p_s < p_e:
                promoters[chrom].append((p_s, p_e, gid))
    total = sum(len(v) for v in promoters.values())
    log(f"       → {total:,} promoter intervals")
    return dict(promoters)


# ---------------------------------------------------------------------------
#  Interval trees
# ---------------------------------------------------------------------------

def build_trees(chrom_data: dict, intron_data: dict, promoter_data: dict,
                requested_features: set) -> dict:
    """
    Build an IntervalTree per chromosome per requested feature.
    Returns {chrom: { 'exon': IntervalTree, 'cds': IntervalTree, … }}
    """
    log("[3/5] Building interval trees …")
    trees: dict = {}
    for chrom in chrom_data.get("__chroms", []):
        trees[chrom] = {}

        # GTF-based features
        for feat_name in ("gene_body", "exon", "cds"):
            if feat_name not in requested_features:
                continue
            gtf_field = FEATURE_SOURCES[feat_name]
            ivs = [(s, e, gid)
                   for s, e, gid, _ in chrom_data.get(chrom, {}).get(gtf_field, [])
                   if s < e]
            trees[chrom][feat_name] = (IntervalTree.from_tuples(ivs) if ivs
                                       else IntervalTree())

        # Computed features
        for feat_name, source in (("intron", intron_data),
                                  ("promoter", promoter_data)):
            if feat_name not in requested_features:
                continue
            ivs = [(s, e, gid) for s, e, gid in source.get(chrom, []) if s < e]
            trees[chrom][feat_name] = (IntervalTree.from_tuples(ivs) if ivs
                                       else IntervalTree())
    return trees


# ---------------------------------------------------------------------------
#  VCF helpers
# ---------------------------------------------------------------------------

def parse_sv_end(pos: int, ref: str, info: str, alt: str) -> int:
    """Extract the end coordinate from a VCF record."""
    # Try END= in INFO first
    for token in info.split(";"):
        if token.startswith("END="):
            try:
                return int(token.split("=")[1])
            except (ValueError, IndexError):
                pass
    # Fallback: use REF length
    return pos + len(ref)


def parse_sv_type(alt: str, info: str) -> str:
    """Extract SVTYPE from INFO or infer from ALT."""
    for token in info.split(";"):
        if token.startswith("SVTYPE="):
            return token.split("=")[1]
    if "<DEL>" in alt or "DEL" in alt:
        return "DEL"
    if "<DUP>" in alt or "DUP" in alt:
        return "DUP"
    if "<INS>" in alt or "INS" in alt:
        return "INS"
    if "<INV>" in alt or "INV" in alt:
        return "INV"
    return "UNK"


# ---------------------------------------------------------------------------
#  Core annotation logic
# ---------------------------------------------------------------------------

def annotate_variant(chrom: str, start: int, end: int,
                     trees: dict, requested_features: list) -> dict:
    """
    Query all requested feature trees for a single variant.
    Returns:
        {
            "features": ["exon", "intron"],        # matched features
            "gene_ids": {"exon": ["GENE1"], ...},  # per-feature gene IDs
            "all_gene_ids": ["GENE1", "GENE2"],    # deduplicated union
        }
    """
    if chrom not in trees:
        return {"features": [], "gene_ids": {}, "all_gene_ids": []}

    matched = []
    gene_ids_by_feat = {}
    all_gids = set()

    if start == end:
        end = start + 1
    if start > end:
        start, end = end, start

    for feat in requested_features:
        tree = trees[chrom].get(feat)
        if tree is None:
            continue
        hits = tree.overlap(start, end)
        if hits:
            matched.append(feat)
            gids = sorted(set(h.data for h in hits))
            gene_ids_by_feat[feat] = gids
            all_gids.update(gids)

    return {
        "features": matched,
        "gene_ids": gene_ids_by_feat,
        "all_gene_ids": sorted(all_gids),
    }


# ---------------------------------------------------------------------------
#  VCF processing  (streaming)
# ---------------------------------------------------------------------------

def process_vcf(vcf_path: str, trees: dict, requested_features: list,
                output_prefix: str, no_filtered: bool, no_annotated: bool,
                chroms_with_data: set):
    """
    Stream through a VCF file.  For each variant:
      - Query overlap against trees
      - Write to per-feature filtered VCFs  (unless --no-filtered)
      - Write to annotated VCF              (unless --no-annotated)
    """
    feature_list = sorted(requested_features)
    is_gzip = vcf_path.endswith(".gz")

    # --- Open output handles ---
    filtered_fhs: dict = {}
    filtered_paths: dict = {}
    annotated_fh = None
    annotated_path = f"{output_prefix}.annotated.vcf"

    if not no_filtered:
        for feat in feature_list:
            fpath = f"{output_prefix}.{feat}.filtered.vcf"
            filtered_paths[feat] = fpath
            filtered_fhs[feat] = open(fpath, "w")

    if not no_annotated:
        annotated_fh = open(annotated_path, "w")

    # --- Stream ---
    stats = {
        "total": 0, "skipped_chrom": 0, "intergenic": 0,
        "written_filtered": {f: 0 for f in feature_list},
    }

    header_lines: list = []
    samples_header: str = ""
    in_header = True
    n_sample_cols = 0

    log(f"[4/5] Streaming VCF: {vcf_path}")

    opener = gzip.open if is_gzip else open
    fh_in = opener(vcf_path, "rt") if is_gzip else open(vcf_path, "r")

    try:
        for line in fh_in:
            line = line.rstrip("\n")

            # ---- Collect headers ----
            if in_header:
                if line.startswith("##"):
                    header_lines.append(line)
                    continue
                if line.startswith("#CHROM"):
                    # Inject INFO header lines
                    header_lines.append(
                        '##INFO=<ID=ANNOTATION,Number=1,Type=String,'
                        'Description="Functional SV Annotation Class">'
                    )
                    header_lines.append(
                        '##INFO=<ID=GENE_ID,Number=1,Type=String,'
                        'Description="Overlapping Feature Gene Identifier(s)">'
                    )
                    header_lines.append(line)
                    samples_header = line
                    n_sample_cols = len(line.split("\t")) - 9
                    # Write headers to all output files
                    if not no_filtered:
                        for fh in filtered_fhs.values():
                            for hl in header_lines:
                                fh.write(hl + "\n")
                    if not no_annotated and annotated_fh:
                        for hl in header_lines:
                            annotated_fh.write(hl + "\n")
                    in_header = False
                    continue

            # ---- Data line ----
            stats["total"] += 1
            if stats["total"] % 5_000_000 == 0:
                parts = ", ".join(f"{f}={stats['written_filtered'][f]:,}"
                                  for f in feature_list)
                log(f"       … scanned {stats['total']:>12,}  [{parts}]  "
                    f"intergenic={stats['intergenic']:,}  "
                    f"skip_chrom={stats['skipped_chrom']:,}")

            fields = line.split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except (ValueError, IndexError):
                continue
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            # Skip translocations
            sv_type = parse_sv_type(alt, info)
            if sv_type == "TRA":
                continue

            end = parse_sv_end(pos, ref, info, alt)

            # Skip chromosomes not in GTF
            if chrom not in chroms_with_data:
                stats["skipped_chrom"] += 1
                ann_str = "ANNOTATION=Intergenic;GENE_ID=None"
                fields[7] = f"{ann_str};{info}"
                if not no_annotated and annotated_fh:
                    annotated_fh.write("\t".join(fields) + "\n")
                continue

            # Annotate
            ann = annotate_variant(chrom, pos, end, trees, feature_list)

            if not ann["features"]:
                # Intergenic variant — not written to any filtered VCF
                stats["intergenic"] += 1
                ann_str = "ANNOTATION=Intergenic;GENE_ID=None"
                fields[7] = f"{ann_str};{info}"
                if not no_annotated and annotated_fh:
                    annotated_fh.write("\t".join(fields) + "\n")
                continue

            # ---- Matched at least one feature ----
            matched_feats = ann["features"]
            all_gids = ",".join(ann["all_gene_ids"]) if ann["all_gene_ids"] else "None"

            # Per-feature filtered VCFs — write only to the matched feature files
            if not no_filtered:
                for feat in matched_feats:
                    gids = ",".join(ann["gene_ids"].get(feat, ["None"]))
                    ann_info = f"ANNOTATION={feat};GENE_ID={gids}"
                    fields[7] = f"{ann_info};{info}"
                    filtered_fhs[feat].write("\t".join(fields) + "\n")
                    stats["written_filtered"][feat] += 1

            # Full annotated VCF (all feature matches combined)
            ann_feat = ",".join(matched_feats)
            ann_info = f"ANNOTATION={ann_feat};GENE_ID={all_gids}"
            fields[7] = f"{ann_info};{info}"
            if not no_annotated and annotated_fh:
                annotated_fh.write("\t".join(fields) + "\n")

    finally:
        fh_in.close()

    # Close output handles
    if not no_filtered:
        for fh in filtered_fhs.values():
            fh.close()
    if not no_annotated and annotated_fh:
        annotated_fh.close()

    return stats, filtered_paths, annotated_path


# ---------------------------------------------------------------------------
#  PAV processing  (streaming — original behaviour)
# ---------------------------------------------------------------------------

def process_pav(pav_path: str, trees: dict, requested_features: list,
                output_prefix: str, chroms_with_data: set) -> dict:
    """Stream through a PAV matrix and write per-feature tab files."""
    feature_list = sorted(requested_features)
    out_fhs: dict = {}
    out_paths: dict = {}

    for feat in feature_list:
        opath = f"{output_prefix}.{feat}.txt"
        out_paths[feat] = opath
        out_fhs[feat] = open(opath, "w")

    stats = {"total": 0, "skipped_chrom": 0,
             "written": {f: 0 for f in feature_list}}
    header_written = False

    log(f"[4/5] Streaming PAV: {pav_path}")

    with open(pav_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if not header_written:
                    header_written = True
                    for feat in feature_list:
                        out_fhs[feat].write(
                            f"{line}\tFeatureType\tOverlappingGeneIDs\n"
                        )
                continue

            stats["total"] += 1
            if stats["total"] % 5_000_000 == 0:
                parts = ", ".join(
                    f"{f}={stats['written'][f]:,}" for f in feature_list
                )
                log(f"       … scanned {stats['total']:>12,}  [{parts}]  "
                    f"skip_chrom={stats['skipped_chrom']:,}")

            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end   = int(fields[2])
            except (ValueError, IndexError):
                continue

            if chrom not in chroms_with_data:
                stats["skipped_chrom"] += 1
                continue

            for feat in feature_list:
                ann = annotate_variant(chrom, start, end, trees, [feat])
                if ann["features"]:
                    gids = ",".join(ann["all_gene_ids"])
                    out_fhs[feat].write(f"{line}\t{feat}\t{gids}\n")
                    stats["written"][feat] += 1

    for fh in out_fhs.values():
        fh.close()
    return stats, out_paths


# ---------------------------------------------------------------------------
#  Utilities
# ---------------------------------------------------------------------------

_quiet = False


def log(msg: str) -> None:
    if not _quiet:
        print(msg, file=sys.stderr)


def report_line(msg: str) -> None:
    print(f"  {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main() -> None:
    global _quiet
    args = parse_args()
    _quiet = args.quiet
    requested_features = set(args.features)

    # --- Detect input format ---
    if args.input_format:
        fmt = args.input_format
    else:
        fname_lower = args.input.lower()
        if fname_lower.endswith((".vcf", ".vcf.gz")):
            fmt = "vcf"
        else:
            fmt = "pav"
    log(f"Detected input format: {fmt.upper()}")

    # --- 1. Parse GTF ---
    chrom_data = load_gtf(args.gtf)
    chroms_with_data = set(chrom_data.get("__chroms", []))

    # --- 2. Compute derived features ---
    intron_data = compute_introns(chrom_data) if "intron" in requested_features else {}
    promoter_data = (compute_promoters(chrom_data, args.promoter_length)
                     if "promoter" in requested_features else {})

    # --- 3. Build trees ---
    trees = build_trees(chrom_data, intron_data, promoter_data, requested_features)

    # --- 4. Process ---
    if fmt == "vcf":
        stats, out_paths, annotated_path = process_vcf(
            args.input, trees, sorted(requested_features),
            args.output_prefix,
            no_filtered=args.no_filtered,
            no_annotated=args.no_annotated,
            chroms_with_data=chroms_with_data,
        )

        # --- 5. Report ---
        sep = "=" * 60
        print(f"\n{sep}", file=sys.stderr)
        print("EXTRACTION COMPLETE  (VCF mode)", file=sys.stderr)
        print(sep, file=sys.stderr)
        report_line(f"Total variants scanned   : {stats['total']:>12,}")
        report_line(f"Skipped (chrom not in GTF): {stats['skipped_chrom']:>12,}")
        report_line(f"Intergenic               : {stats['intergenic']:>12,}")
        if not args.no_filtered:
            print(f"\n  Per-feature filtered VCFs:", file=sys.stderr)
            for feat in sorted(requested_features):
                report_line(
                    f"{feat:<24} : {stats['written_filtered'][feat]:>12,}"
                    f"  →  {out_paths[feat]}"
                )
        if not args.no_annotated:
            report_line(f"\n  Full annotated VCF:")
            report_line(f"{'':24}   {'':>12}  →  {annotated_path}")
        print(sep, file=sys.stderr)

    else:  # PAV
        stats, out_paths = process_pav(
            args.input, trees, sorted(requested_features),
            args.output_prefix, chroms_with_data,
        )

        # --- 5. Report ---
        sep = "=" * 60
        print(f"\n{sep}", file=sys.stderr)
        print("EXTRACTION COMPLETE  (PAV mode)", file=sys.stderr)
        print(sep, file=sys.stderr)
        report_line(f"Total variants scanned   : {stats['total']:>12,}")
        report_line(f"Skipped (chrom not in GTF): {stats['skipped_chrom']:>12,}")
        for feat in sorted(requested_features):
            report_line(
                f"{feat:<24} : {stats['written'][feat]:>12,}"
                f"  →  {out_paths[feat]}"
            )
        print(sep, file=sys.stderr)


if __name__ == "__main__":
    main()
