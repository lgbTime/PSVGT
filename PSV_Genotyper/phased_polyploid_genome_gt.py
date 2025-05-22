import pandas as pd
import argparse

def convert_gt(gt_str):
    """
    polyploid gt to stardard gt
    """
    alleles = set(gt_str.split('|'))
    if '0' in alleles and '1' in alleles:
        return "0/1"
    elif '1' in alleles and '.' not in alleles:
        return "1/1"
    elif '1' in alleles and '.' in alleles:
        return "1/1" ## in hifi reads it will be 1/1
    elif '0' in alleles and '.' in alleles:
        return "./0"
    elif '1' not in alleles and '.' not in alleles:
        return "0/0"
    else:
        return ".".join(alleles)

def main():
    parser = argparse.ArgumentParser(description='Process and combine multiple haplotypes to multi-allele genotype table')
    parser.add_argument('hap_paths', nargs='+', help='Paths to haplotype files (e.g., hap1.tsv hap2.tsv hap3.tsv ...)')
    parser.add_argument('out_multiploid', help='Output file for multi-haplotype phased genotype')
    args = parser.parse_args()
    for path in args.hap_paths:
        try:
            with open(path, 'r'):
                pass
        except FileNotFoundError:
            raise SystemExit(f"file {path} not found")
    haps = []
    for idx, path in enumerate(args.hap_paths, 1):
        df = pd.read_csv(path, sep='\t', header=0, index_col=None)
        df['hap_id'] = f"hap{idx}"
        haps.append(df)
    genotype_map = {
        "1/1": "1",
        "0/0": "0",
        "./.": ".",
        "0/1": "0"  # haplotype genome dont have heterozygous GT
    }
    for df in haps:
        gt_col = df.columns[5]
        df['GT'] = df[gt_col].map(genotype_map).fillna(df[gt_col])
    merged = haps[0].copy()
    gt_columns = [df['GT'] for df in haps]
    merged['GT_polyploid'] = merged.apply(lambda row: "|".join([str(gt.iloc[row.name]) for gt in gt_columns]), axis=1)
    merged['GT_standard'] = merged['GT_polyploid'].apply(convert_gt)
    merged[haps[0].columns[5]] = merged['GT_standard']
    merged.drop(columns=['hap_id', 'GT_polyploid', 'GT_standard'], inplace=False, errors='ignore').to_csv(args.out_multiploid, sep='\t', index=False, header=True)
    print(merged.head())
    merged[haps[0].columns[5]] = merged['GT_polyploid']
    merged.drop(columns=['hap_id', 'GT_polyploid', 'GT_standard'], inplace=False, errors='ignore').to_csv(f"{args.out_multiploid}.polyploid.gt", sep='\t', index=False, header=True)
if __name__ == "__main__":
    main()
