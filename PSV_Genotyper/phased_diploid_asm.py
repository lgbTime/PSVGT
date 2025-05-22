import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process and combine haplotype data')
    parser.add_argument('hap1path', help='Path to the first haplotype file')
    parser.add_argument('hap2path', help='Path to the second haplotype file')
    parser.add_argument('out_diploid', help='output diploid phased genotype file')
    args = parser.parse_args()

    for path in [args.hap1path, args.hap2path]:
        try:
            with open(path, 'r'):
                pass
        except FileNotFoundError:
            raise SystemExit(f"错误：文件 {path} 不存在")

    hap1 = pd.read_csv(args.hap1path, sep='\t', header=0, index_col=None)
    hap2 = pd.read_csv(args.hap2path, sep='\t', header=0, index_col=None)
    gt_h1 = hap1.columns[5]  
    gt_h2 = hap2.columns[5]
    genotype_map = {
        "1/1": "1",
        "0/0": "0",
        "./.": ".",
        "0/1": "0/1"
    }

    hap1['GT'] = hap1[gt_h1].map(genotype_map).fillna(hap1[gt_h1])  # 未知基因型保留原值
    hap2['GT'] = hap2[gt_h2].map(genotype_map).fillna(hap2[gt_h2])
    # ---------------------- 合并相位基因型 ----------------------
    phased = hap1.copy()
    phased['GT'] = hap1['GT'] + "|" + hap2['GT']
    phase_standard = {
        '0/1|1': '0|1',  
        '0/1|0': '1|0',
        '0|0/1': '0|1',
        '1|0/1': '1|0',
	'0/1|0/1':"0/1"
    }
    phased['GT'] = phased['GT'].replace(phase_standard)
    phased[gt_h1] = phased['GT']
    phased.drop(columns=phased.columns[-1], inplace=True)
    phased.to_csv(args.out_diploid, sep='\t', index=False, header=True)

if __name__ == "__main__":
    main()
