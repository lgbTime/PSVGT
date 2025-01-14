import argparse
from os.path import basename
import pandas as pd
def svindeltab2vcf(tab, outvcf):
    with open(tab, 'r') as f:
        lines = f.readlines()
    with open(outvcf, 'w') as outf:
        header = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{basename(tab)}"
        print(header, file=outf)
        for line in lines:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")
            chrom, refbase, sv_start, sv_end, svlen, gt, maps,rate = info[0], "N", info[1], info[2], abs(int(info[3])), info[5],info[6],info[7]
            sv_type = "DEL" if int(info[3]) > 0 else "INS"
            ID = f'{chrom}:{sv_start}-{sv_end}_{svlen}'
            outline = f"{chrom}\t{sv_start}\t{ID}\tN\t<{sv_type}>\t.\tPASS\t{chrom}:{sv_start}-{sv_end};method=SVInDel;{maps};{rate}\tGT\t{gt}"
            print(outline, file=outf)

def main():
    parser = argparse.ArgumentParser(description="Convert SVInDel tab-delimited file to VCF format")
    parser.add_argument('tab', help='Input tab-delimited file with SVInDel information')
    parser.add_argument('outvcf', help='Output VCF file path')
    args = parser.parse_args()
    svindeltab2vcf(args.tab, args.outvcf)
    vcf = pd.read_csv(args.outvcf,header=0,index_col=None,sep="\t")
    vcf.sort_values(by=["#CHROM", "POS", "ID"],inplace=True)
    vcf.to_csv(args.outvcf,header=True,sep="\t",index=None)
if __name__ == "__main__":
    main()

