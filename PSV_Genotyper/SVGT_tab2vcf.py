import argparse
from os.path import basename
import pandas as pd
def svindeltab2vcf(tab, outvcf):
    with open(tab, 'r') as f:
        lines = f.readlines()
    with open(outvcf, 'w') as outf:
        prehead = f"""##fileformat=VCFv4.2
##source=PSVGT1.0
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##INFO=<ID=CHR,Number=1,Type=String,Description="Chromosome">
##INFO=<ID=SV_START,Number=1,Type=Integer,Description="Start position of the structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=maq,Number=1,Type=Integer,Description="The mean mapping quality of SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length difference between the REF and ALT alleles">
##INFO=<ID=method,Number=1,Type=String,Description="Method used to call the variant">
##INFO=<ID=total_map_reads,Number=0,Type=Integer,Description="Total number of mapped reads supporting the variant">
##INFO=<ID=INS_rate,Number=1,Type=Float,Description="Insertion rate (e.g., read depth support)">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (e.g., INS means Insertion)">
##INFO=<ID=total_map_reads_l,Number=1,Type=Integer,Description="Total mapping reads at first breakpoints">
##INFO=<ID=total_map_reads_r,Number=1,Type=Integer,Description="Total mapping reads at second breakpoints">
##INFO=<ID=total_map_reads_bp1,Number=1,Type=Integer,Description="Total mapping reads at first breakpoints">
##INFO=<ID=total_map_reads_bp2,Number=1,Type=Integer,Description="Total mapping reads at second breakpoints">
##INFO=<ID=deles_l_ratio,Number=1,Type=Float,Description="The delestion ratio detected at left breakpoints">
##INFO=<ID=deles_r_ratio,Number=1,Type=Float,Description="The delestion ratio detected at right breakpoints">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{basename(tab)}"""
        print(prehead,file=outf)
        svlist = []
        for line in lines:
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")
            chrom, refbase, sv_start, sv_end, svlen, gt, maps,rate,sv_type = info[0], "N", info[1], info[2], abs(int(info[3])), info[5], info[6], info[7].split(";")[0], info[7].split(";")[1]
            ID = f'{chrom}:{sv_start}-{sv_end}_{svlen}'
            outline = f"{chrom}\t{sv_start}\t{ID}\tN\t<{sv_type}>\t.\tPASS\tEND={sv_end};CHR={chrom};SV_START={sv_start};SVLEN={svlen};method=PSVGT;{maps};{rate};SVTYPE={sv_type}\tGT\t{gt}"
            svlist.append(outline.split("\t"))
        df = pd.DataFrame(svlist)
        df[1] = df[1].astype(int)
        df.sort_values(by=[0, 1, 2],inplace=True)
    df.to_csv(outvcf,header=None,sep="\t",index=None,mode='a')

def main():
    parser = argparse.ArgumentParser(description="Convert SVInDel tab-delimited file to VCF format")
    parser.add_argument('tab', help='Input tab-delimited file with SVInDel information')
    parser.add_argument('outvcf', help='Output VCF file path')
    args = parser.parse_args()
    svindeltab2vcf(args.tab, args.outvcf)
if __name__ == "__main__":
    main()

