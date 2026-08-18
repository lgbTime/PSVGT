import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

def parse_args():
    parser = argparse.ArgumentParser(description='Process GFF and indel files')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file or gtf')
    parser.add_argument('pos_file', type=str, help='Path to the pos file, which must have the 3 columns (#Target_name	Target_start	Target_end)')
    parser.add_argument('output_file', type=str, help='the output prefix name')
    return parser.parse_args()

def load_gff(gff_file):
    gff = pd.read_csv(gff_file, header=None,comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = list(set(gff[0].tolist()))
    return gff, chroms


def extract_gene_cds(gff, gff_file):
    if "gff" in str(gff_file):
        gene = gff[gff[2]=="gene"]
        cds = gff[gff[2]=="CDS"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'ID=([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'Parent=([^;]+)')
    elif "gtf" in str(gff_file):
        gene = gff[gff[2]=="transcript"]
        cds = gff[gff[2]=="exon"]
        cds = cds[cds[3] < cds[4]]
        gene["ID"] = gene[8].str.extract(r'gene_name ([^;]+)')
        cds["ID"] =  cds[8].str.extract(r'gene_name ([^;]+)')

    return gene, cds

def get_3k_promoter(gene):
    pro3k = gene.copy()
    for i in range(pro3k.shape[0]):
        if pro3k.iloc[i,6] == "+":
            pro3k.iloc[i, 4] = pro3k.iloc[i, 3] - 1
            pro3k.iloc[i, 3] = pro3k.iloc[i, 3] - 3000
        else:
            pro3k.iloc[i,3] = pro3k.iloc[i,4] + 1
            pro3k.iloc[i,4] = pro3k.iloc[i,4] + 3000
        pro3k.iloc[i,2] = "promoter"
    return pro3k

def load_indel(indel_file):
    ori_indel = pd.read_csv(indel_file, header=0, sep="\t", index_col=None)
    ori_indel["gene"] = ''
    ori_indel["CDS"] = ''
    ori_indel["3k_promoter"] = ''
    return ori_indel[["#Target_name", "Target_start", "Target_end"]], ori_indel
def overlapper(indel_pos, region, struc, chrom, ori_indel, outeach):
    indel_pos = indel_pos[indel_pos["#Target_name"] == chrom]
    region.columns = ["chr","struc","start","end","ID"]
    region = region[region["chr"] == chrom]
    rtree = IntervalTree()
    
    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        rtree[rStart:rEnd] = region.loc[r, "ID"].replace('"', "")

    for i in indel_pos.index:
        pos1 = indel_pos.loc[i, "Target_start"]
        pos2 = indel_pos.loc[i, "Target_end"]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    print(f'GenePAV\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")

                else:
                    print(f'{struc}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}', file=outeach)
                    IDlst.append(intervalX.data)
                    if struc != "CDS":
                        ori_indel.loc[i, struc] = str(IDlst)[1:-1].replace("'", "")
                    elif struc == "CDS":
                        ori_indel.loc[i, struc] = str(set(IDlst))[1:-1].replace("'", "")
def main():
    args = parse_args()
    gff_file = args.gff_file
    indel_file = args.pos_file
    output_file = args.output_file
    gff,chroms = load_gff(gff_file)
    gene, cds = extract_gene_cds(gff, gff_file)
    pro3k = get_3k_promoter(gene)
    
    gene = gene[[0, 2,3,4,"ID"]]
    gene.index = range(gene.shape[0])
    cds = cds[[0, 2,3,4,"ID"]]
    cds.index = range(cds.shape[0])
    pro3k = pro3k[[0,2,3,4,"ID"]]
    pro3k.index = range(pro3k.shape[0])

    indel, ori_indel = load_indel(indel_file)
    outeach = open(f'{output_file}_each_gene.txt', "w")
    for chrom in tqdm(chroms):
        overlapper(indel, gene, "gene", chrom, ori_indel, outeach)
        overlapper(indel, cds, "CDS", chrom, ori_indel, outeach)
        overlapper(indel, pro3k, "3k_promoter", chrom, ori_indel, outeach)

    ori_indel.to_csv(f"{output_file}_pos.txt", header=1, index=None, sep="\t", quoting=None)
    outeach.close()
if __name__ == '__main__':
    main()
