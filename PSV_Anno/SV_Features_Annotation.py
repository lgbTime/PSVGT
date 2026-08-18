import argparse
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from copy import deepcopy

def del_overlapper(del_pos, region, chrom, ori_indel, outeach):
    chr_del = del_pos[del_pos["#Target_name"] == chrom]
    dele_pos = chr_del[chr_del["Target_end"] - chr_del["Target_start"] >40]
    print(dele_pos.head())
    print(chr_del.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in dele_pos.index:
        pos1 = dele_pos.loc[i, "Target_start"]
        pos2 = dele_pos.loc[i, "Target_end"]
        svsize = dele_pos.loc[i,"Target_size" ]
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if pos1 <= intervalX.begin and pos2 >= intervalX.end:
                    if "cds" in intervalX.data:
                        impact = "AS"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_lost"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLost\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                elif pos1 < intervalX.begin and pos2 <= intervalX.end and pos2 >intervalX.begin:
                    delLeft = pos2 - intervalX.begin + 1
                    if "cds" in intervalX.data and delLeft % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delLeft % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tLeft_Lost:{delLeft}bp\t{impact}')
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    IDlst.append(f"{intervalX.data}:{impact}")
                elif pos1 > intervalX.begin and pos1 <  intervalX.end and pos2 >intervalX.end:
                    delright = intervalX.end - pos1 + 1
                    if "cds" in intervalX.data and delright % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delright % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smalller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRignt_Lost:{delright}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tRight_Lost:{delright}bp\t{impact}')
                elif pos1 >= intervalX.begin and pos2 <= intervalX.end:
                    delcenter = pos2 - pos1 + 1
                    if "cds" in intervalX.data and delcenter % 3 ==0:
                        impact = "AS"
                    elif "cds" in intervalX.data and delcenter % 3 !=0:
                        impact = "frameshift"
                    elif "utr" in intervalX.data:
                        impact = "Exp"
                    elif "intron" in intervalX.data:
                        impact = "intron_smaller"
                    IDlst.append(f"{intervalX.data}:{impact}")
                    ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")
                    print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}', file=outeach)
                    #print(f'SVDel\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\tCenterLost:{delcenter}bp\t{impact}')


def ins_overlapper(indel_pos, region, chrom, ori_indel, outeach):
    chr_indel = indel_pos[indel_pos["#Target_name"] == chrom]
    ins_pos = chr_indel[chr_indel["Target_end"] - chr_indel["Target_start"] == 1]
    print(chr_indel.head())
    print(ins_pos.head())
    region.columns = ["start","end","GeneID"]
    rtree = IntervalTree()

    for r in region.index:
        rStart = region.loc[r, "start"]
        rEnd   = region.loc[r, "end"]
        if rStart == rEnd:
            print(f"{r}\t {rStart}\t{rEnd}\tError gene struc")
        else:
            rtree[rStart:rEnd] = r
    for i in ins_pos.index:
        pos1 = ins_pos.loc[i, "Target_start"]
        pos2 = ins_pos.loc[i, "Target_end"]
        svsize = abs(ins_pos.loc[i, "Target_size"])
        lapper = rtree.overlap(pos1,pos2)
        if lapper:
            IDlst = []
            for intervalX in lapper:
                if "cds" in intervalX.data and svsize % 3 ==0:
                    impact = "AS"
                elif "cds" in intervalX.data and svsize % 3 !=0:
                    impact = "frameshift"
                elif "utr" in intervalX.data:
                    impact = "Exp"
                elif "intron" in intervalX.data:
                    impact = "intron_exten"
                IDlst.append(f"{intervalX.data}:{impact}")
                print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}', file=outeach)
                #print(f'SVIns\t{svsize}\t{chrom}:{pos1}-{pos2}\t{intervalX.data}\t{intervalX.data.split("__")[0]}\t{impact}')
            ori_indel.loc[i, "struc"] = str(set(IDlst))[1:-1].replace("'", "")

import pandas as pd

def extract_gene_cds(gff_file, geneid, cdsid):
    gff = pd.read_csv(gff_file, header=None, comment="#", sep="\t", index_col=None)
    gff[3] = gff[3].astype(int)
    gff[4] = gff[4].astype(int)
    chroms = gff[0].unique()
    if "gff" in str(gff_file):
        gene = gff[gff[2] == "mRNA"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid}=([^;]+)')  # Using f-string for regex
        gene["ID"] = gene["ID"].str.replace('"', '')
        cds["ID"] = cds[8].str.extract(fr'{cdsid}=([^;]+)') 
        cds["ID"] = cds["ID"].str.replace('"', '')
        
    elif "gtf" in str(gff_file):
        gene = gff[gff[2] == "transcript"]
        cds = gff[gff[2] == "CDS"]
        gene["ID"] = gene[8].str.extract(fr'{geneid} ([^;]+)')  
        cds["ID"] = cds[8].str.extract(fr'{cdsid} ([^;]+)')  
    return gene, cds, chroms

def gene_feature(geneid, genedf, cdsdf):
    genes_dict = {}
    gene_intron_regions = {}
    gene_cds_regions = {}
    utr5_region = {}
    utr3_region = {}
    geneIDdf = genedf[genedf["ID"] == geneid]
    geneIDcds = cdsdf[cdsdf["ID"] == geneid]
    cds_regions = []
    ### to build interval tree we dont take region strandness ##
    gene_chr = geneIDdf.iloc[0,0]
    gene_start = geneIDdf.iloc[0, 3]
    gene_end = geneIDdf.iloc[0, 4]
    # Collect CDS regions
    for index in geneIDcds.index:
        cds_region = [geneIDcds.loc[index, 3], geneIDcds.loc[index, 4]]
        cds_regions.append(cds_region)
    gene_cds_regions[geneid] = cds_regions
    # Calculate UTR regions
    first_cds_start = geneIDcds.iloc[0, 3]
    last_cds_end = geneIDcds.iloc[-1, 4]
    gene_strand = geneIDdf.iloc[0, 6]
    utr5 = ""
    utr3 = ""
    if gene_strand == "+":
        if int(gene_start) < int(first_cds_start):
            utr5 = [gene_start, first_cds_start]
        else:
            utr5 = "UTR5 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr3 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5
    elif gene_strand == "-":
        # Calculate UTR regions
        if int(gene_start) < int(first_cds_start):
            utr3 = [gene_start, first_cds_start]
        else:
            utr3 = "UTR3 Info Lost"
        if int(last_cds_end) < int(gene_end):
            utr5 = [last_cds_end, gene_end]
        else:
            utr3 = "UTR3 Info Lost"
    # Store UTR regions
        utr3_region[geneid] = utr3
        utr5_region[geneid] = utr5

    # Calculate intron regions
    intron_regions = []
    if len(cds_regions) > 1:
        for i in range(len(cds_regions) - 1):
            if gene_strand == "+":
                # Intron is between the end of the current CDS and the start of the next CDS
                intron_start = cds_regions[i][1] + 1
                intron_end = cds_regions[i + 1][0] - 1
            else:  # "-" strand
                # Intron is between the end of the next CDS and the start of the current CDS
                intron_start = cds_regions[i + 1][1] + 1
                intron_end = cds_regions[i][0] - 1
            if intron_start <= intron_end:  # Only add valid intron regions
                intron_regions.append([intron_start, intron_end])
    gene_intron_regions[geneid] = intron_regions
    return gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr

############## vcf format to below df ########################
# head ../final.gt.txt 
# POS	#Target_name	Target_start	Target_end	Target_size	Query_sizeERR11436000_1.fastq_dedup	ERR11436001_1.fastq_dedup
# Chr1:17145-17369	Chr1	17145	17369	225	-225	0/0	0/0
# Chr1:132223-132224	Chr1	132223	132224	-127	127	1/1	1/1

def main(gff_file, sv_file, output_file, geneid, cdsid):
    vcf = pd.read_csv(sv_file,sep="\t",header=0,index_col=None)
    use = ['ID', '#CHROM','POS','ALT']
    for acc in vcf.columns[9:]:
        use.append(acc)
    
    sv = vcf[use]
    ## insert target_end
    sv.insert(3,'Target_end',sv['ID'].str.split(":",expand=True)[1])
    sv['Target_end'] = sv['Target_end'].str.split("-|_",expand=True)[1].astype(int)
    ## Target size and query_size
    sv.insert(4,'Target_size',sv['ID'].str.split(":",expand=True)[1])
    sv["Target_size"] = sv["Target_size"].str.split("-|_",expand=True)[2].astype(int)
    sv.insert(5,'Query_size',sv['ID'].str.split(":",expand=True)[1])
    sv["Query_size"] = sv["Query_size"].str.split("-|_",expand=True)[2].astype(int)
    for index in sv.index:
        if sv.loc[index,'ALT'] == "<DEL>":
            sv.loc[index,'Query_size'] = -sv.loc[index,'Query_size']
        elif sv.loc[index,'ALT'] == "<INS>":
            sv.loc[index,'Target_size'] = -sv.loc[index,'Target_size']
    print(sv.head())
    cols = ['POS','#Target_name','Target_start','Target_end','Target_size','Query_size','SV'] 
    for acc in vcf.columns[9:]:
        cols.append(acc)
    sv.columns = cols
    print(sv.head())
    sv.to_csv(f'{sv_file}_tmp.tab',header=True,sep="\t",index=None)

    #print(sv.iloc[0,[3]])
    genedf, cdsdf, chroms = extract_gene_cds(gff_file, geneid, cdsid)
    genedf["ID"] = genedf["ID"].str.replace('"', "")
    cdsdf["ID"] = cdsdf["ID"].str.replace('"', "")
    print(cdsdf.head())
    print(genedf.head())
    # Prepare output DataFrame for SV
    outsv = deepcopy(sv)  # Assuming sv is defined globally or passed as an argument
    outsv["struc"] = ""

    with open(output_file, "w") as outeach:
        print('SV\tSVsize\tChr:Start-End\tGeneFeature\tGeneID\tImpact', file=outeach)
        for chrom in chroms:
            bdict = {}
            gdict = {}
            genedf_chr = genedf[genedf[0] == chrom]
            cdsdf_chr = cdsdf[cdsdf[0] == chrom]
            for gene in list(set(cdsdf_chr["ID"].tolist())):
                gene_cds_regions, gene_intron_regions, utr5_region, utr3_region, gene_chr = gene_feature(gene, genedf_chr, cdsdf_chr)
                for i, cds_region in enumerate(gene_cds_regions[gene]):
                    bdict[f"{gene}__cds{i + 1}"] = cds_region
                for i, intron_region in enumerate(gene_intron_regions[gene]):
                    bdict[f'{gene}__intron{i + 1}'] = intron_region
                    if gene in utr3_region.keys():
                        if utr3_region[gene] != "UTR3 Info Lost":
                            bdict[f'{gene}__utr3'] = utr3_region[gene]
                    if gene in  utr5_region.keys():
                        if utr5_region[gene] != "UTR5 Info Lost":
                            bdict[f'{gene}__utr5'] = utr5_region[gene]
            if bdict:
                gene_body = pd.DataFrame(bdict).T
                gene_body["gene"] = gene_body.index
            
                gene_body["gene"] = gene_body["gene"].str.split("__", expand=True)[0]
                gene_body.dropna(inplace=True)
                print(f"********************* the gene_body format df **********************\n{gene_body.head()}")
                del_overlapper(del_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
                ins_overlapper(indel_pos=sv, region=gene_body, chrom=chrom, ori_indel=outsv, outeach=outeach)
    outsv.to_csv(f"{args.sv}_anno.txt",header=True,index=None,sep="\t")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process gene and structural variant data.")
    parser.add_argument('-g', '--gff', required=True, help='Path to the GFF file (e.g., TAIR10_gene.gff)')
    parser.add_argument('-s', '--sv', required=True, help='Path to the sv vcf file')
    parser.add_argument('-m', '--mRNA', default="ID", help='to extract GeneID, some gff is not in ID=Geneid, here default is ID')
    parser.add_argument('-c', '--cds', default="Parent", help='to extract cds feature coordinated GeneID, some gff is not in Parent=Geneid, here default is Parent')
    parser.add_argument('-o', '--output', default='SVInDels_Lead_Gene_Variants.txt', help='Output file name (default: SVInDels_Lead_Gene_Variants.txt)')
    args = parser.parse_args()
    main(args.gff, args.sv, args.output, args.mRNA, args.cds)

