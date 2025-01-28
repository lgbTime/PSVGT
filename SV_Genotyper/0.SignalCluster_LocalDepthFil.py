import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing          
def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up, 
    since it did not cause time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov =  covinfo[(covinfo['target_chr']==reference) 
                    & (covinfo['target_start']<= start) 
                    & (covinfo['target_end']>= end)
                    ]
        num_maps = cov.shape[0]
    else:
        num_maps = covinfo.count(reference, start, end)
    return num_maps

def svs_clu(df, covinfo, csv, max_diff, svtype, chrom):
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df['SVlen'] = df['SVlen'].astype(int)
    df['maq'] = df['maq'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start', 'SVlen']).reset_index(drop=True)
    print(df.head())
    grouped = df.groupby('SVID')
    clus = []
    for cs_id, cs in grouped:
        clu = cs.iloc[0, :].tolist()
        clu.append(cs.shape[0])
        clus.append(clu)
    clusdf = pd.DataFrame(clus)
    clusdf.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", 'SVType', 'seq', 'clu_size']
    #quantile80 = clusdf['clu_size'].quantile(0.8)
    clusdf = clusdf.sort_values(by=['#Target_name', 'Target_start', 'SVlen']).reset_index(drop=True)
    clusdf['shift_cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(clusdf)):
        old_len = clusdf.loc[i - 1, 'SVlen']
        now_len = clusdf.loc[i, 'SVlen']
        relate_size = round(now_len / old_len, 2)
        ## the shift should base SVType and Len, for clu_size > averageDepth/3 or quantile80? we just need a small shift breakpoints  ##
        if  old_len <= 5000:
            max_diff = 100 +  0.05*old_len
        if old_len > 5000:
            maf_diff = 350 + 0.005*old_len
        if (abs(clusdf.loc[i, 'Target_start'] - clusdf.loc[i - 1, 'Target_start']) <= max_diff or abs(clusdf.loc[i, 'Target_end'] - clusdf.loc[i - 1, 'Target_end']) <= max_diff) and (0.5 < relate_size < 2.0):
            clusdf.loc[i, 'shift_cluster'] = clusdf.loc[i - 1, 'shift_cluster']
        else:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    print(f'The {chrom} raw {svtype} signal: {df.shape} after clustering by SVID and  shifting breakpoints merging to {len(clusdf["shift_cluster"].unique())}')
    shift2clu = []
    shift_clus = clusdf.groupby('shift_cluster')
    for cluID, cs in shift_clus:
        cluster_size = cs['clu_size'].sum()
        Target_name = chrom
        Target_start = cs.loc[cs['clu_size'].idxmax()].Target_start
        Target_end = cs.loc[cs['clu_size'].idxmax()].Target_end
        start_local_map = local_cov(covinfo, Target_name, max(0, Target_start - 150), max(Target_start - 50, 0))
        end_local_map = local_cov(covinfo, Target_name, Target_end + 50, Target_end + 150)
        if cluster_size >= max(start_local_map, end_local_map) * csv or cluster_size >= 6:
            SVlen = cs.loc[cs['clu_size'].idxmax()].SVlen
            SVID = cs.loc[cs['clu_size'].idxmax()].SVID
            seq = cs.loc[cs['clu_size'].idxmax()].seq
            reads_name = "**"
            meanq = cs['maq'].mean()
            shift2clu.append({
                '#Target_name': Target_name,
                'Target_start': Target_start,
                'Target_end': Target_end,
                'SVlen': SVlen,
                'SVID': SVID,
                'SVType': svtype,
                'seq': seq,
                'cluster_size': cluster_size,
                'Query_name': reads_name,
                'maq': meanq
            })
    print(f'The {chrom} {svtype} {len(clusdf["shift_cluster"].unique())} clusters filter to {len(shift2clu)}')
    return shift2clu


def tra_clu(df, covinfo, csv, chrom, max_diff=1000):
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(int)
    df["Target_start2"] = df["Target_start2"].astype(int)
    df['maq'] = df['maq'].astype(int)
    df.sort_values(by=['SVID'], inplace=True)
    grouped = df.groupby('SVID')
    clus = []
    for cluster_id, cs in grouped:
        clu = cs.iloc[0, :].tolist()
        clu.append(cs.shape[0])
        clus.append(clu)
    clusdf = pd.DataFrame(clus)
    clusdf.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq', '#Target_name2', 'clu_size']
    clusdf.sort_values(by=['#Target_name1', '#Target_name2', 'Target_start1', 'Target_start2'], inplace=True)
    clusdf.index = range(len(clusdf))
    clusdf['cluster_bp1'] = -1  # Initializing cluster
    cluster_id = 0
    clusdf.loc[0, 'cluster_bp1'] = cluster_id
    for i in range(1, len(clusdf)):
        if clusdf.loc[i, '#Target_name2'] == clusdf.loc[i - 1, '#Target_name2']:
            if abs(clusdf.loc[i, 'Target_start1'] - clusdf.loc[i - 1, 'Target_start1']) <= max_diff and abs(clusdf.loc[i, 'Target_start2'] - clusdf.loc[i - 1, 'Target_start2']) <= max_diff:
                clusdf.loc[i, 'cluster_bp1'] = clusdf.loc[i - 1, 'cluster_bp1']
            else:
                cluster_id += 1
                clusdf.loc[i, 'cluster_bp1'] = cluster_id
        else:
            cluster_id += 1
            clusdf.loc[i, 'cluster_bp1'] = cluster_id
    print(f'The {chrom} TRA signal shape: {df.shape} after clustering by shift bp1:{int(max_diff)} is {len(clusdf["cluster_bp1"].unique())}')

    bp1_clus = []
    bp1_grouped = clusdf.groupby('cluster_bp1')
    for g, bp1_c in bp1_grouped:
        bp1_size = bp1_c['clu_size'].sum()
        Target_name1, Target_name2 = bp1_c['#Target_name1'].iloc[0], bp1_c['#Target_name2'].iloc[0]
        query_name, seq, SVType, SVlen = "*", "*", "TRA", 0
        Target_start1, Target_start2 = bp1_c.loc[bp1_c['clu_size'].idxmax()].Target_start1, bp1_c.loc[bp1_c['clu_size'].idxmax()].Target_start2
        meanq, SVID = bp1_c['maq'].mean(), bp1_c.loc[bp1_c['clu_size'].idxmax()].SVID
        bp1_size = bp1_c['clu_size'].sum()
        start1_local_map = local_cov(covinfo, Target_name1, max(Target_start1 - 150, 0), max(Target_start1 - 50, 0))
        start2_local_map = local_cov(covinfo, Target_name2, Target_start2 + 50, Target_start2 + 150)
        if bp1_size >= max(start1_local_map, start2_local_map) * csv:
            bp1_clus.append({
                '#Target_name': Target_name1,
                'Target_start': Target_start1,
                'Target_end': Target_start2,
                'SVlen': SVlen,
                'SVID': SVID,
                'SVType': SVType,
                'seq': seq,
                'cluster_size': bp1_size,
                'Query_name': query_name,
                'maq': meanq
            })
    print(f'The TRA signal clusters {len(clusdf["cluster_bp1"].unique())} filter to {len(bp1_clus)}')
    return bp1_clus


def read_file(file_name):
    if not os.path.exists(file_name) or os.path.getsize(file_name) == 0:
        print(f"The file {file_name} does not exist or is empty.")
        return []
    try:
        return pd.read_csv(file_name, sep="\t", header=None, dtype=str,index_col=None)
    except pd.errors.EmptyDataError:
        print(f"The file {file_name} is empty or malformed.")
        return []


def load_and_process_sv_data(args):
    svindel = read_file(args.raw_signal)
    if len(svindel) == 0:
        return [],[]
    print(svindel.head())
    svindel.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType", "seq"]
    chroms = svindel['#Target_name'].unique()
    svindel['SVlen'] = svindel['SVlen'].astype(int)
    svindel = svindel[svindel["SVlen"] <= args.max]
    sv_data = {sv_type: svindel[svindel['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    msv = read_file(f"{args.raw_signal}.suppAlign")
    if len(msv) > 0:
        msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType", "seq"]
        for svtype in sv_data:
            sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
    return sv_data,chroms


def process_svtype(args,sv_data, chroms, svtype):
    """
    cov info parse by covfile or bam file
    """
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    if args.dtype in ['cr']:
        covinfo = pd.read_csv(args.covfile, sep="\t",index_col=None,header=None)
        covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
        covinfo['target_chr'] = covinfo['target_chr'].astype(str)
    else:
        covinfo = pysam.AlignmentFile(args.bam, "rb")
    for chrom in chroms:
        chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
        def process_sv_type(svtype, shift_multiplier, default_max_diff):
            if svtype == "TRA":
                return tra_clu(chrom_data[svtype], covinfo, args.csv, chrom, shift_multiplier * args.shift) if not chrom_data[svtype].empty else []
            else:
                return svs_clu(chrom_data[svtype], covinfo, args.csv, shift_multiplier * args.shift, svtype, chrom) if not chrom_data[svtype].empty else []
        if svtype == "DEL":
            del_clus = process_sv_type("DEL", 1, 200)
        elif svtype == "INS":
            ins_clus = process_sv_type("INS", 1, 100)
        elif svtype == "INV":
            inv_clus = process_sv_type("INV", 2, args.shift)
        elif svtype == "DUP":
            dup_clus = process_sv_type("DUP", 2, args.shift)
        elif svtype == "TRA":
            tra_clus = process_sv_type("TRA", 1, args.shift)
    return tra_clus, del_clus, ins_clus, inv_clus, dup_clus

def candidateSV(args):
    sv_data,chroms = load_and_process_sv_data(args)
    if len(sv_data) !=0:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            tra_clus += result[0]
            del_clus += result[1]
            ins_clus += result[2]
            inv_clus += result[3]
            dup_clus += result[4]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        return [], [], [], [], []


if __name__ == "__main__":
    import os
    parser = argparse.ArgumentParser("signal filtering through support reads ratio", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True, help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=1000, type=int, help="the distance shift of breakpoint to cluster the TRA/INV/DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=6868686, help="the max SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, help="the sequencing type of samples")
    IN.add_argument("-csv", dest="csv", type=float, default=0.2, help="the paramter to filter the sv signal / local_total < 0.2")
    IN.add_argument("--cov", dest="covfile", type=str, help= "Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    #import cProfile
    #cProfile.run('candidateSV(args)')
    if len(tra_clus) >0:
        pd.DataFrame(tra_clus).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if len(dup_clus) >0:
        pd.DataFrame(dup_clus).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if len(inv_clus) >0:
        pd.DataFrame(inv_clus).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if len(ins_clus) >0:
        pd.DataFrame(ins_clus).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if len(del_clus) >0:
        pd.DataFrame(del_clus).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")
