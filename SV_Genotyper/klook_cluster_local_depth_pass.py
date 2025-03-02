import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil


def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cause time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        num_maps = covinfo.count(str(reference), start, end)
    return num_maps


def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        subset_mode = subset.mode().iloc[0]
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)


def candidate_sv(clusdf, opened_bam, support_rate=0.1,add=1):
    """
    to solve the clusters dataframe
    Get: each cluster's sv breakpoints, SV length
    for low depth clus, the cluster must be vary strict and intense, else if will increase false positive a lot. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts()
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.1:
        return []

    top_2clusters = cluster_counts.nlargest(2).index
    if len(top_2clusters) == 2:
        first_count = cluster_counts[top_2clusters[0]]
        second_count = cluster_counts[top_2clusters[1]]
        if second_count >= 0.2 * first_count:
            clusters_to_process = top_2clusters
        else:
            clusters_to_process = [top_2clusters[0]]
    else:
        clusters_to_process = top_2clusters
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        if len(reads_total) < 2:
            continue
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = sv_start + sv_len if svtype != 'INS' else sv_start + 1
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 150), max(sv_start - 50, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 50, sv_end + 150)
        SV_rate = round(sv_eye / np.mean([start_local_map,end_local_map]), 2)
        if sv_eye < (np.mean([start_local_map,end_local_map])*support_rate + add):
            continue
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
    if len(clus) == 2:
        print('found multiple sv allele')
    return clus


def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'first',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index()
    clusdf = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    clusdf.sort_values(by=['Target_start'], inplace=True)
    clusdf.reset_index(drop=True, inplace=True)
    clusdf['shift_cluster'] = -1
    return clusdf


def klook_clusters(clusdf, max_diff_func, len_condition_range):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(clusdf)):
        current_svlen = clusdf.loc[i, 'SVlen']
        current_start = clusdf.loc[i, 'Target_start']
        current_end = clusdf.loc[i, 'Target_end']
        found_cluster = False
        start_index = max(0, i - 4)
        for j in range(i - 1, start_index - 1, -1):
            old_len = clusdf.loc[j, 'SVlen']
            relate_size = round(current_svlen / old_len, 2)
            max_diff = max_diff_func(old_len)
            len_condition = len_condition_range[0] < relate_size < len_condition_range[1]
            pos_condition = (abs(current_start - clusdf.loc[j, 'Target_start']) <= max_diff or
                             abs(current_end - clusdf.loc[j, 'Target_end']) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf


def lowdepth_clu(clusdf, opened_bam, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    def max_diff_func(old_len):
        if old_len <= 100:
            return 25
        elif 100 < old_len <= 500:
            return 50
        elif old_len > 500:
            return 100

    clusdf = klook_clusters(clusdf, max_diff_func, (0.8, 1.2))
    print(clusdf.drop(columns=['seq']))
    candisv = candidate_sv(clusdf, opened_bam, support_rate, add)
    return candisv


def highdepth_clu(clusdf, opened_bam, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)

    def max_diff_func(old_len):
        if old_len <= 100:
            return 150
        elif 100 < old_len <= 500:
            return 300
        elif old_len > 500:
            return 800

    clusdf = klook_clusters(clusdf, max_diff_func, (0.65, 1.5))
    print(clusdf.drop(columns=['seq']))
    candisv = candidate_sv(clusdf, opened_bam, support_rate, add)
    return candisv


def windows_slide(dfs, depth):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if current_svlen < 100:
            window_size = 200 + current_svlen * 0.2
        elif 100 < current_svlen <= 500:
            window_size = 300 + current_svlen * 0.2
        elif 500 < current_svlen <= 1000:
            window_size = 400 + current_svlen * 0.2
        else:
            window_size = 800
        window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 4:
            windows[current_start] = window_df
        i = j
    return windows


def candidate_tra(win, opened_bam):
    """
    Process translocation data and find candidate translocations.
    """
    tra = []
    chr1 = win['#Target_name1'].iloc[0]
    chr1_start = mode_or_median(win['Target_start1'])
    sv_len = 0
    chr2 = win['#Target_name2'].iloc[0]
    chr2_start = mode_or_median(win['Target_start2'])
    maq = mode_or_median(win['maq'])
    readsname = set(win['Query_name'])
    svid = f'{chr2}:{chr2_start}_{chr1}:{chr1_start}'
    sv_eye = len(readsname)
    start_local_map = local_cov(opened_bam, chr1, max(0, chr1_start - 250), max(chr1_start - 150, 0))
    end_local_map = local_cov(opened_bam, chr2, chr2_start + 150, chr2_start + 250)
    if min(start_local_map, end_local_map) == 0 and max(start_local_map, end_local_map) > 0:
        SV_rate = round(sv_eye / max(start_local_map, end_local_map), 2)
    else:
        SV_rate = round(sv_eye / min(start_local_map, end_local_map), 2)
    if SV_rate < 0.2:
        return None
    tra.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tra


def load_and_process_sv_data(args):
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else float(depth_stat.iloc[0, 3])
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"Error reading {supp_align_file}: {e}")
        else:
            if not msv.empty:
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    return sv_data, chroms, depth


def windows_slide4tra(df, shift=800, opened_bam=None):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********\n{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs.sort_values(by=['Target_start1', 'Target_start2'], inplace=True)
        dfs.reset_index(drop=True, inplace=True)
        i = 0
        while i < len(dfs):
            current_chrom = str(dfs.at[i, '#Target_name1'])
            current_start = dfs.at[i, 'Target_start1']
            window_size = shift
            window_end = current_start + window_size
            j = i
            while j < len(dfs) and dfs.at[j, 'Target_start1'] < window_end:
                j += 1
            window_df = dfs[i:j].copy()
            if len(window_df) >= 2:
                windows[current_start] = window_df
            i = j
    return windows

def process_svtype(args, sv_data, chroms, svtype, depth):
    """
        cov info parse by covfile or bam file
    """
    if args.dtype in ['cr', 'sr']:
        try:
            covinfo = pd.read_csv(args.covfile, sep="\t", index_col=None, header=None)
            covinfo.columns = ['query_chr', 'flag', 'target_chr', 'target_start', 'target_end', 'maq', 'cigar']
            covinfo['target_chr'] = covinfo['target_chr'].astype(str)
            bam_path = covinfo
        except FileNotFoundError:
            print(f"Error: Coverage file {args.covfile} not found.")
            return [], [], [], [], []
    else:
        bam_path = args.bam

    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
    try:
        with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}

                def process_sv_type(svtype):
                    if svtype == "TRA":
                        windows = windows_slide4tra(chrom_data[svtype], args.shift, opened_bam)
                        tra_list = []
                        for win in windows.values():
                            tra = candidate_tra(win, opened_bam)
                            if tra:
                                tra_list.extend(tra)
                        return tra_list
                    else:
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] > 45]
                        sv_windows = windows_slide(sv_dfs, depth)
                        candidate_svs = []
                        for key in sv_windows.keys():
                            sv_window = sv_windows[key]
                            if len(sv_window['Query_name']) > 10:
                                candisv = highdepth_clu(sv_window, opened_bam, 0.1, 2)
                            else:
                                candisv = lowdepth_clu(sv_window, opened_bam,0.1, 1)
                            candidate_svs.extend(candisv)
                        return candidate_svs

                if svtype == "DEL":
                    del_clus.extend(process_sv_type("DEL"))
                elif svtype == "INS":
                    ins_clus.extend(process_sv_type("INS"))
                elif svtype == "INV":
                    inv_clus.extend(process_sv_type("INV"))
                elif svtype == "DUP":
                    dup_clus.extend(process_sv_type("DUP"))
                elif svtype == "TRA":
                    tra_clus.extend(process_sv_type("TRA"))
    except Exception as e:
        print(f"Error processing BAM file: {e}")

    return tra_clus, del_clus, ins_clus, inv_clus, dup_clus


def candidateSV(args):
    sv_data, chroms, depth = load_and_process_sv_data(args)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth) for svtype in sv_types])
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
    parser = argparse.ArgumentParser("signal filtering through support reads ratio",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input File ")
    IN.add_argument("-f", dest="raw_signal", required=True,
                    help="the raw sv signal record file from 0 step signalling")
    IN.add_argument("-s", dest="shift", default=800, type=int,
                    help="the distance shift of breakpoint to cluster the TRA/big_INV/big_DUP signal")
    IN.add_argument("-M", dest="max", type=int, default=18888888, help="the max SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("-csv", dest="csv", type=float, default=0.15,
                    help="the paramter to filter the sv signal / local_total < 0.15 if no depth file found")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")
