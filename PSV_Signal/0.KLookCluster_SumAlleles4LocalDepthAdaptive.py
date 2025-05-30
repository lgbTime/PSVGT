import argparse
from time import time
import pandas as pd
import pysam
import multiprocessing
import numpy as np
import os
from math import ceil,floor

def local_cov(covinfo, reference, start, end):
    """
    For Contig samples may be cov file of txt will be better speed up,
    since it did not cost time on open a bam that has long seq mapping
    """
    if isinstance(covinfo, pd.DataFrame):
        cov = covinfo[(covinfo['target_chr'] == reference)
                      & (covinfo['target_start'] <= start)
                      & (covinfo['target_end'] >= end)]
        num_maps = cov.shape[0]
    else:
        try:
            num_maps = covinfo.count(str(reference), start, end)
        except AttributeError:
            print(f"Error: covinfo does not have a 'count' method.")
            num_maps = 0
    return num_maps

def mode_or_median(series, lower_percentile=0.25, upper_percentile=0.75):
    n = len(series)
    lower_value = series.quantile(lower_percentile, interpolation='linear')
    upper_value = series.quantile(upper_percentile, interpolation='linear')
    subset = series[(series >= lower_value) & (series <= upper_value)]
    if not subset.empty:
        subset_median = subset.median()
        subset_mean = subset.mean()
        return ceil(subset_mean)
    else:
        mode_value = series.mode().iloc[0] if not series.mode().empty else np.nan
        median_value = np.nanmedian(series)
        mean_value = np.nanmean(series)
        stats = [mode_value, median_value, mean_value]
        valid_stats = [s for s in stats if not pd.isna(s)]
        return ceil(mean_value)

def candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate=0.1,add=1):
    """
    Solve clusters dataframe
    Get each cluster's sv breakpoints, SV length
    For low depth clus, the cluster must be vary strict and intense, else it will generate a lot of false positive sv. 
    """
    cluster_col = clusdf.columns[-1]
    sv_chrom = clusdf['#Target_name'].iloc[0]
    svtype = clusdf['SVType'].iloc[0]
    clus = []
    cluster_counts = clusdf[cluster_col].value_counts() ## default reversed count
    min_clusters = min(num_hap, len(cluster_counts))
    max_count = cluster_counts.max()
    proportion = max_count / len(clusdf)
    if proportion < 0.1:
        print("proportion is too low, return []")
        return []
    else: 
        print(f"top1 cluster percent is {proportion}")
    
    top_clusters = cluster_counts.head(min_clusters).index
    clusdf = clusdf[clusdf[cluster_col].isin(top_clusters)]
    Total_signal = clusdf.shape[0]
    clusters_to_process =  clusdf[cluster_col].unique()

    msv = []
    for clu in clusters_to_process:
        clu_df = clusdf[clusdf[cluster_col] == clu]
        reads_total = clu_df['Query_name'].unique()
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len = mode_or_median(clu_df['SVlen'])
        sv_end = mode_or_median(clu_df['Target_end'])
        maq = mode_or_median(clu_df['maq'])
        readsname = set(clu_df['Query_name'].tolist())
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        if sv_eye < 2:
            continue
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 150), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 150)
        
        SV_rate = round(sv_eye / min(start_local_map+1,end_local_map+1), 2)
        if sv_eye >= min([start_local_map,end_local_map])*support_rate + add:
            print(f"{svid} sv signal propertion support")
            clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
            msv.append(svid)
        if sv_len > 2000 and svtype=="INS" and sv_eye >= 2 and sv_eye >= min([start_local_map,end_local_map])*support_rate*0.5: ## to capture big INS
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid} big INS,local_map is {min([start_local_map,end_local_map])} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
        if sv_eye >= nreads -1: ## filter 
            if [sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname] not in clus:
                print(f"{sv_eye} reads support {svid},local coverage is {min([start_local_map,end_local_map])} ")
                clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
                msv.append(svid)
            else:
                print(f"{svid}: number sv reads lower than required")
    if len(set(msv)) >= 2:
        outmsv = "\t".join(list(set(msv)))
        print(f'msv:\t{outmsv}')
    return clus

def klook_clusters(clusdf, max_diff_func, len_fold=0.8):
    """
    Assign cluster IDs to each row in the cluster dataframe based on length and position conditions.
    """
    num_signals = len(clusdf)
    if num_signals == 1:
        clusdf['shift_cluster'] = -1
        return clusdf
    
    # Precompute the 'Target_start' and 'Target_end' columns for quicker access
    start_values = clusdf['Target_start'].values
    end_values = clusdf['Target_end'].values
    svlen_values = clusdf['SVlen'].values
    
    cluster_id = 0
    clusdf.loc[0, 'shift_cluster'] = 0
    
    for i in range(1, len(clusdf)):
        current_svlen = svlen_values[i]
        current_start = start_values[i]
        current_end = end_values[i]
        found_cluster = False
        # Adjust the previous cluster range dynamically
        start_index = max(0, i - ceil(num_signals * 0.8))  # previous 80% record

        # Vectorized approach: avoid nested for loops
        for j in range(i - 1, start_index - 1, -1):
            old_len = svlen_values[j]
            relate_size = min(current_svlen / old_len, old_len/current_svlen)
            max_diff = max_diff_func(old_len)
            len_condition = relate_size > len_fold
            pos_condition = (abs(current_start - start_values[j]) <= max_diff or
                             abs(current_end - end_values[j]) <= max_diff)
            if len_condition and pos_condition:
                clusdf.loc[i, 'shift_cluster'] = clusdf.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            clusdf.loc[i, 'shift_cluster'] = cluster_id
    return clusdf

def max_diff_func4DEL(old_len):
    if old_len <= 100:
        return 50 + ceil((old_len-50) * 0.8)
    elif 100 < old_len <= 500:
        return 100 + ceil((old_len-100) * 0.8)
    elif 500 < old_len <= 5000:
        return 450 + ceil((old_len-500) * 0.1)
    elif old_len > 5000:
        return 1000

def max_diff_func4INS(old_len):
    if old_len <= 100:
        return 50 + (old_len-50) * 0.5
    elif 100 < old_len <= 500:
        return 100 + (old_len-100) * 0.5
    elif 500 < old_len <= 5000:
        return 300 + (old_len-500) * 0.1
    elif old_len > 5000:
        return 1000

def onedepth_all_clus(all_signal,svtype, opened_bam, nrate=0.25):
    """
    Process the contig or genome level SV signal.
    Strict condition for merge, single chromosome mode.
    """
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    if all_signal.shape[0] > 1:
        all_signal = merge_and_sort_cr(all_signal)
        print(f"******************signal after merging ****************\n {all_signal}")
        if all_signal.shape[0] > 1:
            all_clus = klook_clusters(all_signal, max_diff_func, 0.8)
        else:
            all_signal['shift_cluster'] = -1
            all_clus = all_signal
    else:
        all_signal['shift_cluster'] = -1
        all_clus = all_signal
    print(all_clus[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    sv_chrom = all_clus["#Target_name"].unique()[0]
    svtype = all_clus["SVType"].unique()[0]
    cluster_col = all_clus.columns[-1]
    clus = []
    msv = []
    for clu in all_clus[cluster_col].unique():
        clu_df = all_clus[all_clus[cluster_col] == clu] 
        print(clu_df[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])   
        sv_start = mode_or_median(clu_df['Target_start'])
        sv_len =   mode_or_median(clu_df['SVlen'])
        sv_end =   mode_or_median(clu_df['Target_end'])
        maq =      mode_or_median(clu_df['maq'])
        readsname = clu_df['Query_name'].tolist()
        svid = f'{sv_chrom}:{sv_start}-{sv_end}_{svtype}={sv_len}'
        print(svid)
        sv_eye = len(readsname)
        start_local_map = local_cov(opened_bam, sv_chrom, max(0, sv_start - 250), max(sv_start - 100, 0))
        end_local_map = local_cov(opened_bam, sv_chrom, sv_end + 100, sv_end + 250)
        depth = max([start_local_map,end_local_map])
        if depth > 0:
            SV_rate = round(sv_eye / depth, 2)
        else:
            SV_rate = 1
        if SV_rate < nrate: ## 0.25 to meet Tetraploid heterozygous
            continue
        print(f"************************cluster done***************************\n{[sv_chrom, sv_start, sv_end, sv_len, svid, svtype, '*', sv_eye, SV_rate, maq, readsname]}")
        clus.append([sv_chrom, sv_start, sv_end, sv_len, svid, svtype, "*", sv_eye, SV_rate, maq, readsname])
        msv.append(svid)
    if len(msv) >=2:
        outmsv = "\t".join(msv)
        print(f'msv:\t{outmsv}')
    return clus

def lowdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=1):
    """
    Process low-depth cluster data.
    Strict condition for clustering.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL    

    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])

    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def highdepth_clu(clusdf, num_hap, svtype, opened_bam, nreads, support_rate=0.1, add=2):
    """
    Process high-depth cluster data.
    """
    clusdf = merge_and_sort(clusdf)
    
    if svtype == "INS":
        max_diff_func= max_diff_func4INS
    else:
        max_diff_func= max_diff_func4DEL
    
    clusdf = klook_clusters(clusdf, max_diff_func, 0.8)
    print(clusdf[['#Target_name','Target_start','Target_end','SVlen','shift_cluster']])
    candisv = candidate_sv(clusdf, num_hap, opened_bam, nreads, support_rate, add)
    return candisv

def merge_and_sort(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    if two signal is from same reads, signals in window will be summed.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates( subset=["Query_name", "Target_start", "Target_end"], keep="first",inplace=True) ## we dont keep the same reads multiple segment
    print(clusdf[["#Target_name","Target_start","Target_end","SVlen","SVType"]])
    clusdf.sort_values(by=['Target_start', 'SVlen'], inplace=True)
    merged_sv = clusdf.groupby('Query_name').agg({
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
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

def merge_and_sort_cr(clusdf):
    """
    Merge and sort the cluster dataframe based on specific columns.
    Optimized version for better performance with large datasets.
    """
    clusdf = clusdf.copy()
    clusdf.drop_duplicates(subset=["Query_name", "Target_start", "Target_end"], keep="first", inplace=True)
    clusdf.sort_values(['Query_name', 'Target_start'], inplace=True)
    clusdf['same_query'] = clusdf['Query_name'] == clusdf['Query_name'].shift(1)
    clusdf['overlap'] = (clusdf['Target_start'] < clusdf['Target_end'].shift(1)) & clusdf['same_query']
    clusdf['group_key'] = (clusdf['Query_name'] != clusdf['Query_name'].shift(1)) | clusdf['overlap']
    
    clusdf['group_key'] = clusdf['group_key'].cumsum()
    merged_sv = clusdf.groupby('group_key').agg({
        'Query_name': 'first',
        '#Target_name': 'first',
        'Target_start': 'first',
        'Target_end': 'max',
        'SVlen': 'sum',
        'SVType': 'first',
        'maq': 'first',
        'seq': 'first'
    }).reset_index(drop=True)
    merged_sv = merged_sv.astype({'Target_start': np.int32, 'Target_end': np.int32, 'SVlen': np.int32, 'maq': np.int16})
    merged_sv.sort_values(by=['Target_start'], inplace=True)
    merged_sv.reset_index(drop=True, inplace=True)
    merged_sv['shift_cluster'] = -1
    return merged_sv


def windows_slide4asm(dfs, svtype, window_size=500):
    """
    fix window gap for assemblies
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        j = i  
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom:
            last_in_window_start = dfs.at[j-1, 'Target_start'] if j > i else current_start
            if svtype == 'DEL':
                dynamic_window_end = last_in_window_start + window_size + 250
            else:
                dynamic_window_end = last_in_window_start + window_size
    
            if dfs.at[j, 'Target_start'] < dynamic_window_end:
                j += 1  
            else:
                break  
        window_df = dfs[i:j].copy()
        if 1 <= len(window_df['Query_name'].unique()):
            window_df.index = range(len(window_df))  
            windows[current_start] = window_df  
        i = j  
    return windows

def windows_slide(dfs, depth, svtype, nreads=2,window_size=500):
    """
    Slide windows over the dataframe and collect valid windows.
    """
    dfs = dfs.sort_values(by=['Target_start', 'SVlen'])
    dfs.index = range(len(dfs))
    windows = {}
    i = 0
    while i < len(dfs):
        current_chrom = str(dfs.at[i, '#Target_name'])
        current_start = dfs.at[i, 'Target_start']
        current_svlen = dfs.at[i, 'SVlen']
        if svtype == 'DEL':
            window_end = current_start + window_size + 250 ## avoid small fragment deletions
        else:
            window_end = current_start + window_size
        j = i
        while j < len(dfs) and str(dfs.at[j, '#Target_name']) == current_chrom and dfs.at[j, 'Target_start'] < window_end:
            j += 1
        window_df = dfs[i:j].copy()
        if 2 <= len(window_df['Query_name'].unique()) < depth * 100:
            window_df.index = range(len(window_df))
            windows[current_start] = window_df
        i = j
    return windows

def windows_slide4tra(df, shift=2000000):
    """
    Slide windows over translocation data and collect valid windows.
    """
    from itertools import product
    df = df.copy()
    df.columns = ["#Target_name1", "Query_name", "Target_start1", "Target_start2", "SVlen", "maq", "SVID", 'SVType', 'seq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(np.int32)
    df["Target_start2"] = df["Target_start2"].astype(np.int32)
    df['maq'] = df['maq'].astype(np.int16)
    chrom1s = df['#Target_name1'].unique()
    chrom2s = df['#Target_name2'].unique()
    chrom_pairs = list(product(chrom1s, chrom2s))
    print(f"************translocation pairs chromsome***********{chrom_pairs}")
    windows = {}
    for chrom1, chrom2 in chrom_pairs:
        dfs = df[(df['#Target_name1'] == chrom1) & (df['#Target_name2'] == chrom2)]
        dfs = dfs.sort_values(by=['Target_start1', 'Target_start2'])
        dfs.index = range(len(dfs))
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
            window_df = window_df.reset_index(drop=True)
            key = f'{chrom1}_{current_start}_{chrom2}'
            windows[key] = window_df
            i = j
    return windows

def klook_clu_tra(win):
    """
    tra clus, the clusdf should be chromosome pair mode,
    each clusdf only allow two chromosome.
    """
    if len(win) == 1:
        win['shift_cluster'] = -1
        return win
    cluster_id = 0
    win.loc[0, 'shift_cluster'] = 0
    for i in range(1, len(win)):
        current_start1 = win.loc[i, 'Target_start1']
        current_start2 = win.loc[i, 'Target_start2']
        found_cluster = False
        start_index = max(0, i - 40)
        for j in range(i - 1, start_index - 1, -1):
            max_diff = 1000
            ## already ensure chrom1 != chrom2
            pos_condition = (abs(current_start1 - win.loc[j, 'Target_start1']) <= max_diff and
                             abs(current_start2 - win.loc[j, 'Target_start2']) <= max_diff)
            if  pos_condition:
                win.loc[i, 'shift_cluster'] = win.loc[j, 'shift_cluster']
                found_cluster = True
                break
        if not found_cluster:
            cluster_id += 1
            win.loc[i, 'shift_cluster'] = cluster_id
    return win

def candidate_tra(window, opened_bam, dtype):
    """
    Process translocation data and find candidate translocations.
    """
    win_clus = klook_clu_tra(window)
    tras = []
    for clu in win_clus['shift_cluster'].unique():
        win = win_clus[win_clus['shift_cluster']==clu]
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
        local_depth = np.mean([start_local_map, end_local_map])
        if local_depth >0:
            SV_rate = round(sv_eye / local_depth, 2)
        else:
            SV_rate = 1
        if dtype in ['ont','hifi', 'pb']:
            ## to one depth ##
            if sv_eye >= (local_depth * 0.1 + 0.7 ):
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
        elif dtype in ['cr', 'sr']:
            if sv_eye > local_depth * 0.25:
                tras.append([chr1, chr1_start, chr2_start, sv_len, svid, "TRA", "*", sv_eye, SV_rate, maq, readsname])
    return tras

def load_and_process_sv_data(args):
    from math import floor
    """
    Load and preprocess structural variation data.
    """
    try:
        sv_indel_data = pd.read_csv(args.raw_signal, sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        print(f"Error: File {args.raw_signal} not found.")
        return {}, [], None
    except pd.errors.EmptyDataError:
        print(f"Warning: File {args.raw_signal} is empty.")
        sv_indel_data = pd.DataFrame()
    try:
        depth_stat = pd.read_csv(f'{args.raw_signal}.depth', sep="\t", header=None, dtype=str, index_col=None)
    except FileNotFoundError:
        depth = None
    else:
        depth = None if depth_stat.empty else ceil(float(depth_stat.iloc[0, 3])+0.3)
        if args.nreads:
            nreads_fil = args.nreads
        else:
            nreads_fil = floor(depth / 10)
    print(f'**************** average depth is {depth} ********************')
    if sv_indel_data.empty:
        return {}, [], depth, nreads_fil
    sv_indel_data.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                             "seq"]
    chroms = sv_indel_data['#Target_name'].unique()
    sv_indel_data = sv_indel_data.copy()
    sv_indel_data['SVlen'] = sv_indel_data['SVlen'].astype(np.int32)
    sv_indel_data['maq'] = sv_indel_data['maq'].astype(np.int16)
    sv_indel_data = sv_indel_data[sv_indel_data["SVlen"] <= args.max]
    sv_indel_data['Target_start'] = sv_indel_data['Target_start'].astype(np.int32)
    sv_indel_data['Target_end'] = sv_indel_data['Target_end'].astype(np.int32)
    sv_data = {sv_type: sv_indel_data[sv_indel_data['SVType'] == sv_type] for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]}
    supp_align_file = f"{args.raw_signal}.suppAlign"
    if os.path.exists(supp_align_file):
        try:
            msv = pd.read_csv(supp_align_file, sep="\t", header=None, dtype=str, index_col=None)
        except Exception as e:
            print(f"***************** empty {supp_align_file} file ******************")
        else:
            if not msv.empty:
                print(f'*************** {args.raw_signal}.suppAlign has {msv.shape[0]} rows ********************')
                msv.columns = ["#Target_name", "Query_name", "Target_start", "Target_end", "SVlen", "maq", "SVID", "SVType",
                               "seq"]
                #msv = msv.drop_duplicates()
                print(f'*************** {args.raw_signal}.suppAlign after duplicates drop has {msv.shape[0]} rows ********************')
                #msv = msv.reset_index(drop=True)
                for svtype in sv_data:
                    sv_data[svtype] = pd.concat([sv_data[svtype], msv[msv['SVType'] == svtype]], axis=0)
                for svtype in sv_data:
                    sv_data[svtype]['SVlen'] = sv_data[svtype]['SVlen'].astype(np.int32)
                    sv_data[svtype]['Target_start'] = sv_data[svtype]['Target_start'].astype(np.int32)
                    sv_data[svtype]['Target_end'] = sv_data[svtype]['Target_end'].astype(np.int32)
                    sv_data[svtype]['maq'] = sv_data[svtype]['maq'].astype(np.int16)
            else:
                print(f"Warning file {supp_align_file} is empty")
    else:
        print(f"Warning file {supp_align_file} not exist")
    print(f'******************** all chromosomes list {chroms} ****************************')
    return sv_data, chroms, depth, nreads_fil

def process_svtype(args, sv_data, chroms, svtype, depth, nreads, minLen):
    """
        cov info parse by covfile or bam file
    """
    print(f"start klook for {args.raw_signal}  SV type: {svtype}") 
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
    if args.dtype in ['pb', 'ont', 'hifi']:
        print(f'data type is {args.dtype}')
        try:
            with pysam.AlignmentFile(bam_path, "rb") as opened_bam:
                for chrom in chroms:
                    chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                    def process_sv_type(svtype):
                        if svtype == "TRA":
                            log = open("log_tra", 'w')
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            windows = windows_slide4tra(chrom_data[svtype], 2000000)
                            print(f'The TRA windows by 2M: \n{windows}', file=log)
                            tra_list = []
                            for win in windows.values():
                                win_clus = klook_clu_tra(win)
                                print(win_clus.iloc[:,[0,2,3,4,5,-1]], file=log)
                                tra = candidate_tra(win_clus, opened_bam, args.dtype)
                                print(tra, file=log)

                                if tra:
                                    tra_list.extend(tra)
                            log.close()
                            return tra_list
                        else:
                            if chrom_data[svtype].empty:
                                print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                                return []
                            sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                            sv_windows = windows_slide(sv_dfs, depth, svtype, nreads,args.window_size)
                            candidate_svs = []
                            if args.dtype == 'hifi':
                                hdepth = 5
                            else:
                                hdepth = 10
                            for sv_window in sv_windows.values():
                                if len(sv_window['Query_name']) > hdepth: 
                                    candisv = highdepth_clu(sv_window,args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
                                else:
                                    candisv = lowdepth_clu(sv_window, args.num_hap, svtype, opened_bam, nreads, args.rate_depth, 1)
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
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    elif args.dtype in ['sr', 'cr']:
        print(f'data type is {args.dtype}')
        try:
            for chrom in chroms:
                chrom_data = {svtype: sv_data[svtype][sv_data[svtype]['#Target_name'] == chrom] for svtype in sv_data}
                def process_sv_type(svtype):
                    if svtype != "TRA":
                        sv_dfs = chrom_data[svtype][chrom_data[svtype]['SVlen'] >= minLen]
                        if sv_dfs.empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        sv_dfs = sv_dfs.sort_values(by=['Target_start', 'SVlen'])
                        sv_dfs.index = range(len(sv_dfs))
                        print(sv_dfs.head(10))
                        print("**************************** Calling one depth all clustering ***********************")
                        candidate_svs = []
                        win_dfs = windows_slide4asm(sv_dfs, svtype, args.window_size)
                        for win_df in win_dfs.values():
                            print(f'**************************window signals*********************\n{win_df}')
                            candidate_sv = onedepth_all_clus(win_df, svtype, bam_path,max(args.rate_depth,0.25))
                            candidate_svs.extend(candidate_sv)
                        return candidate_svs
                    else:
                        if chrom_data[svtype].empty:
                            print(f'******************************** no {svtype} signal found at {chrom} *******************************')
                            return []
                        windows = windows_slide4tra(chrom_data[svtype], 1000000)
                        tra_list = []
                        for value in windows.values():
                            win = value
                            tra = candidate_tra(win, bam_path,args.dtype)
                            if tra:
                                tra_list.extend(tra)
                        print(f'****************************** the tra condidate ***************************')
                        print(tra_list)
                        return tra_list
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
            print(f'Error in processing {args.dtype}: {args.raw_signal}')
            print(f"Error processing {svtype} for chromosome {chrom}: {e}")
            return [],[],[],[],[]
        return tra_clus, del_clus, ins_clus, inv_clus, dup_clus
    else:
        print(f"Error: sequence data dtype error: {args.dtype} is not in [pb,hifi,ont,sr,cr], please check parameter -dtype")
        return [],[],[],[],[]

def candidateSV(args):
    sv_data, chroms, depth, nreads = load_and_process_sv_data(args)
    print(chroms, depth, nreads)
    if sv_data:
        sv_types = ["DEL", "INS", "INV", "DUP", "TRA"]
        with multiprocessing.Pool() as pool:
            results = pool.starmap(process_svtype, [(args, sv_data, chroms, svtype, depth, nreads, args.min) for svtype in sv_types])
        tra_clus, del_clus, ins_clus, inv_clus, dup_clus = [], [], [], [], []
        for result in results:
            if result is not None:
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
    IN.add_argument("-M", dest="max", type=int, default=10000000, help="the max SV length")
    IN.add_argument("-m", dest="min", type=int, default=45, help="the minimum SV length")
    IN.add_argument("-dtype", dest="dtype", type=str, required=True, help="the sequencing type of samples")
    IN.add_argument("--cov", dest="covfile", type=str, help="Coverage File")
    IN.add_argument("--b", dest="bam", type=str, help="the bam file of Individual")
    IN.add_argument("--nreads", dest="nreads", type=int, help="the minimum numbers of reads to support SV, if not provided, we use average_depth / 10 as threshold" )
    IN.add_argument("--rate_depth", dest="rate_depth", type=float, default=0.1, help="the sv supports of local depth ratio to support sv, 0.1 means the percent of local reads shoule support sv")
    IN.add_argument("--window", dest="window_size", type=int, default=500, help="the window size of signal to parse in klook cluster, 500bp suggested")
    IN.add_argument("--num_hap", dest="num_hap", type=int, default=2, help="numbers of haplotypes within local region should be defined by species ploid, 2 for diploid, 4 for Tetraploid")
    args = parser.parse_args()
    start_t = time()
    tra_clus, del_clus, ins_clus, inv_clus, dup_clus = candidateSV(args)
    if tra_clus:
        pd.DataFrame(tra_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_TRA.signal", header=True, sep="\t", index=None)
    if dup_clus:
        pd.DataFrame(dup_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DUP.signal", header=True, sep="\t", index=None)
    if inv_clus:
        pd.DataFrame(inv_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INV.signal", header=True, sep="\t", index=None)
    if ins_clus:
        pd.DataFrame(ins_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_INS.signal", header=True, sep="\t", index=None)
    if del_clus:
        pd.DataFrame(del_clus, columns= ["#Target_name", "Target_start", "Target_end", "SVlen", "SVID", "SVType", "seq", "cluster_size",'sv_rate','maq', 'readsID']).to_csv(f"{args.raw_signal}_DEL.signal", header=True, sep="\t", index=None)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")
