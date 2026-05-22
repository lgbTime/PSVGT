import pandas as pd
import argparse
import os
import math
from time import time
from multiprocessing import Pool
def svs_clu(chrsvdf, svtype, max_diff=30):
    if chrsvdf.empty: return []
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start', 'SVlen'])
    df.index = range(len(df))
    df['cluster'] = -1
    cluster_id = 0
    df.loc[0, 'cluster'] = cluster_id
    
    # Original logic preserved for output identity
    for i in range(1, len(df)):
        old_len = df.loc[i-1, 'SVlen']
        now_len = df.loc[i, 'SVlen']
        relate_size = old_len / now_len if now_len != 0 else 0
        if abs(df.loc[i, 'Target_start'] - df.loc[i-1, 'Target_start']) <= max_diff and \
           abs(df.loc[i, 'Target_end'] - df.loc[i-1, 'Target_end']) <= max_diff and \
           (0.8 < relate_size < 1.2):
            df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id
    
    print(f'Cluster: {svtype} shape: {df.shape} | Clusters: {len(df["cluster"].unique())}')
    
    clu = []
    # Vectorized aggregation would be faster, but keeping this for logic safety
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        max_clu_row = cs.loc[cs['cluster_size'].idxmax()]
        clu.append({
            '#Target_name': max_clu_row['#Target_name'],
            'Target_start': max_clu_row['Target_start'],
            'Target_end': max_clu_row['Target_end'],
            'SVlen': cs['SVlen'].max(),
            'SVID': max_clu_row['SVID'],
            'SVType': max_clu_row['SVType'],
            'seq': "*",
            'maq': int(cs['maq'].mean()),
            'cluster_size': cs['cluster_size'].sum(),
            'sv_rate': cs['sv_rate'].sum()
        })
    return clu

def tra_clu(tradf, max_diff=800):
    if tradf.empty: return []
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','cluster_size','sv_rate', 'maq', 'readsID']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df.sort_values(by=['#Target_name1', '#Target_name2','Target_start1', 'Target_start2'], inplace=True)
    df.index = range(len(df))
    df['cluster'] = -1
    cluster_id = 0
    df.loc[0, 'cluster'] = cluster_id
    
    for i in range(1, len(df)):
        if df.loc[i, '#Target_name1'] == df.loc[i-1, '#Target_name1'] and \
           df.loc[i, '#Target_name2'] == df.loc[i-1, '#Target_name2']:
            if abs(df.loc[i, 'Target_start1'] - df.loc[i-1, 'Target_start1']) <= max_diff and \
               abs(df.loc[i, 'Target_start2'] - df.loc[i-1, 'Target_start2']) <= max_diff:
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id

    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        max_clu_row = cs.loc[cs['cluster_size'].idxmax()]
        clu.append({
            '#Target_name': max_clu_row['#Target_name1'],
            'Target_start': max_clu_row['Target_start1'],
            'Target_end': max_clu_row['Target_start2'],
            'SVlen': 0,
            'SVID': max_clu_row['SVID'],
            'SVType': "TRA",
            'seq': "*",
            'maq': int(cs['maq'].mean()),
            'cluster_size': cs['cluster_size'].sum(),
            'sv_rate': cs['sv_rate'].sum()
        })
    return clu

def drop_ins_from_dup(df):
    """Refined to avoid unnecessary nested loops."""
    ins_mask = df['SVType'] == 'INS'
    dup_mask = df['SVType'] == 'DUP'
    
    if not dup_mask.any() or not ins_mask.any():
        return df

    ins_df = df[ins_mask].sort_values(['#Target_name', 'Target_start'])
    dup_df = df[dup_mask]
    
    drop_svids = set()
    for _, dup in dup_df.iterrows():
        # Only look at INS on the same chromosome within a range
        nearby = ins_df[(ins_df['#Target_name'] == dup['#Target_name']) & 
                        (ins_df['Target_start'] > dup['Target_start'] - 500) & 
                        (ins_df['Target_start'] < dup['Target_end'] + 500)]
        
        for _, ins in nearby.iterrows():
            if (dup['Target_start'] < ins['Target_start'] - 50 < dup['Target_end']):
                ratio = min(ins['SVlen']/dup['SVlen'], dup['SVlen']/ins['SVlen']) if dup['SVlen'] > 0 else 0
                if ratio >= 0.8:
                    drop_svids.add(ins['SVID'])
                    
    print(f'*** Dropped {len(drop_svids)} INS that were likely DUPs ***')
    return df[~df['SVID'].isin(drop_svids)]

# --- 2. Worker Function ---

def process_chromosome_worker(task_bundle):
    """Handles all SV types for one chromosome subset."""
    chrom, chrom_df, shift = task_bundle
    
    ins_df = chrom_df[chrom_df['SVType'] == "INS"]
    del_df = chrom_df[chrom_df['SVType'] == "DEL"]
    tra_df = chrom_df[chrom_df['SVType'] == "TRA"]
    inv_df = chrom_df[chrom_df['SVType'] == "INV"]
    dup_df = chrom_df[chrom_df['SVType'] == "DUP"]

    results = []
    results += svs_clu(del_df, 'DEL', shift)
    results += svs_clu(ins_df, 'INS', shift)
    results += svs_clu(inv_df, 'INV', shift)
    results += svs_clu(dup_df, 'DUP', shift)
    results += tra_clu(tra_df, shift * 2)
    return results

# --- 3. Main Logic ---

def candidateSV(args):
    # Load FAI to get chromosome list
    chrs = pd.read_csv(args.fai, sep="\t", header=None)[0].tolist()
    
    # Collect signal files
    file_lists = []
    for sv_type in ['INS', 'DEL', 'TRA', 'DUP', 'INV']:
        for chrom in chrs:
            fname = f"{args.preffix}_{chrom}.record.txt_{sv_type}.signal"
            if os.path.exists(fname) and os.path.getsize(fname) > 0:
                file_lists.append(fname)
    
    if not file_lists:
        print("No signal files found!")
        return

    # Load and initial filter
    df_list = []
    for f in file_lists:
        tmp = pd.read_csv(f, sep="\t")
        if args.nreads_nrate:
            tmp = tmp[(tmp["cluster_size"] >= args.nreads) & (tmp["sv_rate"] >= args.nrate)]
        else:
            if args.nrate and not args.nreads:
                tmp = tmp[tmp["sv_rate"] >= args.nrate]
            elif args.nreads and not args.nrate:
                tmp = tmp[tmp["cluster_size"] >= args.nreads]
            elif args.nreads and args.nrate:
                tmp = tmp[(tmp["sv_rate"] >= args.nrate) | (tmp["cluster_size"] >= args.nreads)]
        df_list.append(tmp)
        
    sv = pd.concat(df_list, axis=0)
    sv.sort_values(by=["#Target_name", "Target_start", "SVlen", "SVID"], inplace=True)
    
    # Create Multiprocessing Tasks
    unique_chroms = sv['#Target_name'].unique()
    tasks = [(c, sv[sv['#Target_name'] == c].copy(), args.shift) for c in unique_chroms]

    print(f"Parallelizing over {len(unique_chroms)} chromosomes using {os.cpu_count()} cores...")
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(process_chromosome_worker, tasks)

    # Flatten and post-process
    flat_results = [item for sublist in results for item in sublist]
    svs_df = pd.DataFrame(flat_results)
    
    if svs_df.empty:
        print("No SVs remaining after clustering.")
        return

    sv_out = drop_ins_from_dup(svs_df)
    sv_out["SVlen"] = sv_out["SVlen"].astype(int)
    
    # Final filter and save
    final_mask = (sv_out["SVlen"] <= args.max) & (sv_out['maq'] >= args.minimaq)
    output_name = f"{args.preffix}_Clustered_Record.txt"
    sv_out[final_mask].to_csv(output_name, header=True, index=None, sep="\t")
    
    print(f'Done. Output: {output_name}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-preffix", required=True, help="sample preffix")
    parser.add_argument("-fai", required=True, help="reference fai")
    parser.add_argument("-s", dest="shift", default=100, type=int, help="shift distance")
    parser.add_argument("-M", dest="max", default=6868886, type=int, help="max SV length")
    parser.add_argument("--nreads", type=int, help="min support reads")
    parser.add_argument("--nrate", type=float, help="min ratio")
    parser.add_argument("--nn", dest="nreads_nrate", action="store_true", help="strict filter (AND)")
    parser.add_argument("--minimaq", type=int, default=50, help="min mapping quality")
    args = parser.parse_args()
    t0 = time()
    candidateSV(args)
    print(f"Total time: {time() - t0:.2f}s")
