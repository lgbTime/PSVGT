import pandas as pd
import argparse
import os
import math
from time import time
from multiprocessing import Pool

# --- Helper Functions (Must be top-level for Multiprocessing) ---

def most_common(series):
    """
    Returns the most common value from a series. 
    If there are ties, it will return the first mode.
    """
    return series.mode().iloc[0] if not series.mode().empty else series.iloc[0]

def svs_clu(chrsvdf, svtype, chrom, max_diff=50):
    """Identical to original logic to ensure output consistency."""
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start'])
    df.index = range(len(df))
    df['cluster'] = -1
    cluster_id = 0
    df.loc[0, 'cluster'] = cluster_id
    
    for i in range(1, len(df)):
        if df.loc[i, '#Target_name'] == df.loc[i-1, '#Target_name']:
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
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id
            
    print(f'The {chrom} {svtype} data {df.shape} after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        clu.append({
            '#Target_name': most_common(cs['#Target_name']),
            'Target_start': most_common(cs['Target_start']),
            'Target_end': most_common(cs['Target_end']),
            'SVlen': most_common(cs['SVlen']),
            'SVID': most_common(cs['SVID']),
            'SVType': most_common(cs['SVType']),
            'seq': most_common(cs['seq']),
            'maq': int(cs['maq'].mean()),
            'cluster_size_prevalent': most_common(cs['cluster_size']),
            'sv_rate_prevalent': most_common(cs['sv_rate'])
        })
    return clu

def tra_clu(tradf, chrom, max_diff=100):
    """Identical to original logic to ensure output consistency."""
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','maq','cluster_size','sv_rate']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(int)
    df["Target_start2"] = df["Target_start2"].astype(int)
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

    print(f'The {chrom} TRA signal data {tradf.shape} after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        clu.append({
            '#Target_name': most_common(cs['#Target_name1']),
            'Target_start': most_common(cs['Target_start1']),
            'Target_end': most_common(cs['Target_start2']),
            'SVlen': 0,
            'SVID': most_common(cs['SVID']),
            'SVType': "TRA",
            'seq': "*",
            'maq': int(cs['maq'].mean()),
            'cluster_size_prevalent': most_common(cs['cluster_size']),
            'sv_rate_prevalent': most_common(cs['sv_rate'])
        })
    return clu

# --- Multiprocessing Worker ---

def process_chromosome_worker(task_tuple):
    """
    Worker function to process all SV types for a single chromosome.
    Unpacks (chrom_name, chrom_dataframe, shift_value).
    """
    chrom, chrom_df, shift = task_tuple
    
    # Filter by type within the chromosome subset
    inschr = chrom_df[chrom_df['SVType'] == "INS"]
    delchr = chrom_df[chrom_df['SVType'] == "DEL"]
    trachr = chrom_df[chrom_df['SVType'] == "TRA"]
    invchr = chrom_df[chrom_df['SVType'] == "INV"]
    dupchr = chrom_df[chrom_df['SVType'] == "DUP"]

    results = []
    if not delchr.empty: results += svs_clu(delchr, 'DEL', chrom, shift)
    if not inschr.empty: results += svs_clu(inschr, 'INS', chrom, shift)
    if not invchr.empty: results += svs_clu(invchr, 'INV', chrom, shift)
    if not dupchr.empty: results += svs_clu(dupchr, 'DUP', chrom, shift)
    if not trachr.empty: results += tra_clu(trachr, chrom, shift * 2)
    
    return results

# --- Main Script Functions ---

def file_capture(dir_path, suffix):
    return sorted([os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(suffix)])

def read_file(file_name):
    if os.path.exists(file_name) and os.path.getsize(file_name) > 0:
        try:
            return pd.read_csv(file_name, header=0, index_col=None, sep="\t")
        except Exception as e:
            print(f"Error reading {file_name}: {e}")
            return None
    return None

def candidateSV(args):
    file_lists = file_capture(args.sv_dir, "_Clustered_Record.txt")
    if not file_lists:
        print("No files found.")
        return

    print(f'Loading {len(file_lists)} files...')
    sv = pd.concat([read_file(f) for f in file_lists if read_file(f) is not None], axis=0)
    ori = sv.shape
    
    sv.sort_values(by=["#Target_name", "Target_start", "SVlen", "SVID"], inplace=True)
    
    # Prepare tasks: Slice the dataframe by chromosome to minimize cross-process data overhead
    unique_chroms = sv['#Target_name'].unique()
    tasks = []
    for chrom in unique_chroms:
        chrom_subset = sv[sv['#Target_name'] == chrom].copy()
        tasks.append((chrom, chrom_subset, args.shift))

    # EXECUTE MULTIPROCESSING
    print(f"Starting Multiprocessing Pool with {os.cpu_count()} cores...")
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(process_chromosome_worker, tasks)

    # Flatten results
    out = [item for sublist in results for item in sublist]
    
    sv_out = pd.DataFrame(out)
    if sv_out.empty:
        print("No SVs were clustered.")
        return

    sv_out["SVlen"] = sv_out["SVlen"].astype(int)
    sv_out.sort_values(by=['SVID'], inplace=True)
    sv_fil = sv_out.drop_duplicates(subset='SVID', keep='last') 
    
    final_output = f"{args.sv_dir}/PopSV_Candidate_Record.txt"
    sv_fil[sv_fil["SVlen"] <= args.max].to_csv(final_output, header=True, index=None, sep="\t")
    
    print(f'Original shape: {ori}, after clustering: {sv_fil.shape}')
    print(f'Output saved to: {final_output}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", dest="sv_dir", required=True, help="the PSVGT output directory")
    parser.add_argument("-s", dest="shift", default=30, type=int, help="the distance of shifting the breakpoints")
    parser.add_argument("-M", dest="max", default=6868886, type=int, help="the max SV length")
    
    args = parser.parse_args()
    
    start_t = time()
    candidateSV(args)
    print(f"Total processing time: {time() - start_t:.2f} seconds")
