import pandas as pd
import argparse
from time import time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
def mode_or_median(series):
    """
    Returns the most common value from a series. If there are ties, it will return the first mode.
    """
    return series.mode().iloc[0] if not series.mode().empty else math.floor(series.median())

def most_common(series):
    """
    Returns the most common value from a series. If there are ties, it will return the first mode.
    """
    return series.mode().iloc[0] if not series.mode().empty else series.iloc[0]

def svs_clu(chrsvdf, svtype, pop_size,max_diff=100):
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start'])
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1,len(df)):
        if df.loc[i, '#Target_name'] == df.loc[i-1, '#Target_name']:
            old_len = df.loc[i-1, 'SVlen']
            now_len = df.loc[i, 'SVlen']
            relate_size = old_len / now_len
            if abs(df.loc[i, 'Target_start'] - df.loc[i-1, 'Target_start']) <= max_diff and abs(df.loc[i, 'Target_end'] - df.loc[i-1, 'Target_end']) <= max_diff and (0.7 < relate_size < 1.3):
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id
    print(f'The {svtype} data shape: {df.shape}\nThe {svtype} data shape after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        Target_name = most_common(cs['#Target_name'])
        Target_start = most_common(cs['Target_start'])
        Target_end = most_common(cs['Target_end'])
        cluster_size = most_common(cs['cluster_size'])
        SVlen =  most_common(cs['SVlen'])
        SVID =  most_common(cs['SVID'])
        SVType =   most_common(cs['SVType'])
        seq = most_common(cs['seq'])
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start,
            'Target_end': Target_end,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'cluster_size': cluster_size,
                })
    return clu

def tra_clu(tradf, pop_size, max_diff=1000):
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','cluster_size', 'Query_name', 'maq']
    df["#Target_name2"] = df['SVID'].str.split(":", expand=True)[0]
    df["Target_start1"] = df["Target_start1"].astype(int)
    df["Target_start2"] = df["Target_start2"].astype(int)
    df.sort_values(by=['#Target_name1', '#Target_name2','Target_start1', 'Target_start2'],inplace=True)
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1, len(df)):
        if df.loc[i, '#Target_name1'] == df.loc[i-1, '#Target_name1'] and df.loc[i, '#Target_name2'] == df.loc[i-1, '#Target_name2']:
            if abs(df.loc[i, 'Target_start1'] - df.loc[i-1, 'Target_start1']) <= max_diff and abs(df.loc[i, 'Target_start2'] - df.loc[i-1, 'Target_start2']) <= max_diff:
                df.loc[i, 'cluster'] = df.loc[i-1, 'cluster']
            else:
                cluster_id += 1
                df.loc[i, 'cluster'] = cluster_id
        else:
            cluster_id += 1
            df.loc[i, 'cluster'] = cluster_id

    print(f'The TRA signal shape: {df.shape}\nThe TRA signal data shape after clustering by shift:{int(max_diff)} is {len(df["cluster"].unique())}') 
    clu = []
    for c in df['cluster'].unique():
        cs = df[df['cluster'] == c]
        Target_name = most_common(cs['#Target_name1'])
        Target_start1 = most_common(cs['Target_start1'])
        Target_start2 = most_common(cs['Target_start2'])
        cluster_size = most_common(cs['cluster_size'])
        SVlen =  0
        SVID =  most_common(cs['SVID'])
        SVType =  "TRA"
        seq = "*"
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'cluster_size': cluster_size,
            })
    return clu

def file_capture(dir, suffix):
    import os
    captures = []
    all_files = os.listdir(dir)
    for file in all_files:
        if file[-len(suffix):] == suffix:
            captures.append(os.path.join(dir, file))
    return captures
def read_file(file_name):
    import os
    # Check if the file exists and is not empty
    if os.path.exists(file_name) and os.path.getsize(file_name) > 0:
        try:
            svi = pd.read_csv(file_name, header=0, index_col=None, sep="\t")
            return svi
        except pd.errors.EmptyDataError:
            print(f"The file {file_name} is empty or malformed.")
            return None
    else:
        print(f"The file {file_name} does not exist or is empty.")
        return None

def candidateSV(args):
    file_lists = file_capture(args.sv_dir, ".signal")
    file_lists.sort()
    pop_size = len(file_capture(args.sv_dir, ".record.txt"))
    print(f'************************** The population size is {pop_size} ***************************')
    sv = pd.read_csv(file_lists[0],header=0,index_col=None,sep="\t")
    for file_name in file_lists[1:]:
        svi = read_file(file_name)
        #svi = pd.read_csv(file_name,header=0,index_col=None,sep="\t")
        sv = pd.concat([sv,svi],axis=0)
    ori = sv.shape
    sv_types = sv["SVType"].unique()
    sv.sort_values(by=["#Target_name","Target_start","SVlen","SVID"], inplace=True)
    def process_chromosome(chrom):
        inschr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INS")]
        delchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DEL")]
        trachr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="TRA")]
        invchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="INV")]
        dupchr = sv[(sv['#Target_name'] == chrom)&(sv['SVType']=="DUP")]
        if  delchr.empty:
            dels = []
        else:
            dels = svs_clu(delchr,'DEL', pop_size, args.shift)

        if inschr.empty:
            ins = []
        else:
            ins = svs_clu(inschr,'INS', pop_size, args.shift)

        if invchr.empty:
            inv = []
        else:
            inv = svs_clu(invchr, 'INV', pop_size, args.shift)

        if dupchr.empty:
            dup =[]
        else:
            dup = svs_clu(dupchr,'DUP', pop_size, args.shift)

        if  trachr.empty:
            tra = []
        else:
            tra = tra_clu(trachr, pop_size, args.shift*2)
        return tra + dels + ins + inv + dup
    # Use ThreadPoolExecutor to process chromosomes in parallel
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_chromosome, sv['#Target_name'].unique())
    # Collect results from all chromosomes
    out = []
    for result in results:
        out += result
    sv_out = pd.DataFrame(out)
    print(sv_out.head())
    sv_out["SVlen"] = sv_out["SVlen"].astype(int)
    sv_out = sv_out.sort_values(by=['SVID','cluster_size'],inplace=False)
    sv_fil = sv_out.drop_duplicates(subset='SVID', keep='last', inplace=False) ## keep big cluster_size
    sv_fil[sv_fil["SVlen"] <= args.max].to_csv(f"{args.sv_dir}/PopSV_clustered_Record.txt",header=True,index=None,sep="\t")
    now = sv_fil.shape
    print(f'Original data shape: {ori}, after clustering: {now}\nOutput file is {args.sv_dir}/PopSV_clustered_Record.txt')


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-d", dest="sv_dir", required=True, help="the PSVGT output directory")
    IN.add_argument("-s", dest="shift", default=50, type=int, help="the distance of shifting the breakpoints ")
    IN.add_argument("-M", dest="max", default=6868886, type=int, help="the max SV length ")
    args = parser.parse_args()
    start_t = time()
    candidateSV(args)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")

