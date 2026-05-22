import pandas as pd
import argparse
from time import time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
def svs_clu(chrsvdf, svtype,max_diff=30):
    df = chrsvdf.copy()
    df['Target_start'] = df['Target_start'].astype(int)
    df['Target_end'] = df['Target_end'].astype(int)
    df = df.sort_values(by=['#Target_name', 'Target_start','SVlen'])
    df.index = range(len(df))
    df['cluster'] = -1  # Initializing cluster column
    cluster_id = 0
    df.loc[0,'cluster'] = cluster_id
    for i in range(1,len(df)):
        if df.loc[i, '#Target_name'] == df.loc[i-1, '#Target_name']:
            old_len = df.loc[i-1, 'SVlen']
            now_len = df.loc[i, 'SVlen']
            relate_size = old_len / now_len
            if abs(df.loc[i, 'Target_start'] - df.loc[i-1, 'Target_start']) <= max_diff and abs(df.loc[i, 'Target_end'] - df.loc[i-1, 'Target_end']) <= max_diff and (0.8 < relate_size < 1.2):
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
        max_clu_row = cs[cs['cluster_size'] == cs['cluster_size'].max()].iloc[0]
        Target_name = max_clu_row['#Target_name']
        Target_start1 = max_clu_row['Target_start']
        Target_start2 = max_clu_row['Target_end']
        cluster_size = cs['cluster_size'].sum()
        SVlen =  cs['SVlen'].max()
        maq = int(cs['maq'].mean())
        SVID = max_clu_row['SVID']
        SVType = max_clu_row['SVType']
        seq = "*"
        sv_rate = cs['sv_rate'].sum()
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq':maq,
            'cluster_size': cluster_size,
            'sv_rate': sv_rate
         })
    return clu


def tra_clu(tradf,  max_diff=800):
    df = tradf.copy()
    df.columns = ["#Target_name1", "Target_start1","Target_start2", "SVlen","SVID",'SVType','seq','cluster_size','sv_rate', 'maq', 'readsID']
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
        max_clu_row = cs[cs['cluster_size'] == cs['cluster_size'].max()].iloc[0]
        Target_name = max_clu_row['#Target_name1']
        Target_start1 = max_clu_row['Target_start1']
        Target_start2 = max_clu_row['Target_start2']
        cluster_size = cs['cluster_size'].sum()
        maq = int(cs['maq'].mean())
        SVlen = 0
        SVID = max_clu_row['SVID']
        SVType = "TRA"
        seq = "*"
        sv_rate = cs['sv_rate'].sum()
        clu.append({
            '#Target_name': Target_name,
            'Target_start': Target_start1,
            'Target_end': Target_start2,
            'SVlen': SVlen,
            'SVID': SVID,
            'SVType': SVType,
            'seq': seq,
            'maq': maq,
            'cluster_size': cluster_size,
            'sv_rate': sv_rate
         })
    return clu


def drop_ins_from_dup(df):
    out=open("log", 'w')
    ins = df[df['SVType'].isin(['INS', 'DUP'])]
    ins = ins.sort_values(by=['#Target_name', 'Target_start'])
    ins.index = range(len(ins))
    dup_rows = ins[ins['SVType'] == 'DUP']
    if dup_rows.empty:
        return df

    dropINS = []
    for index, dup_row in dup_rows.iterrows():
        start_index = max(0, index - 10)
        end_index = min(len(df), index + 10)
        nearby_ins_rows = ins[(ins.index >= start_index) & (ins.index < end_index)]
        print(nearby_ins_rows, file=out)
        for ins_index, ins_row in nearby_ins_rows.iterrows():
            if (dup_row['Target_start'] < ins_row['Target_start'] - 50 < dup_row['Target_end'] and
                    dup_row['#Target_name'] == ins_row['#Target_name']):
                ratio = min(ins_row['SVlen'] / dup_row['SVlen'], dup_row['SVlen']/ins_row['SVlen'])
                if ratio >= 0.8:
                    dropINS.append(ins_row['SVID'])
    filtered_df = df[~df['SVID'].isin(dropINS)]
    print(f'******************************* drop {len(dropINS)} INS SV since they are DUP **************************************')
    print(f'******************************* drop {dropINS} INS SV since they are DUP **************************************')
    out.close()
    return filtered_df

def read_file(file_name):
    import os
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
    from os.path import basename, exists
    import pandas as pd
    chrs = pd.read_csv(args.fai,sep="\t",index_col=None,header=None)[0].tolist()
    file_lists = []
    for sv in ['INS', 'DEL', 'TRA', 'DUP','INV']:
        for chrom in chrs:
            file_name = f"{args.preffix}_{chrom}.record.txt_{sv}.signal"
            if exists(file_name):
                file_lists.append(file_name)
    print(file_lists)

    sv = read_file(file_lists[0])
    if args.nreads_nrate:
        sv = sv[(sv["cluster_size"] >= args.nreads) & (sv["sv_rate"] >= args.nrate)]
    else:
        if args.nrate and not args.nreads:
            sv = sv[sv["sv_rate"] >= args.nrate]
        elif args.nreads and not args.nrate:
            sv = sv[sv["cluster_size"] >= args.nreads]
        elif args.nreads and args.nrate:
            sv = sv[(sv["sv_rate"] >= args.nrate) | (sv["cluster_size"] >= args.nreads)]

    for file_name in file_lists[1:]:
        svi = read_file(file_name)
        sv = pd.concat([sv,svi],axis=0)
    if args.nreads_nrate:
        sv = sv[(sv["cluster_size"] >= args.nreads) & (sv["sv_rate"] >= args.nrate)]
    else:
        if args.nrate and not args.nreads:
            sv = sv[sv["sv_rate"] >= args.nrate]
        elif args.nreads and not args.nrate:
            sv = sv[sv["cluster_size"] >= args.nreads]
        elif args.nreads and args.nrate:
            sv = sv[(sv["sv_rate"] >= args.nrate) | (sv["cluster_size"] >= args.nreads)]

    ori = sv.shape
    sv_types = sv["SVType"].unique()
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
            dels = svs_clu(delchr,'DEL', args.shift)

        if inschr.empty:
            ins = []
        else:
            ins = svs_clu(inschr,'INS', args.shift)

        if invchr.empty:
            inv = []
        else:
            inv = svs_clu(invchr, 'INV', args.shift)

        if dupchr.empty:
            dup =[]
        else:
            dup = svs_clu(dupchr,'DUP', args.shift)

        if  trachr.empty:
            tra = []
        else:
            tra = tra_clu(trachr, args.shift*2)
        return tra + dels + ins + inv + dup
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_chromosome, sv['#Target_name'].unique())
    out = []
    for result in results:
        out += result
    svs_df = pd.DataFrame(out)
    sv_out = drop_ins_from_dup(svs_df)
    sv_out.loc[:,"SVlen"] = sv_out["SVlen"].astype(int)
    sv_out[(sv_out["SVlen"] <= args.max) & (sv_out['maq']>=args.minimaq)].to_csv(f"{args.preffix}_Clustered_Record.txt",header=True,index=None,sep="\t")
    now = sv_out.shape
    print(f'Original data shape: {ori}, after clustering: {now}\nOutput file is {args.preffix}_Clustered_Record.txt')


if __name__ == "__main__":
    parser = argparse.ArgumentParser("signal filtering through support reads ratio", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    IN = parser.add_argument_group("Input file ")
    IN.add_argument("-preffix", dest="preffix", required=True, help="the sample preffix name of SV raw sginal")
    IN.add_argument("-fai", dest="fai",required=True, help="the reference faidx file")
    IN.add_argument("-s", dest="shift", default=100, type=int, help="the distance of shifting the breakpoints ")
    IN.add_argument("-M", dest="max", default=6868886, type=int, help="the max SV length ")
    IN.add_argument("--nreads", dest="nreads", type=int, help="filter signal by minimum support reads")
    IN.add_argument("--nrate", dest="nrate",   type=float, help="filter signal by minimum ratio of signal")
    IN.add_argument("--nn", dest="nreads_nrate",action="store_true", help="strict filter for signal use both meet threshold of nreads and nrate, if not --nn, keep SV meet nreads or nrate")
    IN.add_argument("--minimaq", dest="minimaq", type=int, default=50,help="minimum mapping quality of SV")
    args = parser.parse_args()
    start_t = time()
    candidateSV(args)
    end_t = time()
    print(f"******************** Time in cluster Cost {end_t - start_t}s *****************************")
