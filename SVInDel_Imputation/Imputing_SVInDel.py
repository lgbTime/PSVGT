import pandas as pd
import numpy as np
from copy import deepcopy
import argparse
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
def calculate_missing_rate(row):
    # Count the number of missing genotypes
    missing_count = sum(genotype in {'.', './.'} for genotype in row)
    total_count = len(row)
    return missing_count / total_count if total_count > 0 else 0

def gtdf_to_haplotype_array_one_hot(df):
    haplotype_dict = {}
    sample_ids = df.columns[1:]  # Assuming the first column contains something like IDs or a placeholder
    for index, row in df.iterrows():
        genotypes = row[1:]  # Skip the first column if it doesn't contain genotype data
        for sample_index, genotype in enumerate(genotypes):
            sample_id = sample_ids[sample_index]
            if sample_id not in haplotype_dict:
                haplotype_dict[sample_id] = []
            if genotype == '0/0':
                haplotype_dict[sample_id].append([1, 0, 0])  # Homozygous reference
            elif genotype == '1/1':
                haplotype_dict[sample_id].append([0, 1, 0])  # Homozygous alternate
            elif genotype in ['0/1', '1/0']:
                haplotype_dict[sample_id].append([0, 0, 1])  # Heterozygous
            elif genotype in ['.', './.']:  # Handle missing genotypes
                haplotype_dict[sample_id].append([0, 0, 0])  # Missing genotype representation
            else:
                haplotype_dict[sample_id].append([0, 0, 0])  # Unknown genotype
    haplotypes_array = np.array([np.array(haplotype_dict[sample_id]).flatten() for sample_id in haplotype_dict])
    return haplotypes_array, sample_ids.tolist()

def blood_hood(gtdf):
    haplotypes_array, sample_ids = gtdf_to_haplotype_array_one_hot(gtdf)
    unique_haplotypes, unique_indices = np.unique(haplotypes_array, axis=0, return_inverse=True)
    n_clusters = len(unique_haplotypes)   ## this values matter a lots
    distance_matrix = pairwise_distances(unique_haplotypes, metric='hamming')
    # Clustering using Agglomerative Clustering
    agg_cluster = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
    clusters = agg_cluster.fit_predict(distance_matrix)
    # Create a DataFrame to map unique haplotypes to their cluster assignments
    cluster_mapping = pd.DataFrame({
        'Unique_Haplotype': [list(hap) for hap in unique_haplotypes],
        'Cluster': clusters
    })
    # Map each sample to its respective cluster based on the unique haplotypes
    sample_clusters = pd.DataFrame({
        'Sample_ID': sample_ids,
        'Cluster': clusters[unique_indices]
    })
    #print("Number of unique clusters:", n_clusters)
    #print("Cluster assignments for each sample:")
    #print(sample_clusters)
    return sample_clusters

def hap_cluster2impute(rowGT: pd.Series, hap_clusters: pd.DataFrame) -> pd.Series:
    """
    Process genotype data and aggregate by haplotype clusters, maintaining original rowGT format.
    Parameters:
    - rowGT: Series containing genotype data for a specific genomic region.
    - hap_clusters: DataFrame containing sample IDs and their corresponding clusters.
    Returns:
    - A Series with aggregated genotype information for each sample, structured like rowGT.
    """
    result_series = pd.Series(index=rowGT.index , dtype=str)
    for sample in rowGT.index:
        result_series[sample] = rowGT[sample]  # Start with the original value
    for cluster_id in hap_clusters['Cluster'].unique():
        cluster_samples = hap_clusters[hap_clusters['Cluster'] == cluster_id]['Sample_ID'].tolist()
        cluster_genotypes = rowGT[cluster_samples]
        non_missing_gt = cluster_genotypes[cluster_genotypes != './.'].unique()
        if len(non_missing_gt) > 0:
            # If there's a valid genotype, use it to fill in missing genotypes
            valid_gt = non_missing_gt[0]
            for sample in cluster_samples:
                if result_series[sample] == './.':
                    result_series[sample] = valid_gt
    return result_series

def imputing(args):
    gt = pd.read_csv(args.SVInDelGT,header=0,sep="\t", index_col=0)
    gt.sort_values(by=["#Target_name", "Target_start"], inplace=True)
    pos = gt.index
    samples = gt.columns[5:]
    gt.index = range(gt.shape[0]) ## avoid same SVInDel break record of Insertion
    out_imputed_gt = deepcopy(gt)
    for chrom in gt["#Target_name"].unique():
        chromgt = gt[gt["#Target_name"] == chrom]
        for i, index in enumerate(chromgt.index):
            rowGT = chromgt.loc[index, :]
            missing = calculate_missing_rate(rowGT)
            if missing <= args.miss and missing > 0 :
                lower_bound = max(0, i - args.Kneigbor)  # Ensure we don't go below 0
                upper_bound = min(len(chromgt) - 1, i + args.Kneigbor)  # Ensure we don't exceed the bounds
                Local_df = chromgt.iloc[lower_bound:upper_bound + 1]  # Include upper bound
                Local_gt = Local_df[samples]
                haplotype_clusters = blood_hood(Local_gt)
                target_gt = Local_gt[samples].loc[index, :]
                imputed_gt = hap_cluster2impute(target_gt, haplotype_clusters)
                out_imputed_gt.loc[index,samples] = imputed_gt
    out_imputed_gt.index = pos
    out_imputed_gt.to_csv(f"{args.SVInDelGT}_imputed_miss{args.miss}_Kneigbor{args.Kneigbor}.txt",sep="\t",header=True,index=True)
    return out_imputed_gt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Haplotypes Cluster to Impute the SVInDel Genotype file")
    parser.add_argument('SVInDelGT', type=str, help='the Genotype File Generated by SVInDel Program.')
    parser.add_argument('--Kneigbor', type=int, help='to Impute the Missing, we use the up down Neigbor SVs to Cluster Haplotype default is 5', default=5)
    parser.add_argument('--miss', type=float,default=0.8, help='missing rate to impute,  0.8 recommended, default is 0.8')
    args = parser.parse_args()
    imputing(args)
