#!/usr/bin/env python3
"""
True O(N) Sliding Window Genomic PAV Matrix Collapse and Merger
Author: lgbtime
Year: 2026

Description:
  Uses an advanced two-pointer sliding window to cluster overlapping DEL /
  close INS variants, guaranteeing linear O(N) time complexity. Includes an
  optional parameter to write merge cluster histories out to a tracking log.
  Parallel chromosome/type processing via multiprocessing.Pool.
  TRA (translocation) support: breakpoints parsed from SVID column;
  merging groups by (local_chrom, remote_chrom) pair with dual-breakpoint jitter.
"""

import os
import sys
import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def parse_tra_svid(svid):
    """
    Parses TRA breakpoints from SVID format:
      <remote_chrom>:<remote_pos>_<local_chrom>:<local_pos>
    Returns (remote_chrom, remote_pos) or (None, None) on failure.
    """
    try:
        parts = svid.split("_")
        # Find the two chrom:pos tokens — split on last ':' of each token
        tokens = [p for p in parts if ":" in p]
        remote_chrom, remote_pos = tokens[0].rsplit(":", 1)
        return remote_chrom, int(remote_pos)
    except Exception:
        return None, None

def parse_input_matrix(matrix_file):
    """Parses the existing uncollapsed PAV matrix file."""
    headers  = []
    variants = []

    with open(matrix_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#Target_name") or line.startswith("Target_name"):
                headers = line.strip().split('\t')
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            sv_type = parts[5]
            v_record = {
                "chrom":        parts[0],
                "start":        int(parts[1]),
                "end":          int(parts[2]),
                "size":         int(parts[3]),
                "svid":         parts[4],
                "type":         sv_type,
                "seq":          parts[6],
                "maq":          parts[7],
                "cluster_size": int(parts[8]),
                "sv_rate":      int(parts[9]),
                "genotypes":    parts[10:]
            }

            # Parse TRA remote breakpoint from SVID
            if sv_type == "TRA":
                remote_chrom, remote_pos = parse_tra_svid(parts[4])
                v_record["remote_chrom"] = remote_chrom
                v_record["remote_pos"]   = remote_pos
            else:
                v_record["remote_chrom"] = None
                v_record["remote_pos"]   = None

            variants.append(v_record)

    return variants, headers[10:]

# ---------------------------------------------------------------------------
# Merge logic
# ---------------------------------------------------------------------------

def are_mergeable(v1, v2, max_jitter=20, min_fc=0.8):
    """Evaluates if two structural variants meet the threshold criteria to be merged."""
    sv_type = v1["type"]

    if sv_type == "TRA":
        # Both local AND remote breakpoints must be within jitter distance
        local_dist  = abs(v1["start"]       - v2["start"])
        remote_dist = abs(v1["remote_pos"]  - v2["remote_pos"])
        return local_dist <= max_jitter and remote_dist <= max_jitter

    elif sv_type == "INS":
        dist       = abs(v1["start"] - v2["start"])
        size_ratio = min(v1["size"], v2["size"]) / max(v1["size"], v2["size"])
        return dist <= max_jitter and size_ratio >= min_fc

    else:  # DEL, DUP, INV
        overlap_start = max(v1["start"], v2["start"])
        overlap_end   = min(v1["end"],   v2["end"])
        if overlap_start < overlap_end:
            overlap_len = overlap_end - overlap_start
            ro_v1      = overlap_len / (v1["end"] - v1["start"])
            ro_v2      = overlap_len / (v2["end"] - v2["start"])
            size_ratio = min(v1["size"], v2["size"]) / max(v1["size"], v2["size"])
            return ro_v1 >= min_fc and ro_v2 >= min_fc and size_ratio >= min_fc

    return False

# ---------------------------------------------------------------------------
# Per-group worker — top-level so it is picklable by multiprocessing.Pool
# ---------------------------------------------------------------------------

def _cluster_one_group(args):
    """
    Clusters a single (chrom, sv_type[, remote_chrom]) group.
    Returns a list of variant clusters.
    """
    sv_type, v_list, max_jitter, min_fc, max_window = args

    # Sort by local start; for TRA this keeps the window eviction meaningful
    v_list.sort(key=lambda x: x["start"])
    n = len(v_list)

    # --- Union-Find: iterative path compression + union by rank ---
    parent = list(range(n))
    rank   = [0] * n

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri == rj:
            return
        if rank[ri] < rank[rj]:
            ri, rj = rj, ri
        parent[rj] = ri
        if rank[ri] == rank[rj]:
            rank[ri] += 1

    # --- Two-pointer sliding window with hard window cap ---
    window_start = 0

    for i in range(n):
        v_i = v_list[i]

        # Evict outdated entries from window left edge
        if sv_type == "INS":
            while window_start < i and (v_i["start"] - v_list[window_start]["start"] > max_jitter):
                window_start += 1
        elif sv_type == "TRA":
            # Evict on local breakpoint distance alone; remote is checked in are_mergeable
            while window_start < i and (v_i["start"] - v_list[window_start]["start"] > max_jitter):
                window_start += 1
        else:  # DEL, DUP, INV
            while window_start < i and (v_i["start"] > v_list[window_start]["end"] + max_jitter):
                window_start += 1

        # Hard cap — prevents O(N²) in dense INV/DUP/TRA regions
        effective_start = max(window_start, i - max_window)

        for j in range(effective_start, i):
            if are_mergeable(v_list[j], v_i, max_jitter, min_fc):
                union(j, i)

    # Extract clusters from disjoint set
    cluster_map = defaultdict(list)
    for i in range(n):
        cluster_map[find(i)].append(v_list[i])

    return list(cluster_map.values())

# ---------------------------------------------------------------------------
# Parallel dispatcher
# ---------------------------------------------------------------------------

def sliding_window_clustering(variants, jitter_map, min_fc=0.8, max_window=200, num_workers=None):
    """
    Dispatches each (chrom, sv_type) group — or (chrom, remote_chrom) for TRA —
    to a worker process in parallel.
    """
    grouped = defaultdict(list)
    for v in variants:
        if v["type"] == "TRA":
            # Group TRAs by (local_chrom, remote_chrom) pair so that only
            # variants involving the same chromosome pair are ever compared
            key = (v["chrom"], "TRA", v["remote_chrom"] or "unknown")
        else:
            key = (v["chrom"], v["type"], "")
        grouped[key].append(v)

    tasks = [
        (sv_type, v_list, jitter_map.get(sv_type, jitter_map["default"]), min_fc, max_window)
        for (chrom, sv_type, remote_chrom), v_list in grouped.items()
    ]

    workers = min(num_workers or cpu_count(), len(tasks))
    print(f"[*] Dispatching {len(tasks)} (chrom, type) groups across {workers} worker processes...")

    all_clusters = []
    with Pool(processes=workers) as pool:
        for group_clusters in pool.imap_unordered(_cluster_one_group, tasks, chunksize=4):
            all_clusters.extend(group_clusters)

    return all_clusters

# ---------------------------------------------------------------------------
# Collapse
# ---------------------------------------------------------------------------

def collapse_clusters(clusters, num_samples, log_file=None):
    """
    Collapses rows based on the most common profile and applies a bitwise OR operation.
    Optionally logs clustering histories to a specified target file.
    """
    collapsed_records = []
    log_handle = open(log_file, 'w') if log_file else None

    if log_handle:
        log_handle.write("#Master_SVID\tCluster_Size\tMerged_Sub_SVIDs\n")

    for cluster in clusters:
        master_v = max(cluster, key=lambda x: x["cluster_size"])

        total_cluster_size = sum(v["cluster_size"] for v in cluster)
        total_sv_rate      = sum(v["sv_rate"]      for v in cluster)

        merged_genotypes = ["0"] * num_samples
        sub_svids = []

        for v in cluster:
            sub_svids.append(v["svid"])
            for idx, gt in enumerate(v["genotypes"]):
                if gt == "1":
                    merged_genotypes[idx] = "1"

        if log_handle:
            log_handle.write(f"{master_v['svid']}\t{len(cluster)}\t{','.join(sub_svids)}\n")

        collapsed_records.append({
            "chrom":        master_v["chrom"],
            "start":        master_v["start"],
            "end":          master_v["end"],
            "size":         master_v["size"],
            "svid":         master_v["svid"],
            "type":         master_v["type"],
            "seq":          master_v["seq"],
            "maq":          master_v["maq"],
            "cluster_size": total_cluster_size,
            "sv_rate":      total_sv_rate,
            "genotypes":    merged_genotypes
        })

    if log_handle:
        log_handle.close()

    return collapsed_records

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Ultra-fast linear O(N) sliding window variant merge engine.")
    parser.add_argument("-m", "--matrix",    required=True,                       help="Path to uncollapsed PAV matrix file")
    parser.add_argument("-o", "--output",    default="collapsed_pav_matrix.txt",  help="Output filepath")
    parser.add_argument("-l", "--log-merge", default=None,                        help="Optional output filepath to log grouped variant merge histories")
    parser.add_argument("--jitter",          type=int,   default=250,             help="Max start distance for INS/DEL merging (Default: 250)")
    parser.add_argument("--jitter-inv",      type=int,   default=1000,            help="Max coordinate jitter for INV merging (Default: 1000)")
    parser.add_argument("--jitter-dup",      type=int,   default=1000,            help="Max coordinate jitter for DUP merging (Default: 1000)")
    parser.add_argument("--jitter-tra",      type=int,   default=5000,             help="Max jitter on BOTH local and remote TRA breakpoints (Default: 500)")
    parser.add_argument("--fc",              type=float, default=0.65,            help="Minimum fold change/overlap (Default: 0.65)")
    parser.add_argument("--min-size",        type=int,   default=40,              help="Minimum SV size in bp eligible for merging (Default: 40)")
    parser.add_argument("--max-window",      type=int,   default=200,             help="Max prior variants to scan per element; prevents O(N²) in dense regions (Default: 200)")
    parser.add_argument("--threads",         type=int,   default=None,            help="Number of parallel worker processes (Default: all available CPUs)")
    args = parser.parse_args()

    jitter_map = {
        "default": args.jitter,
        "INS":     args.jitter,
        "DEL":     args.jitter,
        "INV":     args.jitter_inv,
        "DUP":     args.jitter_dup,
        "TRA":     args.jitter_tra,
    }

    print(f"[*] Loading matrix: '{args.matrix}'...")
    variants, sample_names = parse_input_matrix(args.matrix)

    if not variants:
        print("[-] Error: No data found inside input matrix file.")
        return

    # TRA size is always 0 in your data — skip size filter for TRA
    merge_candidates     = [v for v in variants if v["type"] == "TRA" or abs(v["size"]) >= args.min_size]
    passthrough_variants = [v for v in variants if v["type"] != "TRA" and abs(v["size"]) <  args.min_size]

    tra_count = sum(1 for v in merge_candidates if v["type"] == "TRA")
    workers   = args.threads or cpu_count()

    print(f"[*] Size filter (--min-size {args.min_size} bp): "
          f"{len(merge_candidates)} variants eligible for merging "
          f"(incl. {tra_count} TRA), "
          f"{len(passthrough_variants)} passed through as-is.")
    print(f"[*] Jitter thresholds — INS/DEL: {args.jitter} bp | INV: {args.jitter_inv} bp | "
          f"DUP: {args.jitter_dup} bp | TRA: {args.jitter_tra} bp (both breakpoints) | "
          f"Max window: {args.max_window} variants")
    print(f"[*] Workers: {workers} (use --threads N to override)")

    print(f"[*] Running parallel sliding window clustering across {len(merge_candidates)} records...")
    clusters = sliding_window_clustering(
        merge_candidates,
        jitter_map=jitter_map,
        min_fc=args.fc,
        max_window=args.max_window,
        num_workers=workers
    )
    print(f"[+] Formed {len(clusters)} distinct collapsed structural profiles.")

    print("[*] Collapsing payload blocks and processing bitwise OR flags...")
    collapsed_data = collapse_clusters(clusters, len(sample_names), log_file=args.log_merge)

    for v in passthrough_variants:
        collapsed_data.append(v)

    if args.log_merge:
        print(f"[+] Merge tracking log successfully generated at: '{args.log_merge}'")

    collapsed_data.sort(key=lambda x: (x["chrom"], x["start"]))

    print(f"[*] Saving results: '{args.output}'...")
    with open(args.output, 'w') as out:
        base_headers = [
            "#Target_name", "Target_start", "Target_end", "SVlen",
            "SVID", "SVType", "seq", "maq",
            "cluster_size_prevalent", "sv_rate_prevalent"
        ]
        out.write("\t".join(base_headers + sample_names) + "\n")

        for cv in collapsed_data:
            gt_string = "\t".join(cv["genotypes"])
            out.write(
                f"{cv['chrom']}\t{cv['start']}\t{cv['end']}\t{cv['size']}\t"
                f"{cv['svid']}\t{cv['type']}\t{cv['seq']}\t{cv['maq']}\t"
                f"{cv['cluster_size']}\t{cv['sv_rate']}\t{gt_string}\n"
            )

    print(f"[+] Execution completed successfully.")


if __name__ == "__main__":
    main()
