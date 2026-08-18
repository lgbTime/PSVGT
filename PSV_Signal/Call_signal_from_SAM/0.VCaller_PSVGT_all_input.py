#!/usr/bin/env python
"""
0.VCaller_PSVGT.py - Chromosome-parallel SV caller.

Auto-detects input format by extension:
  .bam / .sam               ? direct variant calling on alignments
  .fa / .fasta / .fq / .fastq ? run minimap2 then call variants from SAM

Output per chromosome:
  <prefix>_<chr>.record.txt           - CIGAR-based DEL/INS
  <prefix>_<chr>.record.txt.cov       - coverage track
  <prefix>_<chr>.record.txt.suppAlign - segment-pair SVs (INS/DEL/DUP/INV/TRA)
  <prefix>_<chr>.record.txt.depth     - per-chromosome depth
  <prefix>_<chr>.record.txt.breakpoint- soft/hard clip breakpoints (LR mode)

Usage:
    python VCaller.py -b alignment.bam -o out -dtype cr -fai ref.fa.fai -msv yes
    python VCaller.py -b alignment.sam -o out -dtype cr -fai ref.fa.fai -msv yes
    python VCaller.py -b reads.fastq -ref genome.fa -fai genome.fa.fai -o out -dtype hifi -msv yes
"""

import os
import sys
import argparse
import subprocess
from time import time
from collections import defaultdict
from multiprocessing import Pool

import sub_VCaller_PSVGT as worker
from sub_Signaling_Module import sam_line_to_dict


# ---------------------------------------------------------------------------
# FAI reader (no pandas dependency)
# ---------------------------------------------------------------------------

def read_fai(fai_path):
    """Read .fai index ? (chromosome_order, {chr: size})."""
    chrom_order = []
    chrom_sizes = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_order.append(parts[0])
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_order, chrom_sizes


# ---------------------------------------------------------------------------
# Helper: detect input format by extension
# ---------------------------------------------------------------------------

def detect_input_type(path):
    """Return 'bam', 'sam', or 'fastx' based on file extension."""
    p = path.lower()
    if p.endswith('.bam'):
        return 'bam'
    if p.endswith('.sam'):
        return 'sam'
    for ext in ('.fasta', '.fa', '.fastq', '.fq',
                '.fasta.gz', '.fa.gz', '.fastq.gz', '.fq.gz'):
        if p.endswith(ext):
            return 'fastx'
    return None


# ---------------------------------------------------------------------------
# Minimap2 pipeline
# ---------------------------------------------------------------------------

MINIMAP2_PRESETS = {
    'cr':   'asm5',
    'sr':   'asm5',
    'pb':   'map-pb',
    'ont':  'map-ont',
    'hifi': 'map-hifi',
}


def run_minimap2(query_path, ref_path, dtype, threads=4):
    """Run minimap2, yield SAM lines (str) from stdout.

    Uses --secondary=no to suppress 0x100 secondary alignments,
    matching the earlier BAM-generation pipeline.
    """
    preset = MINIMAP2_PRESETS.get(dtype, 'map-hifi')
    cmd = ['minimap2', '-ax', preset, '--secondary=no',
           '-t', str(threads), ref_path, query_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)
    for line in proc.stdout:
        yield line
    proc.stdout.close()
    stderr_output = proc.stderr.read()
    proc.wait()
    if proc.returncode != 0:
        sys.stderr.write(stderr_output)
        raise subprocess.CalledProcessError(proc.returncode, cmd)


# ---------------------------------------------------------------------------
# SAM-mode: read + group by chromosome (file path or line iterable)
# ---------------------------------------------------------------------------

def group_sam_by_chromosome(source, chrom_sizes):
    """Read SAM and group records by chromosome.

    source can be a file path (str) or an iterable of SAM lines.
    Returns {chrom: [dict_records]}.
    """
    groups = defaultdict(list)
    if isinstance(source, str):
        fh = open(source)
        own = True
    else:
        fh = source
        own = False
    try:
        for line in fh:
            rec = sam_line_to_dict(line)
            if rec is None:
                continue
            if rec['flag'] == 4 or rec['rname'] == '*':
                continue
            chrom = rec['rname']
            if chrom in chrom_sizes:
                groups[chrom].append(rec)
    finally:
        if own:
            fh.close()
    return dict(groups)


# ---------------------------------------------------------------------------
# Greedy load balancing
# ---------------------------------------------------------------------------

def distribute_by_size(chrom_groups, chrom_sizes, n_workers):
    """Assign chromosomes to workers balancing by record count."""
    sorted_chroms = sorted(chrom_groups.keys(),
                           key=lambda c: -len(chrom_groups[c]))
    worker_batches = [[] for _ in range(n_workers)]
    worker_loads   = [0] * n_workers

    for chrom in sorted_chroms:
        idx = min(range(n_workers), key=lambda i: worker_loads[i])
        worker_batches[idx].append(chrom)
        worker_loads[idx] += len(chrom_groups[chrom])

    return [b for b in worker_batches if b]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Chromosome-parallel SV caller (BAM / SAM / FASTA / FASTQ input)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    io_grp = parser.add_argument_group("Input / Output")
    io_grp.add_argument("-b", dest="input", required=True,
                        help="Input: .bam / .sam alignment, or .fasta/.fastq reads")
    io_grp.add_argument("-o", dest="out", required=True,
                        help="Output prefix")
    io_grp.add_argument("-fai", dest="fai", required=True,
                        help="Reference .fai index file")
    io_grp.add_argument("-ref", dest="ref", default=None,
                        help="Reference genome FASTA (required for FASTA/FASTQ input; "
                             "auto-guessed from .fai if omitted)")

    params = parser.add_argument_group("Parameters")
    params.add_argument("-dtype", dest="dtype", default="cr",
                        choices=["cr", "sr", "pb", "ont", "hifi"],
                        help="Data type (cr/sr = assembly; pb/ont/hifi = long read)")
    params.add_argument("-m", dest="min", type=int, default=40,
                        help="Minimum SV length (bp)")
    params.add_argument("-M", dest="max", type=int, default=200000,
                        help="Maximum SV length (bp)")
    params.add_argument("-maq", dest="maq", type=int, default=50,
                        help="Minimum mapping quality")
    params.add_argument("-msv", dest="msv", default="no",
                        choices=["yes", "no"],
                        help="Detect complex SVs from supplementary alignments")
    params.add_argument("-w", dest="workers", type=int, default=10,
                        help="Number of worker processes")

    args = parser.parse_args()

    chrom_order, chrom_sizes = read_fai(args.fai)
    print(f"  Reference has {len(chrom_order)} chromosomes")
    print(f"  Data type: {args.dtype}  |  msv: {args.msv}  |  workers: {args.workers}")

    input_path = args.input
    input_type = detect_input_type(input_path)

    if input_type is None:
        sys.exit(f"ERROR: unrecognised input format: {input_path}\n"
                 f"  Expected .bam, .sam, .fasta/.fa, or .fastq/.fq")

    n_workers = min(args.workers, len(chrom_order))
    t0 = time()

    # --- FASTA / FASTQ ? minimap2 ? SAM pipeline ----------------------------
    if input_type == 'fastx':
        ref_path = args.ref
        if ref_path is None:
            # Auto-guess: strip .fai extension
            if args.fai.endswith('.fai'):
                candidate = args.fai[:-4]
                if os.path.exists(candidate):
                    ref_path = candidate
            if ref_path is None:
                sys.exit("ERROR: -ref <genome.fasta> is required for FASTA/FASTQ input")

        print(f"  FASTA/FASTQ mode: {input_path}")
        print(f"  Reference: {ref_path}")
        print(f"  Running minimap2 -ax {MINIMAP2_PRESETS.get(args.dtype)} ...")

        # minimap2 ? SAM lines ? group by chromosome
        sam_lines = run_minimap2(input_path, ref_path, args.dtype,
                                 threads=min(args.workers, 8))
        chrom_groups = group_sam_by_chromosome(sam_lines, chrom_sizes)

        n_records = sum(len(v) for v in chrom_groups.values())
        print(f"  {n_records} mapped records across {len(chrom_groups)} chromosomes")

        # Ensure all chromosomes present
        for chrom in chrom_order:
            if chrom not in chrom_groups:
                chrom_groups[chrom] = []

        worker_batches = distribute_by_size(chrom_groups, chrom_sizes, n_workers)
        print(f"  Distributing across {len(worker_batches)} workers:")
        for i, batch in enumerate(worker_batches):
            n = sum(len(chrom_groups[c]) for c in batch)
            print(f"    Worker {i}: {len(batch)} chroms, {n} records  "
                  f"[{', '.join(batch[:4])}{'...' if len(batch) > 4 else ''}]")

        tasks = []
        for chrom in chrom_order:
            tasks.append((
                chrom,
                chrom_sizes.get(chrom, 0),
                chrom_order,
                chrom_groups.get(chrom, []),
                args.min,
                args.max,
                args.maq,
                args.out,
                args.dtype,
                args.msv,
            ))

        with Pool(processes=n_workers) as pool:
            results = pool.starmap(worker.process_chromosome_dicts, tasks)

    # --- SAM mode -----------------------------------------------------------
    elif input_type == 'sam':
        print(f"  SAM mode: reading {input_path} ...")
        chrom_groups = group_sam_by_chromosome(input_path, chrom_sizes)
        n_records = sum(len(v) for v in chrom_groups.values())
        print(f"  {n_records} records across {len(chrom_groups)} chromosomes")

        for chrom in chrom_order:
            if chrom not in chrom_groups:
                chrom_groups[chrom] = []

        worker_batches = distribute_by_size(chrom_groups, chrom_sizes, n_workers)
        print(f"  Distributing across {len(worker_batches)} workers:")
        for i, batch in enumerate(worker_batches):
            n = sum(len(chrom_groups[c]) for c in batch)
            print(f"    Worker {i}: {len(batch)} chroms, {n} records  "
                  f"[{', '.join(batch[:4])}{'...' if len(batch) > 4 else ''}]")

        tasks = []
        for chrom in chrom_order:
            tasks.append((
                chrom,
                chrom_sizes.get(chrom, 0),
                chrom_order,
                chrom_groups.get(chrom, []),
                args.min,
                args.max,
                args.maq,
                args.out,
                args.dtype,
                args.msv,
            ))

        with Pool(processes=n_workers) as pool:
            results = pool.starmap(worker.process_chromosome_dicts, tasks)

    # --- BAM mode -----------------------------------------------------------
    else:
        print(f"  BAM mode: {input_path}")
        tasks = []
        for chrom in chrom_order:
            tasks.append((
                chrom,
                chrom_sizes.get(chrom, 0),
                chrom_order,
                input_path,
                args.min,
                args.max,
                args.maq,
                args.out,
                args.dtype,
                args.msv,
            ))

        with Pool(processes=n_workers) as pool:
            results = pool.starmap(worker.process_chromosome_bam, tasks)

    # Report
    total_sv   = 0
    total_supp = 0
    for chrom, sv_cnt, supp_cnt in results:
        total_sv   += sv_cnt
        total_supp += supp_cnt
        if sv_cnt > 0 or supp_cnt > 0:
            print(f"    {chrom}: {sv_cnt} CIGAR-SVs, {supp_cnt} segment-SVs")

    elapsed = time() - t0
    print(f"\n{'=' * 60}")
    print(f"  SV calling complete  ({elapsed:.1f}s)")
    print(f"  Total CIGAR-based SVs:  {total_sv}")
    print(f"  Total segment-pair SVs: {total_supp}")
    print(f"  Output prefix: {args.out}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
