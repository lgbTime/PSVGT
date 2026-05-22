import os
import subprocess
import argparse
import pandas as pd
import shutil

def run_cmd(cmd):
    print(f"\n>>> Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def count_inv(file_path):
    """Parses the PopSV table and returns the count of INV signals."""
    if not os.path.exists(file_path):
        return 0
    try:
        df = pd.read_csv(file_path, sep='\t')
        inv_count = len(df[df['SVType'] == 'INV'])
        return inv_count
    except Exception as e:
        print(f"Warning: Could not parse {file_path}: {e}")
        return 0

def main():
    parser = argparse.ArgumentParser(description="Recursive PSVGT Inversion Pipeline")
    parser.add_argument("-r", "--ref", required=True, help="Initial reference FASTA file")
    parser.add_argument("-q", "--query", required=True, help="Input query folder (containing genome/reads)")
    parser.add_argument("--psvgt_dir", default="~/bin/PSVGT1.0", help="Base directory for PSVGT scripts")
    args = parser.parse_args()

    # --- Path Configuration ---
    qname= os.listdir(args.query)[0]
    base_dir = os.path.expanduser(args.psvgt_dir)
    nsv_dir = os.path.join(base_dir, "NSV_Genome")
    PSVGT_EXE = os.path.join(base_dir, "PSVGT1.0.py")
    CLUSTER_EXE = os.path.join(base_dir, "PSV_Signal/1.PSV_signal_cluster.py")
    INVERSE_EXE = os.path.join(nsv_dir, "inverse_genome.py") # Assuming it lives here
    SWAP_EXE = os.path.join(nsv_dir, "swap_nested_SV_position.py")
    SWAP_EXE4VCF = os.path.join(nsv_dir, "swap_nested_SV_position4vcf.py")
    VCF_EXE = os.path.join(nsv_dir, "PSVGT2vcf.py")
    CSV_EXE = os.path.join(nsv_dir, "vcf2CSV.py")

    current_ref = os.path.abspath(args.ref)
    query_folder = os.path.abspath(args.query)
    
    round_idx = 0
    prev_inv_count = -1
    inv_lists = []
    output_dirs = []

    # --- Step 1: Iterative Discovery Loop ---
    while True:
        suffix = f"_round{round_idx}" if round_idx > 0 else ""
        out_dir = f"nested_genome{suffix}_PSVGT_out"
        candidate_file = f"{out_dir}/PopSV_Candidate_Record.txt"
        print(f"\n{'='*20} Starting Round {round_idx} {'='*20}")
        
        # Using the defined -q (query) folder here
        run_cmd(f"python {PSVGT_EXE} -cr {query_folder} -o {out_dir} -msv yes -m 45 -r {current_ref} -csv 0.4")
        output_dirs.append(os.path.abspath(out_dir))
        current_inv_count = count_inv(candidate_file)
        print(f"Round {round_idx} Result: Found {current_inv_count} INV signals.")
        
        if current_inv_count == 0:
            print("Termination: 0 INV signals detected.")
            break
        if prev_inv_count != -1 and current_inv_count >= prev_inv_count:
            print(f"Termination: INV count ({current_inv_count}) no longer decreasing.")
            break

        inv_list = os.path.abspath(f"inv_r{round_idx}.list")
        run_cmd(f"grep 'INV' {candidate_file} > {inv_list}")
        
        inv_lists.append(inv_list)
        next_ref = os.path.abspath(f"inv_ref_round{round_idx+1}.fa")
        run_cmd(f"python {INVERSE_EXE} -f {current_ref} -s {inv_list} -o {next_ref}")
        
        current_ref = next_ref
        prev_inv_count = current_inv_count
        round_idx += 1

    if not output_dirs:
        print("No inversions were found in the first round. Exiting.")
        return

    # --- Step 2: Signal Clustering ---
    cluster_dir = os.path.abspath("final_combined_output")
    os.makedirs(cluster_dir, exist_ok=True)
    
    print("\n>>> Preparing Clustering Directory...")
    for i, d in enumerate(output_dirs):
        src = os.path.join(d, f"0_tmp_{qname}_Clustered_Record.txt")
        dst = os.path.join(cluster_dir, f"C{i}_Clustered_Record.txt")
        if os.path.exists(src):
            if os.path.lexists(dst): os.remove(dst)
            os.symlink(src, dst)
    run_cmd(f"cd {cluster_dir} && python {CLUSTER_EXE} -d ./")

    # --- Step 3: Sequential Coordinate Un-flipping ---
    print("\n>>> Correcting Nested SV Coordinates (Un-flipping)...")
    current_swap_file = os.path.join(cluster_dir, "PopSV_Candidate_Record.txt")
    for i, inv_l in enumerate(inv_lists):
        out_swap =os.path.join(cluster_dir, f"PopSV_Candidate_Record.txt_swap{i}")
        vcf = f"nested_genome_round{i+1}_PSVGT_out/2_tmp_{qname}_genotype.vcf"
        run_cmd(f"python {SWAP_EXE4VCF} -i {vcf} -l {inv_l} -o {vcf}.swap{i}")
        run_cmd(f"python {SWAP_EXE} -i {current_swap_file} -l {inv_l} -o {out_swap}")
        run_cmd(f"python {VCF_EXE} -i {out_swap} -o {out_swap}.vcf")
        run_cmd(f"python {CSV_EXE} -i {out_swap}.vcf -o {out_swap}.vcf.csv")
        current_swap_file = out_swap

    # --- Step 4: VCF and Final CSV Generation ---
    vcf_out = "Final_falseGT_Nested_SVs.vcf"
    csv_out = "Final_falseGT_Nested_SVs.csv"

    run_cmd(f"python {VCF_EXE} -i {current_swap_file} -o {vcf_out}")
    run_cmd(f"python {CSV_EXE} -i {vcf_out} -o {csv_out}")
    
    print(f"\n✅ Pipeline Complete. Results in: {os.path.join(cluster_dir, csv_out)}")

if __name__ == "__main__":
    main()
