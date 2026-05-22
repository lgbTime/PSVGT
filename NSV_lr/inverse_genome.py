import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def apply_inversions(fasta_in, sv_file, fasta_out):
    # 1. Load the genome into a dictionary for easy access
    print(f"📖 Loading genome: {fasta_in}")
    genome = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))
    
    # 2. Load SVs (assuming the format you provided)
    # Columns: chr, start, end, size, id, type, ...
    print(f"📋 Loading SV list: {sv_file}")
    svs = pd.read_csv(sv_file, sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'size', 'id', 'type', 'junk1', 'junk2', 'junk3', 'junk4'])
    
    # Filter for only Inversions
    invs = svs[svs['type'] == 'INV'].copy()
    
    # 3. Process by chromosome
    for chrom, group in invs.groupby('chr'):
        if chrom not in genome:
            print(f"⚠️ Warning: {chrom} not found in FASTA. Skipping.")
            continue
            
        print(f"✂️ Processing {len(group)} inversions on {chrom}...")
        
        # Convert chromosome sequence to a mutable list (faster manipulation)
        seq_list = list(str(genome[chrom].seq))
        
        # IMPORTANT: Sort inversions by START position in DESCENDING order
        # This prevents coordinate shifting if you later add DEL/INS logic
        group = group.sort_values('start', ascending=False)
        
        for _, row in group.iterrows():
            # Adjust to 0-based indexing (Standard BED/TSV is usually 0-based start)
            start = int(row['start'])
            end = int(row['end'])
            
            # Extract the segment
            segment = seq_list[start:end]
            
            # Reverse the segment (Inversion)
            # Note: For biological inversions, we usually do Reverse Complement.
            # Use segment[::-1] for simple reverse, or the logic below for RevComp.
            reversed_segment = list(str(Seq("".join(segment)).reverse_complement()))
            
            # Re-insert into the sequence list
            seq_list[start:end] = reversed_segment
        
        # Update the genome object
        genome[chrom].seq = Seq("".join(seq_list))

    # 4. Save the modified genome
    print(f"💾 Saving modified genome to: {fasta_out}")
    with open(fasta_out, "w") as out_handle:
        SeqIO.write(genome.values(), out_handle, "fasta")
    print("🏁 Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify genome with Inversions")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA")
    parser.add_argument("-s", "--sv", required=True, help="Input SV TSV")
    parser.add_argument("-o", "--out", required=True, help="Output FASTA")
    args = parser.parse_args()
    apply_inversions(args.fasta, args.sv, args.out)
