import pandas as pd
import argparse
import sys

def unflip_coordinates(candidate_file, inv_file, output_file):
    inv_cols = ['chrom', 'inv_start', 'inv_end', 'inv_len', 'SVID', 'type', 'seq', 'maq', 'c', 'r']
    try:
        inv_df = pd.read_csv(inv_file, sep='\t', names=inv_cols, header=None)
    except Exception as e:
        print(f"Error loading inversion list: {e}")
        return

    # Load Predicted Candidates (Genotype txt file)
    df = pd.read_csv(candidate_file, sep='\t')
    
    # We will use the original column names from your cat output
    corrected_count = 0

    def process_row(row):
        nonlocal corrected_count
        chrom_str = str(row['#Target_name'])
        old_start = row['Target_start']
        old_end = row['Target_end']
        
        relevant_invs = inv_df[inv_df['chrom'] == chrom_str]
        
        for _, inv in relevant_invs.iterrows():
            # Check if SV is inside OR flanking the inversion (within 1kb buffer)
            # This captures your 1200-1400 case for a 1000-1200 inversion
            if (old_start >= inv['inv_start'] - 10 and old_start <= inv['inv_end'] + 10):
                
                # THE MIRROR TRANSFORMATION
                # To handle the flip and shift correctly for flanking regions:
                # new_point = inv_start - (old_point - inv_end)
                
                t1 = inv['inv_start'] - (old_start - inv['inv_end'])
                t2 = inv['inv_start'] - (old_end - inv['inv_end'])
                
                # Ensure the standard reference has Start < End
                new_start = int(min(t1, t2))
                new_end = int(max(t1, t2))
                
                # Update row
                row['Target_start'] = new_start
                row['Target_end'] = new_end
                corrected_count += 1
                
                # Print log for verification
                print(f'INV={chrom_str}:{inv["inv_start"]}-{inv["inv_end"]} | '
                      f'Mapped {old_start}-{old_end} -> {new_start}-{new_end}')
                row['SVID']= f'{chrom_str}:{new_start}-{new_end}_{row["SVType"]}={row["SVlen"]}'
        return row

    # Apply transformation
    df = df.apply(process_row, axis=1)
    
    # Save output
    df.to_csv(output_file, sep='\t', index=False)
    
    print("-" * 30)
    print(f"✅ Processed {len(df)} lines.")
    print(f"🔄 Corrected {corrected_count} SV positions (including flanking regions).")

def main():
    parser = argparse.ArgumentParser(description="Un-flip SV coordinates (inside and flanking) back to standard reference.")
    parser.add_argument("-i", "--input", required=True, help="Input genotype txt file")
    parser.add_argument("-l", "--list", required=True, help="inv_r2.list (Host inversions)")
    parser.add_argument("-o", "--output", default="Corrected_Genotypes.txt", help="Output filename")
    args = parser.parse_args()
    unflip_coordinates(args.input, args.list, args.output)

if __name__ == "__main__":
    main()
