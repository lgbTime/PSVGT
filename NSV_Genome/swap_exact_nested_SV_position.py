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

    # Load Predicted Candidates
    df = pd.read_csv(candidate_file, sep='\t')
    df.rename(columns={'#Target_name': 'chrom', 'Target_start': 'start', 'Target_end': 'end'}, inplace=True)
    corrected_count = 0

    def process_row(row):
        nonlocal corrected_count
        chrom_str = str(row['chrom'])
        relevant_invs = inv_df[inv_df['chrom'] == chrom_str]
        
        for _, inv in relevant_invs.iterrows():
            # Check if the nested SV is within the host inversion boundaries
            if row['start'] >= inv['inv_start'] and row['end'] <= inv['inv_end']:
                
                # THE MIRROR TRANSFORMATION
                # P_new = Inv_End - (P_old - Inv_Start)
                # This simplifies to: P_new = Inv_End + Inv_Start - P_old
                
                old_start = row['start']
                old_end = row['end']
                
                # Flip the start and end relative to the inversion boundaries
                new_end = inv['inv_end'] - (old_start - inv['inv_start'])
                new_start = inv['inv_end'] - (old_end - inv['inv_start'])
                
                # Update row with the "un-flipped" coordinates
                row['start'] = int(new_start)
                row['end']   = int(new_end)
                corrected_count += 1
                #break # Only one inversion host per nested SV
                if new_end != inv['inv_end'] and new_start != inv['inv_start']:
                    print(f'INV={chrom_str}:{inv["inv_start"]}-{inv["inv_end"]}\t{chrom_str}:{new_start}_{new_end}_{row["SVType"]}={row["SVlen"]}')
        
        return row

    df = df.apply(process_row, axis=1)
    df.rename(columns={'chrom': '#Target_name', 'start': 'Target_start', 'end': 'Target_end'}, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"✅ Processed {len(df)} lines.")
    print(f"🔄 Un-flipped {corrected_count} nested SVs back to Reference Coordinates.")

def main():
    parser = argparse.ArgumentParser(description="Un-flip nested SV coordinates from an inverted reference back to the standard reference.")
    parser.add_argument("-i", "--input", required=True, help="PopSV_Candidate_Record.txt")
    parser.add_argument("-l", "--list", required=True, help="inv_r2.list (Host inversions)")
    parser.add_argument("-o", "--output", default="PopSV_StandardCoords.txt", help="Output filename")
    args = parser.parse_args()
    unflip_coordinates(args.input, args.list, args.output)

if __name__ == "__main__":
    main()
