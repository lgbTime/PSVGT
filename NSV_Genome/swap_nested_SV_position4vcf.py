import pandas as pd
import argparse
import re
import sys

def update_info_field(info_str, new_start, new_end):
    """Updates SV_START and END values inside the VCF INFO string."""
    info_dict = {}
    parts = info_str.split(';')
    for part in parts:
        if '=' in part:
            k, v = part.split('=', 1)
            info_dict[k] = v
        else:
            info_dict[part] = None

    info_dict['SV_START'] = str(new_start)
    info_dict['END'] = str(new_end)
    
    return ";".join([f"{k}={v}" if v is not None else k for k, v in info_dict.items()])

def unflip_vcf(vcf_file, inv_file, output_file):
    # Load Inversion list
    # Expected: chrom, inv_start, inv_end, ...
    inv_cols = ['chrom', 'inv_start', 'inv_end', 'inv_len', 'SVID', 'type', 'seq', 'maq', 'c', 'r']
    try:
        inv_df = pd.read_csv(inv_file, sep='\t', names=inv_cols, header=None)
    except Exception as e:
        print(f"Error loading inversion list: {e}")
        return

    # Extract VCF headers and column names
    header_lines = []
    vcf_cols = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line)
            elif line.startswith('#CHROM'):
                vcf_cols = line.strip().split('\t')
                break
    
    if not vcf_cols:
        print("Error: Could not find #CHROM header line in VCF.")
        return

    df = pd.read_csv(vcf_file, sep='\t', comment='#', names=vcf_cols)
    corrected_count = 0

    def process_vcf_row(row):
        nonlocal corrected_count
        chrom = str(row['#CHROM'])
        old_pos = int(row['POS'])
        
        # Parse END from INFO field
        end_match = re.search(r'END=(\d+)', row['INFO'])
        old_end = int(end_match.group(1)) if end_match else old_pos

        # Filter for inversions on this chromosome
        relevant_invs = inv_df[inv_df['chrom'] == chrom]
        
        for _, inv in relevant_invs.iterrows():
            # Check if SV is inside OR flanking the inversion (within 1kb buffer)
            # This captures your 1200-1400 case for a 1000-1200 inversion
            if (old_pos >= inv['inv_start'] - 10 and old_pos <= inv['inv_end'] + 10):
                # THE MIRROR TRANSFORMATION
                # To handle the flip and shift:
                # new_point = inv_start - (old_point - inv_end)
                
                t1 = inv['inv_start'] - (old_pos - inv['inv_end'])
                t2 = inv['inv_start'] - (old_end - inv['inv_end'])
                
                # Ensure the new VCF record has Start < End
                new_start = int(min(t1, t2))
                new_end = int(max(t1, t2))
                
                # 1. Update POS column
                row['POS'] = new_start
                
                # 2. Update INFO (SV_START and END)
                row['INFO'] = update_info_field(row['INFO'], new_start, new_end)
                
                # 3. Update ID string (format: chrom:start-end_length)
                # We use regex to preserve the SV length suffix
                id_match = re.search(r'_(\d+)$', row['ID'])
                sv_len_suffix = id_match.group(1) if id_match else "len"
                row['ID'] = f"{chrom}:{new_start}-{new_end}_{sv_len_suffix}"
                
                corrected_count += 1
                #break 
        return row

    # Apply processing
    print(f"Processing {len(df)} SV records...")
    df = df.apply(process_vcf_row, axis=1)
    df = df.sort_values(['#CHROM', 'POS'])

    # Write back to VCF
    with open(output_file, 'w') as f:
        f.writelines(header_lines)
        df.to_csv(f, sep='\t', index=False, header=True)

    print("-" * 30)
    print(f"✅ Success! Output written to: {output_file}")
    print(f"🔄 Corrected {corrected_count} records based on inversion coordinates.")

def main():
    parser = argparse.ArgumentParser(description="Un-flip VCF coordinates for nested and flanking SVs.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-l", "--list", required=True, help="Inversion list (TAB separated)")
    parser.add_argument("-o", "--output", default="unflipped_output.vcf", help="Output VCF filename")
    args = parser.parse_args()
    
    unflip_vcf(args.input, args.list, args.output)

if __name__ == "__main__":
    main()
