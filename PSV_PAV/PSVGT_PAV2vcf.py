import argparse
import sys
import os

def convert_pav_to_vcf(input_file, output_file, resolve_duplicates=True):
    """
    Parses a PAV matrix file and converts it into a valid VCF v4.2 file.
    """
    if not os.path.exists(input_file):
        sys.exit(f"Error: Input file '{input_file}' does not exist.")

    print(f"[*] Reading PAV Matrix: {input_file}")
    print(f"[*] Writing VCF Output : {output_file}")
    if resolve_duplicates:
        print("[*] Duplicate coordinate protection: ENABLED (+1bp shift)")

    # Track coordinates to resolve duplicates { chrom: set(pos1, pos2, ...) }
    seen_positions = {}
    record_count = 0
    duplicate_fixes = 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 1. Parse header
        header_line = infile.readline().strip()
        if not header_line:
            sys.exit("Error: The input file is empty.")
            
        if not header_line.startswith("#Target_name"):
            header_line = header_line.lstrip('#')
        
        fields = header_line.split('\t')
        
        # Metadata columns layout verification
        meta_cols = ["Target_name", "Target_start", "Target_end", "SVlen", "SVID", 
                     "SVType", "seq", "maq", "cluster_size_prevalent", "sv_rate_prevalent"]
        
        for i, col in enumerate(meta_cols):
            if i >= len(fields):
                sys.exit(f"Error: Missing expected column metadata at index {i}.")
            actual = fields[i].lower().replace('#', '')
            expected = col.lower()
            if expected not in actual:
                print(f"Warning: Column index {i} expected '{col}', but found '{fields[i]}'")

        # Extract sample names (everything from column index 10 onwards)
        samples = fields[10:]
        if not samples:
            sys.exit("Error: No sample columns detected after the 10 metadata columns.")
        
        # 2. Write VCF Standard Headers
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
        outfile.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
        outfile.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
        outfile.write("##ALT=<ID=TRA,Description=\"Translocation\">\n")
        outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        
        # Write VCF column header line
        vcf_main_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
        outfile.write("\t".join(vcf_main_header) + "\n")
        
        # 3. Process Data Rows
        for line_num, line in enumerate(infile, start=2):
            if not line.strip():
                continue
            
            row = line.strip().split('\t')
            if len(row) < 11:
                print(f"Warning: Skipping malformed row at line {line_num} (insufficient columns).")
                continue
            
            chrom = row[0]
            try:
                pos = int(row[1])
            except ValueError:
                print(f"Warning: Skipping line {line_num} due to invalid integer position: '{row[1]}'.")
                continue
                
            sv_id = row[4]
            sv_type = row[5]
            ref_seq = row[6]
            qual = row[7]
            
            # Sanitize reference allele text markers
            if ref_seq == "*" or ref_seq == '"*"':
                ref_seq = "N"
            
            # Format appropriate structural variant allele types
            if ref_seq != "N" and sv_type == "INS":
                alt_field = ref_seq
            else:
                alt_field = f"<{sv_type}>"
            
            # Apply dynamic deduplication strategy
            original_pos = pos
            if resolve_duplicates:
                if chrom not in seen_positions:
                    seen_positions[chrom] = set()
                
                while pos in seen_positions[chrom]:
                    pos += 1
                    duplicate_fixes += 1
                seen_positions[chrom].add(pos)
            
            # Generate INFO strings fields
            info_field = f"SVTYPE={sv_type};END={row[2]};SVLEN={row[3]}"
            if pos != original_pos:
                info_field += f";ORIGINAL_POS={original_pos}"
                
            format_field = "GT"
            
            # Map PAV elements (0 -> 0/0, 1 -> 1/1)
            genotypes = []
            for gt_val in row[10:]:
                clean_gt = gt_val.strip()
                if clean_gt == "0":
                    genotypes.append("0/0")
                elif clean_gt == "1":
                    genotypes.append("1/1")
                else:
                    genotypes.append("./.")  # Missing or low confidence fallback
            
            vcf_row = [
                chrom, str(pos), sv_id, ref_seq[0], alt_field, 
                qual, "PASS", info_field, format_field
            ] + genotypes
            
            outfile.write("\t".join(vcf_row) + "\n")
            record_count += 1

    print(f"[+] Success! Successfully converted {record_count} variants.")
    if resolve_duplicates:
        print(f"[+] Applied +1bp shifts to {duplicate_fixes} duplicate genomic coordinates.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Presence-Absence Variation (PAV) structural variant matrix into a standard VCF v4.2 file."
    )
    
    # Position arguments or short/long tags
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input PAV matrix text file (tab-separated)."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path where the output VCF file should be written."
    )
    parser.add_argument(
        "--no-offset", 
        action="store_false", 
        dest="resolve_duplicates",
        help="Disable the automatic +1bp offset shift protection for overlapping/duplicated positions."
    )

    args = parser.parse_args()
    
    convert_pav_to_vcf(
        input_file=args.input, 
        output_file=args.output, 
        resolve_duplicates=args.resolve_duplicates
    )

if __name__ == "__main__":
    main()
