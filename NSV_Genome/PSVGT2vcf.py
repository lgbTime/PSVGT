import pandas as pd
import datetime
import argparse
import sys
import re

def get_genotype(chrom):
    """
    Determines genotype based on chromosome name.
    chr1-chr12: 1/1
    chr13-chrX: 0/1
    """
    # Extract the number from the chromosome string (e.g., 'chr10' -> 10, '10' -> 10)
    match = re.search(r'(\d+)', chrom)
    
    if match:
        chrom_num = int(match.group(1))
        if 1 <= chrom_num <= 12:
            return "1/1"
        elif chrom_num >= 13:
            return "0/1"
    
    # Handle X, Y, or other non-numeric chromosomes
    if 'X' in chrom.upper():
        return "0/1"
    
    # Default fallback
    return "./."

def convert_table_to_vcf(input_file, output_vcf):
    try:
        df = pd.read_csv(input_file, sep='\t')
        
        # Natural Sorting
        df['#Target_name'] = df['#Target_name'].astype(str)
        df['chr_sort_num'] = df['#Target_name'].str.extract('(\d+)').fillna(999).astype(int)
        df = df.sort_values(by=['chr_sort_num', '#Target_name', 'Target_start']).drop('chr_sort_num', axis=1)

        with open(output_vcf, 'w') as vcf:
            now = datetime.datetime.now().strftime("%Y%m%d")
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##fileDate={now}\n")
            
            # INFO and FORMAT header definitions
            vcf.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
            vcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
            vcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
            vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            
            # Symbolic alleles
            vcf.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
            vcf.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
            
            # Column headers (Note the addition of FORMAT and SAMPLE)
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            for _, row in df.iterrows():
                chrom = row['#Target_name']
                pos = int(row['Target_start'])
                svid = row['SVID']
                svtype = row['SVType']
                svlen = int(row['SVlen'])
                end = int(row['Target_end'])
                
                # Get dynamic genotype
                gt = get_genotype(chrom)
                
                alt = f"<{svtype}>"
                info = f"SVTYPE={svtype};END={end};SVLEN={svlen}"
                
                # Format: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
                vcf.write(f"{chrom}\t{pos}\t{svid}\tN\t{alt}\t60\tPASS\t{info}\tGT\t{gt}\n")

        print(f"Successfully generated VCF with conditional genotypes: {output_vcf}")

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Convert SV table to VCF with chromosome-specific genotypes.")
    parser.add_argument("-i", "--input", required=True, help="Input variant table (TSV)")
    parser.add_argument("-o", "--output", default="output.vcf", help="Output VCF filename")
    
    args = parser.parse_args()
    convert_table_to_vcf(args.input, args.output)

if __name__ == "__main__":
    main()
