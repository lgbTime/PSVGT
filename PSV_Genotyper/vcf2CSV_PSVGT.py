import argparse
import statistics
import re

def main():
    parser = argparse.ArgumentParser(description='Call consensus nested SV from SV VCF (supports DUP-INV, DUP-DEL, INV-INS, INV-DEL, INV-DUP)')
    parser.add_argument('-i', '--input', required=True, help='Input SV VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output nested SV file')
    parser.add_argument('-g', '--gap', type=int, default=5000, help='Cluster gap tolerance (default: 5000)')
    args = parser.parse_args()

    # Read VCF
    sv_list = []
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chrom = parts[0]
            info = parts[7]

            # Parse INFO
            svtype = re.search(r'SVTYPE=([A-Z]+)', info).group(1)
            sv_start = int(re.search(r'SV_START=(\d+)', info).group(1))
            sv_end = int(re.search(r'END=(\d+)', info).group(1))

            sv_list.append({
                'chr': chrom,
                'start': sv_start,
                'end': sv_end,
                'type': svtype
            })

    if not sv_list:
        print("No SVs found.")
        return

    # Sort
    sv_list.sort(key=lambda x: (x['chr'], x['start']))

    # Cluster
    clusters = []
    current = [sv_list[0]]
    for sv in sv_list[1:]:
        last = current[-1]
        if sv['chr'] == last['chr'] and sv['start'] <= last['end'] + args.gap:
            current.append(sv)
        else:
            clusters.append(current)
            current = [sv]
    clusters.append(current)

    # Allowed nested pairs
    allowed = {
        ('DUP', 'INV'), ('DUP', 'DEL'),
        ('INV', 'INS'), ('INV', 'DEL'), ('INV', 'DUP')
    }

    results = []

    # Process each cluster
    for cl in clusters:
        pairs = []
        for a in cl:
            for b in cl:
                if a is b:
                    continue
                if (a['chr'] == b['chr'] and
                    a['start'] <= b['start'] and
                    a['end'] >= b['end']):
                    if (a['type'], b['type']) in allowed:
                        pairs.append((a, b))
        if not pairs:
            continue

        # Consensus outer
        o_type = pairs[0][0]['type']
        o_starts = [p[0]['start'] for p in pairs]
        o_ends = [p[0]['end'] for p in pairs]
        o_s = int(statistics.median(o_starts))
        o_e = int(statistics.median(o_ends))

        # Consensus inner
        i_type = pairs[0][1]['type']
        i_starts = [p[1]['start'] for p in pairs]
        i_ends = [p[1]['end'] for p in pairs]
        i_s = int(statistics.median(i_starts))
        i_e = int(statistics.median(i_ends))

        chr_name = pairs[0][0]['chr']
        results.append((
            chr_name, o_s, o_e, o_type,
            i_s, i_e, i_type, 'NestedSV'
        ))

    # Write output
    with open(args.output, 'w') as f:
        for line in results:
            f.write('\t'.join(map(str, line)) + '\n')

    print(f"Done! Found {len(results)} nested SV events. Output: {args.output}")

if __name__ == '__main__':
    main()
