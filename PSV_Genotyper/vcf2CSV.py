import argparse
import statistics
import re

def parse_vcf_info(info_str):
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, val = entry.split('=', 1)
            info[key] = val
    return info

def clean_end_position(val):
    if ':' in val:
        val = val.split(':')[-1]
    return int(val)

def main():
    parser = argparse.ArgumentParser(description='Call nested SV from sorted/unsorted VCF (auto-sort enabled)')
    parser.add_argument('-i', '--input', required=True, help='Input SV VCF')
    parser.add_argument('-o', '--output', required=True, help='Output nested SV table')
    parser.add_argument('-g', '--gap', type=int, default=5000, help='Cluster gap tolerance')
    # 修复：这里必须加括号！
    args = parser.parse_args()

    sv_list = []

    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            info_field = parts[7]
            info = parse_vcf_info(info_field)

            svtype = info.get("SVTYPE")
            if not svtype:
                continue

            try:
                end_val = info["END"]
                end = clean_end_position(end_val)
            except:
                end = pos

            if "SV_START" in info:
                start = int(info["SV_START"])
            else:
                start = pos

            sv_list.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "type": svtype
            })

    if not sv_list:
        print("No SVs found.")
        return

    # 自动排序：染色体数字排序 + 坐标排序
    def chr_sort_key(chr_str):
        chr_str = str(chr_str).replace('chr', '')
        if chr_str == "X": return 100
        if chr_str == "Y": return 101
        if chr_str == "M" or chr_str == "MT": return 102
        try:
            return int(chr_str)
        except:
            return 999

    sv_list.sort(key=lambda x: (chr_sort_key(x["chr"]), x["start"]))

    # 聚类
    clusters = []
    current = [sv_list[0]]
    for sv in sv_list[1:]:
        last = current[-1]
        if sv["chr"] == last["chr"] and sv["start"] <= last["end"] + args.gap:
            current.append(sv)
        else:
            clusters.append(current)
            current = [sv]
    clusters.append(current)

    allowed = {
        ("DUP", "INV"), ("DUP", "DEL"),
        ("INV", "INS"), ("INV", "DEL"), ("INV", "DUP")
    }

    results = []
    for cl in clusters:
        pairs = []
        for a in cl:
            for b in cl:
                if a is b:
                    continue
                if (a["chr"] == b["chr"] and
                    a["start"] <= b["start"] and
                    a["end"] >= b["end"]):
                    if (a["type"], b["type"]) in allowed:
                        pairs.append((a, b))
        if not pairs:
            continue

        o_s = int(statistics.median([p[0]["start"] for p in pairs]))
        o_e = int(statistics.median([p[0]["end"] for p in pairs]))
        o_t = pairs[0][0]["type"]

        i_s = int(statistics.median([p[1]["start"] for p in pairs]))
        i_e = int(statistics.median([p[1]["end"] for p in pairs]))
        i_t = pairs[0][1]["type"]

        chr_name = pairs[0][0]["chr"]
        results.append((chr_name, o_s, o_e, o_t, i_s, i_e, i_t, "NestedSV"))

    with open(args.output, 'w') as f:
        for row in results:
            f.write('\t'.join(map(str, row)) + '\n')

    print(f"✅ Done! Found {len(results)} nested SVs. Output: {args.output}")

if __name__ == '__main__':
    main()
