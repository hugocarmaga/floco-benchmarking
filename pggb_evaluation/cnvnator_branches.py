#!/usr/bin/env python3
'''
Script to get CNVnator CN predictions for graph arms.
'''

from collections import defaultdict

def read_chains_bed(chains_bed_file):
    """
    Reads a chains BED file and returns a dictionary with BED record names as keys and tuples of chromosome, start, and end positions as values.
    """
    chains_data = {}
    with open(chains_bed_file, 'r') as f:
        last_chrom = None
        last_start = 0
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            assert last_chrom != chrom or start >= last_end
            last_chrom = chrom
            last_end = end
            if parts[3] not in chains_data:
                chains_data[parts[3]] = [end-start, (chrom, start, end)]
            else:
                chains_data[parts[3]][0] += end-start
                chains_data[parts[3]].append((chrom, start, end))
    return chains_data

def read_cnvnator_output(cnvnator_file, bed_data):
    """
    Reads a CNVnator output file and returns a list of tuples containing the chromosome, start, end, and raw CNVnator depth ratio for CNVs overlapping with the BED regions.
    """
    cnv_data = []
    with open(cnvnator_file, 'r') as f:
        last_chrom = None
        last_end = 0
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[1].split(':')[0]
            start, end = map(int, parts[1].split(':')[1].split('-'))
            start -= 1
            assert last_chrom != chrom or start >= last_end
            last_chrom = chrom
            last_end = end
            if any(bed_chrom == chrom and not (bed_end <= start or bed_start >= end) for bed_chrom, bed_start, bed_end in bed_data):
                ratio = float(parts[3])
                cnv_data.append((chrom, start, end, ratio))
    return cnv_data

def compute_cnvnator_cn(chains_bed_file, cnvnator_file):
    """
    Computes the CNVnator copy number predictions for the graph arms defined in the chains BED file.
    Returns a dictionary with BED record names as keys and their corresponding CNVnator copy number predictions as values.
    """
    chains_data = read_chains_bed(chains_bed_file)
    cnv_predictions = {}

    for bed_name, bed_info in chains_data.items():
        total_length = bed_info[0]
        bed_regions = bed_info[1:]
        cnv_data = read_cnvnator_output(cnvnator_file, bed_regions)
        covered_length = 0
        weighted_cn = 0
        for chrom1, start1, end1, ratio in cnv_data:
            for chrom2, start2, end2 in bed_regions:
                if chrom1 != chrom2:
                    continue
                l = max(0, min(end1, end2) - max(start1, start2))
                covered_length += l
                weighted_cn += l * round(ratio * 2)
        remaining_length = total_length - covered_length
        weighted_cn += 2 * remaining_length
        chain_cn = weighted_cn / total_length * len(bed_regions)
        cnv_predictions[bed_name] = (chain_cn, weighted_cn)
    return cnv_predictions

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Compute CNVnator copy number predictions for graph arms.')
    parser.add_argument('-b', '--bed', required=True, help='Path to the chains BED file.')
    parser.add_argument('-c', '--cnvnator', required=True, help='Path to the CNVnator output file.')
    args = parser.parse_args()

    cnv_predictions = compute_cnvnator_cn(args.bed, args.cnvnator)

    print('chain\tcn\tlength')
    for bed_name, (cn, length) in cnv_predictions.items():
        cn_str = f'{cn:.5f}'.rstrip('0').rstrip('.')
        print(f"{bed_name}\t{cn_str}\t{length}")

if __name__ == "__main__":
    main()
