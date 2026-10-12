# Script to extract haplotype length from CNVnator output files

def read_assembly_fasta(fasta_file):
    """
    Reads a FASTA file and returns the total length of the locus for the sample.
    """
    total_length = 0
    # Check whether file is gzipped or not
    if fasta_file.endswith('.gz'):
        import gzip
        with gzip.open(fasta_file, 'rt') as f:
            for line in f:
                if not line.startswith('>'):
                    total_length += len(line.strip())
    else:
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    total_length += len(line.strip())
    return total_length

def read_bed_file(bed_file):
    """
    Reads a BED file and returns a list of tuples containing the chromosome, start, and end positions.
    """
    bed_data = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                bed_data.append((chrom, start, end))
    return bed_data


def read_cnvnator_output(cnvnator_file, bed_data):
    """
    Reads a CNVnator output file and returns a list of tuples containing the chromosome, start, end, and CN value, for CNVs overlapping with the BED regions.
    """
    cnv_data = []
    with open(cnvnator_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[1].split(':')[0]
            start, end = map(int, parts[1].split(':')[1].split('-'))
            start -= 1
            if any(bed_chrom == chrom and not (bed_end <= start or bed_start >= end) for bed_chrom, bed_start, bed_end in bed_data):
                cn = round(float(parts[3])*2)  # Convert to copy number
                cnv_data.append((chrom, start, end, cn))
    return cnv_data

def compute_cnvnator_haplotype_length(cnv_data, bed_data):
    """
    Computes haplotype length, assuming a baseline copy number of 2.
    CNVs are assumed not to overlap.
    """
    total_length = 0

    for bed_chrom, bed_start, bed_end in bed_data:
        total_length += 2 * (bed_end - bed_start)

        for cnv_chrom, cnv_start, cnv_end, cn in cnv_data:
            if cnv_chrom != bed_chrom:
                continue

            overlap_start = max(bed_start, cnv_start)
            overlap_end = min(bed_end, cnv_end)

            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                total_length += overlap_length * (cn - 2)

    return total_length

def compute_floco_haplotype_length(csv_file):
    """
    Computes the haplotype length from a CSV file containing node names, their lengths and copy numbers.

    Args:
        csv_file (str): Path to the CSV file containing haplotype data.

    Returns:
        int: The total haplotype length.
    """
    total_length = 0
    with open(csv_file, 'r') as f:
        for line in f:
            if line.startswith('Node'):
                continue
            parts = line.strip().split(',')
            total_length += int(parts[1]) * int(parts[3])  # length * copy number
    return total_length

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Compute haplotype length from CNVnator output and BED file.")
    parser.add_argument('--bed', required=True, help="Path to the BED file.")
    parser.add_argument('--cnvnator', required=True, help="Path to the CNVnator output file.")
    parser.add_argument('--floco', required=True, help="Path to the Floco CSV file.")
    parser.add_argument('--assembly', required=True, help="Path to the FASTA file with the assembled haplotypes.")

    args = parser.parse_args()

    bed_data = read_bed_file(args.bed)
    cnv_data = read_cnvnator_output(args.cnvnator, bed_data)
    assembly_length = read_assembly_fasta(args.assembly)
    cnvnator_length = compute_cnvnator_haplotype_length(cnv_data, bed_data)
    floco_length = compute_floco_haplotype_length(args.floco)

    print(f"Total haplotype length from Assembly: {assembly_length}")
    print(f"Total haplotype length from CNVnator (% of assembly): {cnvnator_length} ({cnvnator_length / assembly_length * 100:.2f}%)")
    print(f"Total haplotype length from Floco (% of assembly): {floco_length} ({floco_length / assembly_length * 100:.2f}%)")

if __name__ == "__main__":
    main()
