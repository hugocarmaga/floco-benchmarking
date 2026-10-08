#!/usr/bin/env python3

import sys
from collections import namedtuple, defaultdict
from intervaltree import IntervalTree

class Aln(namedtuple('Aln', 'contig start end branch div')):
    def __len__(self):
        return self.end - self.start


def load_paf(f, min_len):
    alns = []
    for line in f:
        split = line.strip().split('\t')
        contig = split[5]
        start = int(split[7])
        end = int(split[8])
        if end - start < min_len:
            continue
        branch = split[0].rsplit('.', 1)[0]
        div = 1 - int(split[9]) / int(split[10])
        alns.append(Aln(contig, start, end, branch, div))
    return alns


def process_alignments(alns, min_len, min_overlap_frac, out_bed):
    trees = defaultdict(IntervalTree)
    alns.sort(key=lambda aln: aln.div)
    warnings = set()
    def warn(msg):
        if msg not in warnings:
            warnings.add(msg)
            sys.stderr.write(msg + '\n')

    for aln in alns:
        tree = trees[aln.contig]
        new_start = aln.start
        new_end = aln.end

        remove = []
        for prev in tree.overlap(aln.start, aln.end):
            if aln.branch == prev.data:
                overlap_size = min(aln.end, prev.end) - max(aln.start, prev.begin)
                overlap_frac = overlap_size / min(len(aln), len(prev))
                if overlap_frac >= min_overlap_frac:
                    remove.append(prev)
                    new_start = min(new_start, prev.begin)
                    new_end = max(new_start, prev.end)
                else:
                    warn(f'Regions {prev.begin}..{prev.end} and {aln.start}..{aln.end} ({aln.branch}) '
                        f'overlap only by {overlap_frac:.4f}')
        for prev in remove:
            tree.remove(prev)

        for prev in tree.overlap(new_start, new_end):
            if new_start <= prev.begin and prev.end <= new_end:
                warn(f'Region {prev.begin}..{prev.end} ({prev.data}) is enveloped '
                    f'by another {new_start}..{new_end} ({aln.branch})')
                exit(1)
            elif prev.begin < new_start:
                new_start = min(new_end, prev.end)
            elif new_end < prev.end:
                new_end = max(new_start, prev.begin)

        if new_end - new_start < min_len:
            reason = 'fully covered by other branches' if new_start == new_end \
                else f'remainder is too short ({new_end - new_start} bp)'
            warn(f'Skipping {new_start}..{new_end} ({aln.branch}): {reason}')
        else:
            tree.addi(new_start, new_end, aln.branch)

    cn = defaultdict(lambda: [0, 0])
    for contig, tree in trees.items():
        for interval in sorted(tree):
            cn[interval.data][0] += 1
            cn[interval.data][1] += interval.end - interval.begin
            out_bed.write(f'{contig}\t{interval.begin}\t{interval.end}\t{interval.data}\n')
    return cn


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--paf', metavar='FILE', required=True,
        help='Input PAF file.')
    parser.add_argument('-b', '--branches', metavar='FILE', required=True,
        help='CSV file with branches definition.')
    parser.add_argument('-l', '--min-len', metavar='INT', type=int, default=500,
        help='Minimal alignment/output region length [%(default)s].')
    parser.add_argument('-f', '--overlap', metavar='FLOAT', type=float, default=0.9,
        help='Overlap fraction [%(default)s].')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-O', '--out-bed', metavar='FILE', required=False,
        help='Output BED file.')
    args = parser.parse_args()

    with open(args.paf) as f:
        alns = load_paf(f, args.min_len)

    out_bed = open(args.out_bed, 'w') if args.out_bed is not None else None
    cn = process_alignments(alns, args.min_len, args.overlap, out_bed)
    with open(args.branches) as f, open(args.output, 'w') as out:
        out.write('branch\tcn\tlength\n')
        next(f)
        for line in f:
            split = line.strip().split('\t')
            branch = split[0]
            branch_cn, branch_len = cn[branch]
            out.write(f'{branch}\t{branch_cn}\t{branch_len}\n')


if __name__ == '__main__':
    main()
