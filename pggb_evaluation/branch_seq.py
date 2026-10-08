#!/usr/bin/env python3

from collections import defaultdict
from branch_cn import load_branches
import re


def load_graph(f, prefix, exclude):
    exclude = re.compile(exclude) if exclude else None
    node_seqs = dict()
    paths = []
    for line in f:
        split = line.strip().split('\t')
        match split[0]:
            case 'S':
                node_seqs[split[1]] = split[2]
            case 'P':
                hap_name = split[1]
                if hap_name.startswith(prefix) and (exclude is None or not exclude.match(hap_name)):
                    paths.append((split[1], split[2].split(',')))
            case _:
                pass
    assert paths, f'Could not find paths for prefix `{prefix}`'
    return node_seqs, paths


_repl_nt = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }

def rev_comp(seq):
    return ''.join(_repl_nt[nt] for nt in seq[::-1])


def write_branch(prefix, branch_name, branch_occ, seq, out):
    n = branch_occ[branch_name] + 1
    suffix = f'.{n}'
    out.write(f'>{prefix}{branch_name}{suffix}\n{seq}\n')
    branch_occ[branch_name] += 1


def process_path(hap_name, path, node_seqs, node_branch, branch_occ, out):
    prefix = hap_name + '-' if hap_name else ''
    curr_branch = None
    curr_seq = None
    curr_tail = ''
    for node in path:
        strand = node[-1]
        node = node[:-1]
        node_seq = node_seqs[node]
        if strand == '-':
            node_seq = rev_comp(node_seq)

        branch = node_branch.get(node)
        if branch is None:
            curr_tail += node_seq
            continue

        if branch != curr_branch:
            if curr_branch is not None:
                write_branch(prefix, curr_branch, branch_occ, curr_seq, out)
            curr_branch = branch
            curr_seq = node_seq
        else:
            curr_seq += curr_tail + node_seq
        curr_tail = ''

    if curr_branch is not None:
        write_branch(prefix, curr_branch, branch_occ, curr_seq, out)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--graph', metavar='FILE', required=True,
        help='Input GFA graph.')
    parser.add_argument('-b', '--branches', metavar='FILE', required=True,
        help='Input CSV file with branch definitions.')
    parser.add_argument('-p', '--prefix', metavar='STR',
        help='Haplotype prefix.')
    parser.add_argument('-e', '--exclude', metavar='STR',
        help='Optional: regex for removing haplotypes that matched the prefix.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output FASTA file.')
    parser.add_argument('--hap-name', action='store_true',
        help='In the output file, write full haplotype name following by the branch name. '
            'If not specified, prefix and branch names are written.')
    parser.add_argument('--hap-suffix', action='store_true',
        help='When using --hap-name, use separate suffix counters for each haplotype name.')
    args = parser.parse_args()

    with open(args.graph) as f:
        node_seqs, paths = load_graph(f, args.prefix, args.exclude)

    with open(args.branches) as f:
        branches = load_branches(f)
    node_branch = {}
    for branch_name, core_nodes, _ in branches:
        for node in core_nodes:
            assert node not in node_branch
            node_branch[node] = branch_name

    with open(args.output, 'w') as out:
        branch_occ = defaultdict(int)
        for hap_name, path in paths:
            if args.hap_name and args.hap_suffix:
                branch_occ.clear()
            process_path(hap_name if args.hap_name else args.prefix, path, node_seqs, node_branch, branch_occ, out)


if __name__ == '__main__':
    main()
