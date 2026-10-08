#!/usr/bin/env python3

import re
from get_branches import Graph


def load_branches(f):
    branches = []
    next(f)
    for line in f:
        split = line.strip().split('\t')
        name = split[0]
        core_nodes = re.split('[<>]', split[5][1:])
        coll_nodes = split[6].split(',')
        branches.append((name, core_nodes, coll_nodes))
    return branches


def load_cn(f):
    cn = {}
    delim = None
    for line in f:
        if line.startswith('#'):
            continue
        delim = '\t' if '\t' in line else ','
        columns = line.strip().lower().split(delim)
        node_col = columns.index('node')
        cn_col = columns.index('cn') if 'cn' in columns else columns.index('copy_number')
        if 'sum_coverage' in columns:
            cov_col = columns.index('sum_coverage')
        else:
            cov_col = None
        break
    cov = None if cov_col is None else {}

    for line in f:
        split = line.strip().split(delim)
        node = split[node_col]
        cn[node] = int(split[cn_col])
        if cov_col is not None:
            cov[node] = float(split[cov_col])
    return cn, cov


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--graph', metavar='FILE', required=True,
        help='Input graph in the GFA format.')
    parser.add_argument('-c', '--cns', metavar='FILE', required=True,
        help='Input CSV file with columns "node" and "copy_number".')
    parser.add_argument('-b', '--branches', metavar='FILE', required=True,
        help='Input CSV file with branches.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-m', '--mean-cov', metavar='NUM', type=float,
        help='Optional: mean sample read depth.')
    args = parser.parse_args()

    with open(args.graph) as f:
        graph = Graph.load(f)
    with open(args.branches) as f:
        branches = load_branches(f)
    with open(args.cns) as f:
        cn, cov = load_cn(f)

    with open(args.output, 'w') as out:
        out.write(f'branch\tcore_cn')
        n_stars = 1
        if cov is not None:
            out.write('\tcore_cov')
            n_stars += 1
            if args.mean_cov is not None:
                out.write('\tnorm_core_cov')
                n_stars += 1
        out.write('\ttotal_length\n')

        full_length = sum(node_cn * graph.node_lengths[name] for name, node_cn in cn.items())
        out.write(f'graph{"\t*" * n_stars}\t{full_length}\n')

        for branch_name, core_nodes, coll_nodes in branches:
            sum_len = 0
            sum_cov = 0
            baseline_len = 0
            for name in core_nodes:
                length = graph.node_lengths[name]
                baseline_len += length
                sum_len += length * cn.get(name, 0)
                if cov is not None:
                    sum_cov += cov.get(name, 0)

            cn_str = f'{sum_len / baseline_len:.2f}'.rstrip('0').rstrip('.')
            out.write(f'{branch_name}\t{cn_str}')
            if cov is not None:
                aver_cov = sum_cov / baseline_len
                out.write(f'\t{aver_cov:.5f}')
                if args.mean_cov is not None:
                    out.write(f'\t{aver_cov / args.mean_cov:.5f}')

            for name in coll_nodes:
                if name != '*':
                    sum_len += graph.node_lengths[name] * cn.get(name, 0)
            out.write(f'\t{sum_len}\n')


if __name__ == '__main__':
    main()
