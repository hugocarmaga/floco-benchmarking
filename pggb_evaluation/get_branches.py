#!/usr/bin/env python3

import sys
import itertools
from collections import defaultdict


class Edge:
    def __init__(self, name1, side1, name2, side2):
        # assert name1 != name2 or side1 != side2
        if (name1, side1) > (name2, side2):
            name1, name2 = name2, name1
            side1, side2 = side2, side1
        self.name1 = name1
        self.name2 = name2
        self.side1 = side1
        self.side2 = side2

    def __eq__(self, oth):
        return self.name1 == oth.name1 \
            and self.name2 == oth.name2 \
            and self.side1 == oth.side1 \
            and self.side2 == oth.side2

    def __hash__(self):
        return hash((self.name1, self.name2, self.side1, self.side2))

    def is_loop(self):
        return self.name1 == self.name2

    def __str__(self):
        return f'E({self.name1}{"-+"[self.side1]} {self.name2}{"+-"[self.side2]})'

    def __repr__(self):
        return self.__str__()

    def oth_side(self, name, side):
        if name == self.name1 and side == self.side1:
            return self.name2, self.side2
        assert name == self.name2 and side == self.side2
        return self.name1, self.side1


class Node:
    def __init__(self):
        # Left and right edges
        self.edges = [set(), set()]
        self.seq = None

    def set_seq(self, seq):
        self.seq = seq


class Graph:
    def __init__(self):
        self.nodes = defaultdict(Node)
        # Store lengths separately from the nodes as some node entities will be removed.
        self.node_lengths = dict()
        # From edges to sets of collapsed nodes.
        self.collapsed_nodes = defaultdict(set)

    @classmethod
    def load(cls, file):
        this = cls()
        for line in file:
            split = line.strip().split('\t')
            match split[0]:
                case 'S':
                    name = split[1]
                    seq = split[2]
                    seq_len = len(seq)
                    this.nodes[name].set_seq(seq)
                    this.node_lengths[name] = seq_len
                case 'L':
                    assert split[5] == '0M', f'Overlap nodes are not supported ({split[1]} to {split[3]})'
                    # left = false, right = true
                    this.add_edge(Edge(split[1], split[2] == '+', split[3], split[4] == '-'))
                case _:
                    pass
        return this

    def write(self, file):
        all_edges = set()
        for name, node in self.nodes.items():
            file.write(f'S\t{name}\t{node.seq}\n')
            all_edges.update(node.edges[0])
            all_edges.update(node.edges[1])
        for edge in all_edges:
            file.write(f'L\t{edge.name1}\t{"-+"[edge.side1]}\t{edge.name2}\t{"+-"[edge.side2]}\t0M\n')

    def sum_len(self, nodes):
        return sum(self.node_lengths[name] for name in nodes)

    def __len__(self):
        return len(self.nodes)

    def add_edge(self, edge):
        self.nodes[edge.name1].edges[edge.side1].add(edge)
        self.nodes[edge.name2].edges[edge.side2].add(edge)

    def remove_edge(self, edge):
        self.nodes[edge.name1].edges[edge.side1].remove(edge)
        if edge.name1 != edge.name2 or edge.side1 != edge.side2:
            self.nodes[edge.name2].edges[edge.side2].remove(edge)

    def add_transitive_edge(self, left_edge, middle_node, right_edge):
        """
        Add an edge between nodes to the left of `middle_node` and to the right of it.
        """
        if left_edge.is_loop() or right_edge.is_loop():
            return

        if left_edge.name2 == middle_node:
            name1 = left_edge.name1
            side1 = left_edge.side1
        else:
            name1 = left_edge.name2
            side1 = left_edge.side2

        if right_edge.name1 == middle_node:
            name2 = right_edge.name2
            side2 = right_edge.side2
        else:
            name2 = right_edge.name1
            side2 = right_edge.side1
        edge = Edge(name1, side1, name2, side2)
        self.add_edge(edge)
        return edge

    def remove_node(self, name):
        node = self.nodes[name]
        left_edges = list(node.edges[0])
        right_edges = list(node.edges[1])

        for left_edge in left_edges:
            for right_edge in right_edges:
                edge = self.add_transitive_edge(left_edge, name, right_edge)
                collapsed_nodes = self.collapsed_nodes[edge]
                collapsed_nodes.add(name)
                if left_edge in self.collapsed_nodes:
                    collapsed_nodes.update(self.collapsed_nodes[left_edge])
                if right_edge in self.collapsed_nodes:
                    collapsed_nodes.update(self.collapsed_nodes[right_edge])

        for left_edge in left_edges:
            self.remove_edge(left_edge)
        for right_edge in right_edges:
            if not right_edge.is_loop() or right_edge.side1 == right_edge.side2:
                self.remove_edge(right_edge)

        assert not node.edges[0] and not node.edges[1]
        del self.nodes[name]

    def process_snarls(self, f, max_collapsed_length):
        for line in f:
            split = line.strip().split('\t')
            collapsed_nodes = split[3].split(',')
            collapsed_length = self.sum_len(collapsed_nodes)
            if collapsed_length > max_collapsed_length:
                continue
            for name in collapsed_nodes:
                self.remove_node(name)

    def _extend_branch(self, name, seen, forbidden_edges):
        node = self.nodes[name]
        seen.add(name)
        degree_left = sum(edge not in forbidden_edges for edge in node.edges[0])
        degree_right = sum(edge not in forbidden_edges for edge in node.edges[1])
        # Second element = direction, true if right, false if left.
        go_right = degree_left == 0 or degree_right > 0
        core_nodes = [(name, go_right)]
        collapsed_nodes = set()

        while True:
            last_name, last_direction = core_nodes[-1]
            last_node = self.nodes[last_name]
            possible_edges = []
            for edge in last_node.edges[last_direction]:
                if edge in forbidden_edges:
                    continue
                if edge.oth_side(last_name, last_direction)[0] not in seen:
                    possible_edges.append(edge)

            if not possible_edges:
                break
            assert len(possible_edges) == 1
            edge = possible_edges[0]

            name, side = edge.oth_side(last_name, last_direction)
            core_nodes.append((name, not side))
            seen.add(name)
            if edge in self.collapsed_nodes:
                collapsed_nodes.update(self.collapsed_nodes[edge])
        n_pos = sum(direction for _, direction in core_nodes)
        n_neg = len(core_nodes) - n_pos
        if n_neg > n_pos:
            core_nodes = [(name, not direction) for name, direction in core_nodes[::-1]]
        return core_nodes, collapsed_nodes

    def find_branches(self):
        forbidden_edges = set()
        start_list = []
        for name, node in self.nodes.items():
            for side_edges in node.edges:
                if len(side_edges) > 1:
                    forbidden_edges.update(side_edges)
                    start_list.append(name)

        seen = set()
        for name in itertools.chain(start_list, self.nodes.keys()):
            if name not in seen:
                yield self._extend_branch(name, seen, forbidden_edges)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--graph', metavar='FILE', required=True,
        help='Input graph in GFA format.')
    parser.add_argument('-s', '--snarls', metavar='FILE', required=True,
        help='Snarls information, obtained using `vg stats --snarl-contents`.')
    parser.add_argument('-c', '--collapse-len', metavar='INT', type=int, default=300,
        help='Collapse bubbles with sum node size at most INT [%(default)s].')
    parser.add_argument('-l', '--branch-len', metavar='INT', type=int, default=1000,
        help='Output branches with core length at least INT [%(default)s].')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output branches in CSV format.')
    parser.add_argument('-O', '--out-graph', metavar='FILE', required=False,
        help='Optionally: output simplified graph to this GFA file.')
    args = parser.parse_args()

    with open(args.graph) as f:
        graph = Graph.load(f)
    with open(args.snarls) as f:
        graph.process_snarls(f, args.collapse_len)

    with open(args.output, 'w') as out:
        out.write('name\tcore_nodes\tcore_length\tcollapsed_nodes\tcollapsed_length\tcore\tcollapsed\n')
        branches = 0
        total_core_nodes = 0
        total_core_length = 0
        total_nodes = 0
        total_length = 0

        for core_nodes, collapsed_nodes in graph.find_branches():
            core_len = graph.sum_len(node for node, _ in core_nodes)
            if core_len < args.branch_len:
                continue

            assert core_nodes
            collapsed_len = graph.sum_len(collapsed_nodes)
            out.write(f'branch{branches+1}\t{len(core_nodes)}\t{core_len}\t{len(collapsed_nodes)}\t{collapsed_len}\t')
            core_str = ''.join(f'{"<>"[direction]}{name}' for name, direction in core_nodes)
            out.write(f'{core_str}\t{",".join(collapsed_nodes) if collapsed_nodes else "*"}\n')

            branches += 1
            total_core_nodes += len(core_nodes)
            total_core_length += core_len
            total_nodes += len(core_nodes) + len(collapsed_nodes)
            total_length += core_len + collapsed_len

    graph_size = len(graph.node_lengths)
    graph_length = sum(graph.node_lengths.values())
    sys.stderr.write(f'Graph:    {graph_size:7,} nodes,         total length {graph_length:9,} bp\n')
    sys.stderr.write(f'Identified {branches} branches:\n')
    sys.stderr.write('    Core: {:7,} ({:4.1f}%) nodes, total length {:9,} bp ({:4.1f}%)\n'.format(
        total_core_nodes, 100 * total_core_nodes / graph_size,
        total_core_length, 100 * total_core_length / graph_length))
    sys.stderr.write('    All:  {:7,} ({:4.1f}%) nodes, total length {:9,} bp ({:4.1f}%)\n'.format(
        total_nodes, 100 * total_nodes / graph_size,
        total_length, 100 * total_length / graph_length))

    if args.out_graph:
        with open(args.out_graph, 'w') as out:
            graph.write(out)


if __name__ == '__main__':
    main()
