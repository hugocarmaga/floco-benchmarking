# floco-benchmarking

Collection of code used to benchmark Floco

## Processing PGGB graphs

To get graph branches, run
```sh
xz -d graph.gfa.xz
vg stats --snarl-contents -P pggb.gfa > snarls.csv
~/Code/scripts/floco/get_branches.py -g pggb.gfa -s snarls.csv -o branches.csv
```

Scripts `branch_seq.py` allows to extract branch sequence from the graph, `branch_cn` to evaluate branch copy number from a CN predictions for nodes.
Use
```sh
minimap2 -cx asm20 --eqx -N 50 haplotypes/${sample}.fa.gz branches.fa.gz > branches/${sample}.paf
~/Code/scripts/floco/paf_to_chain_cn.py -p branches/${sample}.paf -b branches.csv \
    -o branches/${sample}.csv -O branches/${sample}.bed
```
to obtain branch copy numbers and local coordinates for branches for a given haplotype.