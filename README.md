# floco-benchmarking

Collection of code used to benchmark [Floco](https://github.com/hugocarmaga/floco).

## Making haplotype FASTA files compatible with PGGB

To transform FASTA headers into PSN notation, run:

```sh
zcat haplotypes.fa.gz | awk 'BEGIN{OFS="#";}$0 ~/^>/{split($1,s,"."); print s[1],s[2],$2}$0 !~/^>/{print}' | awk 'BEGIN{FS=OFS="#";}$0 ~/^>/{if ($2 ~/-/) {split($2,id,"-"); split($3,ctg,"-"); print $1,id[1],ctg[1]"-"id[2]"-"ctg[2]} else if ($1 ~/CHM13/ || $1 ~/GRCh38/){print $1,$4,$5} else print} $0 !~/^>/ {print}' | bgzip > haplotypes.PSN.fa.gz
```
## Processing PGGB graphs

In order to get path-derived CN values for SAMPLE from a graph, run:
```sh
awk 'BEGIN{FS="\t";OFS=","; print "Node","Copy_number"}$1=="P" && $2 ~/SAMPLE/{split($3,nodes,","); for (n in nodes){gsub(/[+-]/,"",nodes[n]); cn[nodes[n]]++}}END{for (n in cn){print n,cn[n]}}' graph.gfa | sort -n > copy_numbers-SAMPLE-from-paths.csv
```

To get graph branches, run:
```sh
xz -d graph.gfa.xz
vg stats --snarl-contents -P pggb.gfa > snarls.csv
~/Code/scripts/floco/get_branches.py -g pggb.gfa -s snarls.csv -o branches.csv
```

Scripts `branch_seq.py` allows to extract branch sequence from the graph, `branch_cn` to evaluate branch copy number from a CN predictions for nodes.
Use:
```sh
minimap2 -cx asm20 --eqx -N 50 haplotypes/${sample}.fa.gz branches.fa.gz > branches/${sample}.paf
~/Code/scripts/floco/paf_to_branch_cn.py -p branches/${sample}.paf -b branches.csv \
    -o branches/${sample}.csv -O branches/${sample}.bed
```
to obtain branch copy numbers and local coordinates for branches for a given haplotype.

## Comparing Floco with CNVnator

To compare the total length of the haplotype resulting from the CN predictions from both [Floco](https://github.com/hugocarmaga/floco) and [CNVnator](https://github.com/abyzovlab/CNVnator), using the assembled haplotype as the ground truth, we run:
```sh
python pggb_evaluation/cnvnator_total_len_comparison.py --bed locus-coordinates.bed --assembly assembled_haplotype.fa --cnvnator cnvnator-results.txt --floco floco-results.csv
```
which outputs something like this:
```
Total haplotype length from Assembly: 1717437
Total haplotype length from CNVnator (% of assembly): 1673015 (97.41%)
Total haplotype length from Floco (% of assembly): 1698134 (98.88%)
```

To get CNVnator's prediction for graph branches projected into CHM13-based linear space, we run:
```sh
python pggb_evaluation/cnvnator_branches.py --bed branches.bed --cnvnator cnvnator-results.txt
```
which produces an output file like this:
```
chain   cn          length
chain1  2           4474
chain2  2           141432
chain3  2           53580
chain4  2           23762
chain5  0.65631     1606
chain6  2           26192
chain7  2           56524
chain8  2           172886
chain9  2           127274
chain10 9.93493     34759
chain11 11.20158    58126
chain12 12.99111    34945
```