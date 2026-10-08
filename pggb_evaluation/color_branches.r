#!/usr/bin/env Rscript

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

argv <- commandArgs(trailingOnly = T)

gfa_filename <- argv[1]
branches_filename <- argv[2]
out_filename <- argv[3]

graph <- readLines(gfa_filename)
node_lines <- graph[grepl('^S', graph)]
all_node_names <- str_split_i(node_lines, '\t', 2)

branches <- read.csv(branches_filename, sep = '\t')
core_nodes <- select(branches, name, core) |>
    mutate(
        core = substr(core, 2, str_length(core)),
        core = gsub('<', '>', core),
    ) |>
    separate_rows(core, sep = '>') |>
    rename(branch = name, node = core) |>
    group_by(branch) |>
    mutate(label = c(rep('""', max(0, n() %/% 2 - 1)), branch[1], rep('""', max(0, n() - max(0, n() %/% 2 - 1) - 1))))
coll_nodes <- select(branches, name, collapsed) |>
    separate_rows(collapsed, sep = ',') |>
    rename(branch = name, node = collapsed) |>
    mutate(label = '""')
all_nodes <- rbind(core_nodes, coll_nodes)
all_nodes <- rbind(all_nodes,
    data.frame(node = setdiff(all_node_names, all_nodes$node)) |>
        mutate(branch = 'none', label = '""'))

n_branches <- nrow(branches)
colors <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D",
    "#B6992D", "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#D37295",
    "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6", '#000000')
# branch_colors <- ltc::ltc('hat', 10)
# branch_colors <- branch_colors[c(seq(1, 10, by = 2), seq(2, 10, by = 2))]
n_colors <- length(colors)
branch_colors <- lapply(1:ceiling(n_branches / n_colors), function(i) sample(colors))
branch_colors <- do.call(c, branch_colors) |>
    setNames(branches$name) |>
    c('none' = '#999999')
all_nodes <- mutate(all_nodes, color = branch_colors[branch]) |>
    arrange(branch, node) |>
    select(node, branch, label, color) |>
    write.table(out_filename, sep = '\t', row.names = F, quote = F)
