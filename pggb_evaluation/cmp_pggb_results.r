library(tidyverse)
library(ggplot2)

cowplot::set_null_device('agg')
flanks <- read.csv('data/flanks.csv', sep = '\t')
samples <- readLines('data/samples.txt')
sample_class <- setNames(rep(c('in graph', 'out of graph'), each = 5), samples)

colors <- c('CNVnator' = '#ACC8BE', 'Floco' = '#1A5B5B')

all_lengths <- read.csv('data/length.csv', sep = '\t') |>
    filter(locus != 'AMY' | sample != 'HG01114')
all_lengths2 <- all_lengths |>
    group_by(locus, sample) |>
    mutate(true_length = length[method == 'assembly']) |>
    ungroup() |>
    filter(method != 'assembly') |>
    mutate(
        error = abs(length - true_length),
        method = recode_factor(method, 'cnvnator' = 'CNVnator', 'floco' = 'Floco'),
        sample_class = sample_class[sample],
    )

hap_branch_cn <- read_delim('data/hap_branch_cn.csv.xz', '\t') |>
    mutate(cn = round(cn), length = NULL)
floco_branch_cn <- read_delim('data/floco_branch_cn.csv.xz', '\t') |>
    mutate(cn = round(cn), length = NULL)
cnvnator_branch_cn <- read_delim('data/cnvnator_branch_cn.csv.xz', '\t') |>
    mutate(cn = round(cn), length = NULL)
branch_cn <- full_join(
        hap_branch_cn,
        full_join(floco_branch_cn, cnvnator_branch_cn,
            join_by(locus, sample, branch), suffix = c('.floco', '.cnvnator')),
        join_by(locus, sample, branch)) |>
    mutate(cn.cnvnator = replace_na(cn.cnvnator, 0)) |>
    filter(sample != 'HG01114' | locus != 'AMY') |>
    left_join(mutate(flanks, is_flank = T), join_by(locus, branch)) |>
    mutate(is_flank = replace_na(is_flank, F))

branch_cn_long <- branch_cn |>
    pivot_longer(
        cols = c('cn.floco', 'cn.cnvnator'),
        names_to = 'method', values_to = 'pred_cn') |>
    mutate(method = sub('cn.', '', method) |>
        recode_factor('cnvnator' = 'CNVnator', 'floco' = 'Floco'))

branch_summary <- branch_cn_long |>
    filter(!is_flank) |>
    group_by(locus, sample, method) |>
    summarize(rel_err = mean(abs(cn - pred_cn) / pmax(cn, 1)),
        .groups = 'keep') |>
    mutate(sample_class = sample_class[sample])

absent_sample <- data.frame(locus = '*AMY*', sample = 'HG01114', sample_class = 'in graph')

(g1 <- mutate(all_lengths2, locus = sprintf('*%s*', locus)) |>
ggplot() +
    geom_text(data = absent_sample, aes(sample, y = 3.2, label = 'x'), size = 3) +
    geom_bar(aes(sample, error * 1e-3, fill = method),
        color = 'black', linewidth = 0.2,
        stat = 'identity', position = 'dodge', width = 0.7) +
    annotate('segment', x = 0.65, y = 0, xend = 5.35, yend = 0) +
    scale_fill_manual(NULL, values = colors) +
    facet_grid(locus ~ sample_class, space = 'free', scales = 'free') +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_continuous('Error in predicted length (kb)',
        expand = expansion(mult = c(0, 0.02)),
        breaks = seq(0, 150, 30)) +
    theme_bw() +
    theme(
        text = element_text(family = 'Source Sans 3'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        strip.background = element_rect(color = NA, fill = 'gray90'),
        strip.text.x = element_text(margin = margin(t = 2, b = 2)),
        strip.text.y = ggtext::element_markdown(margin = margin(l = 2, r = 2)),
        legend.position = 'bottom',
        legend.justification = c('right', 'center'),
        legend.key.size = unit(0.7, 'lines'),
        legend.margin = margin(t = -5, b = -2, r = -10),
        panel.spacing.x = unit(0.6, 'lines'),
        panel.spacing.y = unit(0.7, 'lines'),
    ))

(g2 <- mutate(branch_summary, locus = sprintf('*%s*', locus)) |>
ggplot() +
    geom_text(data = absent_sample, aes(sample, y = 0.008, label = 'x'), size = 3) +
    geom_bar(aes(sample, rel_err, fill = method),
        color = 'black', linewidth = 0.2,
        stat = 'identity', position = 'dodge', width = 0.7) +
    annotate('segment', x = 0.65, y = 0, xend = 5.35, yend = 0) +
    facet_grid(locus ~ sample_class, space = 'free', scales = 'free') +
    scale_fill_manual(NULL, values = colors) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_continuous(
        'Mean relative error in branch copy number',
        expand = expansion(mult = c(0, 0.02)),
        breaks = seq(0.0, 1.0, 0.1),
        ) +
    theme_bw() +
    theme(
        text = element_text(family = 'Source Sans 3'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(),
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        strip.background = element_rect(color = NA, fill = 'gray90'),
        strip.text.x = element_text(margin = margin(t = 2, b = 2)),
        strip.text.y = ggtext::element_markdown(margin = margin(l = 2, r = 2)),
        legend.position = 'none',
        legend.key.size = unit(0.7, 'lines'),
        legend.margin = margin(t = -5, b = -2),
        panel.spacing.x = unit(0.6, 'lines'),
        panel.spacing.y = unit(0.7, 'lines'),
    ))

cowplot::plot_grid(g1, g2,
    labels = letters, label_fontfamily = 'Source Sans 3', label_fontface = 'bold')
ggsave('genes_cmp.svg', width = 10, height = 5, device = svglite::svglite,
    scale = 0.85)
