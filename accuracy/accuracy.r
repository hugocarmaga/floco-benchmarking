library(tidyverse)
library(ggplot2)
library(cowplot)

cowplot::set_null_device('agg')

load <- function(fname) {
    read.csv(fname, sep = ',') |> rename(Node = X.Node)
}

# base = baseline = ground truth
calc_metrics <- function(base, pred, node_len = NULL) {
    stopifnot(all(base <= 1))
    if (is.null(node_len)) {
        node_len <- rep(1, length(base))
    }
    
    res <- list()
    res$accuracy <- sum(node_len * (base == pred)) / sum(node_len)
    res$total0 <- sum(node_len * (base == 0))
    res$total1 <- sum(node_len * (base == 1))
    res$tp <- sum(node_len * ((base == 1) & (pred == 1)))
    res$fp <- sum(node_len * ((base == 0) & (pred != 0)))
    res$fn <- sum(node_len * ((base == 1) & (pred != 1)))
    res$tn <- sum(node_len * ((base == 0) & (pred == 0)))
    res$precision <- res$tp / (res$tp + res$fp)
    res$recall <- res$tp / (res$tp + res$fn)
    res$f1 = 2 * res$tp / (2 * res$tp + res$fp + res$fn)
    as.data.frame(res)
}

calc_col_metrics <- function(df, tech) {
    out_df <- data.frame()
    for (weighted in 0:1) {
        curr_df <- Filter(function(x) startsWith(x, 'X'), colnames(df)) |>
            lapply(function(col)
                calc_metrics(df$Ground_Truth, df[[col]],
                        if (weighted) { df$Length  * 1e-9 } else { NULL }) |>
                    add_column(key = col)) %>%
            do.call(rbind, .) |>
            pivot_longer(-key, names_to = 'metric') |>
            separate(key, into = c('cov', 'method'), sep = '_') |>
            mutate(
                tech = tech, weighted = as.logical(weighted),
                cov = gsub('X', '', cov, ignore.case = T) |> as.numeric(),
                )
        out_df <- rbind(out_df, curr_df)
    }
    out_df
}

ont <- load('ont.csv.gz')
hifi <- load('hifi.csv.gz')
chopped <- load('chopped.csv.gz')
res <- rbind(
    calc_col_metrics(ont, 'ONT'),
    calc_col_metrics(hifi, 'HiFi'),
    calc_col_metrics(chopped, 'Assembly')) |>
    mutate(
        tech = factor(tech, levels = c('HiFi', 'ONT', 'Assembly')),
        type = ifelse(weighted, 'bp', 'nodes'),
        method = recode_factor(method, 'ml' = 'Before flow', 'flow' = 'After flow')) |>
    select(!weighted)

colors <- c(
    '#134686',
    '#ED3F27',
    '#FEB21A'
    )
draw_panel <- function(res_subset) {
    title <- sprintf('%s (%s)', tools::toTitleCase(res_subset$metric[1]), res_subset$type[1])
    lower_limit <- case_match(res_subset$metric[1],
        c('accuracy', 'recall') ~ 0.3,
        'precision' ~ 0.85,
    )

    ggplot(res_subset) +
        geom_bar(aes(method, value, fill = tech, alpha = cov, group = tech_cov),
            stat = 'identity', position = position_dodge(), color = 'gray20', linewidth = 0.1) +
        geom_text(aes(method, lower_limit, group = tech_cov, color = tech,
            label = ifelse(value < lower_limit, '▼', ''), family = FONT),
            position = position_dodge(width = 0.9)) +
        geom_vline(xintercept = 1.5, color = 'gray20', linewidth = 0.4, linetype = '13') +
        coord_cartesian(ylim = c(lower_limit, 1)) +
        scale_x_discrete(NULL, expand = expansion(add = 0.5)) +
        scale_y_continuous(NULL, breaks = scales::pretty_breaks(3),
            expand = expansion(mult = c(0.03, 0.0))) +
        scale_fill_manual('Data type', values = colors) +
        scale_color_manual('Data type', values = colors) +
        scale_alpha_manual('Coverage',
            values = seq(0.2, 1, length.out = length(levels(res_subset$cov)))) +
        ggtitle(title) +
        guides(alpha = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            panel.grid = element_blank(),
            plot.title = element_text(size = 10, face = 'bold', hjust = 0.5, margin = margin(t = 2, b = 2)),
            plot.margin = margin(t = 3, r = 3, b = 3, l = 3),
            legend.position = 'bottom',
            legend.key.size = unit(0.8, 'lines'),
            axis.text.x = if (grepl('Recall', title)) {
                    element_text(size = 9.5, color = 'black') } else { element_blank() },
            axis.ticks.x = if (grepl('Recall', title)) { element_line() } else { element_blank() },
            panel.border = element_blank(),
            axis.line = element_line(linewidth = 0.25),
            plot.background = element_rect(color = NA, fill = NA),
        )
}

res2 <- filter(res, metric %in% c('accuracy', 'precision', 'recall')) |>
    filter((tech == 'Assembly' & cov == 20)
        | (tech != 'Assembly' & cov %in% c(1, 10, 20, 30))) |>
    arrange(tech, cov) |>
    mutate(
        cov = factor(cov),
        tech_cov = paste0(tech, cov),
        tech_cov = factor(tech_cov, levels = unique(tech_cov)),
    )
FONT <- 'Source Sans 3'
panels <- list()
for (.metric in unique(res2$metric)) {
    for (.type in unique(res2$type)) {
        panels[[paste0(.metric, '_', .type)]] <-
            draw_panel(filter(res2, metric == .metric & type == .type))
    }
}
no_legend <- theme(legend.position = 'none')

y_lab <- ggplot() + 
    annotate(geom = 'text', x = 0, y = 0,
        label = 'Value',
        family = FONT, angle = 90) +
    coord_cartesian(clip = 'off') +
    theme_void() +
    theme(plot.margin = margin(r = -5))

plot_grid(
    panels$accuracy_nodes + no_legend,
    panels$accuracy_bp + no_legend,
    panels$precision_nodes + no_legend,
    panels$precision_bp + no_legend,
    panels$recall_nodes + no_legend,
    panels$recall_bp + no_legend,
    ncol = 2, align = 'v', rel_heights = c(1, 1, 1.1)
    ) %>%
plot_grid(
    y_lab, .,
    ncol = 2,
    rel_widths = c(0.03, 0.97)) %>%
plot_grid(get_legend(panels$accuracy_nodes), ., ncol = 1, rel_heights = c(0.05, 0.95))
# ggsave('accuracy.png', width = 7, height = 5, dpi = 500, bg = 'white')
ggsave('accuracy.svg', width = 7, height = 5, dpi = 500, bg = 'white')
