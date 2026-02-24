library(tidyverse)
library(ggplot2)
library(patchwork)
library(colorspace)

need_sci_fmt <- function(x) {
    max(abs(log10(x)), na.rm = T) >= 4
}

sci_fmt <- function(x) {
    text <- sprintf('%.1g', x) |>
        str_replace('e(-?)\\+?0*', 'e\\1') |> # e+X -> eX, e0X -> eX
        str_replace('e', ' %*% 10^') |>
        str_replace('^1 %\\*% ', '') |> # 1*10^ -> 10^
        sapply(function(s) parse(text = s))
}

FONT <- 'Source Sans 3'

CN_GROUPS <- 3
CN_LEVELS <- c(as.character(0:(CN_GROUPS - 2)), 'more')
CN_TEXT <- c(as.character(0:(CN_GROUPS - 2)), paste0('≥ ', CN_GROUPS - 1)) |>
    setNames(CN_LEVELS)
CN_COLORS <- c('#D94F01', '#322187', '#23A8C7') |> setNames(CN_LEVELS)
group_cn <- function(cn) {
    ifelse(cn < CN_GROUPS - 1, as.character(cn), 'more')
}

BAR_WIDTH <- 0.8
X_EXPANSION <- 0.5 * BAR_WIDTH + 0.2
Y_EXPANSION <- 0.03

draw_columns <- function(combinations_cropped, col, maxy = NA, ytitle_margin = 0) {
    combinations_cropped$val <- combinations_cropped[[col]]
    
    if (need_sci_fmt(range(combinations_cropped$val))) {
        y_labels = sci_fmt
        y_vjust = 0.4
    } else {
        y_labels = waiver()
        y_vjust = 0.5
    }
    
    yexpansion <- c(0.01, 0.01)
    yexpansion[ifelse(col == 'nodes', 1, 2)] <- Y_EXPANSION
    
    ggplot(combinations_cropped) +
        annotate('segment', x = 1, xend = nrow(combinations_cropped), y = 0, yend = 0,
            color = 'gray20', linewidth = 0.2) +
        geom_bar(aes(code, val, fill = val),
            stat = 'identity', width = BAR_WIDTH, color = 'black', linewidth = 0.2,
            show.legend = F) +
        scale_x_discrete(NULL, expand = expansion(add = X_EXPANSION)) +
        scale_y_continuous(ifelse(col == 'nodes', 'Nodes', 'Length (Gb)'),
            breaks = scales::pretty_breaks(3),
            labels = y_labels,
            expand = expansion(mult = yexpansion),
            limits = c(0, maxy)
        ) +
        scale_fill_continuous_sequential('Grays', begin = 0.1) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            panel.grid = element_blank(),
            legend.position = 'top',
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 10, margin = margin(r = ytitle_margin)),
            axis.text.y = element_text(vjust = y_vjust),
            plot.margin = margin(),
            panel.spacing.x = unit(0, 'lines'),
        )
}

draw_points <- function(combinations_cropped, dataset_labels) {
    datasets <- names(dataset_labels)
    combinations_points <- combinations_cropped[, c(datasets, 'code')] |>
        unique() |>
        pivot_longer(all_of(datasets), names_to = 'dataset', values_to = 'event') |>
        mutate(
            event = factor(event, levels = CN_LEVELS),
            dataset = factor(dataset, levels = rev(datasets)))

    comb_other <- filter(combinations_points, code == 'Other')
    n_datasets <- length(datasets)
    ggplot(combinations_points) +
        geom_line(aes(code, dataset, group = paste0(code, event), color = event),
            show.legend = F) +
        geom_point(aes(code, dataset, color = event),
            size = 3.5) +
        geom_tile(
            aes(code, y = 0.5 * (n_datasets + 1)),
            width = BAR_WIDTH, height = n_datasets - 0.5, fill = 'gray90',
            data = comb_other) +
        geom_text(
            aes(code, y = 0.5 * (n_datasets + 1), label = 'Other'),
            angle = 90, vjust = 0.5, hjust = 0.5, family = FONT, size = 3.5,
            data = comb_other) +
        scale_x_discrete(NULL, expand = expansion(add = X_EXPANSION)) +
        scale_y_discrete(labels = dataset_labels) +
        scale_color_manual('Copy number',
            breaks = CN_LEVELS, labels = CN_TEXT, values = CN_COLORS,
            na.value = 'white', drop = F) +
        guides(color = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            legend.position = 'right',
            legend.box.background = element_rect(fill = NA),
            legend.background = element_rect(fill = NA),
            legend.title = element_text(size = 9, hjust = 0.5, margin = margin(b = -2)),
            legend.text = element_text(margin = margin(l = -1, r = -2)),
            legend.margin = margin(r = 30, b = 5),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = 'black'),
            plot.margin = margin(),
        )
}

draw_totals <- function(totals, col, dataset_labels) {
    totals$val <- totals[[col]]
    
    yexpansion <- c(0.01, 0.01)
    yexpansion[ifelse(col == 'nodes', 1, 2)] <- Y_EXPANSION
    
    ggplot(totals) +
        geom_bar(aes(dataset, val, fill = cn), stat = 'identity',
            #color = 'black', linewidth = 0.2,
            show.legend = T) +
        scale_x_discrete(labels = dataset_labels, limits = rev) +
        scale_y_continuous(ifelse(col == 'nodes', 'Number of nodes', 'Sum length (Gb)'),
            breaks = scales::pretty_breaks(3),
            labels = if (need_sci_fmt(totals$val)) { sci_fmt } else { waiver() },
            expand = expansion(mult = yexpansion),
        ) +
        scale_fill_manual('Copy number', breaks = CN_LEVELS, labels = CN_TEXT,
            values = CN_COLORS, drop = F) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            panel.grid = element_blank(),
            legend.position = 'none',
            legend.background = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank(),
            axis.text.x = if (col == 'nodes') { element_blank() } else {
                element_text(size = 8, hjust = 1, angle = 40, color = 'black') },
            axis.ticks.x = if (col == 'nodes') { element_blank() } else { element_line() },
            plot.margin = margin(),
        )
}

##########################

# callback: function(tag, n_bars, plot, combinations df)
draw_upset <- function(filename_or_data, tag, n_bars, dataset_labels, callback,
        nodes = NULL, title = NULL) {
    
    if (is.character(filename_or_data)) {
        cat('Loading data\n')
        data0 <- read.csv(filename_or_data, sep = ',') |>
            rename('Node' = 'X.Node') |>
            mutate(Length = 1e-9 * Length)
    } else {
        data0 <- filename_or_data
    }

    if (!is.null(nodes)) {
        data0 <- filter(data0, Node %in% nodes)
    }
    dataset_labels <- dataset_labels[names(dataset_labels) %in% colnames(data0)]
    col_ytitle_margin = ifelse(max(str_length(dataset_labels)) > 5, -7, 3)
    datasets <- names(dataset_labels)

    cat('Processing data\n')
    data <- mutate_at(data0, vars(all_of(datasets)), group_cn) |>
        unite('code', all_of(datasets), sep = '-', remove = F)
    total_nodes <- nrow(data)
    sum_length <- sum(data$Length)
    
    combinations <- group_by_at(data, c(datasets, 'code')) |>
        summarize(nodes = length(Node), len = sum(Length), .groups = 'keep') |>
        ungroup() |>
        arrange(-len) |>
        mutate(code = factor(code, levels = c(code, 'Other')))
    
    totals <- pivot_longer(combinations,
            cols = all_of(datasets), names_to = 'dataset', values_to = 'cn') |>
        mutate(
            cn = factor(cn, levels = rev(CN_LEVELS)),
            dataset = factor(dataset, levels = rev(datasets))) |>
        group_by(dataset, cn) |>
        summarize_at(vars(nodes, len), sum) |>
        ungroup()

    n_bars <- sort(n_bars)
    for (curr_bars in n_bars) {
        combinations_cropped <- if (nrow(combinations) > curr_bars + 1) {
            bind_rows(
                combinations[1:curr_bars,],
                summarize_at(combinations[(curr_bars+1):nrow(combinations),],
                    vars(nodes, len), sum) |> add_column(code = 'Other')
            ) |> mutate(code = factor(code, levels = levels(combinations$code)))
        } else {
            combinations
        }
        actual_bars <- nrow(combinations_cropped)

        width2 <- 2.5 / (pmin(curr_bars, 15) + 2.5)
        cat(sprintf('Drawing with %d bars\n', actual_bars))
        g <- draw_columns(combinations_cropped, 'nodes', ytitle_margin = col_ytitle_margin) +
            draw_totals(totals, 'nodes', dataset_labels) +
            draw_columns(combinations_cropped, 'len', ytitle_margin = col_ytitle_margin) +
            draw_totals(totals, 'len', dataset_labels) +
            draw_points(combinations_cropped, dataset_labels) +
            guide_area() +
            plot_layout(
                design = 'AB\nCD\nED\nEF',
                heights = c(0.42, 0.42, 0.00, 0.16),
                widths = c(1 - width2, width2),
                guides = 'collect') +
            plot_annotation(
                title = title,
                theme = theme(
                    plot.margin = margin(1, 1, 1, 1),
                    plot.title = element_text(hjust = 0.5, family = FONT, size = 10),
                ))

        callback(tag, actual_bars, g, combinations_cropped)
        if (actual_bars == nrow(combinations)) {
            break
        }
    }
}

plots <- list()
callback <- function(tag, n_bars, g, combinations_cropped) {
    plots[[tag]] <<- wrap_elements(full = g)

    s <- if ('Other' %in% combinations_cropped$code) {
        sprintf('%d_bars', n_bars) } else { 'all' }
    df <- select(combinations_cropped, -code) |>
        mutate_at(vars(-c(nodes, len)), function(x) CN_TEXT[x]) |>
        write.table(sprintf('csvs/%s-%s.csv', tag, s),
            sep = '\t', quote = F, row.names = F)
}

BARS <- 10
Y_EXPANSION <- 0.04
DATASETS <- c('Chopped' = 'Assembly', 'HiFi' = 'HiFi', 'ONT' = 'ONT')
draw_upset('tech_cmp_HG01258-pang.csv.gz', 'HG01258-pang', BARS, DATASETS, callback,
    title = 'HG01258-pang')
draw_upset('tech_cmp_HG01114-pang.csv.gz', 'HG01114-pang', BARS, DATASETS, callback,
    title = 'HG01114-pang')
draw_upset('tech_cmp_HG01114-asm.csv.gz', 'HG01114-asm', BARS, DATASETS, callback,
    title = 'HG01114-asm')
draw_upset('tech_cmp_Altus.csv.gz', 'Altus', BARS, DATASETS[c('HiFi', 'ONT')], callback,
    title = 'Altus-asm')
plots$`HG01114-asm` + plots$Altus + plots$`HG01258-pang` + plots$`HG01114-pang` +
    plot_layout(design = 'AB\nCD') +
    plot_annotation(tag_levels = 'a') &
    theme(
        plot.tag = element_text(family = FONT, face = 'bold', size = 16),
        plot.tag.position = c(-0.003, 0.985),
        plot.margin = margin(l = 3, t = 2, r = -3, b = -8))
ggsave('upset.svg', width = 10, height = 6, scale = 0.85)

Y_EXPANSION <- 0.03
draw_upset('tech_cmp_HG01258-pang.csv.gz', 'S.HG01258-pang', Inf, DATASETS, callback,
    title = 'HG01258-pang')
draw_upset('tech_cmp_HG01114-pang.csv.gz', 'S.HG01114-pang', Inf ,DATASETS, callback,
    title = 'HG01114-pang')
draw_upset('tech_cmp_HG01114-asm.csv.gz', 'S.HG01114-asm', Inf, DATASETS, callback,
    title = 'HG01114-asm')
draw_upset('tech_cmp_Altus.csv.gz', 'S.Altus', Inf,  DATASETS, callback,
    title = 'Altus-asm')
plots$`S.HG01114-asm` + plots$`S.Altus` + plots$`S.HG01258-pang` + plots$`S.HG01114-pang` +
    plot_layout(design = 'AB\nCD') +
    plot_annotation(tag_levels = 'a') &
    theme(
        plot.tag = element_text(family = FONT, face = 'bold', size = 16),
        plot.tag.position = c(0.007, 0.975),
        plot.margin = margin(l = 1, t = 1, r = -5, b = -10))
ggsave('upset.supp.svg', width = 11, height = 6, scale = 1)

compare_techs <- function(filename_or_data, dataset_labels, max_cn, ..., nodes = NULL) {
    if (is.character(filename_or_data)) {
        cat('Loading data\n')
        df <- read.csv(filename_or_data, sep = ',') |>
            rename('Node' = 'X.Node') |>
            mutate(Length = 1e-9 * Length)
    } else {
        df <- filename_or_data
    }
    if (!is.null(nodes)) {
        df <- filter(df, Node %in% readLines(nodes))
    }
    df <- filter(df, ...)

    datasets <- names(dataset_labels)
    df <- mutate_at(df, vars(all_of(datasets)), function(x) pmin(x, max_cn))
    total_nodes <- nrow(df)
    total_len <- sum(df$Length)
    conc <- function(where) {
        n <- sum(where)
        l <- sum(df$Length[where])
        sprintf('%8s nodes (%.1f%%), %5s Mb (%.1f%%)',
            format(n, big.mark = ','), round(100 * n / total_nodes, 1),
            format(round(1000 * l), big.mark = ','), round(100 * l / total_len, 1))
    }
    
    cat(sprintf('%s %8s nodes,         %5s Mb\n', format('Total', width = 19, justify = 'centre'),
        format(total_nodes, big.mark = ','),
        format(round(1000 * total_len), big.mark = ',')))
    if (length(dataset_labels) == 3) {
        x <- df[[datasets[1]]]
        y <- df[[datasets[2]]]
        z <- df[[datasets[3]]]
        cat(sprintf('%s %s\n', format('All', width = 19, justify = 'centre'),
            conc(x == y & x == z)))
    }
    for (i in 1:(length(datasets) - 1)) {
        for (j in (i+1):length(datasets)) {
            d1 <- datasets[[i]]
            d2 <- datasets[[j]]
            cat(sprintf('%-8s & %-8s %s, r = %.5f\n', dataset_labels[d1], dataset_labels[d2],
                conc(df[[d1]] == df[[d2]]), cor(df[[d1]], df[[d2]])))
        }
    }
}

DATASETS <- c('Chopped' = 'Assembly', 'HiFi' = 'HiFi', 'ONT' = 'ONT')
compare_techs('tech_cmp_HG01114-asm.csv.gz', DATASETS, 4)
compare_techs('tech_cmp_Altus.csv.gz', DATASETS, 4)
compare_techs('tech_cmp_HG01258-pang.csv.gz', DATASETS, 4)
compare_techs('tech_cmp_HG01114-pang.csv.gz', DATASETS, 4)

Y_EXPANSION <- 0.05
chm13_data <- read.csv('tech_cmp_CHM13_ext.csv.gz', sep = ',') |>
    rename('Node' = 'X.Node') |> mutate(Length = 1e-9 * Length)
for (cov in c(30, 1)) {
    for (method in c('flow', 'ml')) {
        chm13_datasets <- setNames(c('Baseline', 'HiFi', 'ONT'),
            c('Ground_Truth', sprintf('%s_%sx_%s', c('HiFi', 'ONT'), cov, method)))
        draw_upset(chm13_data, sprintf('CHM13-%sx_%s', cov, method),
            Inf, chm13_datasets, callback,
            title = sprintf('%s× data, %s network flow',
                cov, ifelse(method == 'flow', 'after', 'before')))
    }
}

plots$`CHM13-30x_ml` + plots$`CHM13-30x_flow` + plots$`CHM13-1x_ml` + plots$`CHM13-1x_flow` +
    plot_layout(design = 'AB\nCD') +
    plot_annotation(tag_levels = 'a') &
    theme(
        plot.tag = element_text(family = FONT, face = 'bold', size = 16),
        plot.tag.position = c(0.007, 0.98),
        plot.margin = margin(l = 1, t = 1, r = -2, b = -10))
ggsave('upset.chm13.svg', width = 10, height = 6, scale = 1)
