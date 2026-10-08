#!/usr/bin/env Rscript

suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument('-i', '--input', metavar = 'FILE', required = T,
    help = 'Input CSV file.')
parser$add_argument('-o', '--output', metavar = 'FILE', required = T,
    help = 'Output CSV file with first column from the input CSV file and the color.')
parser$add_argument('-d', '--delim', metavar='CHAR', default = '\t',
    help = 'Input CSV delimiter [default: tab].')
parser$add_argument('-c', '--column', metavar='STR', required = T,
    help = 'Color the nodes based on this column')
parser$add_argument('-r', '--range', nargs = 2, metavar='NUM', type = 'double', default = c(0.0, 1.0),
    help = 'Truncate the column between these two values [%(default)s].')
parser$add_argument('-p', '--palette', metavar = 'STR', default = 'Viridis',
    help = 'HCL color palette.')
parser$add_argument('-R', '--reverse', action = 'store_true',
    help = 'Reverse palette order.')
args <- parser$parse_args()

palette <- hcl.colors(256, args$palette)
if (args$reverse) {
    palette <- rev(palette)
}

assign_colors <- function(x, colors, minv, maxv) {
    x <- pmin(x, maxv) |> pmax(minv)
    ixs <- round(1 + (length(colors) - 1) * (x - minv) / (maxv - minv))
    colors[ixs]
}

df <- read.csv(args$input, sep = args$delim)
df$color <- assign_colors(df[[args$column]], palette, args$range[1], args$range[2])
df <- df[, c(1, ncol(df))]
write.table(df, args$output, sep = ',', row.names = F, quote = F)
