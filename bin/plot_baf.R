#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)
library(ggplot2)

# options
option_list <- list(make_option("--geno", type = "character"),
                    make_option("--id", type = "character"),
                    make_option("--baf_chrs", type = "character"),
                    make_option("--baf_genes", type = "character"),
                    make_option("--refcds", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# read geno
geno <- readr::read_tsv(opts$geno)

# get baf genes
genes <- ifelse(is.null(opts$baf_genes), NULL,
                unlist(strsplit(opts$baf_genes, ",")))

# plot all chromosomes
p <- alexr::plot_baf(geno, "caveman_snps", genes = genes)
ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_plot.png"),
              width = 5000, height = 1200, res = 300)
print(p)
dev.off()

# plot chromosomes of interest, zoomed in
if (!is.null(opts$baf_chrs)) {
  baf_chrs <- unlist(strsplit(opts$baf_chrs, ","))
  for (chr_i in baf_chrs) {
    p <-
      geno %>%
      dplyr::filter(chr == chr_i) %>%
      alexr::plot_baf(p_source = "caveman_snps", genes = genes,
                      p_alpha = 0.2, p_size = 0.8)
    ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_chr", chr_i, "_plot.png"),
                  width = 5000, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}