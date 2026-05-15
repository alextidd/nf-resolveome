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

# check if any non-zero snps
if (geno %>% dplyr::filter(alt_vaf > 0) %>% nrow() == 0) {
  message("No non-zero BAF SNPs found, skipping BAF plotting")
} else {
  message("Plotting BAF")

  # get baf genes
  genes <- if (is.null(opts$baf_genes)) {
    NULL
  } else {
    unlist(strsplit(opts$baf_genes, ","))
  }

  # plot all chromosomes
  p <- alexr::plot_baf(geno, "caveman_snps", genes = genes)
  pdf(paste0(opts$id, "_caveman_snps_baf_plot.pdf"),
      width = 5000 / 300, height = 1200 / 300)
  print(p)
  dev.off()

  # plot chromosomes of interest, zoomed in
  if (!is.null(opts$baf_chrs)) {
    baf_chrs <- unlist(strsplit(opts$baf_chrs, ","))
    for (chr_i in baf_chrs) {
      p <-
        geno %>%
        dplyr::filter(chr == chr_i) %>%
        alexr::plot_baf(genes = genes, p_alpha = 0.2, p_size = 0.8)
      pdf(paste0(opts$id, "_caveman_snps_baf_chr", chr_i, "_plot.pdf"),
          width = 5000 / 300, height = 1200 / 300)
      print(p)
      dev.off()
    }
  }
}