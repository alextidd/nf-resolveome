#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)
library(ggplot2)

# options
option_list <- list(make_option("--geno", type = "character"),
                    make_option("--id", type = "character"),
                    make_option("--baf_chrs", type = "character"),
                    make_option("--baf_genes", type = "character"),
                    make_option("--refcds", type = "character"),
                    make_option("--min_depth", type = "integer", default = 20),
                    make_option("--raster_dpi", type = "integer", default = 200))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# page size, in inches. this is passed to plot_baf() as well as to the device,
# as the point size and alpha depend on how much of the page the points are
# spread over
p_width <- 16
p_height <- 4

# write a plot to PDF, using cairo so that the rasterised points are embedded
# as compressed images
write_pdf <- function(p, file) {
  if (is.null(p)) return(invisible(NULL))
  cairo_pdf(file, width = p_width, height = p_height)
  print(p)
  dev.off()
}

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
  # (point size and alpha are chosen from how densely the SNPs fall on the
  # page, and the points are rasterised, so the PDF stays small and readable)
  p <-
    alexr::plot_baf(p_dat = geno, p_title = opts$id,
                    refcds = opts$refcds, genes = genes,
                    min_depth = opts$min_depth,
                    raster_dpi = opts$raster_dpi,
                    p_width = p_width, p_height = p_height)
  write_pdf(p, paste0(opts$id, "_snps_baf_plot.pdf"))

  # plot chromosomes of interest, zoomed in
  if (!is.null(opts$baf_chrs)) {
    baf_chrs <- unlist(strsplit(opts$baf_chrs, ","))
    for (chr_i in baf_chrs) {
      p <-
        geno %>%
        dplyr::filter(chr == chr_i) %>%
        alexr::plot_baf(p_title = paste(opts$id, "-", chr_i), genes = genes,
                        refcds = opts$refcds, min_depth = opts$min_depth,
                        raster_dpi = opts$raster_dpi,
                        p_width = p_width, p_height = p_height)
      write_pdf(p, paste0(opts$id, "_snps_baf_", chr_i, "_plot.pdf"))
    }
  }
}
