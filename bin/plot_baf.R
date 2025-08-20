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

# get baf genes coordinates
if (!is.null(opts$baf_genes)) {
  baf_genes <- unlist(strsplit(opts$baf_genes, ","))

  # read refcds
  load(opts$refcds)

  # get gene coordinates
  p_genes <-
    baf_genes %>%
    purrr::set_names() %>%
    purrr::map(function(g) {
      i <- RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x$gene_name == g))]]
      tibble::tibble(chr = i$chr,
                     pos = (min(i$intervals_cds) + max(i$intervals_cds)) / 2)
    }) %>%
    dplyr::bind_rows(.id = "gene")

} else {
  p_genes <- NULL
}

# function: plot BAF
plot_baf <- function(p_dat, p_source, p_genes = NULL,
                     p_alpha = 0.05, p_size = 0.01) {
  # prep data
  p_dat2 <-
    p_dat %>%
    dplyr::filter(total_depth > 20) %>%
    dplyr::mutate(mut_baf = 1 - mut_vaf,
                  chr = factor(sub("^chr", "", chr),
                               levels = c(as.character(1:22), "X", "Y"))) %>%
    # get alternating colours and chr bounds
    dplyr::group_by(chr) %>%
    dplyr::mutate(chr_alternating = dplyr::cur_group_id() %% 2,
                  min_pos = min(pos), max_pos = max(pos)) %>%
    # plot
    tidyr::pivot_longer(cols = c("mut_vaf", "mut_baf"), names_to = "vaf_type")

  # plot
  p <-
    p_dat2 %>%
    ggplot(aes(x = pos, y = value)) +
    # add alternating coloured chromosomes
    geom_rect(aes(ymin = 0, ymax = 1, xmin = min_pos, xmax = max_pos,
                  fill = as.factor(chr_alternating)),
              show.legend = FALSE) +
    # add line at vaf = 0.5
    geom_hline(yintercept = 0.5, colour = "red") +
    # add baf points
    geom_point(size = p_size, alpha = p_alpha) +
    # add baf bands
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("white", "#e3e3e3")) +
    ggh4x::facet_grid2(. ~ chr, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA,
                                      linewidth = 0),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0, linetype = "solid"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = paste(opts$id, "-", p_source))

  # add baf genes
  if (!is.null(p_genes)) {
    p <-
      p +
      geom_vline(
        data = p_genes %>% dplyr::mutate(chr = factor(chr, levels = levels(p_dat2$chr))),
        aes(xintercept = pos, colour = gene))
  }

  # return
  p
}

# plot all chromosomes
p <- plot_baf(geno, "caveman_snps", p_genes = p_genes)
ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_plot.png"), width = 5000, height = 1200, res = 300)
print(p)
dev.off()

# plot chromosomes of interest, zoomed in
if (!is.null(opts$baf_chrs)) {
  baf_chrs <- unlist(strsplit(opts$baf_chrs, ","))
  for (chr_i in baf_chrs) {
    p <-
      geno %>%
      dplyr::filter(chr == chr_i) %>%
      plot_baf(p_source = "caveman_snps",
               p_genes = p_genes %>% dplyr::filter(chr == chr_i),
               p_alpha = 0.2, p_size = 0.8)
    ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_chr", chr_i, "_plot.png"), width = 5000, height = 1200, res = 300)
    print(p)
    dev.off()
  }
}