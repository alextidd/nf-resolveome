#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)
library(ggplot2)
library(dplyr)

# options
option_list <- list(make_option("--regions_bed", type = "character"),
                    make_option("--bait_set_vdj", type = "character"),
                    make_option("--id", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# load bait set
ig_tcr_regions <-
  readr::read_tsv(opts$bait_set_vdj,
                  col_names = c("chr", "start", "end", "gene")) %>%
  mutate(
    chr = as.character(chr),
    type = case_when(substr(gene, 1, 2) == "TR" ~ "TCR",
                     substr(gene, 1, 2) == "IG" ~ "BCR",
                     TRUE ~ NA_character_),
    region = paste0(chr, "_", type),
    segment = substr(gene, 4, 4),
    segment = ifelse(segment %in% c("A", "E", "G", "M"), "C", segment) %>%
      factor(levels = c("V", "D", "J", "C")))

# make plots
dat <-
  readr::read_tsv(gzfile(opts$regions_bed),
                  col_names = c("chr", "start", "end", "gene", "mean_cov"),
                  show_col_types = FALSE) %>%
  mutate(id = opts$id, chr = as.character(chr)) %>%
  left_join(ig_tcr_regions) %>%
  mutate(gene = gene %>% forcats::fct_reorder(start)) %>%
  {split(., .$region)}

purrr::walk2(names(dat), dat, function(chr_i, chr_dat) {
  p <-
    chr_dat %>%
    ggplot(aes(x = gene, y = mean_cov, fill = segment)) +
    geom_col() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    guides(x = guide_axis(angle = -90)) +
    theme_classic() +
    ggtitle(paste0(opts$id, " - chr", chr_i, " genes - mean coverage")) +
    scale_fill_brewer(palette = "Dark2")
  ggsave(paste0(opts$id, "_chr", chr_i, "_mean_cov.png"),
         p, height = 5, width = 20)
})
