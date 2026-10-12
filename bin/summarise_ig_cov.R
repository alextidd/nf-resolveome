#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)
library(dplyr)

# options
option_list <- list(make_option("--ids", type = "character"),
                    make_option("--donor_id", type = "character"),
                    make_option("--seq_type", type = "character"),
                    make_option("--csr_max_norm_cov", type = "double"),
                    make_option("--light_chain_max_norm_cov", type = "double"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

ids <- strsplit(opts$ids, ",")[[1]]

# get the mean coverage of the IG regions per cell
# (the regions bed is named CSR, IGKC and IGLC1-7; regions that fall outside
# the hyb panel are dropped by the intersection, so are missing here)
ig_cov <-
  ids %>%
  purrr::set_names() %>%
  purrr::map(function(id) {
    readr::read_tsv(gzfile(paste0(id, ".ig.regions.bed.gz")),
                    col_names = c("chr", "start", "end", "region", "mean_cov"),
                    col_types = "ciicd")
  }) %>%
  bind_rows(.id = "id") %>%
  mutate(len = end - start)

# get the overall mean coverage of each cell, to normalise against
cell_cov <-
  ids %>%
  purrr::set_names() %>%
  purrr::map(function(id) {
    readr::read_tsv(paste0(id, ".mosdepth.summary.txt"), col_types = readr::cols()) %>%
      dplyr::rename(chr = chrom)
  }) %>%
  bind_rows(.id = "id") %>%
  filter(chr == "total" & opts$seq_type == "dna" |
          chr == "total_region" & opts$seq_type == "dnahyb") %>%
  transmute(id, cell_cov = mean)

# get the length-weighted mean coverage of each set of regions
ig_cov_per_cell <-
  ig_cov %>%
  mutate(set = case_when(region == "CSR" ~ "csr_cov",
                         region == "IGKC" ~ "igk_cov",
                         region %in% paste0("IGLC", 1:7) ~ "igl_cov")) %>%
  filter(!is.na(set)) %>%
  group_by(id, set) %>%
  summarise(cov = sum(mean_cov * len) / sum(len), .groups = "drop") %>%
  # keep every cell and every region, even when nothing was covered
  tidyr::complete(id = ids, set = c("csr_cov", "igk_cov", "igl_cov")) %>%
  tidyr::pivot_wider(names_from = set, values_from = cov) %>%
  left_join(cell_cov, by = "id") %>%
  # normalise to the overall coverage of the cell
  mutate(csr_norm = csr_cov / cell_cov,
         igk_norm = igk_cov / cell_cov,
         igl_norm = igl_cov / cell_cov) %>%
  # call class switch recombination and light chain usage
  # (both calls detect biallelic deletion of the region - a cell that has
  # rearranged only one allele retains ~half coverage)
  mutate(
    csr = csr_norm < opts$csr_max_norm_cov,
    igk_deleted = igk_norm < opts$light_chain_max_norm_cov,
    igl_deleted = igl_norm < opts$light_chain_max_norm_cov,
    light_chain = case_when(
      # a cell with neither light chain locus covered has failed, not switched
      igk_deleted & igl_deleted ~ NA_character_,
      igk_deleted ~ "IGL",
      !igk_deleted ~ "IGK/germline")) %>%
  select(id, cell_cov, csr_cov, igk_cov, igl_cov,
         csr_norm, igk_norm, igl_norm, csr, igk_deleted, light_chain) %>%
  arrange(id)

# save
ig_cov_per_cell %>%
  readr::write_tsv(paste0(opts$donor_id, "_ig_cov.tsv"))
