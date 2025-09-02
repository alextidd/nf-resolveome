#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(
  make_option("--mutations", type = "character"),
  make_option("--bam", type = "character"),
  make_option("--min_bq", type = "numeric"),
  make_option("--mask", type = "numeric"),
  make_option("--min_mq", type = "numeric"),
  make_option("--id", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# genotype mutations
muts <- readr::read_tsv(opts$mutations)
geno <- alexr::genotype_variants(variants = muts,
                                 bam = opts$bam,
                                 min_bq = opts$min_bq,
                                 min_mq = opts$min_mq,
                                 mask = opts$mask)

# write alt calls to out
geno %>%
  dplyr::full_join(muts) %>%
  dplyr::mutate(id = opts$id) %>%
  readr::write_tsv(paste(opts$id, "genotyped_mutations.tsv", sep = "_"))