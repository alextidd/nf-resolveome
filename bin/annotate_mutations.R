#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(make_option("--mutations", type = "character"),
                    make_option("--donor_id", type = "character"),
                    make_option("--refcds", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# load mutations
mutations <- readr::read_tsv(opts$mutations)

# annotate mutations with dndscv
annot_mutations <- alexr::annotate_variants(mutations, refcds = opts$refcds)

# save annotated mutations
annot_mutations %>%
  readr::write_tsv(paste0(opts$donor_id, "_annotated_mutations.tsv"))