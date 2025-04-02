#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(make_option("--mutations", type = "character"),
                    make_option("--donor_id", type = "character"),
                    makr_option("--refcds", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# load mutations
mutations <- readr::read_tsv(opts$mutations)

# annotate mutations with dndscv
dndscv_in <-
  mutations %>%
  dplyr::select(sampleID = donor_id, chr, pos, ref, mut = alt) %>%
  dplyr::distinct()
dndscv_out <-
  dndscv::dndscv(dndscv_in, max_muts_per_gene_per_sample = Inf,
                 max_coding_muts_per_sample = Inf, outp = 1,
                 refdb = params$refcds)

# collapse multiple annotations of the same mutation
# (we want to maintain 1 mutation-x-cell per row)
annots <-
  dndscv_out$annotmuts %>%
  dplyr::rename(donor_id = sampleID) %>%
  dplyr::rename(alt = mut) %>%
  dplyr::group_by(donor_id, chr, pos, ref, alt) %>%
  dplyr::summarise(dplyr::across(everything(), ~ paste(.x, collapse = ",")),
                   .groups = "drop")

# join mutations with dndscv output, collapse multiple annotations
annot_mutations <- dplyr::right_join(annots, mutations)

# save annotated mutations
annot_mutations %>%
  readr::write_tsv(paste0(opts$donor_id, "_annotated_mutations.tsv"))