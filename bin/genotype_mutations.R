#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(
  make_option("--chr", type = "character"),
  make_option("--mutations", type = "character"),
  make_option("--bam", type = "character"),
  make_option("--min_bq", type = "numeric"),
  make_option("--mask", type = "numeric"),
  make_option("--min_mq", type = "numeric"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# get mutations
muts <-
  readr::read_tsv(opts$mutations) %>%
  dplyr::filter(chr == opts$chr) %>%
  dplyr::mutate(
    type = dplyr::case_when(nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
                            nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
                            nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
                            TRUE ~ "complex"))

# check for mutations that are not snv / ins / del
if ("complex" %in% muts$type) {
  message("Complex mutations are not supported!")
  message(paste(muts %>% dplyr::filter(type == "complex") %>% nrow(),
                "complex mutation(s) were found and will be removed."))
  muts <- muts %>% dplyr::filter(type != "complex")
}

# genotype all sites
geno <-
  muts %>%
  dplyr::distinct(chr, pos, ref, alt, type) %>%
  purrr::pmap(function(chr, pos, ref, alt, type) {
    paste(chr, pos, ref, alt, type, "\n") %>% cat()

    # query bam
    calls <- deepSNV::bam2R(opts$bam, chr, pos, pos, q = opts$min_bq,
                            mask = opts$mask, mq = opts$min_mq)

    # count all reads at site
    total_depth <- sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t",
                                  "DEL", "INS", "del", "ins")],
                        na.rm = TRUE)

    # count ref reads at site
    # take first character (in case it is a deletion)
    ref_1 <- substr(ref, 1, 1)
    ref_depth <- sum(calls[, c(ref_1, tolower(ref_1))], na.rm = TRUE)

    # count mutant reads at site
    if (type %in% c("snv", "dnv", "mnv")) {
      # count mutant reads at site
      mut_depth <- sum(calls[, c(alt, tolower(alt))], na.rm = TRUE)
    } else if (type %in% c("ins", "del")) {
      # count ins or del reads at site (don't check sequence)
      mut_depth <- sum(calls[, c(type, toupper(type))], na.rm = TRUE)
    } else {
      stop("alt type not recognised!")
    }

    tibble::tibble(chr = chr, pos = pos, ref = ref, alt = alt,
                   total_depth = total_depth, ref_depth = ref_depth,
                   mut_depth = mut_depth) %>%
      dplyr::mutate(mut_vaf = mut_depth / total_depth)
  }) %>%
  dplyr::bind_rows()

# write alt calls to out
geno %>%
  dplyr::full_join(muts) %>%
  readr::write_tsv("genotyped_mutations.tsv")
