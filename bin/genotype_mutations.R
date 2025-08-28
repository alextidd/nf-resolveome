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

# get mutations
muts <-
  readr::read_tsv(opts$mutations) %>%
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

    # check type
    if (!(type %in% c("snv", "dnv", "mnv", "ins", "del"))) {
      stop(paste("mutation type", type, "not recognized!"))
    }

    # look ahead if deletion
    if (type == "del") {
      pos_i <- pos + 1
    } else {
      pos_i <- pos
    }

    # query bam
    calls <- deepSNV::bam2R(bam, chr, pos_i, pos_i, mask = opts$mask,
                            q = opts$min_bq, mq = opts$min_mq)

    # calculate total depth
		total_depth <-
      sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t", "-", "_")],
          na.rm = TRUE)

    # calculate mut depth
		if (type == "del") {
		  mut_depth <- sum(calls[, c("-", "_")])
		} else if(type == "ins") {
			mut_depth <- sum(calls[, c("INS", "ins")])
		} else {
			mut_i <- unlist(strsplit(alt, ""))[1]
			mut_depth <- calls[1, mut_i] + calls[1, tolower(mut_i)]
		}

    # return
    tibble::tibble(chr = chr, pos = pos, ref = ref, alt = alt,
                   total_depth = total_depth, mut_depth = mut_depth) %>%
      dplyr::mutate(mut_vaf = mut_depth / total_depth)

  }) %>%
  dplyr::bind_rows()

# write alt calls to out
geno %>%
  dplyr::full_join(muts) %>%
  dplyr::mutate(id = opts$id) %>%
  readr::write_tsv(paste(opts$id, "genotyped_mutations.tsv", sep = "_"))