#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# libraries
option_list <- list(make_option("--nr", type = "character"),
                    make_option("--nv", type = "character"),
                    make_option("--id", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# get matrices
nr <- read.table(opts$nr)
nv <- read.table(opts$nv)

# calculate vafs
vafs <- nr / nr

# plot heatmap
pdf(paste0(opts$id, "_mut_VAF_heatmap.pdf"))
pheatmap::pheatmap(vafs, main = opts$id)
dev.off()