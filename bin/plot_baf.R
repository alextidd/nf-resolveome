#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)
library(ggplot2)

# options
option_list <- list(make_option("--geno", type = "character"),
                    make_option("--id", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
print(opts)
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")
# opts <- list(id = "plate3_wellE11_dna_run49882", geno = "out/nf-resolveome/archive_20250313/PD63118/dna/plate3_wellE11_dna_run49882/genotyping/plate3_wellE11_dna_run49882_genotyped_mutations.tsv")

# read geno
geno <- readr::read_tsv(opts$geno) %>% {split(., .$source)}

# function: fit a binomial mixture model
fit_binomial_mixture <- function(alt, cov) {
  data <- data.frame(alt = alt, cov = cov)
  loglik <- function(p1, d) sum(log(0.5 * dbinom(d$alt, d$cov, p1) + 0.5 * dbinom(d$alt, d$cov, 1 - p1)))
  opt <- optimize(function(p) loglik(p, data), interval = c(0, 1), maximum = TRUE)
  out <- list(prob1 = opt$maximum, prob2 = 1 - opt$maximum, mix_prob = 0.5, loglik = opt$objective)
  # return the minimum
  min(out$prob1, out$prob2)
}

# function: plot BAF
plot_baf <- function(p_dat, p_source) {
  p_dat %>%
    dplyr::filter(total_depth > 0) %>%
    # add baf banding from binomial mixture model
    dplyr::mutate(pos_bin = pos %/% 1e6) %>%
    dplyr::group_by(chr, pos_bin) %>%
    dplyr::transmute(pos, mut_vaf, lower_baf_band = fit_binomial_mixture(mut_depth, total_depth)) %>%
    dplyr::mutate(mut_baf = 1 - mut_vaf,
                  chr = factor(sub("^chr", "", chr), levels = c(as.character(1:22), "X", "Y"))) %>%
    # get alternating colours and chr bounds
    dplyr::group_by(chr) %>%
    dplyr::mutate(chr_alternating = dplyr::cur_group_id() %% 2,
                  min_pos = min(pos), max_pos = max(pos)) %>%
    # plot
    tidyr::pivot_longer(cols = c("mut_vaf", "mut_baf"), names_to = "vaf_type") %>%
    ggplot(aes(x = pos, y = value)) +
    # add alternating coloured chromosomes
    geom_rect(aes(ymin = 0, ymax = 1, xmin = min_pos, xmax = max_pos,
                  fill = as.factor(chr_alternating))) +
    # add line at vaf = 0.5
    geom_hline(yintercept = 0.5, colour = "red") +
    # add baf points
    geom_point(size = 0.005, alpha = 0.01) +
    # add baf bands
    geom_point(aes(x = pos, y = lower_baf_band), colour = "green", size = 0.1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("lightgrey", "white")) +
    ggh4x::facet_grid2(. ~ chr, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA,
                                      linewidth = 0),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0, linetype = "solid"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(title = paste(opts$id, "-", p_source),
         subtitle = "BAF distribution across chromosomes")
}

# plot baf
purrr::map2(names(geno), geno, function(p_source, p_dat) {
  p <- plot_baf(p_dat, p_source)
  ragg::agg_png(paste0(opts$id, "_", p_source, "_baf_plot.png"), width = 2500, height = 600, res = 300)
  print(p)
  dev.off()
})