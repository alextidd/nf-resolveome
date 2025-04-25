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

# read geno
geno <- readr::read_tsv(opts$geno)

# function: fit a binomial mixture model
fit_binomial_mixture <- function(alt, cov) {
  data <- data.frame(alt = alt, cov = cov)
  loglik <- function(p1, d) sum(log(0.5 * dbinom(d$alt, d$cov, p1) + 0.5 * dbinom(d$alt, d$cov, 1 - p1)))
  opt <- optimize(function(p) loglik(p, data), interval = c(0, 1), maximum = TRUE)
  out <- list(prob1 = opt$maximum, prob2 = 1 - opt$maximum, mix_prob = 0.5, loglik = opt$objective)
  # return the minimum
  min(out$prob1, out$prob2)
}

# function: compare mixed and single binomial models
fit_best_binomial_model <- function(alt, cov) {
  data <- data.frame(alt = alt, cov = cov)

  # Single binomial model
  p_single <- sum(alt) / sum(cov)
  loglik_single <- sum(dbinom(alt, cov, p_single, log = TRUE))
  aic_single <- -2 * loglik_single + 2  # 1 parameter

  # Binomial mixture model (symmetric, equal mixing proportions)
  loglik_mix <- function(p1, d) {
    sum(log(0.5 * dbinom(d$alt, d$cov, p1) + 0.5 * dbinom(d$alt, d$cov, 1 - p1)))
  }
  opt <- optimize(function(p) loglik_mix(p, data), interval = c(0, 1), maximum = TRUE)
  p_mix1 <- opt$maximum
  p_mix2 <- 1 - p_mix1
  loglik_mixture <- opt$objective
  aic_mixture <- -2 * loglik_mixture + 2  # also 1 parameter (symmetry assumed)

  # Compare models
  if (aic_single <= aic_mixture) {
    return(list(af = p_single, model = "single"))
  } else {
    return(list(af = min(p_mix1, p_mix2), model = "mixture"))
  }
}

# function: plot BAF
plot_baf <- function(p_dat, p_source, incl_binom = TRUE,
                     p_alpha = 0.05, p_size = 0.01) {
  # prep data
  p_dat2 <-
    p_dat %>%
    dplyr::filter(total_depth > 20) %>%
    # compare binom single and mixture models per band, fit best model
    dplyr::mutate(pos_bin = pos %/% 1e6) %>%
    dplyr::group_by(chr, pos_bin) %>%
    dplyr::transmute(
      pos, mut_vaf,
      binom_vaf = fit_best_binomial_model(mut_depth, total_depth)$af,
      binom_model = fit_best_binomial_model(mut_depth, total_depth)$model) %>%
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

  # optionally include best binomial model points
  if (incl_binom) {
    p <-
      p +
      geom_point(aes(x = pos, y = binom_vaf, colour = binom_model),
                 size = 0.1) +
      scale_colour_manual(values = c("mixture" = "blue", "single" = "green")) +
      theme(legend.position = "bottom")
  }

  # return
  p
}

# plot
p <- plot_baf(geno, "caveman_snps")
ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_w_binom_plot.png"), width = 5000, height = 1200, res = 300)
print(p)
dev.off()

p <- plot_baf(geno, "caveman_snps", incl_binom = FALSE)
ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_plot.png"), width = 5000, height = 1200, res = 300)
print(p)
dev.off()

for (chr_i in c("1", "9")) {
  p <- plot_baf(geno %>% dplyr::filter(chr == chr_i), "caveman_snps", p_alpha = 0.2, p_size = 0.8)
  ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_chr", chr_i, "_w_binom_plot.png"), width = 5000, height = 1200, res = 300)
  print(p)
  dev.off()

  p <- plot_baf(geno %>% dplyr::filter(chr == chr_i), "caveman_snps", incl_binom = FALSE, p_alpha = 0.2, p_size = 0.8)
  ragg::agg_png(paste0(opts$id, "_caveman_snps_baf_chr", chr_i, "_plot.png"), width = 5000, height = 1200, res = 300)
  print(p)
  dev.off()
}