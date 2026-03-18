# Neutral Model Pipeline
# Date: 2026-02-24

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(purrr)
library(grDevices)

# Configuration
# Assign your actual phyloseq object here before running
physeq_input <- your_phyloseq_object 
metadata_col <- "island"
presence_th <- 1e-6
det_lim_default <- 1e-6
alpha_level <- 0.05

# Setup Output Directory
out_dir <- paste0("Neutral_Model_Outputs_2026_02_24")
dir.create(out_dir, showWarnings = FALSE)

# Helper Functions
expected_detection_freq <- function(M, p, det_lim) {
  a <- p * M
  b <- (1 - p) / p * a
  # Survival function (1 - CDF)
  pbeta(det_lim, a, b, lower.tail = FALSE) 
}

objective_func <- function(M, p, det_lim, F_obs) {
  sum((F_obs - expected_detection_freq(M, p, det_lim))^2)
}

R2_func <- function(obs, fit) {
  1 - sum((obs - fit)^2) / sum((obs - mean(obs))^2)
}

estimate_M <- function(p, F_obs, det_lim) {
  opt <- optimize(objective_func, interval = c(1e-5, 1e5), p = p, det_lim = det_lim, F_obs = F_obs)
  F_fit <- expected_detection_freq(opt$minimum, p, det_lim)
  list(M = opt$minimum, F_fit = F_fit, R2 = R2_func(F_obs, F_fit))
}

compute_neutral_stats <- function(otu_mat, det_lim, presence_th, alpha) {
  n_samples <- ncol(otu_mat)
  
  k_vec <- rowSums(otu_mat > presence_th)
  F_obs <- k_vec / n_samples
  
  # Relative abundance
  otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), "/")
  otu_rel[is.na(otu_rel)] <- 0
  p_vec <- rowMeans(otu_rel)
  
  keep <- p_vec > 0
  k_vec <- k_vec[keep]
  F_obs <- F_obs[keep]
  p_vec <- p_vec[keep]
  
  fit <- estimate_M(p_vec, F_obs, det_lim)
  
  # Binomial test for classification
  p_vals <- mapply(function(k, p_e) {
    # Constrain probability between 0 and 1 for binom.test
    p_e <- min(max(p_e, 1e-10), 1 - 1e-10)
    binom.test(k, n_samples, p_e)$p.value
  }, k_vec, fit$F_fit)
  
  class <- ifelse(F_obs > fit$F_fit & p_vals < alpha, "Above",
                  ifelse(F_obs < fit$F_fit & p_vals < alpha, "Below", "Neutral"))
  
  tibble(
    Taxon = names(p_vec), 
    MeanFraction = p_vec,
    DetectionFrequency = F_obs, 
    FittedFrequency = fit$F_fit,
    Class = class, 
    p_value = p_vals,
    Residual = F_obs - fit$F_fit,
    M_hat = fit$M,
    R2 = fit$R2
  )
}

plot_neutral_model <- function(stats_df, title_label) {
  label_txt <- sprintf("M = %.2f\nR² = %.3f", unique(stats_df$M_hat)[1], unique(stats_df$R2)[1])
  
  # Pastel color palette for publication
  class_colors <- c("Above" = "#FBB4AE", "Neutral" = "#B3CDE3", "Below" = "#CCEBC5")
  
  ggplot(stats_df, aes(x = MeanFraction, y = DetectionFrequency)) +
    geom_point(aes(fill = Class), shape = 21, color = "grey30", size = 2, alpha = 0.8) +
    geom_line(aes(y = FittedFrequency), color = "grey40", linewidth = 1, linetype = "dashed") +
    scale_fill_manual(values = class_colors) +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotate("text", x = min(stats_df$MeanFraction) * 2, y = 0.9, 
             label = label_txt, hjust = 0, size = 4) +
    labs(title = title_label, x = "Mean Relative Abundance (log10)", y = "Detection Frequency") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# Data Preparation
otu_mat_full <- as(otu_table(physeq_input), "matrix") 
if(taxa_are_rows(physeq_input)) {
  # We want taxa as rows for the math below
} else {
  otu_mat_full <- t(otu_mat_full)
}

meta_df <- as(sample_data(physeq_input), "data.frame") %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(across(all_of(metadata_col), as.character))

islands <- sort(unique(meta_df[[metadata_col]]))

# Run Loop Per Group
all_stats <- list()

for (isl in islands) {
  
  samp_ids <- meta_df %>% 
    filter(.data[[metadata_col]] == isl) %>% 
    pull(sample_id)
  
  sub_otu <- otu_mat_full[, samp_ids, drop = FALSE]
  sub_otu <- sub_otu[rowSums(sub_otu) > 0, colSums(sub_otu) > 0, drop = FALSE]
  
  if (ncol(sub_otu) < 2) {
    warning(paste("Skipping", isl, "due to insufficient samples."))
    next
  }
  
  # Calculate dynamic detection limit or fall back to default
  otu_rel_tmp <- sweep(sub_otu, 2, colSums(sub_otu), "/")
  det_lim_calc <- min(otu_rel_tmp[otu_rel_tmp > 0]) / 10
  current_det_lim <- min(det_lim_calc, det_lim_default)
  
  # Compute stats
  stats_df <- compute_neutral_stats(sub_otu, current_det_lim, presence_th, alpha_level)
  stats_df$Group <- isl
  all_stats[[isl]] <- stats_df
  
  # Generate Plot
  p_plot <- plot_neutral_model(stats_df, isl)
  
  # Save Outputs
  csv_name <- file.path(out_dir, paste0("Neutral_Stats_", isl, "_2026_02_24.csv"))
  write_csv(stats_df, csv_name)
  
  pdf_name <- file.path(out_dir, paste0("Neutral_Plot_", isl, "_2026_02_24.pdf"))
  png_name <- file.path(out_dir, paste0("Neutral_Plot_", isl, "_2026_02_24.png"))
  
  # Use cairo_pdf so text is editable in Inkscape
  ggsave(pdf_name, plot = p_plot, device = cairo_pdf, width = 6, height = 5, units = "in")
  ggsave(png_name, plot = p_plot, width = 6, height = 5, units = "in", dpi = 300)
}

# Overall Model (All Samples Combined)
overall_otu <- otu_mat_full[rowSums(otu_mat_full) > 0, colSums(otu_mat_full) > 0, drop = FALSE]
otu_rel_all <- sweep(overall_otu, 2, colSums(overall_otu), "/")
det_lim_all <- min(otu_rel_all[otu_rel_all > 0]) / 10

overall_stats <- compute_neutral_stats(overall_otu, min(det_lim_all, det_lim_default), presence_th, alpha_level)
overall_stats$Group <- "Overall"

overall_plot <- plot_neutral_model(overall_stats, "Overall Neutral Model")

write_csv(overall_stats, file.path(out_dir, "Neutral_Stats_Overall_2026_02_24.csv"))
ggsave(file.path(out_dir, "Neutral_Plot_Overall_2026_02_24.pdf"), plot = overall_plot, device = cairo_pdf, width = 6, height = 5, units = "in")
ggsave(file.path(out_dir, "Neutral_Plot_Overall_2026_02_24.png"), plot = overall_plot, width = 6, height = 5, units = "in", dpi = 300)

# Combine all stats into one master file
master_stats <- bind_rows(all_stats) %>% bind_rows(overall_stats)
write_csv(master_stats, file.path(out_dir, "Neutral_Stats_Master_2026_02_24.csv"))