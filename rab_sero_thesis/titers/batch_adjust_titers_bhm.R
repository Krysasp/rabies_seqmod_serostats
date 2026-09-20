library(tidyverse)
library(ggplot2)
library(readxl)
library(openxlsx)

source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

batch_adjust_titers_bhm <- function(file_path = "RABV_VNT_DATA.xlsx",
                                    sheet_opt = 3,
                                    sheet_samples = 4,
                                    mab_refs = c("MA1-7079_1:20", "Rab50_1:20"),
                                    output_dir = "batch_adjustment_output"
                                    ) {
  
  `%||%` <- function(a, b) if(!is.null(a)) a else b
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  raw3 <- read_excel(file_path, sheet = sheet_opt)
  raw3 <- as.data.frame(raw3)
  
  if("ana" %in% colnames(raw3)) {
    raw3$ana <- as.character(raw3$ana)
    raw3 <- raw3[!is.na(raw3$ana) & tolower(raw3$ana) == "y", ]
  }
  raw3$experiment <- as.character(raw3$experiment)
  raw3$sample <- as.character(raw3$sample)
  
  for(col in c("log_dil", "RLU_dup1", "RLU_dup2", "RLU_dup3", "avg")) {
    if(col %in% colnames(raw3)) raw3[[col]] <- as.numeric(as.character(raw3[[col]]))
  }
  
  raw3$avg_rlu <- rowMeans(raw3[, c("RLU_dup2", "RLU_dup3")], na.rm = TRUE)
  raw3 <- raw3[is.finite(raw3$avg_rlu) & !is.na(raw3$avg_rlu), ]
  
  mab_data <- raw3[raw3$sample %in% mab_refs, ]
  
  batches <- sort(unique(mab_data$experiment))
  
  
  
  batch_params <- list()
  for(batch in batches) {
    for(mab in mab_refs) {
      sub <- mab_data[mab_data$experiment == batch & mab_data$sample == mab, ]
      sub <- sub[order(sub$log_dil), ]
      if(nrow(sub) < 4) next
      
      fit <- fit_4pl_batch(sub$log_dil, sub$avg_rlu)
      if(fit$success) {
        batch_params[[paste(batch, mab, sep = "|")]] <- c(
          experiment = batch, mAb = mab,
          A = as.numeric(unname(fit$A)), B = as.numeric(unname(fit$B)),
          C = as.numeric(unname(fit$C)), D = as.numeric(unname(fit$D)),
          r_squared = as.numeric(unname(fit$r_squared)), method = fit$method
        )
      }
    }
  }
  
  params_df <- do.call(rbind, lapply(batch_params, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  for(col in c("A", "B", "C", "D", "r_squared")) params_df[[col]] <- as.numeric(as.character(params_df[[col]]))
  
  if(nrow(params_df) == 0) {
    return(NULL)
  }
  
  
  pop_means <- params_df %>%
    group_by(mAb) %>%
    summarise(
      pop_A = mean(A, na.rm = TRUE),
      pop_B = mean(B, na.rm = TRUE),
      pop_C = mean(C, na.rm = TRUE),
      pop_D = mean(D, na.rm = TRUE),
      .groups = "drop"
    )
  
  params_df <- merge(params_df, pop_means, by = "mAb", all.x = TRUE)
  
  
  params_df$batch_effect_A <- params_df$A - params_df$pop_A
  params_df$batch_effect_B <- params_df$B - params_df$pop_B
  params_df$batch_effect_C <- params_df$C - params_df$pop_C
  params_df$batch_effect_D <- params_df$D - params_df$pop_D
  
  
  for(mab in unique(params_df$mAb)) {
    mab_sub <- params_df$mAb == mab
    params_df$shrink_effect_A[mab_sub] <- shrinkage_fun(params_df$A[mab_sub], mean(params_df$A[mab_sub], na.rm = TRUE))
    params_df$shrink_effect_B[mab_sub] <- shrinkage_fun(params_df$B[mab_sub], mean(params_df$B[mab_sub], na.rm = TRUE))
    params_df$shrink_effect_C[mab_sub] <- shrinkage_fun(params_df$C[mab_sub], mean(params_df$C[mab_sub], na.rm = TRUE))
    params_df$shrink_effect_D[mab_sub] <- shrinkage_fun(params_df$D[mab_sub], mean(params_df$D[mab_sub], na.rm = TRUE))
  }
  
  params_df$adjustment_factor <- exp(params_df$shrink_effect_C)
  params_df$batch_effect <- params_df$batch_effect_C
  
  raw4 <- read_excel(file_path, sheet = sheet_samples)
  raw4 <- as.data.frame(raw4)
  
  batch_adjust_map <- params_df %>%
    dplyr::select(experiment, adjustment_factor, batch_effect, shrink_effect_C, pop_C, C) %>%
    dplyr::rename(adj_factor = adjustment_factor, b_effect = batch_effect, 
                  shrink_C = shrink_effect_C, pop_C = pop_C, batch_C = C) %>%
    dplyr::distinct()
  
  raw4$experiment_match <- as.character(raw4$E_plate %||% raw4$rep_plate %||% NA)
  
  sample_adjusted <- raw4 %>%
    dplyr::mutate(
      original_titer = as.numeric(`F4PL/5PLRef2_titer...5`),
      experiment_adj = experiment_match
    ) %>%
    dplyr::left_join(batch_adjust_map, by = c("experiment_adj" = "experiment")) %>%
    dplyr::mutate(
      adjustment_factor = ifelse(is.na(adj_factor), 1, adj_factor),
      batch_effect = ifelse(is.na(b_effect), 0, b_effect),
      adjusted_titer = original_titer * adjustment_factor
    )
  
  output_df <- sample_adjusted %>%
    dplyr::select(
      experiment = experiment_adj,
      sample = comb,
      original_titer = original_titer,
      adjusted_titer = adjusted_titer,
      batch_effect = batch_effect,
      adjustment_factor = adjustment_factor
    ) %>%
    dplyr::filter(!is.na(original_titer)) %>%
    dplyr::distinct()
  
  
  
  plot_dir <- output_dir
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  base_theme <- theme_bw() + theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
  
  for(mab in unique(params_df$mAb)) {
    mab_sub <- params_df[params_df$mAb == mab, ]
    
    x_range <- seq(min(mab_sub$C) - 0.5, max(mab_sub$C) + 0.5, length.out = 200)
    pop_A <- unique(mab_sub$pop_A); pop_B <- unique(mab_sub$pop_B)
    pop_C <- unique(mab_sub$pop_C); pop_D <- unique(mab_sub$pop_D)
    
    curve_df <- data.frame()
    for(batch in unique(mab_sub$experiment)) {
      bsub <- mab_sub[mab_sub$experiment == batch, ]
      if(nrow(bsub) == 0) next
      y_pred <- bsub$A + (bsub$D - bsub$A) / (1 + exp(bsub$B * (x_range - bsub$C)))
      curve_df <- rbind(curve_df, data.frame(
        log_dil = x_range, RLU = y_pred, experiment = batch, type = "Batch-specific"
      ))
    }
    y_pop <- pop_A + (pop_D - pop_A) / (1 + exp(pop_B * (x_range - pop_C)))
    curve_df <- rbind(curve_df, data.frame(
      log_dil = x_range, RLU = y_pop, experiment = "Population", type = "Population"
    ))
    
    p <- ggplot(curve_df, aes(x = log_dil, y = RLU, color = experiment, linetype = type)) +
      geom_line(linewidth = 1) +
      labs(title = NULL,
           x = "Log Dilution", y = "RLU") +
      scale_color_manual(values = c(Population = "#000000", sort(unique(as.character(mab_sub$experiment))))) +
      scale_linetype_manual(values = c(Batch_specific = "solid", Population = "dashed")) +
      base_theme + theme(legend.position = "bottom")
    ggsave(file.path(plot_dir, paste0("4pl_curves_", gsub("[^a-zA-Z0-9]", "_", mab), ".png")),
           p, width = 10, height = 6, dpi = 300)
  }
  
  ba_data <- output_df %>%
    dplyr::mutate(
      above_threshold = original_titer >= 0.5,
      adj_above_threshold = adjusted_titer >= 0.5
    ) %>%
    dplyr::filter(!is.na(original_titer) & !is.na(adjusted_titer))
  
  p_ba <- ggplot(ba_data, aes(x = original_titer, y = adjusted_titer)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_x_log10() + scale_y_log10() +
    labs(title = NULL,
         x = "Original Titer (log10)", y = "Adjusted Titer (log10)") +
    base_theme
  ggsave(file.path(plot_dir, "before_after_adjustment.png"), p_ba, width = 8, height = 6, dpi = 300)
  
  effect_df <- params_df %>%
    dplyr::select(experiment, mAb, batch_effect_A, batch_effect_B, batch_effect_C, batch_effect_D) %>%
    tidyr::pivot_longer(cols = starts_with("batch_effect"),
                        names_to = "parameter", values_to = "effect") %>%
    dplyr::mutate(parameter = gsub("batch_effect_", "", parameter))
  
  p_effect <- ggplot(effect_df, aes(x = experiment, y = effect, fill = parameter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = NULL,
         x = "Experiment (Batch)", y = "Batch Effect (batch - population)") +
    scale_fill_brewer(palette = "Set1") +
    base_theme + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plot_dir, "batch_effect_sizes.png"), p_effect, width = 12, height = 6, dpi = 300)
  
  forest_df <- output_df %>%
    dplyr::group_by(experiment) %>%
    dplyr::summarise(
      mean_orig = mean(original_titer, na.rm = TRUE),
      mean_adj = mean(adjusted_titer, na.rm = TRUE),
      sd_orig = sd(original_titer, na.rm = TRUE),
      sd_adj = sd(adjusted_titer, na.rm = TRUE),
      n = n(),
      se_orig = sd_orig / sqrt(n),
      se_adj = sd_adj / sqrt(n),
      ci_lower_orig = mean_orig - 1.96 * se_orig,
      ci_upper_orig = mean_orig + 1.96 * se_orig,
      ci_lower_adj = mean_adj - 1.96 * se_adj,
      ci_upper_adj = mean_adj + 1.96 * se_adj,
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c(mean_orig, mean_adj, ci_lower_orig, ci_upper_orig, ci_lower_adj, ci_upper_adj),
                        names_to = c("stat", "type"), names_pattern = "(.*)_(.*)") %>%
    tidyr::pivot_wider(names_from = stat, values_from = value)
  
  p_forest <- ggplot(forest_df, aes(x = experiment, y = mean, color = type)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.5), width = 0.2) +
    scale_y_log10() +
    labs(title = NULL,
         x = "Experiment (Batch)", y = "Mean Titer (log10)") +
    scale_color_manual(values = c(orig = "#D55E00", adj = "#009E73"),
                       labels = c("Original", "Adjusted")) +
    base_theme + theme(legend.position = "bottom")
  ggsave(file.path(plot_dir, "forest_plot_adjusted_titers.png"), p_forest, width = 10, height = 6, dpi = 300)
  
  
  export_list <- list(
    Batch_Parameters = params_df %>%
      dplyr::select(experiment, mAb, A, B, C, D, r_squared, method,
                    pop_A, pop_B, pop_C, pop_D,
                    batch_effect_A, batch_effect_B, batch_effect_C, batch_effect_D,
                    shrink_effect_A, shrink_effect_B, shrink_effect_C, shrink_effect_D,
                    adjustment_factor, batch_effect),
    Adjusted_Titers = output_df,
    Batch_Effects = params_df %>%
      dplyr::select(experiment, mAb, C, pop_C, batch_effect_C, shrink_effect_C, adjustment_factor) %>%
      dplyr::rename(batch_C = C, population_C = pop_C),
    Summary = output_df %>%
      dplyr::group_by(experiment) %>%
      dplyr::summarise(
        n_samples = n(),
        mean_original = mean(original_titer, na.rm = TRUE),
        mean_adjusted = mean(adjusted_titer, na.rm = TRUE),
        sd_original = sd(original_titer, na.rm = TRUE),
        sd_adjusted = sd(adjusted_titer, na.rm = TRUE),
        mean_batch_effect = mean(batch_effect, na.rm = TRUE),
        mean_adj_factor = mean(adjustment_factor, na.rm = TRUE),
        .groups = "drop"
      )
  )
  
  output_file <- file.path(output_dir, "batch_adjustment_results.xlsx")
  write.xlsx(export_list, output_file)
  
  
  return(list(
    batch_parameters = params_df,
    adjusted_titers = output_df,
    batch_effects = params_df %>%
      dplyr::select(experiment, mAb, C, pop_C, batch_effect_C, shrink_effect_C, adjustment_factor),
    summary = export_list$Summary,
    output_file = output_file
  ))
}

