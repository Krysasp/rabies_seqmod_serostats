library(tidyverse)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(readxl)
library(openxlsx)

source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

run_sero_neutana <- function(
    neut_file = "neutres",
    meta_file = "input_excel.xlsx",
    sheet_meta = 4,
    ref_choices = "R1",
    elisa_types = "linear",
    output_dir = "sero_neutana",
    exclude_groups = c("unassign", "undesc"),
    low_titer_cutoff = 0.0625,
    bdl = 0.125,
    replicates = "2",
    method = "NULL"
) {
  show_gam <- TRUE
  show_gridlines <- FALSE
  y_scale <- "log2"
  log2_breaks <- c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
  log2_labels <- c("0.0625", "0.125", "0.25", "0.5", "1.0", "2.0", "4.0", "8.0", "16.0", "32.0")
  color_by <- "NULL"
  longitudinal_min_samples_per_comb <- 3
  base_size <- 9
  legend_position <- "NULL"
  legend_ncol <- 2
  longitudinal_font_size <- 11
  scatter_font_size <- 11
  boxplot_x_text_size <- 10
  boxplot_x_text_angle <- 0
  stacked_bar_base_size <- 11
  stacked_bar_title_size <- 12.5
  stacked_bar_axis_title_size <- 12
  stacked_bar_axis_text_size <- 10

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plot_dir <- file.path(output_dir, "plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  .geo_mean <- function(x) {
    x <- x[x > 0 & is.finite(x)]
    if(length(x) == 0) return(NA)
    exp(mean(log(x), na.rm = TRUE))
  }

  .geo_sd_factor <- function(x) {
    x <- x[x > 0 & is.finite(x)]
    if(length(x) < 2) return(NA)
    exp(sd(log(x), na.rm = TRUE))
  }

  data_source <- .resolve_data_source(neut_file, replicates, method)
  merged <- .load_and_merge_data(data_source, neut_file, meta_file, sheet_meta, sample_id_normalize, meta_sample_col, filter_discordant, discordant_filter_type)
  overall_df <- .expand_ref_types(merged, R1_name, R2_name)
  summary_tables <- .compute_summary_tables(overall_df, output_dir, .geo_mean, .geo_sd_factor)
  sero_outcome_df <- .compute_serocontainment_outcomes(overall_df, output_dir, summary_tables$month_group_levels)
  plot_themes <- .setup_plot_themes(show_gridlines, base_size)
  all_plots <- list()
  plot_data_list <- .prepare_all_plot_data(all_data = list(Default = overall_df), ref_choices, elisa_types, ref_name_map, exclude_groups, bdl, low_titer_cutoff, month_group_levels = summary_tables$month_group_levels)
  longitudinal_plots <- .generate_longitudinal_plots(plot_data_list, ref_name_map, elisa_types, exclude_groups, show_gam, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label = paste0(method, "_", replicates, "rep"), longitudinal_min_samples_per_comb, longitudinal_width, longitudinal_height, longitudinal_font_size, base_size, legend_position, legend_ncol, legend_key_size, legend_text_size, legend_title_size, color_by, dot_size, dot_alpha, dot_color, point_size, point_alpha, gm_one_panel_width, gm_one_panel_height, low_titer_cutoff, bdl)
  boxplot_plots <- .generate_boxplots(plot_data_list, ref_name_map, elisa_types, exclude_groups, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label = paste0(method, "_", replicates, "rep"), boxplot_x_text_size, boxplot_x_text_angle, low_titer_cutoff, bdl)
  scatter_plots <- .generate_scatter_plots(plot_data_list, ref_name_map, elisa_types, exclude_groups, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label = paste0(method, "_", replicates, "rep"), gm_width, gm_height, gm_jitter_width, gm_two_panel_width, gm_two_panel_height, scatter_font_size, low_titer_cutoff, bdl, .geo_mean, .geo_sd_factor, month_group_levels = summary_tables$month_group_levels)
  concordance_plot <- .generate_concordance_plot(overall_df, plot_themes, plot_dir, plot_width_one_panel, plot_height_one_panel)
  stacked_bar_plots <- .generate_stacked_bar_plots(overall_df, ref_name_map, elisa_types, exclude_groups, output_dir, plot_dir, method_label = paste0(method, "_", replicates, "rep"), stacked_bar_nrow, stacked_bar_ncol, stacked_bar_two_panel_width, stacked_bar_two_panel_height, stacked_bar_sd_label_size, stacked_bar_gm_label_size, stacked_bar_base_size, stacked_bar_title_size, stacked_bar_axis_title_size, stacked_bar_axis_text_size, bdl, month_group_levels = summary_tables$month_group_levels)
   all_plots <- c(longitudinal_plots, boxplot_plots, scatter_plots, concordance_plot, stacked_bar_plots)

   return(list(
    merged_data = overall_df,
    summary_by_month = summary_tables$summary_by_month,
    summary_by_month_bimonthly = summary_tables$summary_by_month_bimonthly,
    summary_by_group = summary_tables$summary_by_group,
    stacked_bar_data = stacked_bar_plots$stacked_bar_data,
    plots = all_plots
  ))
}

.resolve_data_source <- function(neut_file, replicates, method) {
  excel_sheet_map <- list(
    list(replicates = "2", method = "AICMix", sheet = "Sample_Results_AICMix", label = "AICMix 2-rep"),
    list(replicates = "2", method = "Forced4PL", sheet = "Sample_Results_Forced4PL", label = "Forced 4PL 2-rep"),
    list(replicates = "2", method = "5PL", sheet = "Sample_Results_5PL", label = "5PL 2-rep"),
    list(replicates = "2", method = "Bayes-4PL", sheet = "Sample_Results_Bayes-4PL", label = "Bayes-4PL 2-rep"),
    list(replicates = "3", method = "AICMix", sheet = "Sample_Results_AICMix_3rep", label = "AICMix 3-rep"),
    list(replicates = "3", method = "Forced4PL", sheet = "Sample_Results_Forced4PL_3rep", label = "Forced 4PL 3-rep"),
    list(replicates = "3", method = "5PL", sheet = "Sample_Results_5PL_3rep", label = "5PL 3-rep"),
    list(replicates = "3", method = "Bayes-4PL", sheet = "Sample_Results_Bayes-4PL_3rep", label = "Bayes-4PL 3-rep")
  )
  csv_file_map <- list(
    list(replicates = "2", method = "AICMix", file = "sample_results.csv", label = "AICMix 2-rep"),
    list(replicates = "2", method = "Forced4PL", file = "sample_results_forced4pl.csv", label = "Forced 4PL 2-rep"),
    list(replicates = "2", method = "5PL", file = "sample_results_5pl.csv", label = "5PL 2-rep"),
    list(replicates = "2", method = "Bayes-4PL", file = "sample_results_bayesian.csv", label = "Bayes-4PL 2-rep"),
    list(replicates = "3", method = "AICMix", file = "sample_results_three_replicates.csv", label = "AICMix 3-rep"),
    list(replicates = "3", method = "Forced4PL", file = "sample_results_three_replicates_forced4pl.csv", label = "Forced 4PL 3-rep"),
    list(replicates = "3", method = "5PL", file = "sample_results_three_replicates_5pl.csv", label = "5PL 3-rep"),
    list(replicates = "3", method = "Bayes-4PL", file = "sample_results_three_replicates_bayesian.csv", label = "Bayes-4PL 3-rep")
  )
  selected_sheet <- NULL
  selected_file <- NULL
  selected_label <- NULL
  use_excel_source <- FALSE
  if(grepl("\\.xlsx?$", neut_file, ignore.case = TRUE) && file.exists(neut_file)) {
    for(entry in excel_sheet_map) {
      if(entry$replicates == replicates && entry$method == method) {
        selected_sheet <- entry$sheet
        selected_label <- entry$label
        use_excel_source <- TRUE
        break
      }
    }
  }
  if(!use_excel_source) {
    for(entry in csv_file_map) {
      if(entry$replicates == replicates && entry$method == method) {
        selected_file <- entry$file
        selected_label <- entry$label
        break
      }
    }
  }
  if(is.null(selected_file) && is.null(selected_sheet)) {
    stop("Invalid combination: replicates='", replicates, "', method='", method, "'. Valid options: replicates=c('2','3'), method=c('AICMix','Forced4PL','5PL','Bayes-4PL')")
  }
  list(selected_sheet = selected_sheet, selected_file = selected_file, selected_label = selected_label, use_excel_source = use_excel_source)
}

.load_and_merge_data <- function(neut_file, meta_file, sheet_meta, data_source, sample_id_normalize, meta_sample_col, filter_discordant, discordant_filter_type) {
  meta_raw <- read_excel(meta_file, sheet = sheet_meta)
  meta_df <- as.data.frame(meta_raw)
  if(data_source$use_excel_source) {
    if(!file.exists(neut_file)) stop("Excel results file not found: ", neut_file, ". Run analyze_neutralization_statistics first.")
    neut_raw <- read_excel(neut_file, sheet = data_source$selected_sheet)
    neut_df <- as.data.frame(neut_raw)
  } else {
    if(!file.exists(data_source$selected_file)) stop("Selected data file not found: ", data_source$selected_file, ". Run analyze_neutralization_statistics first.")
    neut_df <- read.csv(data_source$selected_file, stringsAsFactors = FALSE)
  }
  if(filter_discordant) {
    if(discordant_filter_type %in% c("linear", "both")) {
      neut_df <- neut_df[neut_df$Concordance_R1_ELISA != "Discordant" & neut_df$Concordance_R2_ELISA != "Discordant" & (!("Concordance_Ref3_ELISA" %in% colnames(neut_df)) | neut_df$Concordance_Ref3_ELISA != "Discordant"), ]
    }
    if(discordant_filter_type %in% c("nonlin", "both")) {
      neut_df <- neut_df[neut_df$Concordance_R1_Nonlin != "Discordant" & neut_df$Concordance_R2_Nonlin != "Discordant" & (!("Concordance_Ref3_Nonlin" %in% colnames(neut_df)) | neut_df$Concordance_Ref3_Nonlin != "Discordant"), ]
    }
  }
  if(is.null(meta_sample_col)) {
    if("sample_id" %in% colnames(meta_df)) meta_sample_col <- "sample_id"
    else if("sample" %in% colnames(meta_df)) meta_sample_col <- "sample"
    else stop("Column 'sample_id' or 'sample' not found in metadata sheet ", sheet_meta)
  }
  if(!meta_sample_col %in% colnames(meta_df)) stop("Column '", meta_sample_col, "' not found in metadata sheet ", sheet_meta)
  direct_merged <- merge(neut_df, meta_df, by.x = "Sample", by.y = meta_sample_col, all.x = FALSE, all.y = FALSE)
  if(nrow(direct_merged) > 0) {
    merged <- direct_merged
  } else if(is.function(sample_id_normalize)) {
    neut_df[['Sample_norm']] <- sample_id_normalize(neut_df$Sample)
    meta_df[['meta_norm']] <- sample_id_normalize(meta_df[[meta_sample_col]])
    norm_merged <- merge(neut_df, meta_df, by.x = "Sample_norm", by.y = "meta_norm", all.x = FALSE, all.y = FALSE)
    if(nrow(norm_merged) > 0) {
      merged <- norm_merged
      merged$Sample_norm <- NULL
      merged$meta_norm <- NULL
    } else {
      merged <- merge(neut_df, meta_df, by.x = "Sample", by.y = meta_sample_col, all.x = TRUE, all.y = FALSE)
      warning("No matching sample IDs found between neutralization results and metadata. Proceeding with partial data.")
    }
  } else {
    merged <- merge(neut_df, meta_df, by.x = "Sample", by.y = meta_sample_col, all.x = TRUE, all.y = FALSE)
    warning("No matching sample IDs found between neutralization results and metadata. Proceeding with partial data.")
  }
  list(merged = merged, meta_df = meta_df, neut_df = neut_df)
}

.expand_ref_types <- function(merged, R1_name, R2_name) {
  base_cols <- c("Sample", "Experiment_Group", "month", "dose_type", "conc_vax", "set_type", "age_G2", "comb", "groupW", "ELISA_EU_mL", "ELISA_EU_mL_nonlin", "IC50_log_dilution", "Protection_Level", "IC50_Method", "R_squared")
  R1_cols <- c(base_cols, "Titer_R1", "Concordance_R1_ELISA", "Concordance_R1_Nonlin")
  R1_cols <- R1_cols[R1_cols %in% colnames(merged)]
  R1_df <- merged[, R1_cols, drop = FALSE]
  colnames(R1_df)[colnames(R1_df) == "Titer_R1"] <- "Titer"
  colnames(R1_df)[colnames(R1_df) == "Concordance_R1_ELISA"] <- "Concordance_ELISA"
  colnames(R1_df)[colnames(R1_df) == "Concordance_R1_Nonlin"] <- "Concordance_Nonlin"
  colnames(R1_df)[colnames(R1_df) == "ELISA_EU_mL"] <- "ELISA"
  colnames(R1_df)[colnames(R1_df) == "ELISA_EU_mL_nonlin"] <- "ELISA_nonlin"
  R1_df$Ref_Type <- R1_name
  R1_df$Titer_Unit <- "IU/mL"
  R2_cols <- c(base_cols, "Titer_R2", "Concordance_R2_ELISA", "Concordance_R2_Nonlin")
  R2_cols <- R2_cols[R2_cols %in% colnames(merged)]
  R2_df <- merged[, R2_cols, drop = FALSE]
  colnames(R2_df)[colnames(R2_df) == "Titer_R2"] <- "Titer"
  colnames(R2_df)[colnames(R2_df) == "Concordance_R2_ELISA"] <- "Concordance_ELISA"
  colnames(R2_df)[colnames(R2_df) == "Concordance_R2_Nonlin"] <- "Concordance_Nonlin"
  colnames(R2_df)[colnames(R2_df) == "ELISA_EU_mL"] <- "ELISA"
  colnames(R2_df)[colnames(R2_df) == "ELISA_EU_mL_nonlin"] <- "ELISA_nonlin"
  R2_df$Ref_Type <- R2_name
  R2_df$Titer_Unit <- "EU/mL"
  ref3_cols <- c(base_cols, "Titer_Ref3", "Concordance_Ref3_ELISA", "Concordance_Ref3_Nonlin")
  ref3_cols <- ref3_cols[ref3_cols %in% colnames(merged)]
  ref3_df <- merged[, ref3_cols, drop = FALSE]
  colnames(ref3_df)[colnames(ref3_df) == "Titer_Ref3"] <- "Titer"
  colnames(ref3_df)[colnames(ref3_df) == "Concordance_Ref3_ELISA"] <- "Concordance_ELISA"
  colnames(ref3_df)[colnames(ref3_df) == "Concordance_Ref3_Nonlin"] <- "Concordance_Nonlin"
  colnames(ref3_df)[colnames(ref3_df) == "ELISA_EU_mL"] <- "ELISA"
  colnames(ref3_df)[colnames(ref3_df) == "ELISA_EU_mL_nonlin"] <- "ELISA_nonlin"
  ref3_df$Ref_Type <- ifelse("Ref3_Name" %in% colnames(merged) && length(unique(na.omit(merged$Ref3_Name))) > 0 && unique(na.omit(merged$Ref3_Name))[1] != "", unique(na.omit(merged$Ref3_Name))[1], "Ref3")
  ref3_df$Titer_Unit <- "EU/mL"
  mab_cols <- c(base_cols, "Titer_mAb", "Concordance_R1_ELISA", "Concordance_R1_Nonlin")
  mab_cols <- mab_cols[mab_cols %in% colnames(merged)]
  mab_df <- merged[, mab_cols, drop = FALSE]
  colnames(mab_df)[colnames(mab_df) == "Titer_mAb"] <- "Titer"
  colnames(mab_df)[colnames(mab_df) == "Concordance_R1_ELISA"] <- "Concordance_ELISA"
  colnames(mab_df)[colnames(mab_df) == "Concordance_R1_Nonlin"] <- "Concordance_Nonlin"
  colnames(mab_df)[colnames(mab_df) == "ELISA_EU_mL"] <- "ELISA"
  colnames(mab_df)[colnames(mab_df) == "ELISA_EU_mL_nonlin"] <- "ELISA_nonlin"
  mab_df$Ref_Type <- "mAb"
  mab_df$Titer_Unit <- "IU/mL"
  merged_long <- rbind(R1_df, R2_df, ref3_df, mab_df)
  overall_df <- merged_long
  overall_df$Concordance <- overall_df$Concordance_ELISA
  month_group_levels <- MONTH_GROUP_LEVELS
  overall_df$month <- as.numeric(overall_df$month)
  overall_df$month_group <- bucket_months(overall_df$month, month_group_levels)
  list(overall_df = overall_df, month_group_levels = month_group_levels)
}

.compute_summary_tables <- function(overall_df, output_dir, geo_mean, geo_sd_factor) {
  summary_by_month <- overall_df %>%
    group_by(month, Ref_Type) %>%
    summarise(n = n(), n_dogs = n_distinct(comb), GM_ELISA = geo_mean(ELISA), GM_Titer = geo_mean(Titer), GSD_ELISA = geo_sd_factor(ELISA), GSD_Titer = geo_sd_factor(Titer), .groups = "drop") %>%
    mutate(month = as.numeric(month), GM_ELISA_rounded = round(GM_ELISA, 3), GM_Titer_rounded = round(GM_Titer, 3), GSD_ELISA_rounded = round(GSD_ELISA, 3), GSD_Titer_rounded = round(GSD_Titer, 3))
  summary_by_group <- overall_df %>%
    group_by(groupW, Ref_Type) %>%
    summarise(n = n(), n_dogs = n_distinct(comb), GM_ELISA = geo_mean(ELISA), GM_Titer = geo_mean(Titer), GSD_ELISA = geo_sd_factor(ELISA), GSD_Titer = geo_sd_factor(Titer), .groups = "drop") %>%
    mutate(GM_ELISA_rounded = round(GM_ELISA, 3), GM_Titer_rounded = round(GM_Titer, 3), GSD_ELISA_rounded = round(GSD_ELISA, 3), GSD_Titer_rounded = round(GSD_Titer, 3))
  write.csv(overall_df, file.path(output_dir, "merged_data.csv"), row.names = FALSE)
  month_group_levels <- MONTH_GROUP_LEVELS
  summary_by_month_bimonthly <- overall_df %>%
    group_by(month_group, Ref_Type) %>%
    summarise(n = n(), n_dogs = n_distinct(comb), GM_ELISA = geo_mean(ELISA), GM_Titer = geo_mean(Titer), GSD_ELISA = geo_sd_factor(ELISA), GSD_Titer = geo_sd_factor(Titer), .groups = "drop") %>%
    mutate(month_group = factor(month_group, levels = month_group_levels), GM_ELISA_rounded = round(GM_ELISA, 3), GM_Titer_rounded = round(GM_Titer, 3), GSD_ELISA_rounded = round(GSD_ELISA, 3), GSD_Titer_rounded = round(GSD_Titer, 3))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Monthly")
  openxlsx::addWorksheet(wb, "Bimonthly")
  openxlsx::addWorksheet(wb, "Summary_by_Group")
  openxlsx::writeData(wb, "Monthly", summary_by_month)
  openxlsx::writeData(wb, "Bimonthly", summary_by_month_bimonthly)
  openxlsx::writeData(wb, "Summary_by_Group", summary_by_group)
  openxlsx::saveWorkbook(wb, file.path(output_dir, "summary_by_month.xlsx"), overwrite = TRUE)
  list(summary_by_month = summary_by_month, summary_by_month_bimonthly = summary_by_month_bimonthly, summary_by_group = summary_by_group, month_group_levels = month_group_levels)
}

.compute_serocontainment_outcomes <- function(overall_df, output_dir, month_group_levels) {
  sero_outcome_rows <- list()
  for(sz in c(1, 2)) {
    eligible <- if(sz == 1) overall_df else overall_df %>% group_by(comb) %>% filter(n_distinct(Sample) >= 2) %>% ungroup()
    if(nrow(eligible) == 0) next
    long_df <- eligible %>%
      mutate(month_group = bucket_months(month, month_group_levels)) %>%
      pivot_longer(cols = c(ELISA, Titer), names_to = "Assay", values_to = "Value") %>%
      mutate(category = case_when(is.na(Value) ~ NA_character_, Value >= 0.5 ~ "Had >= 0.5 U/mL", Value > 0.125 ~ "Had < 0.5 U/mL but > BDL", TRUE ~ "IgG level is BDL")) %>%
      filter(!is.na(category))
    if(nrow(long_df) == 0) next
    sero_outcome_rows[[as.character(sz)]] <- long_df %>%
      group_by(month_group, Ref_Type, Assay, sample_sz = sz) %>%
      summarise(n_total = n(), n_ge_0_5 = sum(Value >= 0.5, na.rm = TRUE), pct_ge_0_5 = ifelse(n_total > 0, 100 * sum(Value >= 0.5, na.rm = TRUE) / n_total, NA_real_), n_between = sum(Value > 0.125 & Value < 0.5, na.rm = TRUE), pct_between = ifelse(n_total > 0, 100 * sum(Value > 0.125 & Value < 0.5, na.rm = TRUE) / n_total, NA_real_), n_below_bdl = sum(Value <= 0.125, na.rm = TRUE), pct_below_bdl = ifelse(n_total > 0, 100 * sum(Value <= 0.125, na.rm = TRUE) / n_total, NA_real_), .groups = "drop")
  }
  if(length(sero_outcome_rows) > 0 && all(sapply(sero_outcome_rows, function(x) is.data.frame(x) && nrow(x) > 0))) {
    sero_outcome_df <- do.call(rbind, sero_outcome_rows)
    sero_outcome_df$month_group <- factor(sero_outcome_df$month_group, levels = month_group_levels)
    sero_outcome_df <- sero_outcome_df[order(sero_outcome_df$month_group, sero_outcome_df$Ref_Type, sero_outcome_df$Assay, sero_outcome_df$sample_sz), ]
    if(file.exists(file.path(output_dir, "summary_by_month.xlsx"))) {
      wb_out <- openxlsx::loadWorkbook(file.path(output_dir, "summary_by_month.xlsx"))
    } else {
      wb_out <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb_out, "Monthly")
    }
    if(!"Stacked_Bar_Outcome" %in% names(wb_out)) openxlsx::addWorksheet(wb_out, "Stacked_Bar_Outcome")
    openxlsx::writeData(wb_out, "Stacked_Bar_Outcome", sero_outcome_df)
    openxlsx::saveWorkbook(wb_out, file.path(output_dir, "summary_by_month.xlsx"), overwrite = TRUE)
  } else {
    sero_outcome_df <- NULL
  }
  sero_outcome_df
}

.setup_plot_themes <- function(show_gridlines, base_size) {
  theme_publication <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
      theme(text = element_text(family = "sans", color = "black"), axis.title = element_text(size = base_size + 1, face = "bold"), axis.text = element_text(size = base_size, color = "black"), axis.line = element_line(color = "black", linewidth = 0.5), axis.ticks = element_line(color = "black", linewidth = 0.4), panel.grid.major = element_line(color = "grey90", linewidth = 0.3), panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_text(size = base_size, face = "bold"), legend.text = element_text(size = base_size - 1), plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40"))
  }
  if(show_gridlines) {
    base_theme <- theme_bw(base_size = 10)
  } else {
    base_theme <- theme_minimal(base_size = 10) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5), axis.line = element_line(color = "black", linewidth = 0.5), axis.ticks = element_line(color = "black", linewidth = 0.4))
  }
  list(theme_publication = theme_publication, base_theme = base_theme)
}

.prepare_all_plot_data <- function(all_data, ref_choices, elisa_types, ref_name_map, exclude_groups, bdl, low_titer_cutoff, month_group_levels) {
  plot_data_list <- list()
  for(ref_choice in ref_choices) {
    ref_choice_actual <- ref_name_map[ref_choice] %||% ref_choice
    for(elisa_type in elisa_types) {
      elisa_col <- if(elisa_type == "linear") "ELISA" else "ELISA_nonlin"
      conc_col <- if(elisa_type == "linear") "Concordance_ELISA" else "Concordance_Nonlin"
      plot_df <- all_data[["Default"]] %>%
        filter(Ref_Type == ref_choice_actual) %>%
        filter(!is.na(.data[[elisa_col]]) & !is.na(Titer) & !is.na(month)) %>%
        mutate(month = as.numeric(month), ELISA = .data[[elisa_col]], Titer = Titer, Concordance = .data[[conc_col]]) %>%
        filter(!groupW %in% exclude_groups)
      if (!is.null(bdl)) {
        plot_df <- plot_df %>% mutate(ELISA = ifelse(!is.na(ELISA) & ELISA < bdl, bdl, ELISA), Titer = ifelse(!is.na(Titer) & Titer < bdl, bdl, Titer))
      }
      plot_df <- plot_df %>%
        mutate(ELISA_plot = ifelse(ELISA < 0.03125, 0.03125, ELISA), Titer_plot = ifelse(Titer < 0.03125, 0.03125, Titer), Ref_Type = ref_choice_actual, Titer_Unit = ifelse(ref_choice_actual == R1_name, "IU/mL", ifelse(ref_choice_actual == "mAb", "IU/mL", "EU/mL")), Method = method_label, ELISA_Type = elisa_type)
      plot_data_list[[paste(ref_choice_actual, elisa_type)]] <- plot_df
    }
  }
  plot_data_list
}

.build_longitudinal_colors <- function(plot_data, color_by) {
  if(color_by == "comb") {
    unique_combs <- sort(unique(plot_data$comb))
    unique_combs <- unique_combs[!is.na(unique_combs)]
    if(length(unique_combs) > 0) {
      full_palette <- brewer.pal(12, "Paired")
      no_yellow_indices <- setdiff(1:12, c(7, 11))
      filtered_palette <- full_palette[no_yellow_indices]
      palette_13 <- colorRampPalette(filtered_palette)(max(13, length(unique_combs)))
      comb_colors <- setNames(palette_13[1:length(unique_combs)], unique_combs)
      return(list(color_var = "comb", legend_title = "Owner-dog identity", colors = comb_colors))
    }
  } else if(color_by == "conc_vax") {
    vaccine_colors <- c("Nobivac" = "#1c1c84", "Rabisin" = "red", "unsure" = "grey60")
    return(list(color_var = "conc_vax", legend_title = "Vaccine type", colors = vaccine_colors))
  }
  list(color_var = "comb", legend_title = "Owner-dog identity", colors = NULL)
}

.compute_gam_summary <- function(sub_data, k_val, show_gam) {
  unique_months <- length(unique(na.omit(sub_data$month)))
  enough_data <- sum(!is.na(sub_data$ELISA) & !is.na(sub_data$month)) >= 3 && unique_months >= 2
  edf <- NA; r_sq <- NA; dev_expl <- NA
  if(enough_data && show_gam) {
    gam_model <- tryCatch(gam(ELISA ~ s(month, k = k_val), data = sub_data, method = "REML", select = TRUE), error = function(e) NULL)
    if(!is.null(gam_model)) {
      gam_summary <- summary(gam_model)
      edf <- round(gam_summary$s.table[1, "edf"], 2)
      r_sq <- round(gam_summary$r.sq, 3)
      dev_expl <- round(gam_summary$dev.expl * 100, 1)
    }
  }
  list(edf = edf, r_sq = r_sq, dev_expl = dev_expl)
}

.generate_longitudinal_plots <- function(plot_data_list, ref_name_map, elisa_types, exclude_groups, show_gam, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label, longitudinal_min_samples_per_comb, longitudinal_width, longitudinal_height, longitudinal_font_size, base_size, legend_position, legend_ncol, legend_key_size, legend_text_size, legend_title_size, color_by, dot_size, dot_alpha, dot_color, point_size, point_alpha, gm_one_panel_width, gm_one_panel_height, low_titer_cutoff, bdl) {
  all_plots <- list()
  if(y_scale == "log2") {
    y_scale_func <- scale_y_log10
    y_breaks <- log2_breaks
    y_labels <- log2_labels
    y_limits_use <- y_limits
  } else {
    y_scale_func <- scale_y_continuous
    y_breaks <- waiver()
    y_labels <- waiver()
    y_limits_use <- c(0, max(plot_data_list[[1]]$ELISA, plot_data_list[[1]]$Titer, na.rm = TRUE) * 1.1)
  }
  for(ref_choice in names(ref_name_map)) {
    ref_choice_actual <- ref_name_map[ref_choice]
    ref_choice_file <- gsub("[^a-zA-Z0-9_.-]", "_", ref_choice_actual)
    for(elisa_type in elisa_types) {
      plot_df <- plot_data_list[[paste(ref_choice_actual, elisa_type)]]
      if(nrow(plot_df) == 0) next
      if(is.null(legend_key_size)) legend_key_size <- unit(0.3, "lines")
      if(is.null(legend_text_size)) legend_text_size <- base_size - 2
      if(is.null(legend_title_size)) legend_title_size <- base_size - 1
      plot_data <- plot_df %>% arrange(comb, month) %>% mutate(comb = factor(comb))
      eligible_combs <- plot_data %>% group_by(comb) %>% summarise(n_samples = n_distinct(Sample), .groups = "drop") %>% filter(n_samples >= longitudinal_min_samples_per_comb) %>% pull(comb)
      plot_data <- plot_data %>% filter(comb %in% eligible_combs)
      color_info <- .build_longitudinal_colors(plot_data, color_by)
      color_var <- color_info$color_var
      legend_title <- color_info$legend_title
      comb_colors <- color_info$colors
      vaccine_colors <- c("Nobivac" = "#1c1c84", "Rabisin" = "red", "unsure" = "grey60")
      groups <- unique(plot_data$groupW)
      groups <- groups[!is.na(groups) & !groups %in% exclude_groups]
      plot_list <- list()
      plot_counter_long <- 1
      for(g in groups) {
        sub_data <- plot_data %>% filter(groupW == g)
        if(nrow(sub_data) < 3) next
        titer_unit_local <- unique(sub_data$Titer_Unit)[1]
        if(is.na(titer_unit_local)) titer_unit_local <- ""
        unique_months <- length(unique(na.omit(sub_data$month)))
        k_val <- max(3, min(10, unique_months - 1))
        enough_data <- sum(!is.na(sub_data$ELISA) & !is.na(sub_data$month)) >= 3 && unique_months >= 2
        gam_stats <- .compute_gam_summary(sub_data, k_val, show_gam)
        edf <- gam_stats$edf; r_sq <- gam_stats$r_sq; dev_expl <- gam_stats$dev_expl
        plot_sub <- sub_data %>%
          pivot_longer(cols = c(ELISA, Titer), names_to = "Assay", values_to = "Value") %>%
          mutate(Assay = ifelse(Assay == "ELISA", "ELISA", "Neutralization"), Value_plot = ifelse(is.na(Value), NA, ifelse(Value < low_titer_cutoff, low_titer_cutoff, Value)), Value_plot = ifelse(Value_plot < 0.03125, 0.03125, Value_plot))
        p <- ggplot(plot_sub, aes(x = month, y = Value_plot, group = interaction(comb, Assay), color = !!sym(color_var))) +
          geom_point(size = 1.8, alpha = 0.8) + geom_path(alpha = 0.8) + geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.65, color = "black") + y_scale_func(limits = y_limits_use, breaks = y_breaks, labels = y_labels) + scale_x_continuous(breaks = 1:18, limits = c(1, 18.5), expand = expansion(mult = c(0.05, 0))) + facet_wrap(~Assay, ncol = 2, scales = "free_y") + labs(x = "Months since last vaccination", y = "Titer (EU/mL)", color = legend_title) + plot_themes$base_theme + theme(plot.title = element_text(size = longitudinal_font_size + 1, face = "bold", hjust = 0.5, margin = margin(b = 1)), plot.subtitle = element_text(size = longitudinal_font_size - 1, hjust = 0.5, color = "grey40", margin = margin(b = 1)), plot.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), axis.text = element_text(size = longitudinal_font_size, face = "bold"), axis.title = element_text(size = longitudinal_font_size + 1, face = "bold"), strip.text = element_text(size = longitudinal_font_size, face = "bold"), strip.background = element_rect(linewidth = 0), panel.spacing = unit(0.5, "lines"))
        if(legend_position == "inside") {
          p <- p + theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top"), legend.background = element_rect(fill = alpha("white", 0.85), color = NA, linewidth = 0), legend.box.background = element_rect(fill = NA, color = NA), legend.margin = margin(1, 1, 1, 1), legend.spacing.y = unit(0, "lines"), legend.title = element_text(size = legend_title_size, face = "bold"), legend.text = element_text(size = legend_text_size), legend.key.size = legend_key_size, legend.direction = "vertical")
        } else {
          p <- p + theme(legend.position = "bottom", legend.background = element_rect(fill = NA, color = NA), legend.box.background = element_rect(fill = NA, color = NA), legend.margin = margin(1, 1, 1, 1), legend.spacing.y = unit(0, "lines"), legend.title = element_text(size = legend_title_size, face = "bold"), legend.text = element_text(size = legend_text_size), legend.key.size = legend_key_size, legend.direction = "vertical")
        }
        p <- p + guides(color = guide_legend(ncol = legend_ncol, title.position = "top", title.hjust = 0.5, byrow = FALSE))
        if(color_by == "comb" && exists("comb_colors")) p <- p + scale_color_manual(values = comb_colors) else if(color_by == "conc_vax") p <- p + scale_color_manual(values = vaccine_colors)
        if(show_gam && enough_data) {
          p <- p + geom_smooth(data = plot_sub %>% filter(Assay == "ELISA"), aes(x = month, y = Value, group = 1), method = "gam", formula = y ~ s(x, k = k_val), method.args = list(select = TRUE), color = "#242424", size = 0.6, se = TRUE, inherit.aes = FALSE) + geom_smooth(data = plot_sub %>% filter(Assay == "Neutralization"), aes(x = month, y = Value, group = 1), method = "gam", formula = y ~ s(x, k = k_val), method.args = list(select = TRUE), color = "#242424", size = 0.6, se = TRUE, inherit.aes = FALSE)
        }
        plot_list[[plot_counter_long]] <- p
        plot_counter_long <- plot_counter_long + 1
      }
      if(length(plot_list) > 0) {
        n_plots <- length(plot_list)
        ncol_long <- if(!is.null(longitudinal_ncol)) longitudinal_ncol else if(n_plots <= 2) 1 else 2
        nrow_long <- if(!is.null(longitudinal_nrow)) longitudinal_nrow else ceiling(n_plots / ncol_long)
        final_plot <- wrap_plots(plot_list, ncol = ncol_long, nrow = nrow_long)
        file_suffix <- paste0("_", method_label, "_", ref_choice_file, "_", elisa_type)
        ggsave(file.path(plot_dir, paste0("Fig2_Longitudinal", file_suffix, ".png")), final_plot, width = if(!is.null(longitudinal_width)) longitudinal_width else gm_one_panel_width * 2, height = if(!is.null(longitudinal_height)) longitudinal_height else gm_one_panel_height, dpi = 300, limitsize = FALSE)
        all_plots[[paste0("Longitudinal_", method_label, "_", ref_choice_actual, "_", elisa_type)]] <- final_plot
      }
    }
  }
  all_plots
}

.generate_boxplots <- function(plot_data_list, ref_name_map, elisa_types, exclude_groups, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label, boxplot_x_text_size, boxplot_x_text_angle, low_titer_cutoff, bdl) {
  all_plots <- list()
  if(y_scale == "log2") {
    y_scale_func <- scale_y_log10
    y_limits_use <- y_limits
  } else {
    y_scale_func <- scale_y_continuous
    y_limits_use <- c(0, max(plot_data_list[[1]]$ELISA, plot_data_list[[1]]$Titer, na.rm = TRUE) * 1.1)
  }
  for(ref_choice in names(ref_name_map)) {
    ref_choice_actual <- ref_name_map[ref_choice]
    ref_choice_file <- gsub("[^a-zA-Z0-9_.-]", "_", ref_choice_actual)
    for(elisa_type in elisa_types) {
      plot_df <- plot_data_list[[paste(ref_choice_actual, elisa_type)]]
      if(nrow(plot_df) == 0 || length(unique(plot_df$groupW[!is.na(plot_df$groupW)])) == 0) next
      p_box <- plot_df %>%
        pivot_longer(cols = c(ELISA, Titer), names_to = "Assay", values_to = "Value") %>%
        mutate(Assay = ifelse(Assay == "ELISA", "ELISA", "Neutralization"), Value_plot = ifelse(is.na(Value), NA, ifelse(Value < low_titer_cutoff, low_titer_cutoff, Value)), Value_plot = ifelse(Value_plot < 0.03125, 0.03125, Value_plot)) %>%
        ggplot(aes(x = groupW, y = Value_plot, fill = Assay)) + geom_boxplot(alpha = 0.7, outlier.size = 1.5) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) + y_scale_func(limits = y_limits_use) + scale_fill_manual(values = c("ELISA" = "#2E86AB", "Neutralization" = "#F18F01")) + labs(x = "Group", y = "Titer (EU/mL)", fill = "Assay Type") + plot_themes$theme_publication() + theme(axis.text.x = element_text(angle = boxplot_x_text_angle, hjust = ifelse(boxplot_x_text_angle == 0, 0.5, 1), size = boxplot_x_text_size), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
      file_suffix_box <- paste0("_", method_label, "_", ref_choice_file, "_", elisa_type)
      ggsave(file.path(plot_dir, paste0("Fig4_Boxplots_by_Group", file_suffix_box, ".png")), p_box, width = 10, height = 7, dpi = 300)
      all_plots[[paste0("Boxplots_", method_label, "_", ref_choice_actual, "_", elisa_type)]] <- p_box
    }
  }
  all_plots
}

.generate_scatter_plots <- function(plot_data_list, ref_name_map, elisa_types, exclude_groups, y_scale, y_limits, log2_breaks, log2_labels, plot_themes, plot_dir, method_label, gm_width, gm_height, gm_jitter_width, gm_two_panel_width, gm_two_panel_height, scatter_font_size, low_titer_cutoff, bdl, geo_mean, geo_sd_factor, month_group_levels) {
  all_plots <- list()
  if(y_scale == "log2") {
    y_scale_func <- scale_y_log10
    y_breaks <- log2_breaks
    y_labels <- log2_labels
    y_limits_use <- y_limits
  } else {
    y_scale_func <- scale_y_continuous
    y_breaks <- waiver()
    y_labels <- waiver()
    y_limits_use <- c(0, max(plot_data_list[[1]]$ELISA, plot_data_list[[1]]$Titer, na.rm = TRUE) * 1.1)
  }
  for(ref_choice in names(ref_name_map)) {
    ref_choice_actual <- ref_name_map[ref_choice]
    ref_choice_file <- gsub("[^a-zA-Z0-9_.-]", "_", ref_choice_actual)
    for(elisa_type in elisa_types) {
      elisa_col <- if(elisa_type == "linear") "ELISA" else "ELISA_nonlin"
      scatter_plot_df <- plot_data_list[[paste(ref_choice_actual, elisa_type)]] %>%
        mutate(month = as.numeric(month), ELISA = .data[[elisa_col]], Titer = Titer, month_group = bucket_months(month, month_group_levels)) %>%
        filter(!is.na(month_group)) %>% filter(startsWith(Sample, "ibet"))
      if (!is.null(bdl)) {
        scatter_plot_df <- scatter_plot_df %>% mutate(ELISA = ifelse(!is.na(ELISA) & ELISA < bdl, bdl, ELISA), Titer = ifelse(!is.na(Titer) & Titer < bdl, bdl, Titer))
      }
      scatter_plot_df <- scatter_plot_df %>%
        pivot_longer(cols = c(ELISA, Titer), names_to = "Assay", values_to = "Value") %>%
        mutate(Assay = ifelse(Assay == "ELISA", "ELISA", "Neutralization"), Value_plot = ifelse(is.na(Value), NA, ifelse(Value < low_titer_cutoff, low_titer_cutoff, Value)), Value_plot = ifelse(Value_plot < 0.03125, 0.03125, Value_plot)) %>%
        filter(!is.na(Value_plot))
      if(nrow(scatter_plot_df) == 0) next
      n_dogs_scatter <- n_distinct(scatter_plot_df$comb)
      titer_unit_local <- unique(scatter_plot_df$Titer_Unit[scatter_plot_df$Titer_Unit != ""])
      if(length(titer_unit_local) == 0 || all(is.na(titer_unit_local))) titer_unit_local <- ""
      x_levels <- month_group_levels
      vline_positions <- seq(0.5, length(x_levels) - 0.5, by = 1)
      seg_offsets <- 0.4
      gm_scatter_df <- scatter_plot_df %>% group_by(month_group, Assay) %>% summarise(n = n(), geo_mean_v = geo_mean(Value_plot), .groups = "drop") %>% mutate(x_num = as.numeric(month_group))
      scatter_base <- list(geom_point(color = "#323232", size = 2.2, alpha = 0.6, position = position_jitter(width = 0.3, height = 0)), geom_vline(xintercept = vline_positions, linetype = "dashed", linewidth = 0.5, color = "#191919"), geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.8, color = "black"), y_scale_func(limits = y_limits_use, breaks = y_breaks, labels = y_labels), scale_x_discrete(drop = FALSE), labs(x = "months since last vaccination", title = NULL, subtitle = NULL), plot_themes$base_theme, theme(plot.title = element_text(size = scatter_font_size + 1, face = "bold", hjust = 0.5, margin = margin(b = 1)), plot.subtitle = element_text(size = scatter_font_size - 1, hjust = 0.5, color = "grey40", margin = margin(b = 1)), axis.text = element_text(size = scatter_font_size, face = "bold"), axis.title = element_text(size = scatter_font_size + 1, face = "bold"), strip.text = element_text(size = scatter_font_size, face = "bold"), strip.background = element_rect(linewidth = 0), panel.spacing = unit(0.5, "lines")))
      gm_segment_elisa <- geom_segment(data = gm_scatter_df %>% filter(Assay == "ELISA"), aes(x = x_num - seg_offsets, xend = x_num + seg_offsets, y = geo_mean_v, yend = geo_mean_v), color = "black", linewidth = 1.2, inherit.aes = FALSE)
      gm_segment_neut <- geom_segment(data = gm_scatter_df %>% filter(Assay == "Neutralization"), aes(x = x_num - seg_offsets, xend = x_num + seg_offsets, y = geo_mean_v, yend = geo_mean_v), color = "black", linewidth = 1.2, inherit.aes = FALSE)
      p_scatter_elisa <- ggplot(scatter_plot_df %>% filter(Assay == "ELISA"), aes(x = month_group, y = Value_plot)) + scatter_base + gm_segment_elisa + labs(y = "anti-RABV-G IgG (EU/mL)")
      p_scatter_neut <- ggplot(scatter_plot_df %>% filter(Assay == "Neutralization"), aes(x = month_group, y = Value_plot)) + scatter_base + gm_segment_neut + labs(y = "PBNA titer (EU/mL)")
      p_scatter <- (p_scatter_elisa | p_scatter_neut) + plot_layout(guides = "collect")
      file_suffix_scatter <- paste0("_", method_label, "_", ref_choice_file, "_", elisa_type)
      ggsave(file.path(plot_dir, paste0("Fig3_Bimonthly_Scatter", file_suffix_scatter, ".png")), p_scatter, width = gm_two_panel_width, height = gm_two_panel_height, dpi = 300, limitsize = FALSE)
      all_plots[[paste0("Bimonthly_Scatter_", method_label, "_", ref_choice_actual, "_", elisa_type)]] <- p_scatter
    }
  }
  all_plots
}

.generate_concordance_plot <- function(overall_df, plot_themes, plot_dir, plot_width_one_panel, plot_height_one_panel) {
  all_plots <- list()
  if("Concordance" %in% colnames(overall_df)) {
    concordance_data <- overall_df %>% group_by(Ref_Type, Concordance) %>% summarise(Count = n(), .groups = "drop") %>% filter(Concordance != "Unknown" & !is.na(Concordance))
    if(nrow(concordance_data) > 0) {
      concordance_data <- concordance_data %>% group_by(Ref_Type) %>% mutate(Percentage = Count / sum(Count) * 100) %>% ungroup()
      p_conc <- ggplot(concordance_data, aes(x = Concordance, y = Count, fill = Concordance)) + geom_bar(stat = "identity", alpha = 0.85) + geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")), vjust = -0.3, size = 3.5) + ylim(0, 110) + scale_fill_manual(values = c("Concordant" = "#228B22", "Discordant" = "#DC143C")) + facet_wrap(~Ref_Type) + labs(subtitle = NULL, x = "Concordance Status", y = "Number of Samples") + plot_themes$theme_publication() + theme(legend.position = "none", strip.text = element_text(size = 10, face = "bold"))
      ggsave(file.path(plot_dir, "Fig5_Concordance.png"), p_conc, width = plot_width_one_panel, height = plot_height_one_panel, dpi = 300)
      all_plots[["Concordance"]] <- p_conc
    }
  }
  all_plots
}

.assign_sample_categories <- function(data_prepped, bdl_val = 0.125) {
  data_prepped$elisa_category <- ifelse(is.na(data_prepped$ELISA), NA_character_, ifelse(data_prepped$ELISA >= 0.5, "Had >= 0.5 U/mL", ifelse(data_prepped$ELISA > bdl_val, "Had < 0.5 U/mL but > BDL", "IgG level is BDL")))
  data_prepped$titer_category <- ifelse(is.na(data_prepped$Titer), NA_character_, ifelse(data_prepped$Titer >= 0.5, "Had >= 0.5 U/mL", ifelse(data_prepped$Titer > bdl_val, "Had < 0.5 U/mL but > BDL", "IgG level is BDL")))
  data_prepped[!is.na(data_prepped$elisa_category) & !is.na(data_prepped$titer_category), , drop = FALSE]
}

.generate_stacked_bar_plots <- function(overall_df, ref_name_map, elisa_types, exclude_groups, output_dir, plot_dir, method_label, stacked_bar_nrow, stacked_bar_ncol, stacked_bar_two_panel_width, stacked_bar_two_panel_height, stacked_bar_sd_label_size, stacked_bar_gm_label_size, stacked_bar_base_size, stacked_bar_title_size, stacked_bar_axis_title_size, stacked_bar_axis_text_size, bdl, month_group_levels) {
  all_plots <- list()
  all_stacked_bar_data <- list()
  .build_stacked_bar_panel <- function(data, sample_sz = 1, ref_name = "", method_label = "", bdl_val = 0.125) {
    if(nrow(data) == 0) return(NULL)
    data_prepped <- data
    data_prepped$month <- as.numeric(data_prepped$month)
    data_prepped$ELISA <- as.numeric(data_prepped$ELISA)
    data_prepped$Titer <- as.numeric(data_prepped$Titer)
    data_prepped$month_group <- bucket_months(data_prepped$month, month_group_levels)
    data_prepped <- data_prepped[!is.na(data_prepped$month_group), , drop = FALSE]
    if(nrow(data_prepped) == 0) return(NULL)
    data_prepped <- .assign_sample_categories(data_prepped, bdl_val)
    if(nrow(data_prepped) == 0) return(NULL)
    x_levels <- month_group_levels
    data_prepped$month_group <- factor(data_prepped$month_group, levels = x_levels)
    data_prepped$elisa_category <- factor(data_prepped$elisa_category, levels = c("IgG level is BDL", "Had < 0.5 U/mL but > BDL", "Had >= 0.5 U/mL"))
    data_prepped$titer_category <- factor(data_prepped$titer_category, levels = c("IgG level is BDL", "Had < 0.5 U/mL but > BDL", "Had >= 0.5 U/mL"))
    bar_list <- list()
    for(mg in levels(data_prepped$month_group)) {
      sub_mg <- data_prepped[data_prepped$month_group == mg, , drop = FALSE]
      for(cat in levels(data_prepped$elisa_category)) {
        n_in <- sum(sub_mg$elisa_category == cat, na.rm = TRUE)
        bar_list[[length(bar_list) + 1]] <- data.frame(month_group = factor(mg, levels = x_levels), adj_category = factor(cat, levels = levels(data_prepped$elisa_category)), n_in = n_in, Assay = "ELISA", stringsAsFactors = FALSE)
      }
      for(cat in levels(data_prepped$titer_category)) {
        n_in <- sum(sub_mg$titer_category == cat, na.rm = TRUE)
        bar_list[[length(bar_list) + 1]] <- data.frame(month_group = factor(mg, levels = x_levels), adj_category = factor(cat, levels = levels(data_prepped$titer_category)), n_in = n_in, Assay = "Neutralization", stringsAsFactors = FALSE)
      }
    }
    bar_data <- do.call(rbind, bar_list)
    rownames(bar_data) <- NULL
    bar_data <- bar_data[bar_data$n_in > 0, , drop = FALSE]
    if(nrow(bar_data) == 0) return(NULL)
    bar_data$ref_name <- ref_name
    bar_data$method_label <- method_label
    bar_data$sample_sz <- sample_sz
    bar_data$elisa_type <- elisa_type
    all_stacked_bar_data[[length(all_stacked_bar_data) + 1]] <- bar_data
    total_per_group <- tapply(bar_data$n_in, list(bar_data$month_group, bar_data$Assay), sum)
    bar_data$percent <- NA_real_
    for(i in seq_len(nrow(bar_data))) {
      mg_i <- as.character(bar_data$month_group[i])
      assay_i <- as.character(bar_data$Assay[i])
      total_i <- total_per_group[mg_i, assay_i]
      if(!is.na(total_i) && total_i > 0) bar_data$percent[i] <- 100 * bar_data$n_in[i] / total_i
    }
    bar_data <- as.data.frame(bar_data)
    geo_data <- as.data.frame(do.call(rbind, lapply(levels(data_prepped$month_group), function(mg) {
      sub <- data_prepped[data_prepped$month_group == mg & data_prepped$elisa_category != "IgG level is BDL", , drop = FALSE]
      do.call(rbind, lapply(c("Had >= 0.5 U/mL", "Had < 0.5 U/mL but > BDL"), function(cat) {
        v <- sub$ELISA[sub$elisa_category == cat]
        gm <- geo_mean(v)
        if(is.na(gm)) return(NULL)
        data.frame(month_group = factor(mg, levels = x_levels), adj_category = factor(cat, levels = c("Had < 0.5 U/mL but > BDL", "Had >= 0.5 U/mL")), gm = gm, Assay = "ELISA", stringsAsFactors = FALSE)
      }))
    })), stringsAsFactors = FALSE)
    if(nrow(geo_data) > 0) geo_data$month_group <- factor(geo_data$month_group, levels = x_levels)
    geo_data_titer <- as.data.frame(do.call(rbind, lapply(levels(data_prepped$month_group), function(mg) {
      sub <- data_prepped[data_prepped$month_group == mg & data_prepped$titer_category != "IgG level is BDL", , drop = FALSE]
      do.call(rbind, lapply(c("Had >= 0.5 U/mL", "Had < 0.5 U/mL but > BDL"), function(cat) {
        v <- sub$Titer[sub$titer_category == cat]
        gm <- geo_mean(v)
        if(is.na(gm)) return(NULL)
        data.frame(month_group = factor(mg, levels = x_levels), adj_category = factor(cat, levels = c("Had < 0.5 U/mL but > BDL", "Had >= 0.5 U/mL")), gm = gm, Assay = "Neutralization", stringsAsFactors = FALSE)
      }))
    })), stringsAsFactors = FALSE)
    if(nrow(geo_data_titer) > 0) geo_data_titer$month_group <- factor(geo_data_titer$month_group, levels = x_levels)
    geo_data <- rbind(geo_data, geo_data_titer)
    data_mg <- data_prepped
    data_mg$month_group <- factor(data_mg$month_group, levels = x_levels)
    dog_counts_list <- list()
    for(mg in levels(data_mg$month_group)) {
      parts <- strsplit(mg, "-")[[1]]
      months_in_group <- as.integer(parts[1]):as.integer(parts[length(parts)])
      sub <- data[data$month %in% months_in_group, , drop = FALSE]
      dog_counts_list[[length(dog_counts_list) + 1]] <- data.frame(month_group = factor(mg, levels = x_levels), D = length(unique(sub$comb)), S = sum(!is.na(sub$ELISA) | !is.na(sub$Titer)), stringsAsFactors = FALSE)
    }
    dog_counts <- do.call(rbind, dog_counts_list)
    rownames(dog_counts) <- NULL
    dog_counts$month_group <- factor(dog_counts$month_group, levels = x_levels)
    y_pos_df <- aggregate(percent ~ month_group, data = bar_data, FUN = sum)
    names(y_pos_df)[2] <- "y_pos"
    sd_labels <- dog_counts
    sd_labels$label <- paste0("S: ", sd_labels$S, "/D: ", sd_labels$D)
    sd_labels <- merge(sd_labels, y_pos_df, by = "month_group", all.x = TRUE)
    sd_labels$y_pos[is.na(sd_labels$y_pos)] <- 0
    sd_labels$y_pos <- sd_labels$y_pos + 5
    fill_colors <- c("IgG level is BDL" = "#4c4c4c", "Had < 0.5 U/mL but > BDL" = "#FFCC33", "Had >= 0.5 U/mL" = "#4FA861")
    line_colors <- c("Had >= 0.5 U/mL" = "#e41a1c", "Had < 0.5 U/mL but > BDL" = "#377eb8")
    plot_title <- paste0(">= ", sample_sz, " blood sample(s) (", ref_name, ")")
    gm_ratio <- 25
    p <- ggplot() + geom_bar(data = bar_data, aes(x = month_group, y = percent, fill = adj_category), stat = "identity", position = "stack", width = 0.8) + geom_text(data = sd_labels, aes(x = month_group, y = y_pos, label = label), size = stacked_bar_sd_label_size, fontface = "bold", color = "black") + scale_fill_manual(values = fill_colors, name = "IgG Ab Response Range\n(S: sera datapoints, D: dogs)") + scale_y_continuous(name = "Sample proportion (%)", limits = c(0, 100), sec.axis = sec_axis(~ . * 0.04, name = "Geometric Mean Titer (U/mL)", breaks = seq(0, 4, 0.5))) + facet_wrap(~Assay, scales = "free_y")
    if(nrow(geo_data) > 0) {
      p <- p + ggnewscale::new_scale_color() +
        {if(any(geo_data$adj_category == "Had >= 0.5 U/mL")) {
          geo_hi <- geo_data[geo_data$adj_category == "Had >= 0.5 U/mL", , drop = FALSE]
          list(geom_line(data = geo_hi, aes(x = month_group, y = gm * gm_ratio, group = 1, color = "Had >= 0.5 U/mL"), linewidth = 1.2), geom_point(data = geo_hi, aes(x = month_group, y = gm * gm_ratio, color = "Had >= 0.5 U/mL"), size = 2.5), geom_text(data = geo_hi, aes(x = month_group, y = gm * gm_ratio - 3, label = sprintf("%.2f", gm)), color = "#6c1010", size = stacked_bar_gm_label_size, fontface = "bold"))
        }} +
        {if(any(geo_data$adj_category == "Had < 0.5 U/mL but > BDL")) {
          geo_lo <- geo_data[geo_data$adj_category == "Had < 0.5 U/mL but > BDL", , drop = FALSE]
          list(geom_line(data = geo_lo, aes(x = month_group, y = gm * gm_ratio, group = 1, color = "Had < 0.5 U/mL but > BDL"), linetype = "dashed", linewidth = 1.2), geom_point(data = geo_lo, aes(x = month_group, y = gm * gm_ratio, color = "Had < 0.5 U/mL but > BDL"), size = 2.5), geom_text(data = geo_lo, aes(x = month_group, y = gm * gm_ratio - 3, label = sprintf("%.2f", gm)), color = "#003747", size = stacked_bar_gm_label_size, fontface = "bold"))
        }} +
        scale_color_manual(values = line_colors, name = "Geometric Mean Lines")
    }
    p + labs(x = "Months since last vaccination", title = plot_title) + theme_classic(base_size = stacked_bar_base_size) + theme(axis.text = element_text(size = stacked_bar_axis_text_size, face = "bold"), axis.title = element_text(size = stacked_bar_axis_title_size, face = "bold"), legend.title = element_text(face = "bold", size = 9), legend.position = "none", plot.title = element_text(face = "bold", size = stacked_bar_title_size, hjust = 0.5))
  }
  for(elisa_type in elisa_types) {
    elisa_col <- if(elisa_type == "linear") "ELISA" else "ELISA_nonlin"
    panels_ge1 <- list()
    panels_ge2 <- list()
    for(ref_choice in names(ref_name_map)) {
      ref_choice_actual <- ref_name_map[ref_choice]
      ref_choice_file <- gsub("[^a-zA-Z0-9_.-]", "_", ref_choice_actual)
      ref_data <- overall_df %>% filter(Ref_Type == ref_choice_actual) %>% filter(!is.na(.data[[elisa_col]]) & !is.na(Titer) & !is.na(month)) %>% mutate(month = as.numeric(month), ELISA = .data[[elisa_col]], Titer = Titer) %>% filter(!groupW %in% exclude_groups)
      if(!is.null(bdl)) {
        ref_data <- ref_data %>% mutate(ELISA = ifelse(!is.na(ELISA) & ELISA < bdl, bdl, ELISA), Titer = ifelse(!is.na(Titer) & Titer < bdl, bdl, Titer))
      }
      eligible_ge1 <- ref_data
      eligible_ge2 <- ref_data %>% group_by(comb) %>% filter(n_distinct(Sample) >= 2) %>% ungroup()
      p1 <- .build_stacked_bar_panel(eligible_ge1, sample_sz = 1, ref_name = ref_choice_actual, method_label = method_label, bdl_val = bdl)
      p2 <- .build_stacked_bar_panel(eligible_ge2, sample_sz = 2, ref_name = ref_choice_actual, method_label = method_label, bdl_val = bdl)
      panels_ge1[[ref_choice_actual]] <- p1
      panels_ge2[[ref_choice_actual]] <- p2
    }
    valid_ge1 <- panels_ge1[!sapply(panels_ge1, is.null)]
    if(length(valid_ge1) > 0) {
      combined_ge1 <- wrap_plots(valid_ge1) + plot_layout(ncol = stacked_bar_ncol, nrow = stacked_bar_nrow, guides = "collect")
      file_suffix_sb <- paste0("_", method_label, "_", elisa_type, "_ge1")
      ggsave(file.path(plot_dir, paste0("Fig3_Stacked_Bars", file_suffix_sb, ".png")), combined_ge1, width = stacked_bar_two_panel_width, height = stacked_bar_two_panel_height, dpi = 300)
      all_plots[[paste0("Stacked_Bars_ge1_", method_label, "_", elisa_type)]] <- combined_ge1
    }
    valid_ge2 <- panels_ge2[!sapply(panels_ge2, is.null)]
    if(length(valid_ge2) > 0) {
      combined_ge2 <- wrap_plots(valid_ge2) + plot_layout(ncol = stacked_bar_ncol, nrow = stacked_bar_nrow, guides = "collect")
      file_suffix_sb <- paste0("_", method_label, "_", elisa_type, "_ge2")
      ggsave(file.path(plot_dir, paste0("Fig3_Stacked_Bars", file_suffix_sb, ".png")), combined_ge2, width = stacked_bar_two_panel_width, height = stacked_bar_two_panel_height, dpi = 300)
      all_plots[[paste0("Stacked_Bars_ge2_", method_label, "_", elisa_type)]] <- combined_ge2
    }
  }
  all_bar_data_df <- NULL
  if(length(all_stacked_bar_data) > 0) {
    all_bar_data_df <- do.call(rbind, all_stacked_bar_data)
    wb_bar <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_bar, "Stacked_Bar_Data")
    openxlsx::writeData(wb_bar, "Stacked_Bar_Data", all_bar_data_df)
    openxlsx::saveWorkbook(wb_bar, file.path(output_dir, "stacked_bar_data.xlsx"), overwrite = TRUE)
  }
  list(plots = all_plots, stacked_bar_data = all_bar_data_df)
}


