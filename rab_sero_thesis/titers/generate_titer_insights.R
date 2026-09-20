library(tidyverse)
library(ggplot2)
library(viridis)
source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

generate_titer_insights <- function(
    stats_results = NULL,
    file_path = "RABV_VNT_DATA.xlsx",
    sheet = 3,
    output_dir = "titer_insights",
    reference_standard = "eRIG_1:600",
    reference_concentration = 0.5,
    ed50_tolerance = 0.5,
    sample_prefix_filter = NULL,
    generate_plots = TRUE
) {
  use_dual_reference <- FALSE
  second_reference_standard <- NULL
  second_reference_concentration <- 0.5
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
  plot_dir <- file.path(output_dir, "plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive=TRUE)
  
  use_dual_reference <- !is.null(second_reference_standard)
  standard_unit <- if(reference_standard=="eRIG_1:600") "IU/mL" else "EU/mL"
  second_unit <- if(use_dual_reference && second_reference_standard=="eRIG_1:600") "IU/mL" else "EU/mL"
  
  if(is.null(stats_results)) {
    csv_dirs <- unique(c(output_dir, "neutralization_statistics"))
    cs_main <- cs_f4 <- cs_bf <- NULL
    for(d in csv_dirs) {
      candidate_main <- file.path(d, "sample_results.csv")
      if(file.exists(candidate_main)) {
        cs_main <- candidate_main
        cs_f4   <- file.path(d, "sample_results_forced4pl.csv")
        cs_bf   <- file.path(d, "sample_results_5pl.csv")
        break
      }
    }
    
    if(!is.null(cs_main) && file.exists(cs_main)) {
      combined_summary_csv <- read.csv(cs_main, stringsAsFactors = FALSE)
      cs_f4_csv  <- if(!is.null(cs_f4) && file.exists(cs_f4)) read.csv(cs_f4, stringsAsFactors = FALSE) else NULL
      cs_bf_csv  <- if(!is.null(cs_bf) && file.exists(cs_bf)) read.csv(cs_bf, stringsAsFactors = FALSE) else NULL
        stats_results <- list(
          combined_summary = combined_summary_csv,
          combined_summary_forced4pl = cs_f4_csv,
          combined_summary_sk = cs_bf_csv
        )
    } else {
      stats_results <- tryCatch({
        analyze_neutralization_statistics(
          file_path = file_path, sheet = sheet, output_dir = output_dir,
          reference_standard = reference_standard, reference_concentration = reference_concentration,
          second_reference_standard = second_reference_standard, second_reference_concentration = second_reference_concentration,
          sample_prefix_filter = NULL,
          generate_publication_graphs = FALSE
        )
      }, error = function(e) {
        NULL
      })
      
      if(is.null(stats_results) || is.null(stats_results$combined_summary)) {
        if(!is.null(cs_main) && file.exists(cs_main)) {
          combined_summary_csv <- read.csv(cs_main, stringsAsFactors = FALSE)
          cs_f4_csv  <- if(!is.null(cs_f4) && file.exists(cs_f4)) read.csv(cs_f4, stringsAsFactors = FALSE) else NULL
          cs_bf_csv  <- if(!is.null(cs_bf) && file.exists(cs_bf)) read.csv(cs_bf, stringsAsFactors = FALSE) else NULL
            stats_results <- list(
              combined_summary = combined_summary_csv,
              combined_summary_forced4pl = cs_f4_csv,
              combined_summary_sk = cs_bf_csv
            )
        }
      }
    }
  }
  
  combined_summary <- stats_results$combined_summary
  if(!is.null(stats_results$combined_summary_forced4pl) &&
     "Titer_Ref1" %in% colnames(stats_results$combined_summary_forced4pl)) {
    cs_f4 <- stats_results$combined_summary_forced4pl
    cs_f4 <- cs_f4[!is.na(cs_f4$Sample), c("Sample", "Titer_Ref1")]
    colnames(cs_f4)[colnames(cs_f4) == "Titer_Ref1"] <- "Titer_Forced4PL"
    combined_summary <- merge(combined_summary, cs_f4, by = "Sample", all.x = TRUE)
  }
  if(!is.null(stats_results$combined_summary_sk) &&
     "Titer_Ref1" %in% colnames(stats_results$combined_summary_sk)) {
    cs_bf <- stats_results$combined_summary_sk
    cs_bf <- cs_bf[!is.na(cs_bf$Sample), c("Sample", "Titer_Ref1")]
    colnames(cs_bf)[colnames(cs_bf) == "Titer_Ref1"] <- "Titer_5PL"
    combined_summary <- merge(combined_summary, cs_bf, by = "Sample", all.x = TRUE)
  }
  if(is.null(combined_summary) || nrow(combined_summary) == 0) {
    return(invisible(NULL))
  }
  
  combined_summary <- ref_filter(combined_summary, sample_prefix_filter = sample_prefix_filter)
  
  if(nrow(combined_summary) == 0 && !is.null(stats_results$combined_summary) &&
     nrow(stats_results$combined_summary) > 0 && !is.null(sample_prefix_filter)) {
    combined_summary <- stats_results$combined_summary
    if(!is.null(stats_results$combined_summary_forced4pl) &&
       "Titer_Ref1" %in% colnames(stats_results$combined_summary_forced4pl)) {
      cs_f4 <- stats_results$combined_summary_forced4pl
      cs_f4 <- cs_f4[!is.na(cs_f4$Sample), c("Sample", "Titer_Ref1")]
      colnames(cs_f4)[colnames(cs_f4) == "Titer_Ref1"] <- "Titer_Forced4PL"
      combined_summary <- merge(combined_summary, cs_f4, by = "Sample", all.x = TRUE)
    }
    if(!is.null(stats_results$combined_summary_sk) &&
       "Titer_Ref1" %in% colnames(stats_results$combined_summary_sk)) {
      cs_bf <- stats_results$combined_summary_sk
      cs_bf <- cs_bf[!is.na(cs_bf$Sample), c("Sample", "Titer_Ref1")]
      colnames(cs_bf)[colnames(cs_bf) == "Titer_Ref1"] <- "Titer_5PL"
      combined_summary <- merge(combined_summary, cs_bf, by = "Sample", all.x = TRUE)
    }
  }
  
  if(nrow(combined_summary) == 0) {
    return(invisible(NULL))
  }
  
   method_colors <- c("AICMix"="#1B9E77", "Forced 4PL"="#D55E00", "Forced 5PL"="#6A3D9A", "Bayes-4PL"="#E41A1C", "Reed-Muench"="#88CCEE")
  conc_colors <- c("Concordant"="#228B22", "Discordant"="#DC143C", "Unknown"="#808080")
  
            "(prefix filter: ", ifelse(is.null(sample_prefix_filter), "none", sample_prefix_filter), ")")
  
  insight_theme <- function(base_size = 12) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "sans", color = "black"),
        axis.title = element_text(size = base_size + 2, face = "bold"),
        axis.text = element_text(size = base_size, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size - 1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank()
      )
  }
  
  save_plot <- function(p, name, w = 10, h = 6) {
    if(!is.null(p)) {
      ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
    }
  }
  
  has_elisa_nonlin <- "ELISA_EU_mL_nonlin" %in% colnames(combined_summary)
  
   method_titer_cols <- list(
     "AICMix" = "Titer_Ref1",
     "Forced 4PL" = "Titer_Forced4PL",
     "Forced 5PL" = "Titer_5PL"
   )
  
  if("RM_IC50_log" %in% colnames(combined_summary) && !all(is.na(combined_summary$RM_IC50_log))) {
    method_titer_cols[["Reed-Muench"]] <- "RM_Titer_Ref1"
  }
  
  method_titer_cols <- method_titer_cols[sapply(method_titer_cols, function(c) c %in% colnames(combined_summary))]
  
  titer_matrix <- combined_summary[, "Sample", drop = FALSE]
  for(meth_name in names(method_titer_cols)) {
    titer_col <- method_titer_cols[[meth_name]]
    if(titer_col %in% colnames(combined_summary)) {
      titer_matrix[[meth_name]] <- combined_summary[[titer_col]]
    }
  }
  
  active_methods <- sapply(names(method_titer_cols), function(col) {
    col %in% colnames(titer_matrix) && any(!is.na(titer_matrix[[col]]))
  })
  method_titer_cols <- method_titer_cols[active_methods]
  
  if(length(method_titer_cols) == 0) {
    titer_matrix <- NULL
  } else {
    titer_matrix <- titer_matrix[, c("Sample", names(method_titer_cols)), drop = FALSE]
    titer_matrix <- titer_matrix[complete.cases(titer_matrix[, names(method_titer_cols), drop = FALSE]), ]
  }
  
  if(!is.null(titer_matrix) && nrow(titer_matrix) >= 2) {
    cor_matrix <- cor(log10(titer_matrix[, names(method_titer_cols), drop = FALSE]), use = "complete.obs")
    cor_df <- as.data.frame(as.table(cor_matrix))
    colnames(cor_df) <- c("Method1", "Method2", "Correlation")
    cor_df$Method1 <- factor(cor_df$Method1, levels = names(method_titer_cols))
    cor_df$Method2 <- factor(cor_df$Method2, levels = names(method_titer_cols))
    cor_df$Label <- ifelse(cor_df$Method1 == cor_df$Method2, "1.00", sprintf("%.2f", cor_df$Correlation))
    
    p1 <- ggplot(cor_df, aes(x = Method1, y = Method2, fill = Correlation)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Label), color = "black", size = 4, fontface = "bold") +
      scale_fill_gradient2(low = "#DC143C", mid = "#FFD700", high = "#228B22", midpoint = 0.9,
                           limits = c(0.5, 1), oob = squish) +
      labs(subx = "", y = "", fill = "Correlation") +
      insight_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")
    save_plot(p1, "Plot01_Method_Agreement_Matrix", w = 8, h = 6)
  }
  
  elisa_titer_data <- combined_summary[!is.na(combined_summary$ELISA_EU_mL) & !is.na(combined_summary$Titer_Ref1) & !combined_summary$Is_Ref1 & !combined_summary$Is_Ref2 & !combined_summary$Is_Ref3 & !combined_summary$Is_mAb, ]
  if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(elisa_titer_data)) {
    prefix_pattern <- paste0("^", sample_prefix_filter)
    elisa_titer_data <- elisa_titer_data[grepl(prefix_pattern, elisa_titer_data$Sample, ignore.case = FALSE), ]
  }
  
  if(nrow(elisa_titer_data) >= 3) {
    elisa_titer_data$log_elisa <- log10(elisa_titer_data$ELISA_EU_mL)
    elisa_titer_data$log_titer <- log10(elisa_titer_data$Titer_Ref1)
    
    lims <- range(c(elisa_titer_data$log_elisa, elisa_titer_data$log_titer), na.rm = TRUE)
    lims <- lims + c(-0.5, 0.5)
    
    p2 <- ggplot() +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[1], ymax = lims[2] + ed50_tolerance,
               fill = "#228B22", alpha = 0.15) +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[2] + ed50_tolerance, ymax = lims[2] + 1.0,
               fill = "#FFD700", alpha = 0.15) +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[2] + 1.0, ymax = lims[2],
               fill = "#DC143C", alpha = 0.15) +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[2] - ed50_tolerance, ymax = lims[2],
               fill = "#228B22", alpha = 0.15) +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[2] - 1.0, ymax = lims[2] - ed50_tolerance,
               fill = "#FFD700", alpha = 0.15) +
      annotate("rect", xmin = lims[1], xmax = lims[2], ymin = lims[2], ymax = lims[1],
               fill = "#DC143C", alpha = 0.15) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      geom_point(data = elisa_titer_data, aes(x = log_elisa, y = log_titer, color = Concordance_Ref1_ELISA),
                 size = 3, alpha = 0.8) +
      scale_color_manual(values = conc_colors, name = "Agreement") +
      scale_x_continuous(limits = lims) +
      scale_y_continuous(limits = lims) +
      labs(x = "log10(ELISA EU/mL)", y = paste0("log10(Titer ", standard_unit, ")")) +
      insight_theme() +
      theme(legend.position = "bottom")
    save_plot(p2, "Plot02_ELISA_Titer_Agreement_Zones", w = 8, h = 7)
  }
  
  if(!is.null(titer_matrix) && nrow(titer_matrix) >= 2) {
    scorecard_data <- data.frame(Method = names(method_titer_cols), stringsAsFactors = FALSE)
    scorecard_data$Concordance_Rate <- NA
    scorecard_data$Mean_Abs_Bias <- NA
    scorecard_data$Median_Abs_Bias <- NA
    
    for(i in seq_len(nrow(scorecard_data))) {
      meth <- scorecard_data$Method[i]
      titer_col <- method_titer_cols[[meth]]
      if(titer_col %in% colnames(combined_summary)) {
        sub <- combined_summary[!is.na(combined_summary[[titer_col]]) & !combined_summary$Is_Ref1 & !combined_summary$Is_Ref2 & !combined_summary$Is_Ref3 & !combined_summary$Is_mAb, ]
        if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(sub)) {
          prefix_pattern <- paste0("^", sample_prefix_filter)
          sub <- sub[grepl(prefix_pattern, sub$Sample, ignore.case = FALSE), ]
        }
        if(nrow(sub) > 0) {
          conc_col <- if(meth == "Reed-Muench") "Concordance_Ref1_ELISA_RM" else "Concordance_Ref1_ELISA"
          if(conc_col %in% colnames(sub)) {
            valid_conc <- sub[[conc_col]] %in% c("Concordant", "Discordant")
            scorecard_data$Concordance_Rate[i] <- mean(sub[[conc_col]][valid_conc] == "Concordant", na.rm = TRUE) * 100
          }
          
          if(meth != "Reed-Muench") {
            log_diff <- log10(sub[[titer_col]]) - log10(sub$ELISA_EU_mL)
            scorecard_data$Mean_Abs_Bias[i] <- mean(abs(log_diff), na.rm = TRUE)
            scorecard_data$Median_Abs_Bias[i] <- median(abs(log_diff), na.rm = TRUE)
          }
        }
      }
    }
    
    scorecard_data <- scorecard_data[!is.na(scorecard_data$Concordance_Rate) | !is.na(scorecard_data$Mean_Abs_Bias), ]
    
    if(nrow(scorecard_data) > 0) {
      p3a <- ggplot(scorecard_data, aes(x = reorder(Method, Concordance_Rate), y = Concordance_Rate, fill = Method)) +
        geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
        geom_text(aes(label = ifelse(!is.na(Concordance_Rate), paste0(round(Concordance_Rate, 0), "%"), "N/A")),
                  vjust = -0.5, size = 3.5, fontface = "bold") +
        scale_fill_manual(values = method_colors) +
        ylim(0, 105) +
        labs(x = "", y = "Concordant (%)") +
        insight_theme() +
        theme(legend.position = "none")
      save_plot(p3a, "Plot03a_Concordance_Rate", w = 7, h = 5)
      
      bias_data <- scorecard_data[!is.na(scorecard_data$Mean_Abs_Bias), ]
      if(nrow(bias_data) > 0) {
        p3b <- ggplot(bias_data, aes(x = reorder(Method, Mean_Abs_Bias), y = Mean_Abs_Bias, fill = Method)) +
          geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
          geom_text(aes(label = paste0(round(Mean_Abs_Bias, 2), " log")), vjust = -0.5, size = 3.5, fontface = "bold") +
          scale_fill_manual(values = method_colors) +
          labs(x = "", y = "Mean |log10(Titer) - log10(ELISA)|") +
          insight_theme() +
          theme(legend.position = "none")
        save_plot(p3b, "Plot03b_Mean_Bias", w = 7, h = 5)
      }
    }
  }
  
  if(use_dual_reference && "Titer_Ref2" %in% colnames(combined_summary)) {
    ref_data <- combined_summary[!is.na(combined_summary$Titer_Ref1) & !is.na(combined_summary$Titer_Ref2) & !combined_summary$Is_Ref1 & !combined_summary$Is_Ref2 & !combined_summary$Is_Ref3 & !combined_summary$Is_mAb, ]
    if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(ref_data)) {
      prefix_pattern <- paste0("^", sample_prefix_filter)
      ref_data <- ref_data[grepl(prefix_pattern, ref_data$Sample, ignore.case = FALSE), ]
    }
    
    if(nrow(ref_data) >= 2) {
      ref_data$Sample <- factor(ref_data$Sample, levels = ref_data$Sample[order(ref_data$Titer_Ref1)])
      ref_data$log_titer_ref1 <- log10(ref_data$Titer_Ref1)
      ref_data$log_titer_ref2 <- log10(ref_data$Titer_Ref2)
      ref_data$max_log <- pmax(ref_data$log_titer_ref1, ref_data$log_titer_ref2, na.rm = TRUE)
      ref_data$min_log <- pmin(ref_data$log_titer_ref1, ref_data$log_titer_ref2, na.rm = TRUE)
      ref_data$log_change <- ref_data$log_titer_ref2 - ref_data$log_titer_ref1
      
      p4 <- ggplot(ref_data) +
        geom_segment(aes(x = min_log, xend = max_log, y = Sample, yend = Sample), color = "grey70", linewidth = 1.2) +
        geom_point(aes(x = log_titer_ref1, y = Sample), color = method_colors["AICMix"], size = 3) +
        geom_point(aes(x = log_titer_ref2, y = Sample), color = method_colors["Forced 4PL"], size = 3) +
        geom_text(aes(x = log_titer_ref1, y = Sample, label = paste0("  ", reference_standard)), 
                  hjust = -0.1, size = 2.8, color = method_colors["AICMix"]) +
        geom_text(aes(x = log_titer_ref2, y = Sample, label = paste0(second_reference_standard, "  ")), 
                  hjust = 1.1, size = 2.8, color = method_colors["Forced 4PL"]) +
        labs(x = paste0("log10(Titer) (", standard_unit, ")"), y = "") +
        insight_theme() +
        theme(legend.position = "none")
      save_plot(p4, "Plot04_Reference_Dumbbell", w = 10, h = 6)
    }
  }
  
  robust_data <- titer_matrix
  if(!is.null(robust_data) && nrow(robust_data) >= 2) {
    robust_data$Min_Titer <- apply(robust_data[, names(method_titer_cols), drop = FALSE], 1, min, na.rm = TRUE)
    robust_data$Max_Titer <- apply(robust_data[, names(method_titer_cols), drop = FALSE], 1, max, na.rm = TRUE)
    robust_data$Titer_Range_log <- log10(robust_data$Max_Titer) - log10(robust_data$Min_Titer)
    robust_data$Mean_Titer_log <- apply(log10(robust_data[, names(method_titer_cols), drop = FALSE]), 1, mean, na.rm = TRUE)
    robust_data$Sample <- factor(robust_data$Sample, levels = robust_data$Sample[order(robust_data$Titer_Range_log)])
    
    p5 <- ggplot(robust_data, aes(x = Titer_Range_log, y = Sample)) +
      geom_vline(xintercept = 0.3, linetype = "dashed", color = "#228B22", linewidth = 0.8) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "#DC143C", linewidth = 0.8) +
      geom_segment(aes(x = 0, xend = Titer_Range_log, y = Sample, yend = Sample), color = "grey60", linewidth = 0.8) +
      geom_point(aes(color = ifelse(Titer_Range_log <= 0.3, "Robust", ifelse(Titer_Range_log <= 0.5, "Moderate", "Sensitive"))), size = 3) +
      scale_color_manual(values = c("Robust" = "#228B22", "Moderate" = "#FFD700", "Sensitive" = "#DC143C"), name = "Robustness") +
      scale_x_continuous(limits = c(0, NA)) +
      labs(subx = "Titer Range (log10 scale)", y = "") +
      insight_theme() +
      theme(legend.position = "bottom")
    save_plot(p5, "Plot05_Sample_Robustness", w = 9, h = 6)
  }
  
  conc_ladder_data <- combined_summary[!is.na(combined_summary$Concordance_Ref1_ELISA) & combined_summary$Concordance_Ref1_ELISA %in% c("Concordant", "Discordant"), ]
  if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(conc_ladder_data)) {
    prefix_pattern <- paste0("^", sample_prefix_filter)
    conc_ladder_data <- conc_ladder_data[grepl(prefix_pattern, conc_ladder_data$Sample, ignore.case = FALSE), ]
  }
  
  if(nrow(conc_ladder_data) > 0) {
    conc_ladder_data$Agreement_Level <- "All methods agree"
    
    for(i in seq_len(nrow(conc_ladder_data))) {
      conc_values <- c(conc_ladder_data$Concordance_Ref1_ELISA[i])
      if(has_elisa_nonlin && !is.na(conc_ladder_data$Concordance_Ref1_Nonlin[i])) {
        conc_values <- c(conc_values, conc_ladder_data$Concordance_Ref1_Nonlin[i])
      }
      if(use_dual_reference) {
        if(!is.na(conc_ladder_data$Concordance_Ref2_ELISA[i])) conc_values <- c(conc_values, conc_ladder_data$Concordance_Ref2_ELISA[i])
        if(has_elisa_nonlin && !is.na(conc_ladder_data$Concordance_Ref2_Nonlin[i])) conc_values <- c(conc_values, conc_ladder_data$Concordance_Ref2_Nonlin[i])
      }
      
      n_conc <- sum(conc_values == "Concordant", na.rm = TRUE)
      n_total <- sum(!is.na(conc_values))
      
      if(n_total > 0 && n_conc == n_total) {
        conc_ladder_data$Agreement_Level[i] <- "All agree"
      } else if(n_total > 0 && n_conc >= n_total / 2) {
        conc_ladder_data$Agreement_Level[i] <- "Most agree"
      } else {
        conc_ladder_data$Agreement_Level[i] <- "Methods disagree"
      }
    }
    
    ladder_summary <- as.data.frame(table(conc_ladder_data$Agreement_Level))
    colnames(ladder_summary) <- c("Level", "Count")
    ladder_summary$Level <- factor(ladder_summary$Level, levels = c("All agree", "Most agree", "Methods disagree"))
    ladder_summary$Pct <- ladder_summary$Count / sum(ladder_summary$Count) * 100
    
    p6 <- ggplot(ladder_summary, aes(x = "", y = Pct, fill = Level)) +
      geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
      geom_text(aes(label = paste0(Level, "\n", round(Pct, 0), "% (", Count, " samples)")),
                position = position_stack(vjust = 0.5), color = "white", size = 4, fontface = "bold") +
      scale_fill_manual(values = c("All agree" = "#228B22", "Most agree" = "#FFD700", "Methods disagree" = "#DC143C")) +
      coord_polar("y", start = 0) +
      labs(subfill = "") +
      insight_theme() +
      theme(axis.text.x = element_blank(), axis.ticks = element_blank())
    save_plot(p6, "Plot06_Concordance_Ladder", w = 7, h = 6)
  }
  
  if(!is.null(titer_matrix) && nrow(titer_matrix) >= 2) {
    long_titer <- titer_matrix %>%
      pivot_longer(cols = names(method_titer_cols), names_to = "Method", values_to = "Titer") %>%
      filter(!is.na(Titer))
    long_titer$log_Titer <- log10(long_titer$Titer)
    long_titer$Method <- factor(long_titer$Method, levels = names(method_titer_cols))
    
    p7 <- ggplot(long_titer, aes(x = Method, y = log_Titer, color = Method)) +
      geom_jitter(width = 0.25, size = 2.5, alpha = 0.7) +
      geom_boxplot(alpha = 0.2, outlier.shape = NA, linewidth = 0.6) +
      scale_color_manual(values = method_colors) +
      labs(subx = "", y = paste0("log10(Titer) (", standard_unit, ")")) +
      insight_theme() +
      theme(legend.position = "none")
    save_plot(p7, "Plot07_Titer_by_Method", w = 8, h = 6)
  }
  
  if(has_elisa_nonlin) {
    nonlin_data <- combined_summary[!is.na(combined_summary$ELISA_EU_mL) & !is.na(combined_summary$ELISA_EU_mL_nonlin) & !is.na(combined_summary$Titer_Ref1) & !combined_summary$Is_Ref1 & !combined_summary$Is_Ref2 & !combined_summary$Is_Ref3 & !combined_summary$Is_mAb, ]
    if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(nonlin_data)) {
      prefix_pattern <- paste0("^", sample_prefix_filter)
      nonlin_data <- nonlin_data[grepl(prefix_pattern, nonlin_data$Sample, ignore.case = FALSE), ]
    }
    
    if(nrow(nonlin_data) >= 3) {
      nonlin_data$log_elisa_lin <- log10(nonlin_data$ELISA_EU_mL)
      nonlin_data$log_elisa_nonlin <- log10(nonlin_data$ELISA_EU_mL_nonlin)
      
      p8 <- ggplot(nonlin_data, aes(x = log_elisa_lin, y = log_elisa_nonlin)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
        geom_point(size = 3, alpha = 0.7, color = "#0072B2") +
        geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80", alpha = 0.3, linewidth = 0.8) +
        scale_x_continuous(limits = c(0, NA)) +
        scale_y_continuous(limits = c(0, NA)) +
        labs(subx = "log10(Linear ELISA EU/mL)", y = "log10(Non-linear ELISA EU/mL)") +
        insight_theme()
      save_plot(p8, "Plot08_Linear_vs_Nonlinear_ELISA", w = 7, h = 6)
    }
  }
  
  if(!is.null(titer_matrix) && nrow(titer_matrix) >= 3) {
    rec_data <- titer_matrix
    rec_data$Recommended <- "AICMix"
    
    for(i in seq_len(nrow(rec_data))) {
      row_vals <- unlist(rec_data[i, names(method_titer_cols), drop = FALSE])
      row_vals <- row_vals[!is.na(row_vals)]
      if(length(row_vals) < 2) next
      
      row_log <- log10(row_vals)
      cv <- sd(row_log, na.rm = TRUE) / mean(row_log, na.rm = TRUE) * 100
      
      if(cv < 10) {
        rec_data$Recommended[i] <- "Robust (any method)"
      } else if(cv < 25) {
        rec_data$Recommended[i] <- "AICMix"
      } else {
        rec_data$Recommended[i] <- "Review needed"
      }
    }
    
    rec_summary <- as.data.frame(table(rec_data$Recommended))
    colnames(rec_summary) <- c("Recommendation", "Count")
    rec_summary$Recommendation <- factor(rec_summary$Recommendation, 
                                         levels = c("Robust (any method)", "AICMix", "Review needed"))
    
    p9 <- ggplot(rec_summary, aes(x = reorder(Recommendation, Count), y = Count, fill = Recommendation)) +
      geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
      geom_text(aes(label = paste0(Count, " samples")), vjust = -0.5, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = c("Robust (any method)" = "#228B22", "AICMix" = "#56B4E9", "Review needed" = "#DC143C")) +
      labs(subx = "", y = "Number of Samples") +
      insight_theme() +
      theme(legend.position = "none")
    save_plot(p9, "Plot09_Recommendation_Summary", w = 8, h = 5)
  }
  
  if(!is.null(titer_matrix) && nrow(titer_matrix) >= 3) {
    heatmap_data <- titer_matrix
    heatmap_data$Sample <- factor(heatmap_data$Sample, levels = heatmap_data$Sample[order(apply(log10(heatmap_data[, names(method_titer_cols), drop = FALSE]), 1, mean, na.rm = TRUE))])
    
    long_heatmap <- heatmap_data %>%
      pivot_longer(cols = names(method_titer_cols), names_to = "Method", values_to = "Titer") %>%
      filter(!is.na(Titer)) %>%
      mutate(log_Titer = log10(Titer))
    
    if(nrow(long_heatmap) > 0) {
      p10 <- ggplot(long_heatmap, aes(x = Method, y = Sample, fill = log_Titer)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = ifelse(!is.na(log_Titer), sprintf("%.2f", log_Titer), "")), 
                  color = "black", size = 2.5) +
        scale_fill_viridis(option = "D", name = "log10(Titer)") +
        labs(subx = "", y = "") +
        insight_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(p10, "Plot10_Titer_Heatmap", w = 8, h = 6)
    }
  }
  
  return(invisible(list(plot_dir = plot_dir)))
}