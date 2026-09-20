library(tidyverse)
library(ggplot2)
library(scales)
library(readxl)

source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

predict_elisacat <- function(
    neut_titer = NULL,
    file_path = NULL,
    sheet = "sample4PL",
    sample_id_column = "specify",
    titer_column = "specify",
    elisa_column = "specify",
    tolerance = 0.5,
    elisa_categories = c("<0.125 EU/mL", "0.125 to <0.5 EU/mL",
                         "0.5 to <2.0 EU/mL", "2.0 to <4.0 EU/mL", ">=4.0 EU/mL"),
    cat_bounds = list(
      c(0, 0.125),
      c(0.125, 0.5),
      c(0.5, 2.0),
      c(2.0, 4.0),
      c(4.0, Inf)
    ),
    output_dir = "predict_elisacat_results",
    generate_plots = TRUE
  ) {
    base_font_size <- 11
    dpi <- 300

    .assign_elisacat <- function(elisa_value) {
    if(is.na(elisa_value)) return(NA)
    if(elisa_value < 0.125) return("<0.125 EU/mL")
    if(elisa_value < 0.5) return("0.125 to <0.5 EU/mL")
    if(elisa_value < 2.0) return("0.5 to <2.0 EU/mL")
    if(elisa_value < 4.0) return("2.0 to <4.0 EU/mL")
    return(">=4.0 EU/mL")
  }
  
  .theme_pub <- function(bsize = base_font_size) {
    theme_minimal(bsize = bsize) +
      theme(
        text = element_text(family = "sans", color = "black"),
        axis.title = element_text(size = bsize + 2, face = "bold"),
        axis.text = element_text(size = bsize, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        plot.title = element_text(size = bsize + 3, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = bsize, hjust = 0.5, color = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = bsize, face = "bold"),
        legend.text = element_text(size = bsize - 1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank()
      )
  }
  
  .save_plot <- function(p, name, w = 8, h = 6) {
    if(!is.null(p) && generate_plots) {
      ggsave(file.path(plot_dir, paste0(name, ".png")), p,
             width = w, height = h, dpi = dpi, bg = "white")
    }
  }
  
  .compute_prediction <- function(titer_val) {
    if(is.na(titer_val) || titer_val <= 0) {
      return(list(
        neut_titer = titer_val,
        tolerance = tolerance,
        tolerance_lower = NA,
        tolerance_upper = NA,
        category_probabilities = setNames(rep(NA, length(elisa_categories)), elisa_categories),
        most_probable_category = NA,
        max_probability = NA,
        num_probable_categories = 0,
        prediction = NA,
        confidence_interval = NA
      ))
    }
    
    log_titer <- log10(titer_val)
    tol_lower <- 10^(log_titer - tolerance)
    tol_upper <- 10^(log_titer + tolerance)
    
    overlaps <- sapply(cat_bounds, function(b) {
      a_int <- max(tol_lower, b[1])
      b_int <- min(tol_upper, b[2])
      pmax(0, b_int - a_int)
    })
    
    total_overlap <- sum(overlaps, na.rm = TRUE)
    if(total_overlap <= 0) {
      return(list(
        neut_titer = titer_val,
        tolerance = tolerance,
        tolerance_lower = tol_lower,
        tolerance_upper = tol_upper,
        category_probabilities = setNames(rep(NA, length(elisa_categories)), elisa_categories),
        most_probable_category = NA,
        max_probability = NA,
        num_probable_categories = 0,
        prediction = NA,
        confidence_interval = NA
      ))
    }
    
    probs <- overlaps / total_overlap
    probs[is.nan(probs) | is.infinite(probs)] <- NA
    max_prob <- max(probs, na.rm = TRUE)
    most_prob_cat <- names(probs)[which.max(probs)]
    n_prob_cats <- sum(probs > 0.01, na.rm = TRUE)
    
    prediction <- if(max_prob >= 0.5) most_prob_cat else "Uncertain"
    confidence_interval <- if(max_prob >= 0.5) {
      paste0(round(max_prob * 100, 1), "% (", most_prob_cat, ")")
    } else {
      paste0("No category >50%; highest = ", round(max_prob * 100, 1), "% (", most_prob_cat, ")")
    }
    
    list(
      neut_titer = titer_val,
      tolerance = tolerance,
      tolerance_lower = tol_lower,
      tolerance_upper = tol_upper,
      category_probabilities = setNames(round(probs, 6), elisa_categories),
      most_probable_category = most_prob_cat,
      max_probability = round(max_prob, 6),
      num_probable_categories = n_prob_cats,
      prediction = prediction,
      confidence_interval = confidence_interval
    )
  }
  
  .load_and_prepare_data <- function() {
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_dir <<- file.path(output_dir, "plots")
    if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    
    data <- read_excel(file_path, sheet = sheet)
    data <- as.data.frame(data)
    colnames(data) <- gsub(" ", "_", colnames(data))
    colnames(data) <- gsub("/", "_", colnames(data))
    colnames(data) <- gsub("\\.", "_", colnames(data))
    
    data <- data[!is.na(data[[titer_column]]) & data[[titer_column]] > 0, ]
    data$Reference_elisacat <- sapply(data[[elisa_column]], .assign_elisacat)
    data$Sample_ID <- data[[sample_id_column]]
    
    data
  }
  
  .compute_all_predictions <- function(data) {
    results <- list()
    for(i in 1:nrow(data)) {
      titer <- data[[titer_column]][i]
      pred <- .compute_prediction(titer)
      row <- data[i, ]
      results[[i]] <- data.frame(
        Sample_ID = row$Sample_ID,
        Measured_Titer = titer,
        Tolerance_Lower = pred$tolerance_lower,
        Tolerance_Upper = pred$tolerance_upper,
        Predicted_Category = pred$most_probable_category,
        Confidence = pred$max_probability,
        N_Probable_Categories = pred$num_probable_categories,
        Prediction = pred$prediction,
        Reference_elisacat = row$Reference_elisacat,
        Reference_ELISA_Value = row[[elisa_column]],
        stringsAsFactors = FALSE
      )
    }
    
    results_df <- do.call(rbind, results)
    rownames(results_df) <- NULL
    
    prob_col_names <- c("ProbCLT0.125", "Prob_Cat_0.125_to_0.5", "Prob_Cat_0.5_to_2.0",
                        "Prob_Cat_2.0_to_4.0", "Prob_Cat_GTE_4.0")
    for(j in seq_along(prob_col_names)) {
      results_df[[prob_col_names[j]]] <- NA_real_
    }
    
    for(i in 1:nrow(results_df)) {
      titer <- results_df$Measured_Titer[i]
      if(is.na(titer) || titer <= 0) next
      pred <- .compute_prediction(titer)
      for(j in seq_along(elisa_categories)) {
        results_df[i, prob_col_names[j]] <- pred$category_probabilities[elisa_categories[j]]
      }
    }
    
    results_df$Concordance_Status <- ifelse(
      is.na(results_df$Predicted_Category) | is.na(results_df$Reference_elisacat) | is.na(results_df$Reference_ELISA_Value),
      NA,
      ifelse(results_df$Reference_ELISA_Value >= results_df$Tolerance_Lower & results_df$Reference_ELISA_Value <= results_df$Tolerance_Upper,
             "Concordant", "Discordant")
    )
    results_df$Prob_Cat_Concord <- NA_character_
    for(i in 1:nrow(results_df)) {
      ref_val <- results_df$Reference_ELISA_Value[i]
      tol_lower <- results_df$Tolerance_Lower[i]
      tol_upper <- results_df$Tolerance_Upper[i]
      ref_cat <- results_df$Reference_elisacat[i]
      
      if(is.na(ref_val) || is.na(tol_lower) || is.na(tol_upper) || is.na(ref_cat)) {
        results_df$Prob_Cat_Concord[i] <- NA
      } else if(ref_val >= tol_lower && ref_val <= tol_upper) {
        prob <- results_df[i, prob_col_names[match(ref_cat, elisa_categories)]]
        if(!is.na(prob) && prob > 0) {
          results_df$Prob_Cat_Concord[i] <- paste0(ref_cat, ", ", round(prob * 100, 1), "%")
        } else {
          results_df$Prob_Cat_Concord[i] <- "Discordant"
        }
      } else {
        results_df$Prob_Cat_Concord[i] <- "Discordant"
      }
    }
    results_df$Concordance_Probability <- ifelse(
      is.na(results_df$Predicted_Category) | is.na(results_df$Reference_elisacat),
      NA,
      ifelse(results_df$Predicted_Category == results_df$Reference_elisacat,
             results_df$Confidence, 0)
    )
    
    results_df
  }
  
  .compute_overall_summary <- function(results_df) {
    data.frame(
      Metric = c("N_Samples", "N_With_Reference_ELISA", "N_Correctly_Predicted",
                  "Prob_Cat_Concord_Rate_percent", "Mean_Confidence",
                 "N_MultiCategory", "Pct_MultiCategory", "N_Uncertain",
                 "Pct_Uncertain", "Mean_Titer_EU_mL", "Tolerance_log10"),
      Value = c(
        nrow(results_df),
        sum(!is.na(results_df$Reference_elisacat)),
        sum(results_df$Predicted_Category == results_df$Reference_elisacat, na.rm = TRUE),
        round(mean(results_df$Predicted_Category == results_df$Reference_elisacat, na.rm = TRUE) * 100, 1),
        round(mean(results_df$Confidence, na.rm = TRUE), 4),
        sum(results_df$N_Probable_Categories > 1, na.rm = TRUE),
        round(mean(results_df$N_Probable_Categories > 1, na.rm = TRUE) * 100, 1),
        sum(results_df$Prediction == "Uncertain", na.rm = TRUE),
        round(mean(results_df$Prediction == "Uncertain", na.rm = TRUE) * 100, 1),
        round(mean(results_df$Measured_Titer, na.rm = TRUE), 4),
        tolerance
      ),
      stringsAsFactors = FALSE
    )
  }
  
  .compute_category_breakdown <- function(results_df) {
    by_predicted_cat <- results_df %>%
      group_by(Predicted_Category) %>%
      summarise(
        N_Samples = n(),
        Mean_Confidence = mean(Confidence, na.rm = TRUE),
        Min_Confidence = min(Confidence, na.rm = TRUE),
        Max_Confidence = max(Confidence, na.rm = TRUE),
        N_MultiCategory = sum(N_Probable_Categories > 1, na.rm = TRUE),
        .groups = "drop"
      )
    
    by_reference_cat <- results_df %>%
      group_by(Reference_elisacat) %>%
      summarise(
        N_Samples = n(),
        N_Correctly_Predicted = sum(Predicted_Category == Reference_elisacat, na.rm = TRUE),
        Prediction_Rate_percent = round(mean(Predicted_Category == Reference_elisacat, na.rm = TRUE) * 100, 1),
        Mean_Confidence = mean(Confidence, na.rm = TRUE),
        .groups = "drop"
      )
    
    list(by_predicted_cat = by_predicted_cat, by_reference_cat = by_reference_cat)
  }
  
  .generate_confidence_by_titer_plot <- function(results_df) {
    category_colors <- c(
      "<0.125 EU/mL" = "#56B4E9",
      "0.125 to <0.5 EU/mL" = "#0072B2",
      "0.5 to <2.0 EU/mL" = "#E69F00",
      "2.0 to <4.0 EU/mL" = "#D55E00",
      ">=4.0 EU/mL" = "#CC79A7"
    )
    
    p1 <- ggplot(results_df, aes(x = Measured_Titer, y = Confidence, color = Predicted_Category)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.8) +
      scale_x_log10(labels = scales::comma) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_color_manual(values = category_colors, name = "Predicted Category") +
      labs(x = "Measured Neutralization Titer (EU/mL)", y = "Classification Confidence") +
      .theme_pub()
    .save_plot(p1, "Fig01_Confidence_by_Titer", w = 10, h = 6)
  }
  
  .generate_confidence_distribution_plot <- function(results_df) {
    p2 <- ggplot(results_df, aes(x = Confidence)) +
      geom_histogram(fill = "#0072B2", color = "black", alpha = 0.7, bins = 20) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "#D55E00", linewidth = 1) +
      annotate("text", x = 0.5, y = Inf, label = "Decision threshold (50%)",
               vjust = 1.5, hjust = 0, size = 3.5, color = "#D55E00") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1),
                         labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Maximum Category Probability", y = "Frequency") +
      .theme_pub()
    .save_plot(p2, "Fig02_Confidence_Distribution", w = 8, h = 6)
  }
  
  .generate_multi_category_plot <- function(results_df) {
    results_df$Assignment_Type <- factor(
      ifelse(results_df$N_Probable_Categories == 1, "Single",
             ifelse(results_df$N_Probable_Categories == 2, "Two", "Three or More")),
      levels = c("Single", "Two", "Three or More")
    )
    p3 <- ggplot(results_df, aes(x = Assignment_Type, fill = Assignment_Type)) +
      geom_bar(color = "black", linewidth = 0.3, width = 0.6) +
      geom_text(stat = "count", aes(label = after_stat(count)),
                vjust = -0.5, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = c("Single" = "#228B22", "Two" = "#D55E00",
                                  "Three or More" = "#DC143C")) +
      labs(x = "Assignment Type", y = "Number of Samples") +
      .theme_pub() +
      theme(legend.position = "none")
    .save_plot(p3, "Fig03_MultiCategory_Breakdown", w = 8, h = 6)
  }
  
  .generate_heatmap_plot <- function(results_df) {
    if(nrow(results_df) <= 150) {
      p4_data_long <- results_df %>%
        tidyr::pivot_longer(cols = all_of(prob_col_names),
                            names_to = "Category",
                            values_to = "Probability") %>%
        dplyr::mutate(Category = gsub("Prob_Cat_", "", Category),
                      Category = gsub("_", " ", Category),
                      Category = gsub("LT 0\\.125", "<0.125", Category),
                      Category = gsub("0\\.125 to 0\\.5", "0.125 to <0.5", Category),
                      Category = gsub("0\\.5 to 2\\.0", "0.5 to <2.0", Category),
                      Category = gsub("2\\.0 to 4\\.0", "2.0 to <4.0", Category),
                      Category = gsub("GTE 4\\.0", ">=4.0", Category))
      
      p4 <- ggplot(p4_data_long, aes(x = Category, y = Sample_ID, fill = Probability)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(low = "#56B4E9", mid = "#FFD700", high = "#DC143C",
                             midpoint = 0.5, limits = c(0, 1),
                             breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                             labels = c("0%", "25%", "50%", "75%", "100%")) +
         labs(x = "ELISA Titer Category", y = "Sample", fill = "Probability") +
        .theme_pub() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      .save_plot(p4, "Fig04_Category_Heatmap", w = 10, h = min(10, nrow(results_df) * 0.15))
    }
  }
  
  .generate_predicted_vs_reference_plot <- function(results_df) {
    if(!all(is.na(results_df$Reference_elisacat))) {
      cat_to_numeric <- function(cat) {
        if(is.na(cat)) return(NA)
        match(cat, elisa_categories)
      }
      results_df$Predicted_Num <- sapply(results_df$Predicted_Category, cat_to_numeric)
      results_df$Reference_Num <- sapply(results_df$Reference_elisacat, cat_to_numeric)
      
      p5 <- ggplot(results_df, aes(x = Reference_Num, y = Predicted_Num, color = Concordance_Status)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
        geom_point(size = 3, alpha = 0.7) +
        scale_x_continuous(breaks = 1:5, labels = elisa_categories) +
        scale_y_continuous(breaks = 1:5, labels = elisa_categories) +
        scale_color_manual(values = c("Concordant" = "#228B22", "Discordant" = "#DC143C")) +
         labs(x = "Reference ELISA Category", y = "Predicted Category", color = "Concordance") +
        .theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      .save_plot(p5, "Fig05_Predicted_vs_Reference", w = 10, h = 8)
    }
  }
  
  .generate_scorecard_plot <- function(results_df) {
    scorecard_data <- data.frame(
      Metric = c("Total Samples", "Concordance Rate", "Mean Confidence",
                 "Multi-Category %", "Uncertain %", "Tolerance"),
      Plot_Value = c(
        nrow(results_df),
        round(mean(results_df$Predicted_Category == results_df$Reference_elisacat, na.rm = TRUE) * 100, 1),
        round(mean(results_df$Confidence, na.rm = TRUE) * 100, 1),
        round(mean(results_df$N_Probable_Categories > 1, na.rm = TRUE) * 100, 1),
        round(mean(results_df$Prediction == "Uncertain", na.rm = TRUE) * 100, 1),
        tolerance
      ),
      Label = c(
        paste0(nrow(results_df), " samples"),
        paste0(round(mean(results_df$Predicted_Category == results_df$Reference_elisacat, na.rm = TRUE) * 100, 1), "%"),
        paste0(round(mean(results_df$Confidence, na.rm = TRUE) * 100, 1), "%"),
        paste0(round(mean(results_df$N_Probable_Categories > 1, na.rm = TRUE) * 100, 1), "%"),
        paste0(round(mean(results_df$Prediction == "Uncertain", na.rm = TRUE) * 100, 1), "%"),
        paste0(tolerance, " log10")
      ),
      stringsAsFactors = FALSE
    )
    
    p6 <- ggplot(scorecard_data, aes(x = reorder(Metric, Plot_Value), y = Plot_Value)) +
      geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.8, width = 0.6) +
      geom_text(aes(label = Label), vjust = -0.5, size = 3.5, fontface = "bold") +
      coord_flip() +
      labs(x = "", y = "Value") +
      .theme_pub() +
      theme(legend.position = "none")
    .save_plot(p6, "Fig06_Summary_Scorecard", w = 8, h = 5)
  }
  
  .generate_all_plots <- function(results_df) {
    if(!generate_plots || nrow(results_df) == 0) return(invisible(NULL))
    
    .generate_confidence_by_titer_plot(results_df)
    .generate_confidence_distribution_plot(results_df)
    .generate_multi_category_plot(results_df)
    .generate_heatmap_plot(results_df)
    .generate_predicted_vs_reference_plot(results_df)
    .generate_scorecard_plot(results_df)
    
    invisible(list(
      prediction_results = results_df,
      output_dir = output_dir
    ))
  }
  
  .validate_inputs()
  
  if(!is.null(neut_titer)) {
    if(length(neut_titer) == 1) {
      return(.compute_prediction(neut_titer))
    } else {
      results <- lapply(neut_titer, .compute_prediction)
      return(results)
    }
  }
  
  data <- .load_and_prepare_data()
  results_df <- .compute_all_predictions(data)
  
  overall_summary <- .compute_overall_summary(results_df)
  category_breakdown <- .compute_category_breakdown(results_df)
  
  write.csv(results_df, file.path(output_dir, "prediction_results.csv"), row.names = FALSE)
  write.csv(overall_summary, file.path(output_dir, "overall_summary.csv"), row.names = FALSE)
  write.csv(category_breakdown$by_predicted_cat, file.path(output_dir, "by_predicted_category.csv"), row.names = FALSE)
  write.csv(category_breakdown$by_reference_cat, file.path(output_dir, "by_reference_category.csv"), row.names = FALSE)
  
  excel_sheets <- list(
    Prediction_Results = results_df,
    Overall_Summary = overall_summary,
    By_Predicted_Category = category_breakdown$by_predicted_cat,
    By_Reference_Category = category_breakdown$by_reference_cat
  )
  excel_sheets <- Filter(function(x) is.data.frame(x) && ncol(x) > 0, excel_sheets)
  write_xlsx(excel_sheets, file.path(output_dir, "predict_elisacat_results.xlsx"))
  
  .generate_all_plots(results_df)
  
  invisible(list(
    prediction_results = results_df,
    overall_summary = overall_summary,
    by_predicted_category = category_breakdown$by_predicted_cat,
    by_reference_category = category_breakdown$by_reference_cat,
    output_dir = output_dir
  ))
}
