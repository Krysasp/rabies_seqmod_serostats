library(tidyverse)
library(ggplot2)
library(grid)
library(readxl)
library(writexl)
library(scales)
library(splines)

source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

normalize_luminescence <- function(input_file,
                                    output_file = NULL,
                                    reference_experiment = NULL,
                                    curve_method = "log_log",
                                    plot_results = TRUE) {
  plot_params <- .prepare_plot_params(list())
  df <- .load_and_validate_data(input_file)
  fixed_gain_cols <- .get_fixed_gain_cols(df)
  auto_cols <- .get_auto_cols(df)
  experiments <- unique(df$experiment)
  reference_experiment <- .resolve_reference_experiment(experiments, reference_experiment)
  ref_data <- df[df$experiment == reference_experiment, ]
  output_df <- df

  for(auto_col in auto_cols) {
    norm_col_name <- paste0(auto_col, "_rev")
    output_df[[norm_col_name]] <- NA
    ref_indices <- which(output_df$experiment == reference_experiment)
    output_df[ref_indices, norm_col_name] <- output_df[ref_indices, auto_col]

    ref_auto <- c(); ref_gain <- c()
    for(i in seq_len(nrow(ref_data))) {
      auto_val <- ref_data[[auto_col]][i]
      if(!is.na(auto_val) && auto_val > 0) {
        for(gain_col in fixed_gain_cols) {
          gain_val <- ref_data[[gain_col]][i]
          if(!is.na(gain_val) && gain_val > 0) {
            ref_auto <- c(ref_auto, auto_val)
            ref_gain <- c(ref_gain, gain_val)
          }
        }
      }
    }
    if(length(ref_auto) < 3) next

    curve_model <- .build_reference_curve(ref_auto, ref_gain, curve_method)

    for(exp in experiments[experiments != reference_experiment]) {
      exp_indices <- which(output_df$experiment == exp)
      exp_data <- output_df[exp_indices, ]
      if(nrow(exp_data) > 0) {
        for(i in seq_len(nrow(exp_data))) {
          row_idx <- exp_indices[i]
          gain_vals <- c()
          for(gain_col in fixed_gain_cols) {
            gain_val <- exp_data[[gain_col]][i]
            if(!is.na(gain_val) && gain_val > 0) gain_vals <- c(gain_vals, gain_val)
          }
          if(length(gain_vals) > 0) {
            avg_gain <- mean(gain_vals, na.rm = TRUE)
            predicted_auto <- .predict_from_curve(curve_model, avg_gain, curve_method)
            if(!is.na(predicted_auto) && predicted_auto > 0) {
              output_df[row_idx, norm_col_name] <- predicted_auto
            }
          }
        }
      }
    }

    if(plot_results && plot_params$show_curve_plot) {
      p_curve <- .generate_curve_plot(ref_gain, ref_auto, curve_model, curve_method, auto_col, plot_params)
      assign(paste0("p_curve_", gsub(" ", "_", auto_col)), p_curve, envir = environment())
      if(plot_params$save_plots) .save_curve_plot(p_curve, auto_col, plot_params)
    }
  }

  if(plot_results) {
    for(auto_col in auto_cols) {
      norm_col <- paste0(auto_col, "_rev")
      if(sum(!is.na(output_df[[auto_col]])) > 0 && sum(!is.na(output_df[[norm_col]])) > 0) {
        p <- .generate_distribution_plot(output_df, auto_col, norm_col, reference_experiment, curve_method, plot_params)
        print(p)
        if(plot_params$save_plots) .save_distribution_plot(p, auto_col, plot_params)
      }
    }
  }

  if(plot_results && plot_params$show_combined_plot && length(auto_cols) > 0) {
    for(auto_col in auto_cols) {
      norm_col <- paste0(auto_col, "_rev")
      if(sum(!is.na(output_df[[auto_col]])) > 0 && sum(!is.na(output_df[[norm_col]])) > 0) {
        p_norm <- .generate_distribution_plot(output_df, auto_col, norm_col, reference_experiment, curve_method, plot_params)
        curve_plot_name <- paste0("p_curve_", gsub(" ", "_", auto_col))
        curve_plot <- get(curve_plot_name, envir = environment())
        if(!is.null(curve_plot)) {
          combined_plot <- p_norm | curve_plot
          print(combined_plot)
          if(plot_params$save_plots) .save_combined_plot(combined_plot, auto_col, plot_params)
        }
      }
    }
  }

  if(!is.null(output_file)) {
    write_xlsx(output_df, output_file)
  }

  .print_normalization_summary(output_df, reference_experiment, curve_method, fixed_gain_cols, auto_cols, experiments)
  return(invisible(output_df))
}

.prepare_plot_params <- function(plot_params) {
  default_plot_params <- list(
    title_size = 14, subtitle_size = 12, axis_title_size = 12, axis_text_size = 10,
    legend_title_size = 12, legend_text_size = 10, strip_text_size = 11,
    original_color = "#E69F00", normalized_color = "#56B4E9", curve_color = "#D55E00",
    point_size = 1.5, point_alpha = 0.6, box_alpha = 0.3, box_width = 0.7, jitter_width = 0.2,
    scale_type = "log10", y_axis_limits = NULL, use_integer_labels = TRUE,
    base_family = "Arial", legend_position = "right", axis_line_width = 0.5, panel_grid = FALSE,
    save_plots = FALSE, plot_width = 8, plot_height = 6, plot_dpi = 300, plot_format = "png", plot_dir = "plots",
    show_curve_plot = TRUE, curve_confidence_interval = 0.95,
    show_combined_plot = TRUE, combined_plot_width = 16, combined_plot_height = 6
  )
  if(length(plot_params) > 0) {
    for(param in names(plot_params)) {
      if(param %in% names(default_plot_params)) {
        default_plot_params[[param]] <- plot_params[[param]]
      }
    }
  }
  default_plot_params
}

.load_and_validate_data <- function(input_file) {
  df <- read_excel(input_file)
  if(!"experiment" %in% colnames(df)) stop("Input file must contain an 'experiment' column")
  df
}

.get_fixed_gain_cols <- function(df) {
  cols <- colnames(df)[grepl("gain", colnames(df), ignore.case = TRUE)]
  if(length(cols) == 0) stop("No fixed gain columns found (columns containing 'gain')")
  cols
}

.get_auto_cols <- function(df) {
  cols <- colnames(df)[grepl("^Auto", colnames(df))]
  if(length(cols) == 0) stop("No Auto columns found (columns starting with 'Auto')")
  cols
}

.resolve_reference_experiment <- function(experiments, reference_experiment) {
  if(is.null(reference_experiment)) {
    reference_experiment <- experiments[1]
  }
  if(!reference_experiment %in% experiments) {
    stop(paste("Reference experiment", reference_experiment, "not found in data"))
  }
  reference_experiment
}



.generate_curve_plot <- function(ref_gain, ref_auto, curve_model, curve_method, auto_col, plot_params) {
  gain_range <- seq(min(ref_gain, na.rm = TRUE), max(ref_gain, na.rm = TRUE), length.out = 100)
  pred_auto <- sapply(gain_range, function(g) .predict_from_curve(curve_model, g, curve_method))

  curve_data <- data.frame(gain = ref_gain, auto = ref_auto, experiment = reference_experiment)
  pred_data <- data.frame(gain = gain_range, auto = pred_auto)
  curve_type <- switch(curve_method,
    linear = "Linear",
    log_log = "Log-Log",
    polynomial = "Polynomial (degree 2)",
    spline = "Natural Spline",
    curve_method
  )

  ggplot() +
    geom_point(data = curve_data, aes(x = gain, y = auto), color = plot_params$original_color, size = plot_params$point_size, alpha = 0.7) +
    geom_line(data = pred_data, aes(x = gain, y = auto), color = plot_params$curve_color, size = 1.2) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = comma, limits = c(1000, NA)) +
    labs(title = NULL,
         subtitle = NULL,
         x = "Gain Value (average of gain columns)", y = "RLU (predicted)") +
    theme_bw() +
    theme(
      axis.text = element_text(size = plot_params$axis_text_size),
      axis.title = element_text(size = plot_params$axis_title_size, face = "bold"),
      plot.title = element_text(size = plot_params$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = plot_params$subtitle_size, hjust = 0.5),
      panel.grid = element_blank()
    )
}

.generate_distribution_plot <- function(output_df, auto_col, norm_col, reference_experiment, curve_method, plot_params) {
  plot_df <- data.frame(
    experiment = output_df$experiment,
    original = output_df[[auto_col]],
    normalized = output_df[[norm_col]]
  )
  plot_df <- plot_df[!is.na(plot_df$original) | !is.na(plot_df$normalized), ]

  original_data <- data.frame(experiment = plot_df$experiment, value = plot_df$original, type = "Original")
  original_data <- original_data[!is.na(original_data$value), ]
  normalized_data <- data.frame(experiment = plot_df$experiment, value = plot_df$normalized, type = "Normalized")
  normalized_data <- normalized_data[!is.na(normalized_data$value), ]
  plot_data <- rbind(original_data, normalized_data)

  if(nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = experiment, y = value, fill = type, color = type)) +
    geom_boxplot(alpha = plot_params$box_alpha, width = 0.35, outlier.shape = NA, position = position_dodge(width = 0.7)) +
    geom_jitter(alpha = plot_params$point_alpha, size = plot_params$point_size, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15, jitter.height = 0)) +
    scale_fill_manual(values = setNames(c(plot_params$original_color, plot_params$normalized_color), c("Original", "Normalized"))) +
    scale_color_manual(values = setNames(c(plot_params$original_color, plot_params$normalized_color), c("Original", "Normalized")))

  y_limits <- plot_params$y_axis_limits
  if(is.null(y_limits)) y_limits <- c(1000, 2000000)

  if(plot_params$scale_type == "linear") {
    p <- if(plot_params$use_integer_labels) p + scale_y_continuous(labels = comma, limits = y_limits) else p + scale_y_continuous(limits = y_limits)
  } else if(plot_params$scale_type == "log2") {
    p <- if(plot_params$use_integer_labels) p + scale_y_log2(labels = comma, limits = y_limits) else p + scale_y_log2(limits = y_limits)
  } else {
    p <- if(plot_params$use_integer_labels) p + scale_y_log10(labels = comma, limits = y_limits) else p + scale_y_log10(limits = y_limits)
  }

  p + labs(
    title = NULL,
    subtitle = NULL,
    x = "Experiment", y = paste0("RLU (", plot_params$scale_type, " scale)"),
    fill = "Data type", color = "Data type"
  ) +
  theme_bw(base_family = plot_params$base_family) +
  theme(
    axis.line = element_line(color = "black", size = plot_params$axis_line_width),
    axis.ticks = element_line(color = "black", size = plot_params$axis_line_width),
    axis.text = element_text(size = plot_params$axis_text_size, color = "black"),
    axis.title = element_text(size = plot_params$axis_title_size, color = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = plot_params$axis_line_width),
    panel.grid = if(plot_params$panel_grid) element_line(color = "gray90", size = 0.3) else element_blank(),
    panel.grid.minor = if(plot_params$panel_grid) element_line(color = "gray95", size = 0.15) else element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = plot_params$legend_position,
    legend.title = element_text(size = plot_params$legend_title_size, face = "bold"),
    legend.text = element_text(size = plot_params$legend_text_size),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing = unit(0.1, "cm"),
    plot.title = element_text(size = plot_params$title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = plot_params$subtitle_size, hjust = 0.5),
    strip.text = element_text(size = plot_params$strip_text_size, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black", size = 0.5)
  )
}

.save_curve_plot <- function(p_curve, auto_col, plot_params) {
  if(!dir.exists(plot_params$plot_dir)) dir.create(plot_params$plot_dir, recursive = TRUE)
  curve_filename <- paste0(auto_col, "_standard_curve.", plot_params$plot_format)
  curve_path <- file.path(plot_params$plot_dir, curve_filename)
  ggsave(curve_path, plot = p_curve, width = plot_params$plot_width, height = plot_params$plot_height, dpi = plot_params$plot_dpi, units = "in")
}

.save_distribution_plot <- function(p, auto_col, plot_params) {
  if(!dir.exists(plot_params$plot_dir)) dir.create(plot_params$plot_dir, recursive = TRUE)
  plot_filename <- paste0(auto_col, "_normalization.", plot_params$plot_format)
  plot_path <- file.path(plot_params$plot_dir, plot_filename)
  ggsave(plot_path, plot = p, width = plot_params$plot_width, height = plot_params$plot_height, dpi = plot_params$plot_dpi, units = "in")
}

.save_combined_plot <- function(combined_plot, auto_col, plot_params) {
  if(!dir.exists(plot_params$plot_dir)) dir.create(plot_params$plot_dir, recursive = TRUE)
  combined_filename <- paste0(auto_col, "_combined.", plot_params$plot_format)
  combined_path <- file.path(plot_params$plot_dir, combined_filename)
  ggsave(combined_path, plot = combined_plot, width = plot_params$combined_plot_width, height = plot_params$combined_plot_height, dpi = plot_params$plot_dpi, units = "in")
}

.print_normalization_summary <- function(output_df, reference_experiment, curve_method, fixed_gain_cols, auto_cols, experiments) {
  cat("\n=== Normalization Summary ===\n")
  cat("Reference experiment:", reference_experiment, "\n")
  cat("Curve method:", curve_method, "\n")
  cat("Fixed gain columns used:", paste(fixed_gain_cols, collapse = ", "), "\n")
  cat("Auto columns normalized:", paste(auto_cols, collapse = ", "), "\n")
  cat("Experiments:", paste(experiments, collapse = ", "), "\n")

  cat("\n=== Normalization Statistics ===\n")
  for(auto_col in auto_cols) {
    norm_col <- paste0(auto_col, "_rev")
    ref_vals <- output_df[output_df$experiment == reference_experiment, norm_col]
    ref_vals <- ref_vals[!is.na(ref_vals)]
    if(length(ref_vals) > 0) {
      ref_mean <- mean(ref_vals)
      ref_sd <- sd(ref_vals)
      cat("\n", auto_col, ":\n")
      cat("  Reference (", reference_experiment, ") mean:", round(ref_mean, 0), "±", round(ref_sd, 0), "\n")
      for(exp in experiments[experiments != reference_experiment]) {
        exp_vals <- output_df[output_df$experiment == exp, norm_col]
        exp_vals <- exp_vals[!is.na(exp_vals)]
        if(length(exp_vals) > 0) {
          exp_mean <- mean(exp_vals)
          exp_sd <- sd(exp_vals)
          fold_diff <- exp_mean / ref_mean
          cat("  ", exp, " mean:", round(exp_mean, 0), "±", round(exp_sd, 0), " (", round(fold_diff, 3), "x of reference)\n")
        }
      }
    }
  }
}
