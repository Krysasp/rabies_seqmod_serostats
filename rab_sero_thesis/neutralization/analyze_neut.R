library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(RColorBrewer)
library(readxl)
source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

analyze_neut <- function(
    file_path = "your_excel.xlsx",
    sheet = 3,
    out_dir = "neutres",
    exclude_samples = NULL,
    exclude_groups = NULL,
    titer_threshold = 0.5,
    bdl = 0.125,
    sample_results_file = NULL,
    concordance_type = "linear",
    titer_replicates = "NULL",
    color_by_groupW = FALSE,
    groupW_sheet = 4,
    groupW_name = 'ELISA_Titer',
    sample_colors = NULL,
    sample_linetypes = NULL
  ) {
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  

  .norm <- function(x) if(is.null(x)) NULL else trimws(as.character(x))
  .key <- function(x) {
    if(is.null(x)) return(NULL)
    gsub("[^a-z0-9]", "", tolower(as.character(x)))
  }
  sample_colors_n <- if(is.null(sample_colors)) NULL else setNames(sample_colors, .norm(names(sample_colors)))
  sample_linetypes_n <- if(is.null(sample_linetypes)) NULL else setNames(sample_linetypes, .norm(names(sample_linetypes)))
  .scol <- function(s) { ks <- .norm(s); if(is.null(sample_colors_n) || !ks %in% names(sample_colors_n)) NULL else sample_colors_n[[ks]] }
  .slty <- function(s) { ks <- .norm(s); if(is.null(sample_linetypes_n) || !ks %in% names(sample_linetypes_n)) NULL else sample_linetypes_n[[ks]] }
  .has_scol <- function(s) !is.null(.scol(s))
  .has_slty <- function(s) !is.null(.slty(s))
  
  sample_concordance <- NULL
  sample_titer_3rep <- NULL
  if(!is.null(sample_results_file) && file.exists(sample_results_file)) {
    sr_data <- read.csv(sample_results_file, stringsAsFactors = FALSE)
    
    if(concordance_type == "nonlin") {
      conc_col <- "Concordance_R1_Nonlin"
    } else {
      conc_col <- "Concordance_R1_ELISA"
    }
    if(concordance_type == "nonlin" && "Concordance_R1_Nonlin" %in% colnames(sr_data)) {
      sample_concordance <- setNames(sr_data$Concordance_R1_Nonlin, sr_data$Sample)
    } else if("Concordance_R1_ELISA" %in% colnames(sr_data)) {
      sample_concordance <- setNames(sr_data$Concordance_R1_ELISA, sr_data$Sample)
    }
    
    titer_file <- gsub("sample_results", "sample_results_three_replicates", sample_results_file)
    if(titer_replicates == "3" && file.exists(titer_file)) {
      titer_data <- read.csv(titer_file, stringsAsFactors = FALSE)
      sample_titer_3rep <- setNames(titer_data$Titer_R1, titer_data$Sample)
    }
  }
  
  sample_groupW <- NULL
  if(color_by_groupW || !is.null(exclude_groups)) {
    if(file.exists(file_path)) {
      meta_raw <- read_excel(file_path, sheet = groupW_sheet)
      meta_df <- as.data.frame(meta_raw)
      
      sample_col_meta <- NULL
      groupW_col_meta <- NULL
      cn_lower <- tolower(colnames(meta_df))
      gw_key <- gsub("[^a-z0-9]", "", tolower(as.character(groupW_name)))
      for(i in seq_along(colnames(meta_df))) {
        cn <- colnames(meta_df)[i]
        if(cn_lower[i] == "sample" || cn_lower[i] == "sample_id") sample_col_meta <- cn
        if(cn_lower[i] == tolower(groupW_name) ||
           gsub("[^a-z0-9]", "", cn_lower[i]) == gw_key) {
          groupW_col_meta <- cn
        }
      }
      if(is.null(groupW_col_meta)) {
        warning("Could not find a group column matching groupW_name='", groupW_name,
                "' in sheet ", groupW_sheet, ". groupW coloring/filtering disabled.")
      }
      
      if(!is.null(sample_col_meta) && !is.null(groupW_col_meta)) {
        meta_samples <- trimws(as.character(meta_df[[sample_col_meta]]))
        meta_groups  <- trimws(as.character(meta_df[[groupW_col_meta]]))
        meta_groups[meta_groups == ""] <- NA
        sample_groupW <- setNames(meta_groups, .key(meta_samples))
        sample_groupW <- sample_groupW[!is.na(names(sample_groupW)) & names(sample_groupW) != ""]
      } else {
        warning("Could not find 'sample' and 'groupW' columns in sheet ", groupW_sheet, ". groupW coloring/filtering disabled.")
      }
    }
  }
  
  all_groupW_colors <- NULL
  if(!is.null(sample_groupW)) {
    gw_all_vals <- sort(unique(na.omit(unlist(sample_groupW))))
    if(length(gw_all_vals) > 0) {
      gw_pal <- if(color_palette == "grayscale") {
        if(length(gw_all_vals) <= 2) c("#999999", "#333333") else grey(seq(0.9, 0.2, length.out = length(gw_all_vals)))
      } else if(color_palette %in% c("okabe_colors", "okabe_ito")) {
        okabe_colors <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000",
                          "#88CCEE", "#CC6677", "#117733", "#332288", "#AA4499", "#E06C75", "#5E4FA2")
        if(length(gw_all_vals) <= length(okabe_colors)) okabe_colors[1:length(gw_all_vals)] else colorRampPalette(okabe_colors)(length(gw_all_vals))
      } else if(color_palette == "viridis") {
        viridis(length(gw_all_vals), option = "D")
      } else {
        brewer.pal(min(length(gw_all_vals), 8), color_palette)
      }
      all_groupW_colors <- setNames(gw_pal[1:length(gw_all_vals)], gw_all_vals)
    }
  }
}
  
  get_color_scale <- function(sample_names) {
    n <- length(sample_names)
    
    explicit_colors <- if(!is.null(sample_colors_n)) {
      unlist(lapply(sample_names, function(s) {
        sc <- .scol(s)
        if(is.null(sc)) NULL else sc
      }))
    } else NULL
    
    generate_palette_colors <- function(k) {
      if(k <= 0) return(character(0))
      if(color_palette == "grayscale") {
        if(use_grayscale_variants) {
          grey(seq(0.95, 0.15, length.out = max(k, 2)))
        } else {
          grey(seq(0.9, 0.2, length.out = max(k, 2)))
        }
      } else if(color_palette == "okabe_colors" || color_palette == "okabe_ito") {
        if(k <= length(OKABE_ITO_COLORS)) {
          OKABE_ITO_COLORS[1:k]
        } else {
          colorRampPalette(OKABE_ITO_COLORS)(k)
        }
      } else if(color_palette == "viridis") {
        viridis(k, option = "D")
      } else {
        brewer.pal(min(k, 8), color_palette)
      }
    }
    
    n_explicit <- if(is.null(explicit_colors)) 0 else length(explicit_colors)
    n_nonexplicit <- n - n_explicit
    
    if(n_nonexplicit > 0) {
      needed <- n_nonexplicit + max(n_explicit, 2)
      cand <- generate_palette_colors(needed)
      used <- character(0)
      nonexplicit_colors <- character(n_nonexplicit)
      ci <- 1
      for(i in seq_len(n_nonexplicit)) {
        picked <- NA
        while(ci <= length(cand)) {
          c <- cand[ci]; ci <- ci + 1
          if(!(c %in% explicit_colors) && !(c %in% used)) {
            picked <- c
            break
          }
        }
        if(is.na(picked)) picked <- paste0("#", sprintf("%06X", sample(0:16777215, 1)))
        nonexplicit_colors[i] <- picked
        used <- c(used, picked)
      }
    } else {
      nonexplicit_colors <- character(0)
    }
    
    base_colors <- rep(NA, n)
    names(base_colors) <- sample_names
    ni <- 1
    for(s in sample_names) {
      sc <- .scol(s)
      if(!is.null(sc)) {
        base_colors[s] <- sc
      } else {
        base_colors[s] <- nonexplicit_colors[ni]
        ni <- ni + 1
      }
    }
    
    if(color_by_concordance && !is.null(sample_concordance)) {
      for(s in sample_names) {
        if(s %in% names(sample_concordance)) {
          conc_status <- sample_concordance[s]
          if(conc_status == "Concordant") {
            base_colors[s] <- "#009E73"
          } else if(conc_status == "Discordant") {
            base_colors[s] <- "#D55E00"
          } else {
            base_colors[s] <- "#999999"
          }
        }
      }
    }
    
    if(color_by_groupW && !is.null(sample_groupW)) {
      groupW_vals <- unique(na.omit(unlist(sample_groupW)))
      groupW_vals <- sort(groupW_vals)
      n_groups <- length(groupW_vals)
      if(n_groups > 0) {
        group_palette <- if(color_palette == "grayscale") {
          if(n_groups <= 2) c("#999999", "#333333") else grey(seq(0.9, 0.2, length.out = n_groups))
        } else if(color_palette == "okabe_colors" || color_palette == "okabe_ito") {
          okabe_colors <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000",
                            "#88CCEE", "#CC6677", "#117733", "#332288", "#AA4499", "#E06C75", "#5E4FA2")
          if(n_groups <= length(okabe_colors)) okabe_colors[1:n_groups] else colorRampPalette(okabe_colors)(n_groups)
        } else if(color_palette == "viridis") {
          viridis(n_groups, option = "D")
        } else {
          brewer.pal(min(n_groups, 8), color_palette)
        }
        group_color_map <- setNames(group_palette, groupW_vals)
        for(s in sample_names) {
          if(.key(s) %in% names(sample_groupW)) {
            gw <- sample_groupW[.key(s)]
            if(!is.na(gw) && gw %in% names(group_color_map)) {
              base_colors[s] <- group_color_map[[gw]]
            }
          }
        }
      }
    }
    
    return(base_colors)
  }
  
  build_color_mapping <- function(samples, sample_groupW, sample_concordance,
                                  color_by_groupW, color_by_concordance,
                                  color_palette, use_grayscale_variants,
                                  exclude_groups = NULL,
                                  base_colors = NULL,
                                  sample_colors = NULL,
                                  groupW_name = "Group") {
    color_aes <- "sample"
    color_values <- if(is.null(base_colors)) get_color_scale(samples) else base_colors
    color_labels <- setNames(samples, samples)
    color_legend_name <- "sample \nEU/mL"
    
    if(!color_by_groupW && !color_by_concordance) {
      return(list(
        color_aes = color_aes,
        color_values = color_values,
        color_labels = color_labels,
        color_legend_name = color_legend_name
      ))
    }
    
    gw_vals <- if(color_by_groupW && !is.null(sample_groupW)) {
      sapply(samples, function(s) {
        if(.key(s) %in% names(sample_groupW)) {
          val <- sample_groupW[.key(s)]
          if(!is.null(exclude_groups) && val %in% exclude_groups) NA else val
        } else NA
      })
    } else rep(NA_character_, length(samples))
    
    conc_vals <- if(color_by_concordance && !is.null(sample_concordance)) {
      sapply(samples, function(s) {
        if(s %in% names(sample_concordance)) sample_concordance[s] else NA
      })
    } else rep(NA_character_, length(samples))
    
    if(color_by_groupW && color_by_concordance) {
      color_aes <- "combined_label"
      combined <- sapply(seq_along(samples), function(i) {
        parts <- na.omit(c(conc_vals[i], gw_vals[i]))
        if(length(parts) == 0) NA else paste(parts, collapse=", ")
      })
      names(combined) <- samples
      
      present_labels <- unique(na.omit(combined))
      
      all_meta_samples <- unique(c(names(sample_groupW), names(sample_concordance)))
      all_combined <- sapply(all_meta_samples, function(s) {
        parts <- na.omit(c(
          if(!is.null(sample_concordance) && s %in% names(sample_concordance)) sample_concordance[s] else NULL,
          if(!is.null(sample_groupW) && s %in% names(sample_groupW)) sample_groupW[s] else NULL
        ))
        if(length(parts) == 0) NA else paste(parts, collapse=", ")
      })
      all_labels <- sort(unique(na.omit(all_combined)))
      if(!is.null(exclude_groups) && length(all_labels) > 0) {
        all_labels <- all_labels[!sapply(all_labels, function(lbl) {
          any(exclude_groups %in% strsplit(lbl, ", ")[[1]])
        })]
      }
      if(length(all_labels) == 0) all_labels <- present_labels
      
      unique_labels <- intersect(all_labels, present_labels)
      if(length(unique_labels) == 0) {
        return(list(
          color_aes = "sample", color_values = color_values,
          color_labels = setNames(samples, samples), color_legend_name = "sample \nEU/mL"
        ))
      }
      n_labels <- length(unique_labels)
      
      pal <- generate_palette(n_labels, color_palette, use_grayscale_variants)
      
      color_values <- setNames(pal[1:n_labels], unique_labels)
      color_labels <- setNames(unique_labels, unique_labels)
      color_legend_name <- "Concordance, Pattern"
      
    } else if(color_by_concordance) {
      color_aes <- "Concordance"
      unique_labels <- unique(na.omit(conc_vals))
      keep_labels <- intersect(unique_labels, c("Concordant", "Discordant"))
      if(length(keep_labels) == 0) {
        return(list(
          color_aes = "sample", color_values = color_values,
          color_labels = setNames(samples, samples), color_legend_name = "sample \nEU/mL"
        ))
      }
      
      pal <- c(Concordant = "#009E73", Discordant = "#D55E00")
      pal <- pal[keep_labels]
      
      color_values <- pal
      color_labels <- setNames(names(pal), names(pal))
      color_legend_name <- "Concordance"
      
    } else if(color_by_groupW) {
      color_aes <- "groupW"
      present_labels <- unique(na.omit(gw_vals))
      if(!is.null(exclude_groups)) present_labels <- setdiff(present_labels, exclude_groups)
      if(length(present_labels) == 0) {
        return(list(
          color_aes = "sample", color_values = color_values,
          color_labels = setNames(samples, samples), color_legend_name = "sample \nEU/mL"
        ))
      }
      if(!is.null(all_groupW_colors) && length(all_groupW_colors) > 0) {
        unique_labels <- intersect(names(all_groupW_colors), present_labels)
        color_values <- all_groupW_colors[unique_labels]
      } else {
        unique_labels <- present_labels
        pal <- generate_palette(length(unique_labels), color_palette, use_grayscale_variants)
        color_values <- setNames(pal[1:length(unique_labels)], unique_labels)
      }
      color_labels <- setNames(unique_labels, unique_labels)
      color_legend_name <- groupW_name
    }
    
    if(!is.null(sample_colors_n) && length(sample_colors_n) > 0) {
      for(s in samples) {
        sc <- .scol(s)
        if(is.null(sc)) next
        lbl <- NULL
        if(color_aes == "groupW" && .key(s) %in% names(sample_groupW)) {
          lbl <- sample_groupW[.key(s)]
        } else if(color_aes == "Concordance" && .norm(s) %in% names(sample_concordance)) {
          lbl <- sample_concordance[.norm(s)]
        } else if(color_aes == "combined_label") {
          parts <- na.omit(c(if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)],
                             if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)]))
          if(length(parts) > 0) lbl <- paste(parts, collapse=", ")
        }
        if(!is.null(lbl) && lbl %in% names(color_values)) {
          color_values[lbl] <- sc
        }
      }
    }
    
    return(list(
      color_aes = color_aes,
      color_values = color_values,
      color_labels = color_labels,
      color_legend_name = color_legend_name
    ))
  }
  
  generate_palette <- function(n, color_palette, use_grayscale_variants) {
    if(color_palette == "grayscale") {
      if(use_grayscale_variants) grey(seq(0.95, 0.15, length.out = max(n, 2)))
      else grey(seq(0.9, 0.2, length.out = max(n, 2)))
    } else if(color_palette %in% c("okabe_colors", "okabe_ito")) {
      if(n <= length(OKABE_ITO_COLORS)) OKABE_ITO_COLORS[1:n]
      else colorRampPalette(OKABE_ITO_COLORS)(n)
    } else if(color_palette == "viridis") {
      viridis(max(n, 2), option = "D")
    } else {
      brewer.pal(min(n, 8), color_palette)
    }
  }
  
  get_linetype_scale <- function(sample_names) {
    n <- length(sample_names)
    
    if(!is.null(sample_linetypes_n)) {
      default_linetypes <- rep("solid", n)
      names(default_linetypes) <- sample_names
      for(s in sample_names) {
        sl <- .slty(s)
        if(!is.null(sl)) {
          default_linetypes[s] <- sl
        }
      }
    } else if(use_various_linetypes) {
      available_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
      if(n <= length(available_linetypes)) {
        default_linetypes <- available_linetypes[1:n]
      } else {
        default_linetypes <- rep(available_linetypes, length.out = n)
      }
      names(default_linetypes) <- sample_names
    } else {
      default_linetypes <- rep("solid", n)
      names(default_linetypes) <- sample_names
    }
    
    if(highlight_high_titer && !is.null(sample_titer_3rep)) {
      for(s in sample_names) {
        if(s %in% names(sample_titer_3rep)) {
          titer_val <- sample_titer_3rep[s]
          if(!is.na(titer_val) && titer_val > titer_threshold) {
            default_linetypes[s] <- "dotted"
          }
        }
      }
    }
    
    return(default_linetypes)
  }
  
  
  simple_theme_neut <- function() {
    theme_minimal(base_size = neut_base_font) + theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title = element_text(size = neut_base_font + 2, face = "bold"),
      axis.text = element_text(size = neut_base_font, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = neut_base_font + 3, face = "bold", hjust = 0.5),
      legend.position = legend_pos,
      legend.title = element_text(size = neut_base_font, face = "bold"),
      legend.text = element_text(size = neut_base_font, hjust = legend_text_align),
      legend.text.align = 0,
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(legend_key_width, "cm"),
      legend.key.height = unit(legend_key_height, "cm"),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0.1, 0.3, 0.1, 0.3), "cm")
    )
  }
  
  simple_theme_inhib <- function() {
    theme_minimal(base_size = inhib_base_font) + theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title = element_text(size = inhib_base_font + 2, face = "bold"),
      axis.text = element_text(size = inhib_base_font, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = inhib_base_font + 3, face = "bold", hjust = 0.5),
      legend.position = legend_pos,
      legend.title = element_text(size = inhib_base_font, face = "bold"),
      legend.text = element_text(size = inhib_base_font, hjust = legend_text_align),
      legend.text.align = 0,
      legend.key.size = unit(0.7, "cm"),
      legend.key.width = unit(legend_key_width, "cm"),
      legend.key.height = unit(legend_key_height, "cm"),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0.1, 0.3, 0.1, 0.3), "cm")
    )
  }
  
  
  calculate_error_bars <- function(data) {
    data$ymin <- NA
    data$ymax <- NA
    data$sd_val <- NA
    data$n_replicates <- 0
    
    for(i in 1:nrow(data)) {
      rlu_vals <- c()
      if(!is.na(data$RLU_dup1[i]) && data$RLU_dup1[i] != "") rlu_vals <- c(rlu_vals, as.numeric(data$RLU_dup1[i]))
      if(!is.na(data$RLU_dup2[i]) && data$RLU_dup2[i] != "") rlu_vals <- c(rlu_vals, as.numeric(data$RLU_dup2[i]))
      if(!is.na(data$RLU_dup3[i]) && data$RLU_dup3[i] != "") rlu_vals <- c(rlu_vals, as.numeric(data$RLU_dup3[i]))
      
      data$n_replicates[i] <- length(rlu_vals)
      
      if(length(rlu_vals) >= 2) {
        data$ymin[i] <- min(rlu_vals, na.rm = TRUE)
        data$ymax[i] <- max(rlu_vals, na.rm = TRUE)
        data$sd_val[i] <- sd(rlu_vals, na.rm = TRUE)
      }
    }
    return(data)
  }
  
  
  
  .make_legend_plot <- function(df, color_map, linetype_map, title,
                                font = neut_base_font, key_w = legend_key_width,
                                key_h = legend_key_height, color_labels = NULL,
                                symbol_size_mult = 1.0, ncol = NULL,
                                fixed_size = FALSE) {
    cats <- as.character(unique(df$cat))
    n_cats <- length(cats)
    if(is.null(ncol)) {
      ncol <- legend_ncol
    }
    ncol <- max(1, min(ncol, n_cats))
    leg_data <- data.frame(
      cat = rep(cats, each = 2),
      x = rep(c(1, 2), times = n_cats),
      y = rep(seq_len(n_cats), each = 2)
    )
    base_lwd <- max(0.6, font / 8) * symbol_size_mult
    leg <- ggplot(leg_data, aes(x = x, y = y)) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = font, face = "bold"),
            legend.text = element_text(size = font),
            legend.key.width = unit(key_w, "cm"),
            legend.key.height = unit(key_h, "cm"),
            legend.key = element_blank(),
            legend.box = "vertical",
            legend.box.spacing = unit(0, "cm"),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0),
            legend.spacing.x = unit(0, "cm"),
            legend.spacing.y = unit(0, "cm"),
            plot.margin = margin(0, 0, 0, 0))
    
    if(!is.null(color_map)) {
      lt_override <- if(!is.null(linetype_map)) as.character(linetype_map[cats]) else rep("solid", n_cats)
      col_override <- as.character(color_map[cats])
      lab_override <- if(!is.null(color_labels)) color_labels[cats] else cats
      
      leg <- leg + 
        geom_line(aes(color = cat), linewidth = base_lwd, show.legend = TRUE) +
        scale_color_manual(
          values = setNames(col_override, cats),
          labels = lab_override,
          name = title,
          guide = guide_legend(
            override.aes = list(linetype = lt_override, linewidth = base_lwd),
            ncol = ncol, byrow = TRUE
          )
        )
    }
    if(!is.null(linetype_map) && is.null(color_map)) {
      leg <- leg + 
        geom_line(aes(linetype = cat), linewidth = base_lwd, color = "black", show.legend = TRUE) +
        scale_linetype_manual(values = setNames(as.character(linetype_map[cats]), cats),
                              name = title,
                              guide = guide_legend(
                                override.aes = list(linewidth = base_lwd),
                                ncol = ncol, byrow = TRUE
                              ))
    }
    attr(leg, "n_cats") <- n_cats
    attr(leg, "legend_ncol") <- ncol
    attr(leg, "fixed_size") <- fixed_size
    leg
  }
  
  build_explicit_sample_legend <- function(samples, title = "Sample (specified)",
                                           font = neut_base_font, key_w = legend_key_width,
                                           key_h = legend_key_height,
                                           symbol_size_mult = 1.0,
                                           fixed_size = FALSE) {
    explicit_samples <- samples[sapply(samples, function(s) .has_scol(s) || .has_slty(s))]
    if(length(explicit_samples) == 0) return(NULL)
    
    leg_df <- data.frame(
      cat = explicit_samples,
      col = sapply(explicit_samples, function(s) {
        c <- .scol(s); if(is.null(c)) "black" else c
      }),
      lty = sapply(explicit_samples, function(s) {
        l <- .slty(s); if(is.null(l)) "solid" else l
      }),
      stringsAsFactors = FALSE
    )
    leg_df$cat <- factor(leg_df$cat, levels = explicit_samples)
    
    .make_legend_plot(leg_df,
                      color_map = setNames(as.character(leg_df$col), explicit_samples),
                      linetype_map = setNames(as.character(leg_df$lty), explicit_samples),
                      title = title, font = font, key_w = key_w, key_h = key_h,
                      symbol_size_mult = symbol_size_mult,
                      fixed_size = FALSE)
  }
  
  build_concordance_legend <- function(font = neut_base_font, key_w = legend_key_width,
                                       key_h = legend_key_height,
                                       symbol_size_mult = 1.0) {
    if(is.null(sample_concordance)) return(NULL)
    vals <- unique(na.omit(sample_concordance))
    vals <- intersect(vals, c("Concordant", "Discordant"))
    if(length(vals) == 0) return(NULL)
    col_map <- c(Concordant = "#009E73", Discordant = "#D55E00")
    col_map <- col_map[vals]
    df <- data.frame(cat = factor(vals, levels = vals))
    .make_legend_plot(df, color_map = col_map, linetype_map = NULL, title = "Concordance",
                      font = font, key_w = key_w, key_h = key_h,
                      symbol_size_mult = symbol_size_mult)
  }
  
  build_highlight_legend <- function(font = neut_base_font, key_w = legend_key_width,
                                     key_h = legend_key_height,
                                     symbol_size_mult = 1.0) {
    if(!highlight_high_titer || is.null(sample_titer_3rep)) return(NULL)
    df <- data.frame(cat = factor(c("High titer", "Other"), levels = c("High titer", "Other")))
    .make_legend_plot(df,
                      color_map = NULL,
                      linetype_map = c("High titer" = "dotted", "Other" = "solid"),
                      title = "Highlight", font = font, key_w = key_w, key_h = key_h,
                      symbol_size_mult = symbol_size_mult)
  }
  
  build_color_legend <- function(color_values, color_labels,
                                 title = "Pattern", font = neut_base_font,
                                 key_w = legend_key_width, key_h = legend_key_height,
                                 symbol_size_mult = 1.0) {
    if(is.null(color_values) || length(color_values) == 0) return(NULL)
    cats <- names(color_values)
    if(length(cats) == 0) return(NULL)
    df <- data.frame(cat = factor(cats, levels = cats))
    lab_map <- if(!is.null(color_labels) && length(color_labels) > 0) color_labels else setNames(cats, cats)
    .make_legend_plot(df,
                      color_map = setNames(as.character(color_values), cats),
                      linetype_map = NULL,
                      title = title, font = font, key_w = key_w, key_h = key_h,
                      color_labels = lab_map,
                      symbol_size_mult = symbol_size_mult)
  }
  
  is_sigmoidal_data <- function(x, y) {
    if(length(x) < 4) return(FALSE)
    y_range <- diff(range(y, na.rm = TRUE))
    x_range <- diff(range(x, na.rm = TRUE))
    if(y_range <= 0 || x_range <= 0) return(FALSE)
    
    y_median <- median(abs(y), na.rm = TRUE)
    if(y_range / max(y_median, 1) < 0.2) return(FALSE)
    
    ord <- order(x)
    x_s <- x[ord]
    y_s <- y[ord]
    
    dy <- diff(y_s)
    n_pos <- sum(dy > 0)
    n_neg <- sum(dy < 0)
    if(min(n_pos, n_neg) / max(1, length(dy)) > 0.25) return(FALSE)
    
    y_lo <- quantile(y_s, 0.2)
    y_hi <- quantile(y_s, 0.8)
    if(!any(y_s <= y_lo) || !any(y_s >= y_hi)) return(FALSE)
    
    if(x_range < 1.5) return(FALSE)
    
    return(TRUE)
  }
  
  append_legends <- function(p, leg_list, n_samples = 1, pattern_legend = NULL) {
    leg_list <- Filter(Negate(is.null), leg_list)
    if(!is.null(pattern_legend)) leg_list <- c(list(pattern_legend), leg_list)
    leg_list <- Filter(Negate(is.null), leg_list)
    if(length(leg_list) == 0) return(p)
    
    n_keys <- max(n_samples, sum(sapply(leg_list, function(g) {
      if(inherits(g, "ggplot")) {
        n_c <- attr(g, "n_cats")
        if(!is.null(n_c)) return(n_c)
        nrow_guide <- tryCatch({
          lg <- cowplot::get_legend(g)
          sum(sapply(lg$grobs, function(x) max(1, length(x$grobs))))
        }, error = function(e) n_samples)
        max(1, nrow_guide)
      } else {
        max(1, length(g$grobs))
      }
    })))
    scale_f <- if(n_keys <= 4) 1.0 else if(n_keys <= 8) 0.85 else if(n_keys <= 14) 0.7 else 0.55
    leg_font <- max(6, round(neut_base_font * scale_f))
    leg_kw <- max(0.5, legend_key_width * scale_f)
    leg_kh <- max(0.4, legend_key_height * scale_f)
    
    rebuilt <- list()
    for(leg in leg_list) {
      if(inherits(leg, "ggplot")) {
        is_fixed <- !is.null(attr(leg, "fixed_size")) && attr(leg, "fixed_size")
        if(is_fixed) {
          rebuilt <- c(rebuilt, list(leg))
        } else {
          rebuilt <- c(rebuilt, list(leg + theme(
            legend.title = element_text(size = leg_font, face = "bold"),
            legend.text = element_text(size = leg_font),
            legend.key.width = unit(leg_kw, "cm"),
            legend.key.height = unit(leg_kh, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.spacing.x = unit(0, "cm"),
            legend.spacing.y = unit(0, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(0, 0, 0, 0))))
        }
      } else {
        rebuilt <- c(rebuilt, list(leg))
      }
    }
    leg_grobs <- lapply(rebuilt, function(l) if(inherits(l, "ggplot")) cowplot::get_legend(l) else l)
    not_null <- !sapply(leg_grobs, is.null)
    leg_grobs <- leg_grobs[not_null]
    if(length(leg_grobs) == 0) return(p)
    
    p <- p + theme(legend.position = "none")
    
    actual_heights <- sapply(leg_grobs, function(g) {
      tryCatch({
        gt <- if(inherits(g, "ggplot")) {
          ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
        } else if(inherits(g, "gtable")) {
          g
        } else {
          return(1)
        }
        h <- grid::grobHeight(gt)
        as.numeric(grid::convertHeight(h, "cm", valueOnly = TRUE))
      }, error = function(e) 1)
    })
    actual_heights <- pmax(actual_heights, 0.01)
    
    bottom <- if(length(leg_grobs) == 1) {
      leg_grobs[[1]]
    } else {
      cowplot::plot_grid(plotlist = leg_grobs, ncol = 1,
                         rel_heights = actual_heights,
                         align = "v")
    }
    
    total_leg_h <- sum(actual_heights)
    max_leg_h <- 3.0
    leg_ratio <- min(total_leg_h, max_leg_h)
    plot_ratio <- max(2, 10 - leg_ratio)
    cowplot::plot_grid(p, bottom, ncol = 1, align = "v", rel_heights = c(plot_ratio, leg_ratio))
  }
  
  
  apply_neut_y_scale <- function(p, y_values) {
    y_values <- y_values[is.finite(y_values) & !is.na(y_values)]
    if(length(y_values) == 0) {
      y_min_use <- if(!is.null(neut_y_min)) neut_y_min else 0
      y_max_use <- if(!is.null(neut_y_max)) neut_y_max else 1000
      y_max_limit <- if(!is.null(neut_y_max)) neut_y_max * 1.1 else y_max_use
      return(p + scale_y_continuous(breaks = c(y_min_use, y_max_use), labels = scales::comma, limits = c(y_min_use, y_max_limit)))
    }
    
    y_max_actual <- max(y_values, na.rm = T)
    y_min_actual <- min(y_values, na.rm = T)
    
    y_max_use <- if(!is.null(neut_y_max)) {
      neut_y_max
    } else if(fix_y_max) {
      y_max_fixed <- 10^ceiling(log10(y_max_actual))
      if(!is.finite(y_max_fixed)) y_max_fixed <- 1000
      y_max_fixed
    } else {
      y_max_actual
    }
    
    y_max_limit <- if(!is.null(neut_y_max)) {
      neut_y_max * 1.1
    } else {
      y_max_use
    }
    
    y_min_use <- if(!is.null(neut_y_min)) {
      neut_y_min
    } else {
      y_min_actual
    }
    
    if(y_axis_scale == "log10") {
      y_max <- y_max_limit
      y_min <- min(y_values[y_values > 0], na.rm = T)
      if(!is.finite(y_min) || y_min <= 0) y_min <- 1
      max_break <- if(!is.null(neut_y_max)) neut_y_max else y_max_limit
      breaks <- 10^(floor(log10(y_min)):ceiling(log10(max_break)))
      p <- p + scale_y_log10(breaks = breaks, labels = scales::comma, limits = c(y_min, y_max))
      p <- p + annotation_logticks(sides = "l")
    } else if(y_axis_scale == "log2") {
      y_max <- y_max_limit
      y_min <- min(y_values[y_values > 0], na.rm = T)
      if(!is.finite(y_min) || y_min <= 0) y_min <- 1
      max_break <- if(!is.null(neut_y_max)) neut_y_max else y_max_limit
      breaks <- 2^(floor(log2(y_min)):ceiling(log2(max_break)))
      p <- p + scale_y_continuous(trans = log2_trans(), breaks = breaks, labels = scales::comma, limits = c(y_min, y_max))
    } else {
      y_max <- y_max_limit
      y_min <- y_min_use
      if(!is.finite(y_min)) y_min <- 0
      if(!is.finite(y_max)) y_max <- 1000
      breaks <- pretty(c(y_min, y_max_use), n = 6)
      p <- p + scale_y_continuous(breaks = breaks, labels = scales::comma, limits = c(y_min, y_max))
    }
    return(p)
  }
  
  
  .reorder_data_columns <- function(data) {
    cn <- tolower(colnames(data))
    exp_idx <- which(cn == "experiment")
    type_idx <- which(cn %in% c("type", "test_type"))
    if(length(exp_idx) == 0) exp_idx <- 1
    if(length(type_idx) == 0) type_idx <- 3
    exp_idx <- exp_idx[1]; type_idx <- type_idx[1]
    if(exp_idx != 1 || type_idx != 3) {
      rest <- setdiff(seq_len(ncol(data)), c(exp_idx, type_idx))
      col2 <- if(length(rest) >= 1) rest[1] else integer(0)
      new_order <- c(exp_idx, col2, type_idx, setdiff(rest, col2))
      new_order <- new_order[!is.na(new_order) & new_order > 0]
      if(length(new_order) == ncol(data)) {
        data <- data[, new_order]
        message("Reordered columns (was experiment col ", exp_idx, ", type col ", type_idx, ")")
      }
    }
    data
  }
  
  .compute_derived_columns <- function(data) {
    if("RLU_dup3" %in% colnames(data)) {
      data$RLU_dup3 <- as.numeric(data$RLU_dup3)
    } else {
      data$RLU_dup3 <- NA
    }
    
    if("RLU_dup1" %in% colnames(data)) {
      data$RLU_dup1 <- as.numeric(as.character(data$RLU_dup1))
    }
    if("RLU_dup2" %in% colnames(data)) {
      data$RLU_dup2 <- as.numeric(as.character(data$RLU_dup2))
    }
    
    if("ELISA_EU/mL" %in% colnames(data)) {
      data$`ELISA_EU/mL` <- as.numeric(as.character(data$`ELISA_EU/mL`))
    }
    
    if("use" %in% colnames(data)) {
      keep_use <- is.na(data$use) | !grepl("^no$", data$use, ignore.case = TRUE)
      data <- data[keep_use, ]
    }
    
    if("log_dil" %in% colnames(data)) {
      data$log_dil <- as.numeric(data$log_dil)
    }
    
    if("avg" %in% colnames(data)) {
      data$avg <- as.numeric(data$avg)
    }
    
    .rlu_cols <- c("RLU_dup1", "RLU_dup2", "RLU_dup3")[c("RLU_dup1", "RLU_dup2", "RLU_dup3") %in% colnames(data)]
    if(length(.rlu_cols) >= 1) {
      data$avg <- rowMeans(data[, .rlu_cols, drop = FALSE], na.rm = TRUE)
      .na_avg <- is.na(data$avg) & !is.na(data$RLU_dup1)
      if(any(.na_avg)) {
        data$avg[.na_avg] <- data$RLU_dup1[.na_avg]
      }
    } else if("avg" %in% colnames(data)) {
      data$avg <- as.numeric(data$avg)
    }
    
    if(!"log_dil" %in% colnames(data) && "vir_dil" %in% colnames(data)) {
      data$log_dil <- NA
      for(i in 1:nrow(data)) {
        if(!is.na(data$vir_dil[i])) {
          val <- as.character(data$vir_dil[i])
          if(grepl("1:", val)) {
            num <- as.numeric(gsub("1:", "", val))
            if(!is.na(num) && num > 0) {
              data$log_dil[i] <- log10(1/num)
            }
          } else if(grepl("%", val)) {
            num <- as.numeric(gsub("%", "", val))/100
            if(!is.na(num) && num > 0) {
              data$log_dil[i] <- log10(num)
            }
          }
        }
      }
    }
    
    data
  }
  
  .filter_and_split_data <- function(data) {
    if("ana" %in% colnames(data)) {
      keep_ana <- !is.na(data$ana) & data$ana == "y"
      data <- data[keep_ana, ]
    }
    
    if(ncol(data) >= 3) {
      test_data <- data[data[[3]] == "test", ]
      opt_data <- data[data[[3]] == "optimize", ]
    } else {
      test_data <- data
      opt_data <- data.frame()
    }
    list(test_data = test_data, opt_data = opt_data)
  }
  
  prepare_data <- function(fpath, sheet) {
    data <- read_xlsx(fpath, sheet = sheet)
    data <- .reorder_data_columns(data)
    data <- as.data.frame(data)
    data <- .compute_derived_columns(data)
    data <- calculate_error_bars(data)
    filtered <- .filter_and_split_data(data)
    test_data <- filtered$test_data
    opt_data <- filtered$opt_data
    
    ctrl_pattern <- "VC|CC|control|ctrl|standard|positive|negative|\\bPC\\b|\\bNC\\b|vector|cell|virus|mock|blank"
    is_ctrl_name <- grepl(ctrl_pattern, as.character(data$sample), ignore.case = TRUE)
    controls <- data[is_ctrl_name, , drop = FALSE]
    if("use" %in% colnames(data)) {
      use_ctrl <- grepl("control|standard|ctrl|positive|negative", as.character(data$use), ignore.case = TRUE)
      controls <- data[is_ctrl_name | use_ctrl, , drop = FALSE]
    }
    if(nrow(controls) == 0) {
      controls <- data.frame()
    }
    
    samples_test <- test_data[!grepl(ctrl_pattern, as.character(test_data$sample), ignore.case = TRUE), ]
    if("use" %in% colnames(test_data)) {
      samples_test <- samples_test[!grepl("control|standard|ctrl|positive|negative", as.character(samples_test$use), ignore.case = TRUE), ]
    }
    
    if(!is.null(exclude_samples)) {
      samples_test <- samples_test[!samples_test$sample %in% exclude_samples, ]
    }
    
    
    samples_test$inhibition <- rep(NA_real_, nrow(samples_test))
    
    samples_test$.exp_norm <- .norm(samples_test[[1]])
    controls$.exp_norm <- .norm(controls[[1]])
    all_experiments <- unique(samples_test$.exp_norm)
    all_experiments <- all_experiments[!is.na(all_experiments) & all_experiments != ""]
    
    for(exp in all_experiments) {
      exp_val <- as.character(exp)
      
      exp_data <- samples_test[samples_test$.exp_norm == exp, ]
      
      exp_controls <- controls[controls$.exp_norm == exp, ]
      
      if(nrow(exp_controls) < 2 && nchar(exp_val) == 10) {
        base_exp <- substr(exp_val, 1, 8)
        exp_controls <- controls[controls$.exp_norm == base_exp, ]
      }
      
      if(nrow(exp_controls) < 2) {
        exp_controls <- controls
      }
      
      common_controls <- FALSE
      if(nrow(exp_controls) >= 2) {
        s_low <- tolower(as.character(exp_controls$sample))
        vc_idx <- grepl("vc|virus|vector|positive|\\bpc\\b", s_low) & !grepl("cc|cell|negative|\\bnc\\b", s_low)
        if(!any(vc_idx)) vc_idx <- grepl("vc", s_low)
        cc_idx <- grepl("cc|cell|negative|\\bnc\\b", s_low) & !grepl("vc|virus|vector|positive|\\bpc\\b", s_low)
        if(!any(cc_idx)) cc_idx <- grepl("cc", s_low)
        
        vc_val <- mean(exp_controls$avg[vc_idx], na.rm = TRUE)
        cc_val <- mean(exp_controls$avg[cc_idx], na.rm = TRUE)
        
        
        if(!is.na(vc_val) && !is.na(cc_val) && vc_val != cc_val) {
          common_controls <- TRUE
          inhib_vals <- (1 - (exp_data$avg - cc_val) / (vc_val - cc_val)) * 100
          inhib_vals <- pmax(pmin(inhib_vals, 100), 0)
          
          exp_idx <- which(samples_test$.exp_norm == exp)
          samples_test$inhibition[exp_idx] <- inhib_vals
          
                          "to", round(max(inhib_vals, na.rm = TRUE), 2)))

      } else {
                        "- found", nrow(exp_controls), "control rows; will use data-derived fallback"))
      }
      
      if(!common_controls) {
        samp_avg <- exp_data$avg[!is.na(exp_data$avg)]
        if(length(samp_avg) >= 2) {
          vc_val <- max(samp_avg, na.rm = TRUE)
          cc_val <- min(samp_avg, na.rm = TRUE)
          if(is.na(cc_val) || cc_val < 0) cc_val <- 0
          if(vc_val != cc_val) {
            inhib_vals <- (1 - (exp_data$avg - cc_val) / (vc_val - cc_val)) * 100
            inhib_vals <- pmax(pmin(inhib_vals, 100), 0)
            exp_idx <- which(samples_test$.exp_norm == exp)
             samples_test$inhibition[exp_idx] <- inhib_vals
             message("Inhibition calculated for ", length(exp_idx), " rows; VC~", round(vc_val, 1), " CC~", round(cc_val, 1))


      }
    }
    
    valid_inhib <- sum(!is.na(samples_test$inhibition))
    total_samples <- nrow(samples_test)
    
    return(list(test = samples_test, opt = opt_data, controls = controls))
  }
  
  
  get_ic50_values <- function(data, sample_names) {
    ic50_map <- list()
    
    if("ELISA_EU/mL" %in% colnames(data)) {
      
      for(s in sample_names) {
        sample_rows <- data[data$sample == s, ]
        
        ic50_vals <- sample_rows$`ELISA_EU/mL`
        ic50_vals <- ic50_vals[!is.na(ic50_vals)]
        
        if(length(ic50_vals) > 0) {
          ic50_val <- ic50_vals[1]
          ic50_map[[s]] <- ic50_val

      }

    
    return(ic50_map)
  }
  
  
  create_legend_labels <- function(sample_names, ic50_map) {
    legend_labels <- sapply(sample_names, function(s) {
      if(s %in% names(ic50_map) && !is.na(ic50_map[[s]])) {
        ic50_val <- ic50_map[[s]]
        if(ic50_val == round(ic50_val)) {
          ic50_formatted <- as.character(round(ic50_val))
        } else {
          ic50_formatted <- format(ic50_val, scientific = FALSE, trim = TRUE, digits = 4)
        }
        label <- paste0(s, "\n", ic50_formatted)
      } else {
        label <- s
      }
      return(label)
    })
    names(legend_labels) <- sample_names
    return(legend_labels)
  }
  
  
  plot_combined_neut <- function(data, add_title = TRUE, add_subtitle = TRUE, add_legend_strip = TRUE) {
    plots <- list()
    if(nrow(data$test) == 0) return(plots)
    
    exp_groups <- identify_experiment_groups(data$test[[1]])
    
    if(length(exp_groups) == 0) {
      return(plots)
    }
    
    
    for(group_name in names(exp_groups)) {
      group_info <- exp_groups[[group_name]]
      
      group_data <- data$test[data$test[[1]] %in% group_info$experiments, ]
      
      if(nrow(group_data) == 0) {
        next
      }
      
      samples_unique <- unique(group_data$sample)
      samples_unique <- samples_unique[!grepl("VC|CC", samples_unique, ignore.case = TRUE)]
      
      if(length(samples_unique) == 0) {
        next
      }
      
      if(!is.null(exclude_groups) && !is.null(sample_groupW)) {
        keep_samples <- sapply(samples_unique, function(s) {
          if(.key(s) %in% names(sample_groupW)) {
            !sample_groupW[.key(s)] %in% exclude_groups
          } else {
            TRUE
          }
        })
        samples_unique <- samples_unique[keep_samples]
        if(length(samples_unique) == 0) {
          next
        }
        group_data <- group_data[group_data$sample %in% samples_unique, ]
      }
        if(length(samples_unique) == 0) {
          next
        }
        group_data <- group_data[group_data$sample %in% samples_unique, ]
      }
      
                      paste(samples_unique, collapse=", ")))
      if(group_info$is_grouped) {
      }
      
      colors <- get_color_scale(samples_unique)
      linetypes <- get_linetype_scale(samples_unique)
      
      color_map <- build_color_mapping(samples_unique, sample_groupW, sample_concordance,
                                       color_by_groupW, color_by_concordance,
                                       color_palette, use_grayscale_variants,
                                       exclude_groups,
                                       base_colors = colors,
                                       sample_colors = sample_colors,
                                       groupW_name = "Group")
      
      show_sample_legend <- !color_map$color_aes %in% c("Concordance", "groupW", "combined_label")
      
      group_data$color_aes_col <- NA
      if(color_map$color_aes == "Concordance" && !is.null(sample_concordance)) {
        conc_vals <- sample_concordance[match(group_data$sample, names(sample_concordance))]
        group_data$Concordance <- ifelse(is.na(conc_vals), NA, conc_vals)
        group_data$color_aes_col <- group_data$Concordance
      } else if(color_map$color_aes == "groupW" && !is.null(sample_groupW)) {
        gw_vals <- sample_groupW[match(.key(group_data$sample), names(sample_groupW))]
        group_data$groupW <- ifelse(is.na(gw_vals), NA, gw_vals)
        group_data$color_aes_col <- group_data$groupW
      } else if(color_map$color_aes == "combined_label") {
        conc <- if(!is.null(sample_concordance)) sample_concordance[match(group_data$sample, names(sample_concordance))] else NA
        gw <- if(!is.null(sample_groupW)) sample_groupW[match(.key(group_data$sample), names(sample_groupW))] else NA
        group_data$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                            paste(conc, gw, sep=", "),
                                            ifelse(!is.na(conc), conc,
                                                   ifelse(!is.na(gw), gw, NA)))
        group_data$color_aes_col <- group_data$combined_label
      }
      
      ic50_map <- get_ic50_values(group_data, samples_unique)
      
      legend_labels <- create_legend_labels(samples_unique, ic50_map)
      
      if(!color_by_groupW && !color_by_concordance) {
        color_map$color_labels <- legend_labels
      }
      
      for(s in names(legend_labels)) {
      }
      
      if(group_info$is_grouped) {
        title_text <- paste("Neutralization -", group_name, "(with base", group_info$base, ")")
      } else {
        title_text <- paste("Neutralization -", group_name)
      }
      
      p <- ggplot() + 
        labs(x = "log[fold dilution]", y = "RLU") + 
        simple_theme_neut()
      
      color_aes_name <- color_map$color_aes
      
      fit_info <- list()
      has_fitted_curves <- FALSE
      
      for(s in samples_unique) {
        samp_data <- group_data[group_data$sample == s, ]
        samp_data <- samp_data[!is.na(samp_data$log_dil) & !is.na(samp_data$avg), ]
        samp_data$log_dil <- as.numeric(samp_data$log_dil)
        samp_data <- samp_data[is.finite(samp_data$log_dil), ]
        samp_data <- samp_data %>% arrange(log_dil)
        
        use_explicit <- .has_scol(s) || !is.null(.slty(s))
        
        if(nrow(samp_data) < 2) {
          if(nrow(samp_data) > 0) {
            if(use_explicit) {
              p <- p + geom_point(data=samp_data, aes(x=log_dil, y=avg), shape=15,
                                  color=if(.has_scol(s)) .scol(s) else "black", size=1.8)
            } else {
              p <- p + geom_point(data=samp_data, aes(x=log_dil, y=avg, color=.data[[color_aes_name]]), shape=15, size=1.8)
            }
            if(show_error_bars && "ymin" %in% colnames(samp_data) && "ymax" %in% colnames(samp_data)) {
              eb_data <- samp_data[!is.na(samp_data$ymin) & !is.na(samp_data$ymax) &
                                     samp_data$ymin < samp_data$ymax, ]
              if(nrow(eb_data) > 0) {
                if(use_explicit) {
                  p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=ymin, ymax=ymax),
                                         color=if(.has_scol(s)) .scol(s) else "black",
                                         width=error_bar_width, linewidth=0.6, na.rm=TRUE)
                } else {
                  p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=ymin, ymax=ymax, color=.data[[color_aes_name]]),
                                         width=error_bar_width, linewidth=0.6, na.rm=TRUE)
                }
              }
            }
          }
          next
        }
        
        line_x <- samp_data$log_dil
        line_y <- samp_data$avg
        
        if(neut_use_sigmoid && nrow(samp_data) >= 4 && is_sigmoidal_data(samp_data$log_dil, samp_data$avg)) {
          y_min_obs <- min(samp_data$avg, na.rm = TRUE)
          y_max_obs <- max(samp_data$avg, na.rm = TRUE)
          fit <- fit_sigmoid_robust(samp_data, "log_dil", "avg",
                                    y_min = y_min_obs, y_max = y_max_obs)
          
          if(!is.null(fit)) {
            fit_type <- fit$type
            if(is.null(fit_type)) fit_type <- "unknown"
            
            x_min_fit <- min(samp_data$log_dil, na.rm=T)
            x_max_fit <- max(samp_data$log_dil, na.rm=T)
            
            pred_df <- predict_sigmoid_smooth(fit, c(x_min_fit, x_max_fit), n_points = 200,
                                              pred_min = y_min_obs, pred_max = y_max_obs)
            
            if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
              pred_df$sample <- s
              pred_df$.data_col <- NA
              if(color_aes_name == "Concordance") {
                pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                pred_df$.data_col <- pred_df$Concordance
              } else if(color_aes_name == "groupW") {
                pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                pred_df$.data_col <- pred_df$groupW
              } else if(color_aes_name == "combined_label") {
                conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                 paste(conc, gw, sep=", "),
                                                 ifelse(!is.na(conc), conc,
                                                        ifelse(!is.na(gw), gw, NA)))
                pred_df$.data_col <- pred_df$combined_label
              } else {
                pred_df$.data_col <- s
              }
              
              if(use_explicit) {
                p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                   color=if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width * 1.2)
              } else {
                line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
              }
              has_fitted_curves <- TRUE
              fit_info[[s]] <- fit_type
            } else {
              if(nrow(samp_data) >= 3) {
                tryCatch({
                  lo <- loess(avg ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
                  xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
                  ys <- predict(lo, newdata = data.frame(log_dil = xs))
                  if(use_explicit) {
                    p <- p + geom_line(data = data.frame(x = xs, y = ys),
                                       aes(x = x, y = y),
                                       color = if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width, na.rm = TRUE)
                  } else {
                    loess_df <- data.frame(x = xs, y = ys, sample = s)
                    loess_df$.data_col <- s
                    if(color_aes_name == "Concordance") {
                      loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                      loess_df$.data_col <- loess_df$Concordance
                    } else if(color_aes_name == "groupW") {
                      loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      loess_df$.data_col <- loess_df$groupW
                    } else if(color_aes_name == "combined_label") {
                      conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                      gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                        paste(conc, gw, sep=", "),
                                                        ifelse(!is.na(conc), conc,
                                                               ifelse(!is.na(gw), gw, NA)))
                      loess_df$.data_col <- loess_df$combined_label
                    }
                    line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                    p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
                  }
                  fit_info[[s]] <- "loess_fallback"
                }, error = function(e) {
                  if(use_explicit) {
                    p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                       aes(x = x, y = y),
                                       color = if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width, na.rm = TRUE)
                  } else {
                    line_aes <- aes(x=log_dil, y=avg, color=.data[[color_aes_name]], linetype=sample)
                    p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
                  }
                  fit_info[[s]] <- "direct"
                })
              } else {
                if(use_explicit) {
                  p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                     aes(x = x, y = y),
                                     color = if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width, na.rm = TRUE)
                } else {
                  line_aes <- aes(x=log_dil, y=avg, color=.data[[color_aes_name]], linetype=sample)
                  p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
                }
                fit_info[[s]] <- "direct"
              }
            }
          } else {
            poly_fit <- NULL
            if(nrow(samp_data) >= 4) {
              tryCatch({
                poly_fit <- fit_polynomial_robust(samp_data, "log_dil", "avg",
                                                  y_min = y_min_obs, y_max = y_max_obs)
              }, error = function(e) {})
            }
            
            if(!is.null(poly_fit)) {
              x_min_fit <- min(samp_data$log_dil, na.rm=T)
              x_max_fit <- max(samp_data$log_dil, na.rm=T)
              pred_df <- predict_sigmoid_smooth(poly_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                pred_min = y_min_obs, pred_max = y_max_obs)
              if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                pred_df$sample <- s
                pred_df$.data_col <- s
                if(color_aes_name == "Concordance") {
                  pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                  pred_df$.data_col <- pred_df$Concordance
                } else if(color_aes_name == "groupW") {
                  pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$.data_col <- pred_df$groupW
                } else if(color_aes_name == "combined_label") {
                  conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                  gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                   paste(conc, gw, sep=", "),
                                                   ifelse(!is.na(conc), conc,
                                                          ifelse(!is.na(gw), gw, NA)))
                  pred_df$.data_col <- pred_df$combined_label
                }
                if(use_explicit) {
                  p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width * 1.2)
                } else {
                  line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                  p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                }
                fit_info[[s]] <- poly_fit$type
              } else {
                poly_fit <- NULL
              }
            }
            
            if(is.null(poly_fit) && nrow(samp_data) >= 4) {
              nonlinear_fit <- NULL
              tryCatch({
                nonlinear_fit <- fit_nonlinear_robust(samp_data, "log_dil", "avg",
                                                      y_min = y_min_obs, y_max = y_max_obs)
              }, error = function(e) {})
              
              if(!is.null(nonlinear_fit)) {
                x_min_fit <- min(samp_data$log_dil, na.rm=T)
                x_max_fit <- max(samp_data$log_dil, na.rm=T)
                pred_df <- predict_sigmoid_smooth(nonlinear_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                  pred_min = y_min_obs, pred_max = y_max_obs)
                if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                  pred_df$sample <- s
                  pred_df$.data_col <- s
                  if(color_aes_name == "Concordance") {
                    pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                    pred_df$.data_col <- pred_df$Concordance
                  } else if(color_aes_name == "groupW") {
                    pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    pred_df$.data_col <- pred_df$groupW
                  } else if(color_aes_name == "combined_label") {
                    conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                    gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                     paste(conc, gw, sep=", "),
                                                     ifelse(!is.na(conc), conc,
                                                            ifelse(!is.na(gw), gw, NA)))
                    pred_df$.data_col <- pred_df$combined_label
                  }
                  if(use_explicit) {
                    p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                       color=if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width * 1.2)
                  } else {
                    line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                    p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                  }
                  fit_info[[s]] <- nonlinear_fit$type
                }
              }
            }
            
            if(is.null(poly_fit) && is.null(nonlinear_fit) && nrow(samp_data) >= 3) {
              tryCatch({
                lo <- loess(avg ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
                xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
                ys <- predict(lo, newdata = data.frame(log_dil = xs))
                if(use_explicit) {
                  p <- p + geom_line(data = data.frame(x = xs, y = ys),
                                     aes(x = x, y = y),
                                     color = if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width, na.rm = TRUE)
                } else {
                  loess_df <- data.frame(x = xs, y = ys, sample = s)
                  loess_df$.data_col <- s
                  if(color_aes_name == "Concordance") {
                    loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                    loess_df$.data_col <- loess_df$Concordance
                  } else if(color_aes_name == "groupW") {
                    loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    loess_df$.data_col <- loess_df$groupW
                  } else if(color_aes_name == "combined_label") {
                    conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                    gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                      paste(conc, gw, sep=", "),
                                                      ifelse(!is.na(conc), conc,
                                                             ifelse(!is.na(gw), gw, NA)))
                    loess_df$.data_col <- loess_df$combined_label
                  }
                  line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                  p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
                }
                fit_info[[s]] <- "loess_fallback"
              }, error = function(e) {
                if(use_explicit) {
                  p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                     aes(x = x, y = y),
                                     color = if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width, na.rm = TRUE)
                } else {
                  line_aes <- aes(x=log_dil, y=avg, color=.data[[color_aes_name]], linetype=sample)
                  p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
                }
                fit_info[[s]] <- "direct"
              })
            } else if(is.null(poly_fit)) {
              if(use_explicit) {
                p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                   aes(x = x, y = y),
                                   color = if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width, na.rm = TRUE)
              } else {
                line_aes <- aes(x=log_dil, y=avg, color=.data[[color_aes_name]], linetype=sample)
                p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
              }
              fit_info[[s]] <- "direct"
            }
          }
        } else {
          curve_added <- FALSE
          
          if(nrow(samp_data) >= 3 && neut_use_loess) {
            tryCatch({
              lo <- loess(avg ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
              xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
              ys <- predict(lo, newdata = data.frame(log_dil = xs))
              loess_df <- data.frame(x = xs, y = ys, sample = s)
              loess_df$.data_col <- s
              if(color_aes_name == "Concordance") {
                loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                loess_df$.data_col <- loess_df$Concordance
              } else if(color_aes_name == "groupW") {
                loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                loess_df$.data_col <- loess_df$groupW
              } else if(color_aes_name == "combined_label") {
                conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                  paste(conc, gw, sep=", "),
                                                  ifelse(!is.na(conc), conc,
                                                         ifelse(!is.na(gw), gw, NA)))
                loess_df$.data_col <- loess_df$combined_label
              }
              if(use_explicit) {
                p <- p + geom_line(data = loess_df,
                                   aes(x = x, y = y),
                                   color = if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width)
              } else {
                line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
              }
              fit_info[[s]] <- "loess"
              curve_added <- TRUE
            }, error = function(e) {})
          }
          
          if(nrow(samp_data) >= 4 && (neut_use_sigmoid || neut_use_loess)) {
            poly_fit <- NULL
            tryCatch({
              poly_fit <- fit_polynomial_robust(samp_data, "log_dil", "avg",
                                                y_min = min(samp_data$avg, na.rm=TRUE),
                                                y_max = max(samp_data$avg, na.rm=TRUE))
            }, error = function(e) {})
            
            if(!is.null(poly_fit)) {
              x_min_fit <- min(samp_data$log_dil, na.rm=T)
              x_max_fit <- max(samp_data$log_dil, na.rm=T)
              pred_df <- predict_sigmoid_smooth(poly_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                pred_min = min(samp_data$avg, na.rm=TRUE),
                                                pred_max = max(samp_data$avg, na.rm=TRUE))
              if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                pred_df$sample <- s
                pred_df$.data_col <- s
                if(color_aes_name == "Concordance") {
                  pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                  pred_df$.data_col <- pred_df$Concordance
                } else if(color_aes_name == "groupW") {
                  pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$.data_col <- pred_df$groupW
                } else if(color_aes_name == "combined_label") {
                  conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                  gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                   paste(conc, gw, sep=", "),
                                                   ifelse(!is.na(conc), conc,
                                                          ifelse(!is.na(gw), gw, NA)))
                  pred_df$.data_col <- pred_df$combined_label
                }
                if(use_explicit) {
                  p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width * 1.2)
                } else {
                  line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                  p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                }
                fit_info[[s]] <- poly_fit$type
                curve_added <- TRUE
              }
            }
          }
          
          if(!curve_added && nrow(samp_data) >= 4 && (neut_use_sigmoid || neut_use_loess)) {
            nonlinear_fit <- NULL
            tryCatch({
              nonlinear_fit <- fit_nonlinear_robust(samp_data, "log_dil", "avg",
                                                    y_min = min(samp_data$avg, na.rm=TRUE),
                                                    y_max = max(samp_data$avg, na.rm=TRUE))
            }, error = function(e) {})
            
            if(!is.null(nonlinear_fit)) {
              x_min_fit <- min(samp_data$log_dil, na.rm=T)
              x_max_fit <- max(samp_data$log_dil, na.rm=T)
              pred_df <- predict_sigmoid_smooth(nonlinear_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                pred_min = min(samp_data$avg, na.rm=TRUE),
                                                pred_max = max(samp_data$avg, na.rm=TRUE))
              if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                pred_df$sample <- s
                pred_df$.data_col <- s
                if(color_aes_name == "Concordance") {
                  pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                  pred_df$.data_col <- pred_df$Concordance
                } else if(color_aes_name == "groupW") {
                  pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$.data_col <- pred_df$groupW
                } else if(color_aes_name == "combined_label") {
                  conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                  gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                   paste(conc, gw, sep=", "),
                                                   ifelse(!is.na(conc), conc,
                                                          ifelse(!is.na(gw), gw, NA)))
                  pred_df$.data_col <- pred_df$combined_label
                }
                if(use_explicit) {
                  p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width * 1.2)
                } else {
                  line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                  p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                }
                fit_info[[s]] <- nonlinear_fit$type
                curve_added <- TRUE
              }
            }
          }
          
          if(!curve_added && nrow(samp_data) >= 4 && (neut_use_sigmoid || neut_use_loess)) {
            spline_fit <- NULL
            tryCatch({
              spline_fit <- fit_spline_robust(samp_data, "log_dil", "avg",
                                              y_min = y_min_obs, y_max = y_max_obs)
            }, error = function(e) {})
            
            if(!is.null(spline_fit)) {
              x_min_fit <- min(samp_data$log_dil, na.rm=T)
              x_max_fit <- max(samp_data$log_dil, na.rm=T)
              pred_df <- predict_sigmoid_smooth(spline_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                pred_min = y_min_obs, pred_max = y_max_obs)
              if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                pred_df$sample <- s
                pred_df$.data_col <- s
                if(color_aes_name == "Concordance") {
                  pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                  pred_df$.data_col <- pred_df$Concordance
                } else if(color_aes_name == "groupW") {
                  pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$.data_col <- pred_df$groupW
                } else if(color_aes_name == "combined_label") {
                  conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                  gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                   paste(conc, gw, sep=", "),
                                                   ifelse(!is.na(conc), conc,
                                                          ifelse(!is.na(gw), gw, NA)))
                  pred_df$.data_col <- pred_df$combined_label
                }
                if(use_explicit) {
                  p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width * 1.2)
                } else {
                  line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                  p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                }
                fit_info[[s]] <- "spline"
                curve_added <- TRUE
              }
            }
          }
          
          if(!curve_added && nrow(samp_data) >= 3 && (neut_use_sigmoid || neut_use_loess)) {
            tryCatch({
              lo <- loess(avg ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
              xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
              ys <- predict(lo, newdata = data.frame(log_dil = xs))
              if(use_explicit) {
                p <- p + geom_line(data = data.frame(x = xs, y = ys),
                                   aes(x = x, y = y),
                                   color = if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width, na.rm = TRUE)
              } else {
                loess_df <- data.frame(x = xs, y = ys, sample = s)
                loess_df$.data_col <- s
                if(color_aes_name == "Concordance") {
                  loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                  loess_df$.data_col <- loess_df$Concordance
                } else if(color_aes_name == "groupW") {
                  loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  loess_df$.data_col <- loess_df$groupW
                } else if(color_aes_name == "combined_label") {
                  conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                  gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                  loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                    paste(conc, gw, sep=", "),
                                                    ifelse(!is.na(conc), conc,
                                                           ifelse(!is.na(gw), gw, NA)))
                  loess_df$.data_col <- loess_df$combined_label
                }
                line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
              }
              fit_info[[s]] <- "loess"
              curve_added <- TRUE
            }, error = function(e) {})
          }
          if(!curve_added) {
            if(use_explicit) {
              p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                 aes(x = x, y = y),
                                 color = if(.has_scol(s)) .scol(s) else "black",
                                 linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                 linewidth = line_width, na.rm = TRUE)
            } else {
              line_aes <- aes(x=log_dil, y=avg, color=.data[[color_aes_name]], linetype=sample)
              p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
            }
            fit_info[[s]] <- "direct"
          }
        }
        
        if(use_explicit) {
          p <- p + geom_point(data=samp_data, aes(x=log_dil, y=avg), shape=15,
                              color=if(.has_scol(s)) .scol(s) else "black", size=1.8)
        } else {
          p <- p + geom_point(data=samp_data, aes(x=log_dil, y=avg, color=.data[[color_aes_name]]), shape=15, size=1.8)
        }
        
        if(show_error_bars && "ymin" %in% colnames(samp_data) && "ymax" %in% colnames(samp_data)) {
          eb_data <- samp_data[!is.na(samp_data$ymin) & !is.na(samp_data$ymax) &
                                 samp_data$ymin < samp_data$ymax, ]
          if(nrow(eb_data) > 0) {
            if(use_explicit) {
              p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=ymin, ymax=ymax),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     width=error_bar_width, linewidth=0.6, na.rm=TRUE)
            } else {
              p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=ymin, ymax=ymax, color=.data[[color_aes_name]]),
                                     width=error_bar_width, linewidth=0.6, na.rm=TRUE)
            }
          }
        }
      }
      
      all_x <- group_data$log_dil[!is.na(group_data$log_dil) & is.finite(group_data$log_dil)]
      if(length(all_x) > 0) {
        x_min <- floor(min(all_x, na.rm = T))
        x_max <- ceiling(max(all_x, na.rm = T))
      } else {
        x_min <- 0
        x_max <- 5
      }
      if(!is.finite(x_min) || is.na(x_min)) x_min <- 0
      if(!is.finite(x_max) || is.na(x_max)) x_max <- 5
      x_breaks <- seq(x_min, x_max, by = 1)
      
      has_mapped_color <- any(sapply(samples_unique, function(s) !.has_scol(s)))
      has_mapped_linetype <- any(sapply(samples_unique, function(s) is.null(.slty(s))))
      has_any_explicit <- !is.null(sample_colors_n) && length(sample_colors_n) > 0 ||
        !is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0
      
      p <- p + scale_x_continuous(breaks = x_breaks, labels = sprintf("%.0f", x_breaks))
      if(has_mapped_color) {
        p <- p + scale_color_manual(values = color_map$color_values,
                                    name = color_map$color_legend_name,
                                    labels = color_map$color_labels)
      }
      if(has_mapped_linetype) {
        p <- p + scale_linetype_manual(values = linetypes,
                                       name = "Sample",
                                       labels = legend_labels,
                                       guide = if(has_any_explicit) "none" else "legend")
      }
      subtitle_parts <- list()
      if(!is.null(sample_colors_n) && length(sample_colors_n) > 0) {
        color_samples <- intersect(names(sample_colors_n), .norm(samples_unique))
        if(length(color_samples) > 0) {
          subtitle_parts[[length(subtitle_parts) + 1]] <- paste0("Fixed colors: ", paste(color_samples, collapse=", "))
        }
      }
      if(!is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0) {
        lty_samples <- intersect(names(sample_linetypes_n), .norm(samples_unique))
        if(length(lty_samples) > 0) {
          subtitle_parts[[length(subtitle_parts) + 1]] <- paste0("Fixed linetypes: ", paste(lty_samples, collapse=", "))
        }
      }
      if(neut_use_sigmoid && has_fitted_curves) {
        method_text <- paste(unique(unlist(fit_info)), collapse = ", ")
        subtitle_parts[[length(subtitle_parts) + 1]] <- paste("Fitted curves:", method_text)
      } else if(neut_use_sigmoid && !has_fitted_curves) {
        subtitle_parts[[length(subtitle_parts) + 1]] <- "Sigmoid fitting failed - showing point connections"
      }
      if(add_subtitle && length(subtitle_parts) > 0) {
        p <- p + labs(subtitle = NULL)
      }
      
      all_y <- group_data$avg[!is.na(group_data$avg) & is.finite(group_data$avg)]
      p <- apply_neut_y_scale(p, all_y)
      
      if(show_legend && length(samples_unique) > 0) {
        n_color <- length(color_map$color_labels)
        n_cols_color <- min(legend_ncol, n_color)
        n_rows_color <- ceiling((n_color + 1) / n_cols_color)
        
        color_guide <- guide_legend(
          nrow = n_rows_color, 
          ncol = n_cols_color, 
          byrow = TRUE,
          override.aes = list(linewidth = line_width * 1.5)
        )
        if(n_rows_color * n_cols_color <= n_color) {
          color_guide <- guide_legend(
            ncol = n_cols_color,
            byrow = TRUE,
            override.aes = list(linewidth = line_width * 1.5)
          )
        }
        
        guide_list <- list(color = color_guide)
        
        if(!color_by_groupW && !color_by_concordance && !has_any_explicit) {
          n_linetype <- length(legend_labels)
          n_cols_linetype <- min(legend_ncol, n_linetype)
          n_rows_linetype <- ceiling((n_linetype + 1) / n_cols_linetype)
          
          linetype_guide <- guide_legend(
            nrow = n_rows_linetype, 
            ncol = n_cols_linetype, 
            byrow = TRUE,
            override.aes = list(linewidth = line_width * 1.5)
          )
          if(n_rows_linetype * n_cols_linetype <= n_linetype) {
            linetype_guide <- guide_legend(
              ncol = n_cols_linetype,
              byrow = TRUE,
              override.aes = list(linewidth = line_width * 1.5)
            )
          }
          guide_list$linetype <- linetype_guide
        }
        
        p <- p + guides(guide_list)
        
        if(color_by_groupW || color_by_concordance) {
          p <- p + guides(linetype = "none")
        }
      }
      
      if(!show_legend) p <- p + theme(legend.position = "none")
      
      if(show_legend && length(samples_unique) > 0) {
        pat_legend <- NULL
        if(!color_by_groupW && !color_by_concordance) {
          explicit_sample_names <- c()
          if(!is.null(sample_colors_n) && length(sample_colors_n) > 0) {
            explicit_sample_names <- union(explicit_sample_names, names(sample_colors_n))
          }
          if(!is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0) {
            explicit_sample_names <- union(explicit_sample_names, names(sample_linetypes_n))
          }
          pat_color_values <- color_map$color_values[!(names(color_map$color_values) %in% explicit_sample_names)]
          pat_color_labels <- color_map$color_labels[!(names(color_map$color_labels) %in% explicit_sample_names)]
          
          if(length(pat_color_values) > 0) {
            pat_legend <- build_color_legend(
              color_values = pat_color_values,
              color_labels = pat_color_labels,
              title = color_map$color_legend_name,
              symbol_size_mult = legend_symbol_size_multiplier
            )
          }
        } else if(color_by_groupW || color_by_concordance) {
          pat_legend <- build_color_legend(
            color_values = color_map$color_values,
            color_labels = color_map$color_labels,
            title = color_map$color_legend_name,
            symbol_size_mult = legend_symbol_size_multiplier
          )
        }
        
        has_any_explicit <- !is.null(sample_colors_n) && length(sample_colors_n) > 0 ||
          !is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0
        
        if(add_legend_strip) {
          p <- append_legends(p, list(
            if(has_any_explicit) build_explicit_sample_legend(samples_unique, symbol_size_mult = legend_symbol_size_multiplier),
            if(color_by_concordance) build_concordance_legend(symbol_size_mult = legend_symbol_size_multiplier),
            if(highlight_high_titer) build_highlight_legend(symbol_size_mult = legend_symbol_size_multiplier)
          ), n_samples = length(samples_unique), pattern_legend = pat_legend)
        }
      }
      
      plot_name <- paste0("neutralization_", group_name)
      plots[[plot_name]] <- p
    }
    
    return(plots)
  }
  
  
  plot_combined_inhib <- function(data, add_title = TRUE, add_subtitle = TRUE, add_legend_strip = TRUE) {
    plots <- list()
    if(nrow(data$test) == 0 || is.null(data$test$inhibition)) return(plots)
    
    exp_groups <- identify_experiment_groups(data$test[[1]])
    
    if(length(exp_groups) == 0) {
      return(plots)
    }
    
    
    for(group_name in names(exp_groups)) {
      group_info <- exp_groups[[group_name]]
      
      group_data <- data$test[data$test[[1]] %in% group_info$experiments, ]
      
      if(nrow(group_data) == 0) {
        next
      }
      
      samples_unique <- unique(group_data$sample)
      samples_unique <- samples_unique[!grepl("VC|CC", samples_unique, ignore.case = TRUE)]
      
      if(length(samples_unique) == 0) {
        next
      }
      
      if(!is.null(exclude_groups) && !is.null(sample_groupW)) {
        keep_samples <- sapply(samples_unique, function(s) {
          if(.key(s) %in% names(sample_groupW)) {
            !sample_groupW[.key(s)] %in% exclude_groups
          } else {
            TRUE
          }
        })
        samples_unique <- samples_unique[keep_samples]
        if(length(samples_unique) == 0) {
          next
        }
        group_data <- group_data[group_data$sample %in% samples_unique, ]
      }
        if(length(samples_unique) == 0) {
          next
        }
        group_data <- group_data[group_data$sample %in% samples_unique, ]
      }
      
                      paste(samples_unique, collapse=", ")))
      if(group_info$is_grouped) {
      }
      
      valid_inhib_rows <- sum(!is.na(group_data$inhibition) & 
                                !grepl("VC|CC", group_data$sample, ignore.case = TRUE))
      
      group_data_valid <- group_data[!is.na(group_data$inhibition) & 
                                       !grepl("VC|CC", group_data$sample, ignore.case = TRUE), ]
      
      if(nrow(group_data_valid) == 0) {
        next
      }
      
      samples_with_data <- unique(group_data_valid$sample)
      
      colors <- get_color_scale(samples_with_data)
      linetypes <- get_linetype_scale(samples_with_data)
      
      color_map <- build_color_mapping(samples_with_data, sample_groupW, sample_concordance,
                                       color_by_groupW, color_by_concordance,
                                       color_palette, use_grayscale_variants,
                                       exclude_groups,
                                       base_colors = colors,
                                       sample_colors = sample_colors,
                                       groupW_name = "Group")
      
      group_data_valid$inhib_ymin <- NA
      group_data_valid$inhib_ymax <- NA
      ctrl_rows <- group_data[grepl("VC|CC|control|ctrl|standard|positive|negative|virus|cell|mock|blank",
                                    as.character(group_data$sample), ignore.case = TRUE), , drop = FALSE]
      if(nrow(ctrl_rows) >= 2) {
        s_low <- tolower(as.character(ctrl_rows$sample))
        vc_idx <- grepl("vc|virus|vector|positive|\\bpc\\b", s_low) & !grepl("cc|cell|negative|\\bnc\\b", s_low)
        if(!any(vc_idx)) vc_idx <- grepl("vc", s_low)
        cc_idx <- grepl("cc|cell|negative|\\bnc\\b", s_low) & !grepl("vc|virus|vector|positive|\\bpc\\b", s_low)
        if(!any(cc_idx)) cc_idx <- grepl("cc", s_low)
        vc_val <- mean(ctrl_rows$avg[vc_idx], na.rm = TRUE)
        cc_val <- mean(ctrl_rows$avg[cc_idx], na.rm = TRUE)
        if(!is.na(vc_val) && !is.na(cc_val) && vc_val != cc_val) {
          rlu_to_inhib <- function(rlu) pmax(pmin((1 - (rlu - cc_val) / (vc_val - cc_val)) * 100, 100), 0)
          for(i in 1:nrow(group_data_valid)) {
            rlu_vals <- c(group_data_valid$RLU_dup1[i], group_data_valid$RLU_dup2[i], group_data_valid$RLU_dup3[i])
            rlu_vals <- rlu_vals[!is.na(rlu_vals)]
            if(length(rlu_vals) >= 2) {
              group_data_valid$inhib_ymin[i] <- min(rlu_to_inhib(rlu_vals), na.rm = TRUE)
              group_data_valid$inhib_ymax[i] <- max(rlu_to_inhib(rlu_vals), na.rm = TRUE)
            }
          }
        }
      }
      show_sample_legend <- !color_map$color_aes %in% c("Concordance", "groupW", "combined_label")
      
      group_data_valid$color_aes_col <- NA
      if(color_map$color_aes == "Concordance" && !is.null(sample_concordance)) {
        conc_vals <- sample_concordance[match(group_data_valid$sample, names(sample_concordance))]
        group_data_valid$Concordance <- ifelse(is.na(conc_vals), NA, conc_vals)
        group_data_valid$color_aes_col <- group_data_valid$Concordance
      } else if(color_map$color_aes == "groupW" && !is.null(sample_groupW)) {
        gw_vals <- sample_groupW[match(.key(group_data_valid$sample), names(sample_groupW))]
        group_data_valid$groupW <- ifelse(is.na(gw_vals), NA, gw_vals)
        group_data_valid$color_aes_col <- group_data_valid$groupW
      } else if(color_map$color_aes == "combined_label") {
        conc <- if(!is.null(sample_concordance)) sample_concordance[match(group_data_valid$sample, names(sample_concordance))] else NA
        gw <- if(!is.null(sample_groupW)) sample_groupW[match(.key(group_data_valid$sample), names(sample_groupW))] else NA
        group_data_valid$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                  paste(conc, gw, sep=", "),
                                                  ifelse(!is.na(conc), conc,
                                                         ifelse(!is.na(gw), gw, NA)))
        group_data_valid$color_aes_col <- group_data_valid$combined_label
      }
      
      ic50_map <- get_ic50_values(group_data_valid, samples_with_data)
      
      legend_labels <- create_legend_labels(samples_with_data, ic50_map)
      
      for(s in names(legend_labels)) {
      }
      
      if(group_info$is_grouped) {
        title_text <- paste("Inhibition -", group_name, "(with base", group_info$base, ")")
      } else {
        title_text <- paste("Inhibition -", group_name)
      }
      
      p <- ggplot() + 
        labs(x = "log[fold dilution]", y = "Inhibition (%)") +
        simple_theme_inhib()
      
      fit_info <- list()
      has_fitted_curves <- FALSE
      
      color_aes_name <- color_map$color_aes
      
      for(s in samples_with_data) {
        
        use_explicit <- .has_scol(s) || !is.null(.slty(s))
        
        samp_data <- group_data_valid[group_data_valid$sample == s, ]
        samp_data <- samp_data[!is.na(samp_data$inhibition) & !is.na(samp_data$log_dil), ]
        samp_data$log_dil <- as.numeric(samp_data$log_dil)
        samp_data <- samp_data[is.finite(samp_data$log_dil), ]
        samp_data <- samp_data %>% arrange(log_dil)
        
        
        if(nrow(samp_data) < 2) {
          if(nrow(samp_data) > 0) {
            if(use_explicit) {
              p <- p + geom_point(data=samp_data,
                                  aes(log_dil, inhibition),
                                  color=if(.has_scol(s)) .scol(s) else "black",
                                  shape=15, size=2.5)
            } else {
              p <- p + geom_point(data=samp_data,
                                  aes(log_dil, inhibition, color = .data[[color_aes_name]]),
                                  shape=15, size=2.5)
            }
            if(show_error_bars && "inhib_ymin" %in% colnames(samp_data) && "inhib_ymax" %in% colnames(samp_data)) {
              eb_data <- samp_data[!is.na(samp_data$inhib_ymin) & !is.na(samp_data$inhib_ymax) &
                                     samp_data$inhib_ymin < samp_data$inhib_ymax, ]
              if(nrow(eb_data) > 0) {
                if(use_explicit) {
                  p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=inhib_ymin, ymax=inhib_ymax),
                                         color=if(.has_scol(s)) .scol(s) else "black",
                                         width=error_bar_width, linewidth=0.6, na.rm=TRUE)
                } else {
                  p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=inhib_ymin, ymax=inhib_ymax, color=.data[[color_aes_name]]),
                                         width=error_bar_width, linewidth=0.6, na.rm=TRUE)
                }
              }
            }
          }
          next
        }
        
        line_x <- samp_data$log_dil
        line_y <- samp_data$inhibition
        if(inhib_use_loess && !inhib_use_sigmoid && nrow(samp_data) >= 3) {
          tryCatch({
            lo <- loess(inhibition ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
            xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
            ys <- predict(lo, newdata = data.frame(log_dil = xs))
            line_x <- xs
            line_y <- ys
          }, error = function(e) {})
        }
        
        if(inhib_use_sigmoid && !inhib_use_loess && nrow(samp_data) >= 4) {
          fit <- fit_sigmoid_robust(samp_data, "log_dil", "inhibition")
          
          if(!is.null(fit)) {
            fit_type <- fit$type
            if(is.null(fit_type)) fit_type <- "unknown"
            
            x_min_fit <- min(samp_data$log_dil, na.rm=T)
            x_max_fit <- max(samp_data$log_dil, na.rm=T)
            
            pred_df <- predict_sigmoid_smooth(fit, c(x_min_fit, x_max_fit), n_points = 200)
            
            if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
              pred_df$sample <- s
              pred_df$.data_col <- s
              if(color_aes_name == "Concordance") {
                pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                pred_df$.data_col <- pred_df$Concordance
              } else if(color_aes_name == "groupW") {
                pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                pred_df$.data_col <- pred_df$groupW
              } else if(color_aes_name == "combined_label") {
                conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                 paste(conc, gw, sep=", "),
                                                 ifelse(!is.na(conc), conc,
                                                        ifelse(!is.na(gw), gw, NA)))
                pred_df$.data_col <- pred_df$combined_label
              }
              
              if(use_explicit) {
                p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                   color=if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width * 1.2)
              } else {
                line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
              }
              has_fitted_curves <- TRUE
              fit_info[[s]] <- fit_type
            } else {
              poly_fit <- NULL
              if(nrow(samp_data) >= 4) {
                tryCatch({
                  poly_fit <- fit_polynomial_robust(samp_data, "log_dil", "inhibition",
                                                    y_min = 0, y_max = 100)
                }, error = function(e) {})
              }
              
              if(!is.null(poly_fit)) {
                x_min_fit <- min(samp_data$log_dil, na.rm=T)
                x_max_fit <- max(samp_data$log_dil, na.rm=T)
                pred_df <- predict_sigmoid_smooth(poly_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                  pred_min = 0, pred_max = 100)
                if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                  pred_df$sample <- s
                  pred_df$.data_col <- s
                  if(color_aes_name == "Concordance") {
                    pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                    pred_df$.data_col <- pred_df$Concordance
                  } else if(color_aes_name == "groupW") {
                    pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    pred_df$.data_col <- pred_df$groupW
                  } else if(color_aes_name == "combined_label") {
                    conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                    gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                    pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                     paste(conc, gw, sep=", "),
                                                     ifelse(!is.na(conc), conc,
                                                            ifelse(!is.na(gw), gw, NA)))
                    pred_df$.data_col <- pred_df$combined_label
                  }
                  if(use_explicit) {
                    p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                       color=if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width * 1.2)
                  } else {
                    line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                    p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                  }
                  fit_info[[s]] <- poly_fit$type
                } else {
                  poly_fit <- NULL
                }
              }
              
              if(is.null(poly_fit)) {
                nonlinear_fit <- NULL
                tryCatch({
                  nonlinear_fit <- fit_nonlinear_robust(samp_data, "log_dil", "inhibition",
                                                        y_min = 0, y_max = 100)
                }, error = function(e) {})
                
                if(!is.null(nonlinear_fit)) {
                  x_min_fit <- min(samp_data$log_dil, na.rm=T)
                  x_max_fit <- max(samp_data$log_dil, na.rm=T)
                  pred_df <- predict_sigmoid_smooth(nonlinear_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                    pred_min = 0, pred_max = 100)
                  if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                    pred_df$sample <- s
                    pred_df$.data_col <- s
                    if(color_aes_name == "Concordance") {
                      pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                      pred_df$.data_col <- pred_df$Concordance
                    } else if(color_aes_name == "groupW") {
                      pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      pred_df$.data_col <- pred_df$groupW
                    } else if(color_aes_name == "combined_label") {
                      conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                      gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                       paste(conc, gw, sep=", "),
                                                       ifelse(!is.na(conc), conc,
                                                              ifelse(!is.na(gw), gw, NA)))
                      pred_df$.data_col <- pred_df$combined_label
                    }
                    if(use_explicit) {
                      p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                         color=if(.has_scol(s)) .scol(s) else "black",
                                         linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                         linewidth = line_width * 1.2)
                    } else {
                      line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                      p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                    }
                    fit_info[[s]] <- nonlinear_fit$type
                  } else {
                    nonlinear_fit <- NULL
                  }
                }
              }
              
              if(is.null(poly_fit) && is.null(nonlinear_fit)) {
                spline_fit <- NULL
                tryCatch({
                  spline_fit <- fit_spline_robust(samp_data, "log_dil", "inhibition",
                                                  y_min = 0, y_max = 100)
                }, error = function(e) {})
                
                if(!is.null(spline_fit)) {
                  x_min_fit <- min(samp_data$log_dil, na.rm=T)
                  x_max_fit <- max(samp_data$log_dil, na.rm=T)
                  pred_df <- predict_sigmoid_smooth(spline_fit, c(x_min_fit, x_max_fit), n_points = 200,
                                                    pred_min = 0, pred_max = 100)
                  if(!is.null(pred_df) && nrow(pred_df) > 0 && any(!is.na(pred_df$y) & is.finite(pred_df$y))) {
                    pred_df$sample <- s
                    pred_df$.data_col <- s
                    if(color_aes_name == "Concordance") {
                      pred_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                      pred_df$.data_col <- pred_df$Concordance
                    } else if(color_aes_name == "groupW") {
                      pred_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      pred_df$.data_col <- pred_df$groupW
                    } else if(color_aes_name == "combined_label") {
                      conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                      gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      pred_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                       paste(conc, gw, sep=", "),
                                                       ifelse(!is.na(conc), conc,
                                                              ifelse(!is.na(gw), gw, NA)))
                      pred_df$.data_col <- pred_df$combined_label
                    }
                    if(use_explicit) {
                      p <- p + geom_line(data=pred_df, aes(x = x, y = y),
                                         color=if(.has_scol(s)) .scol(s) else "black",
                                         linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                         linewidth = line_width * 1.2)
                    } else {
                      line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                      p <- p + geom_line(data=pred_df, line_aes, linewidth = line_width * 1.2)
                    }
                    fit_info[[s]] <- "spline"
                  }
                }
              }
              
              if(is.null(poly_fit) && is.null(nonlinear_fit) && is.null(spline_fit) && nrow(samp_data) >= 3) {
                tryCatch({
                  lo <- loess(inhibition ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
                  xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
                  ys <- predict(lo, newdata = data.frame(log_dil = xs))
                  if(use_explicit) {
                    p <- p + geom_line(data = data.frame(x = xs, y = ys),
                                       aes(x = x, y = y),
                                       color = if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width, na.rm = TRUE)
                  } else {
                    loess_df <- data.frame(x = xs, y = ys, sample = s)
                    loess_df$.data_col <- s
                    if(color_aes_name == "Concordance") {
                      loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                      loess_df$.data_col <- loess_df$Concordance
                    } else if(color_aes_name == "groupW") {
                      loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      loess_df$.data_col <- loess_df$groupW
                    } else if(color_aes_name == "combined_label") {
                      conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                      gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                      loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                        paste(conc, gw, sep=", "),
                                                        ifelse(!is.na(conc), conc,
                                                               ifelse(!is.na(gw), gw, NA)))
                      loess_df$.data_col <- loess_df$combined_label
                    }
                    line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                    p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
                  }
                  fit_info[[s]] <- "loess"
                }, error = function(e) {
                  if(use_explicit) {
                    p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                       aes(x = x, y = y),
                                       color = if(.has_scol(s)) .scol(s) else "black",
                                       linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                       linewidth = line_width, na.rm = TRUE)
                  } else {
                    line_aes <- aes(log_dil, inhibition, color = .data[[color_aes_name]], linetype = sample)
                    p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
                  }
                  fit_info[[s]] <- "direct"
                })
              } else if(is.null(poly_fit) && is.null(nonlinear_fit) && is.null(spline_fit)) {
                if(use_explicit) {
                  p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                     aes(x = x, y = y),
                                     color = if(.has_scol(s)) .scol(s) else "black",
                                     linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                     linewidth = line_width, na.rm = TRUE)
                } else {
                  line_aes <- aes(log_dil, inhibition, color = .data[[color_aes_name]], linetype = sample)
                  p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
                }
                fit_info[[s]] <- "direct"
              }
            }
          } else {
            if(use_explicit) {
              p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                 aes(x = x, y = y),
                                 color = if(.has_scol(s)) .scol(s) else "black",
                                 linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                 linewidth = line_width, na.rm = TRUE)
            } else {
              line_aes <- aes(log_dil, inhibition, color = .data[[color_aes_name]], linetype = sample)
              p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
            }
            fit_info[[s]] <- "direct"
          }
        } else {
          curve_added <- FALSE
          
          if(nrow(samp_data) >= 3 && inhib_use_loess) {
            tryCatch({
              lo <- loess(inhibition ~ log_dil, data = samp_data, span = loess_span, degree = loess_degree)
              xs <- seq(min(samp_data$log_dil), max(samp_data$log_dil), length.out = 100)
              ys <- predict(lo, newdata = data.frame(log_dil = xs))
              loess_df <- data.frame(x = xs, y = ys, sample = s)
              loess_df$.data_col <- s
              if(color_aes_name == "Concordance") {
                loess_df$Concordance <- if(s %in% names(sample_concordance)) sample_concordance[s] else NA
                loess_df$.data_col <- loess_df$Concordance
              } else if(color_aes_name == "groupW") {
                loess_df$groupW <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                loess_df$.data_col <- loess_df$groupW
              } else if(color_aes_name == "combined_label") {
                conc <- if(.norm(s) %in% names(sample_concordance)) sample_concordance[.norm(s)] else NA
                gw <- if(.key(s) %in% names(sample_groupW)) sample_groupW[.key(s)] else NA
                loess_df$combined_label <- ifelse(!is.na(conc) & !is.na(gw),
                                                  paste(conc, gw, sep=", "),
                                                  ifelse(!is.na(conc), conc,
                                                         ifelse(!is.na(gw), gw, NA)))
                loess_df$.data_col <- loess_df$combined_label
              }
              if(use_explicit) {
                p <- p + geom_line(data = loess_df,
                                   aes(x = x, y = y),
                                   color = if(.has_scol(s)) .scol(s) else "black",
                                   linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                   linewidth = line_width)
              } else {
                line_aes <- aes(x = x, y = y, color = .data$.data_col, linetype = sample)
                p <- p + geom_line(data=loess_df, line_aes, linewidth = line_width)
              }
              fit_info[[s]] <- "loess"
              curve_added <- TRUE
            }, error = function(e) {})
          }
          
          if(!curve_added) {
            if(use_explicit) {
              p <- p + geom_line(data = data.frame(x = line_x, y = line_y),
                                 aes(x = x, y = y),
                                 color = if(.has_scol(s)) .scol(s) else "black",
                                 linetype = if(!is.null(.slty(s))) .slty(s) else "solid",
                                 linewidth = line_width, na.rm = TRUE)
            } else {
              line_aes <- aes(log_dil, inhibition, color = .data[[color_aes_name]], linetype = sample)
              p <- p + geom_line(data=samp_data, line_aes, linewidth = line_width)
            }
            fit_info[[s]] <- "direct"
          }
        }
        
        if(use_explicit) {
          p <- p + geom_point(data=samp_data,
                              aes(log_dil, inhibition),
                              color=if(.has_scol(s)) .scol(s) else "black",
                              shape=15, size=2.5)
        } else {
          p <- p + geom_point(data=samp_data,
                              aes(log_dil, inhibition, color = .data[[color_aes_name]]),
                              shape=15, size=2.5)
        }
        if(show_error_bars && "inhib_ymin" %in% colnames(samp_data) && "inhib_ymax" %in% colnames(samp_data)) {
          eb_data <- samp_data[!is.na(samp_data$inhib_ymin) & !is.na(samp_data$inhib_ymax) &
                                 samp_data$inhib_ymin < samp_data$inhib_ymax, ]
          if(nrow(eb_data) > 0) {
            if(use_explicit) {
              p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=inhib_ymin, ymax=inhib_ymax),
                                     color=if(.has_scol(s)) .scol(s) else "black",
                                     width=error_bar_width, linewidth=0.6, na.rm=TRUE)
            } else {
              p <- p + geom_errorbar(data=eb_data, aes(x=log_dil, ymin=inhib_ymin, ymax=inhib_ymax, color=.data[[color_aes_name]]),
                                     width=error_bar_width, linewidth=0.6, na.rm=TRUE)
            }
          }
        }
      }
      
      all_x <- group_data_valid$log_dil[!is.na(group_data_valid$log_dil) & is.finite(group_data_valid$log_dil)]
      if(length(all_x) > 0) {
        x_min <- floor(min(all_x, na.rm = T))
        x_max <- ceiling(max(all_x, na.rm = T))
      } else {
        x_min <- 0
        x_max <- 5
      }
      if(!is.finite(x_min) || is.na(x_min)) x_min <- 0
      if(!is.finite(x_max) || is.na(x_max)) x_max <- 5
      x_breaks <- seq(x_min, x_max, by = 1)
      
      has_mapped_color <- any(sapply(samples_with_data, function(s) !.has_scol(s)))
      has_mapped_linetype <- any(sapply(samples_with_data, function(s) is.null(.slty(s))))
      has_any_explicit <- !is.null(sample_colors_n) && length(sample_colors_n) > 0 ||
        !is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0
      
      p <- p + geom_hline(yintercept=50, linetype="dashed", color="gray50", linewidth = line_width * 0.6)
      if(has_mapped_color) {
        p <- p + scale_color_manual(values = color_map$color_values,
                                    name = color_map$color_legend_name,
                                    labels = color_map$color_labels)
      }
      if(has_mapped_linetype) {
        p <- p + scale_linetype_manual(values = linetypes,
                                       name = "Sample",
                                       labels = legend_labels,
                                       guide = if(has_any_explicit) "none" else "legend")
      }
      p <- p +
        scale_x_continuous(breaks = x_breaks, labels = sprintf("%.0f", x_breaks)) +
        scale_y_continuous(breaks = seq(0, 100, by = 20),
                           limits = c(inhib_y_min, inhib_y_max))
      
      subtitle_parts <- list()
      if(!is.null(sample_colors_n) && length(sample_colors_n) > 0) {
        color_samples <- intersect(names(sample_colors_n), .norm(samples_with_data))
        if(length(color_samples) > 0) {
          subtitle_parts[[length(subtitle_parts) + 1]] <- paste0("Fixed colors: ", paste(color_samples, collapse=", "))
        }
      }
      if(!is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0) {
        lty_samples <- intersect(names(sample_linetypes_n), .norm(samples_with_data))
        if(length(lty_samples) > 0) {
          subtitle_parts[[length(subtitle_parts) + 1]] <- paste0("Fixed linetypes: ", paste(lty_samples, collapse=", "))
        }
      }
      if(inhib_use_sigmoid && has_fitted_curves) {
        method_text <- paste(unique(unlist(fit_info)), collapse = ", ")
        subtitle_parts[[length(subtitle_parts) + 1]] <- paste("Fitted curves:", method_text)
      } else if(inhib_use_sigmoid && !has_fitted_curves) {
        subtitle_parts[[length(subtitle_parts) + 1]] <- "Sigmoid fitting failed - showing point connections"
      }
      
      if(add_subtitle && length(subtitle_parts) > 0) {
        p <- p + labs(subtitle = NULL)
      }
      
      if(show_legend && length(color_map$color_labels) > 0) {
        n_color <- length(color_map$color_labels)
        n_cols_color <- min(legend_ncol, n_color)
        n_rows_color <- ceiling((n_color + 1) / n_cols_color)
        
        color_guide <- guide_legend(
          nrow = n_rows_color, 
          ncol = n_cols_color, 
          byrow = TRUE,
          override.aes = list(linewidth = line_width * 1.5)
        )
        if(n_rows_color * n_cols_color <= n_color) {
          color_guide <- guide_legend(
            ncol = n_cols_color,
            byrow = TRUE,
            override.aes = list(linewidth = line_width * 1.5)
          )
        }
        
        guide_list <- list(color = color_guide)
        
        if(!color_by_groupW && !color_by_concordance && !has_any_explicit) {
          n_linetype <- length(legend_labels)
          n_cols_linetype <- min(legend_ncol, n_linetype)
          n_rows_linetype <- ceiling((n_linetype + 1) / n_cols_linetype)
          
          linetype_guide <- guide_legend(
            nrow = n_rows_linetype, 
            ncol = n_cols_linetype, 
            byrow = TRUE,
            override.aes = list(linewidth = line_width * 1.5)
          )
          if(n_rows_linetype * n_cols_linetype <= n_linetype) {
            linetype_guide <- guide_legend(
              ncol = n_cols_linetype,
              byrow = TRUE,
              override.aes = list(linewidth = line_width * 1.5)
            )
          }
          guide_list$linetype <- linetype_guide
        }
        
        p <- p + guides(guide_list)
        
        if(color_by_groupW || color_by_concordance) {
          p <- p + guides(linetype = "none")
        }
      }
      
      if(!show_legend) p <- p + theme(legend.position = "none")
      
      if(show_legend && length(samples_with_data) > 0) {
        pat_legend <- NULL
        if(!color_by_groupW && !color_by_concordance) {
          explicit_sample_names <- c()
          if(!is.null(sample_colors_n) && length(sample_colors_n) > 0) {
            explicit_sample_names <- union(explicit_sample_names, names(sample_colors_n))
          }
          if(!is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0) {
            explicit_sample_names <- union(explicit_sample_names, names(sample_linetypes_n))
          }
          pat_color_values <- color_map$color_values[!(names(color_map$color_values) %in% explicit_sample_names)]
          pat_color_labels <- color_map$color_labels[!(names(color_map$color_labels) %in% explicit_sample_names)]
          
          if(length(pat_color_values) > 0) {
            pat_legend <- build_color_legend(
              color_values = pat_color_values,
              color_labels = pat_color_labels,
              title = color_map$color_legend_name,
              symbol_size_mult = legend_symbol_size_multiplier
            )
          }
        } else if(color_by_groupW || color_by_concordance) {
          pat_legend <- build_color_legend(
            color_values = color_map$color_values,
            color_labels = color_map$color_labels,
            title = color_map$color_legend_name,
            symbol_size_mult = legend_symbol_size_multiplier
          )
        }
        
        has_any_explicit <- !is.null(sample_colors_n) && length(sample_colors_n) > 0 ||
          !is.null(sample_linetypes_n) && length(sample_linetypes_n) > 0
        
        if(add_legend_strip) {
          p <- append_legends(p, list(
            if(has_any_explicit) build_explicit_sample_legend(samples_with_data, symbol_size_mult = legend_symbol_size_multiplier),
            if(color_by_concordance) build_concordance_legend(symbol_size_mult = legend_symbol_size_multiplier),
            if(highlight_high_titer) build_highlight_legend(symbol_size_mult = legend_symbol_size_multiplier)
          ), n_samples = length(samples_with_data), pattern_legend = pat_legend)
        }
      }
      
      plot_name <- paste0("inhibition_", group_name)
      plots[[plot_name]] <- p
    }
    
    return(plots)
  }
  
  
  build_combined_neut_inhib <- function(neut_plots, inhib_plots, out_dir) {
    if(length(neut_plots) == 0) return(list())
    
    get_shared_legend <- function(p) {
      g <- ggplotGrob(p + theme(legend.position = "bottom"))
      idx <- which(g$layout$name == "guide-box")
      if(length(idx) == 0) return(NULL)
      gtable::gtable_filter(g, "guide-box", trim = TRUE)
    }
    
    placeholder_grob <- function(txt = "No inhibition data") {
      grid::textGrob(txt, gp = grid::gpar(fontsize = 10, col = "grey40"))
    }
    
    group_names <- sub("^neutralization_", "", names(neut_plots))
    base_of <- function(g) {
      if(grepl("^[0-9]{10}$", g)) substr(g, 1, 8) else g
    }
    bases <- sort(unique(sapply(group_names, base_of)))
    
    panel_w <- max(neut_fig_width, inhib_fig_width)
    panel_h <- max(neut_fig_height, inhib_fig_height)
    legend_h <- 1.4
    
    combined_out <- list()
    
    for(base in bases) {
      subs <- sort(group_names[sapply(group_names, base_of) == base])
      if(length(subs) == 0) next
      
      panel_list <- list()
      for(sg in subs) {
        np <- neut_plots[[paste0("neutralization_", sg)]]
        ip <- inhib_plots[[paste0("inhibition_", sg)]]
        np_g <- if(!is.null(np)) np + theme(legend.position = "none") else placeholder_grob()
        ip_g <- if(!is.null(ip)) ip + theme(legend.position = "none") else placeholder_grob()
        panel_list <- c(panel_list, list(np_g, ip_g))
      }
      
      panel_grid <- arrangeGrob(grobs = panel_list, ncol = 2)
      
      ref_panel <- NULL
      for(sg in subs) {
        cand <- neut_plots[[paste0("neutralization_", sg)]]
        if(!is.null(cand)) { ref_panel <- cand; break }
        cand <- inhib_plots[[paste0("inhibition_", sg)]]
        if(!is.null(cand)) { ref_panel <- cand; break }
      }
      shared_leg <- if(!is.null(ref_panel)) get_shared_legend(ref_panel) else NULL
      
      if(!is.null(shared_leg)) {
        combined <- arrangeGrob(panel_grid, shared_leg, nrow = 2,
                                heights = unit(c(length(subs) * panel_h, legend_h), "in"))
      } else {
        combined <- panel_grid
      }
      
      n_sub <- length(subs)
      total_w <- panel_w * 2
      total_h <- n_sub * panel_h + ifelse(is.null(shared_leg), 0, legend_h)
      
      fname <- file.path(out_dir, paste0("combined_neut_inhib_", base, ".png"))
      ggsave(fname, combined, width = total_w, height = total_h, dpi = 300, units = "in")
      
      combined_out[[paste0("combined_neut_inhib_", base)]] <- combined
    }
    
    return(combined_out)
  }
  
  
  
  data <- prepare_data(file_path, sheet)
  
  all_plots <- list()
  
  if(combine_neut_inhib) {
    neut_plots <- plot_combined_neut(data, add_title = FALSE, add_subtitle = FALSE, add_legend_strip = FALSE)
    inhib_plots <- plot_combined_inhib(data, add_title = FALSE, add_subtitle = FALSE, add_legend_strip = FALSE)
    combined_plots <- build_combined_neut_inhib(neut_plots, inhib_plots, out_dir)
    all_plots <- c(all_plots, combined_plots)
  } else {
    neut_plots <- plot_combined_neut(data)
    all_plots <- c(all_plots, neut_plots)
    
    inhib_plots <- plot_combined_inhib(data)
    all_plots <- c(all_plots, inhib_plots)
    
    if(length(inhib_plots) == 0 && nrow(data$test) > 0) {
    }
  }
  
  if(length(all_plots) == 0) {

  
  for(name in names(all_plots)) {
    if(grepl("^combined_neut_inhib_", name)) next
    filename <- file.path(out_dir, paste0(name, ".png"))
    if(grepl("^neutralization_", name)) {
      ggsave(filename, all_plots[[name]], width = neut_fig_width, height = neut_fig_height, dpi = 300)
    } else if(grepl("^inhibition_", name)) {
      ggsave(filename, all_plots[[name]], width = inhib_fig_width, height = inhib_fig_height, dpi = 300)
    } else {
      ggsave(filename, all_plots[[name]], width = fig_width, height = fig_height, dpi = 300)
    }
  }
  
  analysis_statistics_sheets <- list()
  
  if(exists("combined_summary_dup2_forced4pl") && !is.null(combined_summary_dup2_forced4pl) && nrow(combined_summary_dup2_forced4pl) > 0) {
    df_summ <- combined_summary_dup2_forced4pl
    
    analysis_statistics_sheets$Sample_Results_Forced4PL <- df_summ %>%
      select(
        Experiment_Group, Sample, Is_R1, Is_R2, Is_Ref3, Is_mAb, IC50_Method,
        IC50_log_dilution, IC50_dilution, IC50_SE, IC50_95CI_Lower, IC50_95CI_Upper,
        Hill_Slope, Lower_Asymptote, Upper_Asymptote, Dynamic_Range,
        R_squared, RMSE, AIC,
        Mean_CV_percent, Median_CV_percent, Min_CV_percent, Max_CV_percent,
        R1_Name, R1_IC50_log, R1_IC50_dilution, Titer_R1, Titer_R1_SE,
        Titer_R1_95CI_Lower, Titer_R1_95CI_Upper, Titer_R1_Unit, Potency_Ratio_R1,
        R2_Name, R2_IC50_log, R2_IC50_dilution, Titer_R2, Titer_R2_SE,
        Titer_R2_95CI_Lower, Titer_R2_95CI_Upper, Titer_R2_Unit, Potency_Ratio_R2,
        Ref3_Name, Ref3_IC50_log, Ref3_IC50_dilution, Titer_Ref3, Titer_Ref3_SE,
        Titer_Ref3_95CI_Lower, Titer_Ref3_95CI_Upper, Titer_Ref3_Unit, Potency_Ratio_Ref3,
        Titer_mAb, Titer_mAb_SE, Titer_mAb_95CI_Lower, Titer_mAb_95CI_Upper, Titer_mAb_Unit, Potency_Ratio_mAb,
        Protection_Level, ELISA_EU_mL, ELISA_EU_mL_nonlin,
        Concordance_R1_ELISA, Concordance_R2_ELISA, Concordance_Ref3_ELISA, Concordance_R1_Nonlin, Concordance_R2_Nonlin, Concordance_Ref3_Nonlin,
         Mutual_Concordance_ELISA, Mutual_Concordance_Nonlin,
         ED50_Pct_Infection, RLU_Config, RM_IC50_log, RM_Method
      )
    
    analysis_statistics_sheets$Experiment_Summary_Forced4PL <- df_summ %>%
      group_by(Experiment_Group) %>%
      summarise(
        N_Samples = n(),
        N_R1 = sum(Is_R1, na.rm = TRUE),
        N_R2 = sum(Is_R2, na.rm = TRUE),
        N_Ref3 = sum(Is_Ref3, na.rm = TRUE),
        N_mAb = sum(Is_mAb, na.rm = TRUE),
        Mean_IC50_log = mean(IC50_log_dilution, na.rm = TRUE),
        SD_IC50_log = sd(IC50_log_dilution, na.rm = TRUE),
        Median_IC50_log = median(IC50_log_dilution, na.rm = TRUE),
        Min_IC50_log = min(IC50_log_dilution, na.rm = TRUE),
        Max_IC50_log = max(IC50_log_dilution, na.rm = TRUE),
        Mean_R_squared = mean(R_squared, na.rm = TRUE),
        Min_R_squared = min(R_squared, na.rm = TRUE),
        Mean_RMSE = mean(RMSE, na.rm = TRUE),
        Mean_AIC = mean(AIC, na.rm = TRUE),
        Mean_CV_percent = mean(Mean_CV_percent, na.rm = TRUE),
        Median_CV_percent = median(Median_CV_percent, na.rm = TRUE),
        Max_CV_percent = max(Max_CV_percent, na.rm = TRUE),
         ED50_Pct_Infection = mean(ED50_Pct_Infection, na.rm = TRUE),
        N_Above_Threshold = sum(Titer_R1 >= 0.5, na.rm = TRUE),
        N_Below_Threshold = sum(Titer_R1 < 0.5, na.rm = TRUE),
        Mean_Titer_R1 = mean(Titer_R1, na.rm = TRUE),
        Median_Titer_R1 = median(Titer_R1, na.rm = TRUE),
        Min_Titer_R1 = min(Titer_R1, na.rm = TRUE),
        Max_Titer_R1 = max(Titer_R1, na.rm = TRUE),
        Mean_Potency_Ratio_R1 = mean(Potency_Ratio_R1, na.rm = TRUE),
        Median_Potency_Ratio_R1 = median(Potency_Ratio_R1, na.rm = TRUE),
        Mean_Titer_mAb = mean(Titer_mAb, na.rm = TRUE),
        Median_Titer_mAb = median(Titer_mAb, na.rm = TRUE),
        N_Concordant = sum(Concordance_R1_ELISA == "Concordant", na.rm = TRUE),
        N_Discordant = sum(Concordance_R1_ELISA == "Discordant", na.rm = TRUE),
        N_Enhanced_4PL = sum(IC50_Method == "Enhanced_4PL_Good_Fit", na.rm = TRUE),
        N_Extrapolated = sum(grepl("Extrapolated", IC50_Method), na.rm = TRUE),
        .groups = "drop"
      )
    
    if(exists("method_stats_forced")) {
      analysis_statistics_sheets$Method_Statistics <- data.frame(
        IC50_Method = names(method_stats_forced),
        Count = as.numeric(unlist(method_stats_forced)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if(exists("concordance_summary") && !is.null(concordance_summary) && nrow(concordance_summary) > 0) {
    analysis_statistics_sheets$Concordance_Summary <- concordance_summary
  }
  
  if(exists("sensitivity_specificity") && !is.null(sensitivity_specificity) && nrow(sensitivity_specificity) > 0) {
    analysis_statistics_sheets$Sensitivity_Specificity <- sensitivity_specificity
  }
  
  if(exists("comparison_aicmix_forced") && !is.null(comparison_aicmix_forced) && nrow(comparison_aicmix_forced) > 0) {
    analysis_statistics_sheets$Method_Comparison_All_Five <- comparison_aicmix_forced
  }
  
  if(length(analysis_statistics_sheets) > 0) {
    tryCatch({
      write_xlsx(analysis_statistics_sheets, file.path(out_dir, "forced4pl_analysis_statistics.xlsx"))
    }, error = function(e) {
    })
  }
  
  
  return(list(plots = all_plots, data = data))
}
}

}
}


