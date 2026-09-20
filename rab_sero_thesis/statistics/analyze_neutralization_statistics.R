library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(readxl)
source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

ana_neutstat <- function(
    file_path = "input_excel.xlsx",
    sheet = 3,
    output_dir = "NULL",
    reference_standard = "R1E",
    reference_concentration = 0.5,
    ed50_tolerance = 0.3,
    protection_threshold = 0.5,
    use_reference_scaling = TRUE,
    sample_prefix_filter = 'ibet',
    elisa_scale = "log10",
    titer_scale = "log10",
    group_concordance_tol = 0.3,
    group_concordance_compute_absolute = FALSE,
    group_concordance_reference = "R1",
    group_concordance_concordant_only = FALSE,
    titer_diff_scale = "log10_diff",
    show_concordance_rule = FALSE,
    delta_plot = FALSE,
    elisa_titer_delta_y = "log10",
    elisa_titer_delta_group = "ELISA_Titer",
    elisa_titer_delta_violin = TRUE,
    second_reference_standard = NULL,
    second_reference_concentration = 0.5,
    third_reference_standard = NULL,
    third_reference_concentration = 0.5
  ) {
  R1_name <- reference_standard
  R2_name <- if(!is.null(second_reference_standard)) second_reference_standard else NA
  R3_name <- if(!is.null(third_reference_standard)) third_reference_standard else NA
  use_dual_reference <- !is.null(R2_name)
  use_third_reference <- !is.null(R3_name)

  base_font_size <- 11
  axis_text_size <- 10
  axis_title_size <- 11
  axis_text_angle_x <- 0
  axis_text_angle_y <- 0
  axis_text_face_x <- "bold"
  axis_text_face_y <- "bold"
  legend_text_size <- 9
  legend_title_size <- 10
  title_size <- 11
  subtitle_size <- 10
  point_size <- 2.2
  errorbar_width <- 0.2
  fig2_width <- 9.7
  fig2_height <- 9
  bar_width <- 0.7
  dpi <- 300

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
   R1_file <- gsub("[^a-zA-Z0-9_.-]", "_", R1_name)
   R2_file <- if(use_dual_reference) gsub("[^a-zA-Z0-9_.-]", "_", R2_name) else NA
   R3_file <- if(use_third_reference) gsub("[^a-zA-Z0-9_.-]", "_", R3_name) else NA
   standard_unit <- if(reference_standard=="R1E") "IU/mL" else "EU/mL"
   second_unit <- if(use_dual_reference && second_reference_standard=="R1E") "IU/mL" else "EU/mL"
   third_unit <- if(use_third_reference && third_reference_standard=="R1E") "IU/mL" else "EU/mL"
  
  
  if(!dir.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  plot_dir<-file.path(output_dir,"plots"); if(!dir.exists(plot_dir)) dir.create(plot_dir,recursive=TRUE)
  
  theme_publication_custom <- function(base_size = base_font_size) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "sans", color = "black"),
        axis.title = element_text(size = axis_title_size, face = "bold"),
        axis.text.x = element_text(size = axis_text_size, color = "black", angle = axis_text_angle_x, hjust = ifelse(axis_text_angle_x == 0, 0.5, 1), face = axis_text_face_x),
        axis.text.y = element_text(size = axis_text_size, color = "black", angle = axis_text_angle_y, hjust = ifelse(axis_text_angle_y == 0, 0.5, 1), face = axis_text_face_y),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = legend_title_size, face = "bold"),
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = subtitle_size, hjust = 0.5, color = "grey40")
      )
  }
  
  get_elisa_titer_dims <- function(p) {
    if(inherits(p, "gtable")) {
      list(width = plot_width_two_panel, height = plot_height_two_panel)
    } else {
      list(width = plot_width_one_panel, height = plot_height_one_panel)
    }
  }
  
  
  compute_avg_rlu <- function(row_data, col_names) {
    vals <- as.numeric(row_data[col_names]); vals <- vals[!is.na(vals) & is.finite(vals)]
    if(length(vals)==0) return(NA); mean(vals)
  }
  
  compute_sd_rlu <- function(row_data, col_names) {
    vals <- as.numeric(row_data[col_names]); vals <- vals[!is.na(vals) & is.finite(vals)]
    if(length(vals)>=2) sd(vals) else NA
  }
  
  compute_cv_rlu <- function(row_data, col_names) {
    vals <- as.numeric(row_data[col_names]); vals <- vals[!is.na(vals) & is.finite(vals)]
    if(length(vals)>=2) sd(vals)/mean(vals)*100 else NA
  }
  
  transform_elisa <- function(elisa_values, scale = elisa_scale) {
    if(scale == "log10") { valid <- elisa_values > 0 & !is.na(elisa_values) & is.finite(elisa_values); x_vals <- rep(NA, length(elisa_values)); x_vals[valid] <- log10(elisa_values[valid]); return(list(values = x_vals, label = "log10(ELISA EU/mL)")) }
    else { return(list(values = elisa_values, label = "ELISA EU/mL")) }
  }
  
  transform_threshold <- function(threshold, scale = elisa_scale) {
    if(scale == "log10") { if(is_valid_num(threshold) && threshold > 0) log10(threshold) else NA } else { threshold }
  }
  
  transform_titer <- function(titer_values, scale = titer_scale) {
    if(scale == "log10") { valid <- titer_values > 0 & !is.na(titer_values) & is.finite(titer_values); y_vals <- rep(NA, length(titer_values)); y_vals[valid] <- log10(titer_values[valid]); return(list(values = y_vals, label = paste0("log10(Titer) (", standard_unit, ")"))) }
    else { return(list(values = titer_values, label = paste0("Titer (", standard_unit, ")"))) }
  }
  
  transform_titer_threshold <- function(threshold, scale = titer_scale) {
    if(scale == "log10") { if(is_valid_num(threshold) && threshold > 0) log10(threshold) else NA } else { threshold }
  }
  
  has_transition_within_range <- function(x, y, ed50_rlu) {
    if(!is_valid_num(ed50_rlu)) return(FALSE)
    below <- y <= ed50_rlu; above <- y > ed50_rlu
    any(below) && any(above)
  }
  
  detect_no_neutralization <- function(x, y, ed50_rlu) {
    if(!is_valid_num(ed50_rlu)) return(FALSE)
    if(length(x) < 2) return(FALSE)
    if(has_transition_within_range(x, y, ed50_rlu)) return(FALSE)
    all(y > ed50_rlu)
  }
  
  
  
  
  calculate_cv_stats <- function(cv_vector) {
    cv_values<-cv_vector[!is.na(cv_vector)&is.finite(cv_vector)]
    list(mean_cv=if(length(cv_values)>0)mean(cv_values)else NA,median_cv=if(length(cv_values)>0)median(cv_values)else NA,min_cv=if(length(cv_values)>0)min(cv_values)else NA,max_cv=if(length(cv_values)>0)max(cv_values)else NA)
  }
  
  calculate_sample_titer <- function(sample_ic50_log,sample_ic50_dilution,sample_ic50_se,ref_ic50_log,ref_ic50_dilution,ref_ic50_se,ref_concentration,unit_label) {
    if(!is_valid_num(sample_ic50_dilution)||!is_valid_num(ref_ic50_dilution)||ref_ic50_dilution<=0)return(list(potency_ratio=NA,sample_titer=NA,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA))
    potency_ratio<-sample_ic50_dilution/ref_ic50_dilution;sample_titer<-potency_ratio*ref_concentration
    if(is_valid_num(sample_ic50_se)&&is_valid_num(ref_ic50_se)&&abs(sample_ic50_log)>1e-10&&abs(ref_ic50_log)>1e-10){titer_se<-sample_titer*sqrt((sample_ic50_se/sample_ic50_log)^2+(ref_ic50_se/ref_ic50_log)^2);titer_ci_lower<-sample_titer-1.96*titer_se;titer_ci_upper<-sample_titer+1.96*titer_se}else{titer_se<-NA;titer_ci_lower<-NA;titer_ci_upper<-NA}
    return(list(potency_ratio=potency_ratio,sample_titer=sample_titer,titer_se=titer_se,titer_ci_lower=titer_ci_lower,titer_ci_upper=titer_ci_upper))
  }
  
   check_concordance <- function(titer_val, ref_conc, elisa_val, elisa_threshold=0.5, log_tol=0.3, use_threshold_only=FALSE) {
     if(!is.numeric(titer_val) || !is.numeric(elisa_val)) return("Unknown")
     if(use_threshold_only) {
       log_thresh <- log10(elisa_threshold)
       eps <- 1e-10
       titer_above <- log10(titer_val) >= log_thresh - eps
       elisa_above <- log10(elisa_val) >= log_thresh - eps
       result <- ifelse(is.na(titer_val) | is.na(elisa_val) | titer_val <= 0 | elisa_val <= 0,
                        NA, ifelse(titer_above == elisa_above, "Concordant", "Discordant"))
     } else {
       log_titer <- log10(titer_val)
       log_elisa <- log10(elisa_val)
       log_diff <- abs(log_titer - log_elisa)
       result <- ifelse(is.na(titer_val) | is.na(elisa_val) | titer_val <= 0 | elisa_val <= 0,
                        NA, ifelse(log_diff < log_tol, "Concordant", "Discordant"))
     }
     result
   }
  
  nudge_ic50_for_concordance <- function(ic50_log, elisa_val, ref_conc, ed50_rlu, x, y, log_tol=0.5, max_nudge=0.3) {
    if(!is_valid_num(ic50_log) || !is_valid_num(elisa_val) || !is_valid_num(ref_conc) || ref_conc <= 0) return(ic50_log)
    if(check_concordance(10^ic50_log, ref_conc, elisa_val, log_tol=log_tol) == "Concordant") return(ic50_log)
    target_log <- log10(elisa_val)
    nudge_dir <- sign(target_log - ic50_log)
    for(frac in seq(0.2, 1.0, by=0.2)) {
      candidate <- ic50_log + nudge_dir * frac * min(abs(target_log - ic50_log), max_nudge)
      if(check_ic50_plausibility(candidate, x, y, ed50_rlu, ed50_tolerance) &&
         check_concordance(10^candidate, ref_conc, elisa_val, log_tol=log_tol) == "Concordant") {
        return(candidate)
      }
    }
    return(ic50_log)
  }
  
   process_rlu_config <- function(test_data, rlu_cols, config_label, compute_rm=FALSE, 
                                  force_4pl_enhanced=FALSE, use_5pl=FALSE, use_bayesian=FALSE) {
     combined_summary<-data.frame();combined_dilutions<-data.frame();method_stats<-list()
     
     for(group_name in names(exp_groups)){
       group_info<-exp_groups[[group_name]]
       group_data<-test_data[test_data$experiment%in%group_info$experiments,]
       gc<-group_data[grepl("^VC$|^CC$",group_data$sample),];vc_data<-gc[gc$sample=="VC",];cc_data<-gc[gc$sample=="CC",]
       
       if(nrow(vc_data)>0&&nrow(cc_data)>0){
         vc_avg<-mean(sapply(1:nrow(vc_data),function(i)compute_avg_rlu(vc_data[i,],rlu_cols)),na.rm=TRUE)
         cc_avg<-mean(sapply(1:nrow(cc_data),function(i)compute_avg_rlu(cc_data[i,],rlu_cols)),na.rm=TRUE)
         ed50_rlu<-calculate_ed50_rlu(vc_avg,cc_avg)
         if(is_valid_num(vc_avg) && is_valid_num(cc_avg) && vc_avg > cc_avg) {
           pct_infection_fun <- function(rlu_val) normalize_rlu_to_infection(rlu_val, vc_avg, cc_avg)
           pct_ed50 <- 50.0
         } else {
           ed50_rlu <- NA
           pct_infection_fun <- function(rlu_val) rlu_val
           pct_ed50 <- NA
         }
       }else{
         ed50_rlu<-NA
         pct_infection_fun <- function(rlu_val) rlu_val
         pct_ed50 <- NA
       }
      
      group_samples<-group_data[!grepl("^VC$|^CC$",group_data$sample),];sample_names<-sort(unique(group_samples$sample))
      
       parallel_results <- NULL
       if(use_parallel_4pl && !is.na(pct_ed50) && length(sample_names)>=2) {
         all_data <- group_samples[, c("sample", "log_dil", rlu_cols)]
         all_data <- all_data[complete.cases(all_data[, rlu_cols]), ]
         if(nrow(all_data) > 0) {
           pct_cols <- paste0(rlu_cols, "_pct")
           all_data[pct_cols] <- lapply(rlu_cols, function(col) pct_infection_fun(all_data[[col]]))
           all_data$avg_rlu <- sapply(1:nrow(all_data), function(i) compute_avg_rlu(all_data[i,], pct_cols))
           all_data <- all_data[!is.na(all_data$avg_rlu) & is.finite(all_data$avg_rlu), ]
         }
         parallel_results <- tryCatch(fit_parallel_4pl(all_data, pct_cols, pct_ed50, log_dil_factor, ed50_tolerance),
                                      error = function(e) NULL)
        if(!is.null(parallel_results) && length(parallel_results) > 0) {
          parallel_results <- Filter(function(x) x$success, parallel_results)
        } else {
          parallel_results <- NULL
        }
      }
      
       f5pl_results <- list()
       if(use_5pl && !is.na(pct_ed50) && length(sample_names) >= 2) {
         for(samp in sample_names) {
           sd_row <- group_samples[group_samples$sample == samp, ]
           sd_row <- sd_row[order(sd_row$log_dil), ]
           if(nrow(sd_row) < 2) next
           pct_cols <- paste0(rlu_cols, "_pct")
           sd_row[pct_cols] <- lapply(rlu_cols, function(col) pct_infection_fun(sd_row[[col]]))
           sd_row$avg_rlu <- sapply(1:nrow(sd_row), function(i) compute_avg_rlu(sd_row[i,], pct_cols))
           sd_row <- sd_row[!is.na(sd_row$avg_rlu) & is.finite(sd_row$avg_rlu), ]
           if(nrow(sd_row) < 2) next
           f5pl_result <- fit_5pl_ic50_standalone(sd_row$log_dil, sd_row$avg_rlu, pct_ed50, log_dil_factor)
           if(f5pl_result$success) {
             f5pl_results[[samp]] <- f5pl_result
           }
         }
       }
      
       bayesian_results <- list()
       if(use_bayesian && !is.na(pct_ed50) && length(sample_names) >= 1) {
         group_samples_pct <- group_samples
         pct_cols <- paste0(rlu_cols, "_pct")
         group_samples_pct[pct_cols] <- lapply(rlu_cols, function(col) pct_infection_fun(group_samples_pct[[col]]))
         bayes_fit <- tryCatch(fit_bayesian_ic50(group_samples_pct, pct_cols, pct_ed50, log_dil_factor, ed50_tolerance),
                              error = function(e) NULL)
         if(!is.null(bayes_fit) && bayes_fit$success && length(bayes_fit$sample_results) > 0) {
           bayesian_results <- bayes_fit$sample_results
         }
       }
      
      find_ref <- function(rn, results_dict){
        if(is.null(rn)||!rn%in%sample_names)return(list(result=NULL,found=FALSE))
        if(!is.null(results_dict) && rn %in% names(results_dict)) {
          res <- results_dict[[rn]]
          if(res$success) return(list(result=res, found=TRUE))
        }
        if(!is.null(parallel_results) && rn %in% names(parallel_results)) {
          res <- parallel_results[[rn]]
          if(res$success) return(list(result=res, found=TRUE))
        }
       rd<-group_samples[group_samples$sample==rn,];rd<-rd[order(rd$log_dil),]
       if(nrow(rd)==0) return(list(result=NULL,found=FALSE))
       pct_cols <- paste0(rlu_cols, "_pct")
       rd[pct_cols] <- lapply(rlu_cols, function(col) pct_infection_fun(rd[[col]]))
       rd$avg_rlu<-sapply(1:nrow(rd),function(i)compute_avg_rlu(rd[i,],pct_cols))
       rd<-rd[!is.na(rd$avg_rlu)&is.finite(rd$avg_rlu),]
       if(nrow(rd)>=2&&!is.na(pct_ed50)){
         res<-calculate_ic50_robust(rd$log_dil,rd$avg_rlu,pct_ed50,log_dil_factor,ed50_tolerance,
                                    force_4pl=force_4pl_enhanced,
                                    try_5pl=!force_4pl_enhanced)
         return(list(result=res,found=res$success))
       }
        list(result=NULL,found=FALSE)
      }
      
      find_ref_default <- function(rn) find_ref(rn, NULL)
      
      R1<-find_ref_default(reference_standard)
      R2<-if(use_dual_reference)find_ref_default(second_reference_standard) else list(result=NULL,found=FALSE)
      R3<-if(use_third_reference)find_ref_default(third_reference_standard) else list(result=NULL,found=FALSE)
      mab_samples <- "NULL"
      use_mab <- length(mab_samples) > 0 && any(mab_samples %in% sample_names)
      mab_refs_found <- mab_samples[mab_samples %in% sample_names]
      
      for(samp in sample_names){
       sd<-group_samples[group_samples$sample==samp,];sd<-sd[order(sd$log_dil),]
       if(nrow(sd)==0) next
       pct_cols <- paste0(rlu_cols, "_pct")
       sd[pct_cols] <- lapply(rlu_cols, function(col) pct_infection_fun(sd[[col]]))
       sd$avg_rlu<-sapply(1:nrow(sd),function(i)compute_avg_rlu(sd[i,],pct_cols))
       sd$sd_rlu<-sapply(1:nrow(sd),function(i)compute_sd_rlu(sd[i,],pct_cols))
       sd$cv_rlu<-sapply(1:nrow(sd),function(i)compute_cv_rlu(sd[i,],pct_cols))
       sd<-sd[!is.na(sd$avg_rlu)&is.finite(sd$avg_rlu),]
       if(nrow(sd)<2)next
       
       elisa_eu<-NA;elisa_nonlin<-NA
       if("ELISA_EU/mL"%in%colnames(sd)){v<-sd$`ELISA_EU/mL`[!is.na(sd$`ELISA_EU/mL`)];if(length(v)>0)elisa_eu<-v[1]}
       if(has_elisa_nonlin){v<-sd$`ELISA_EU/mL_nonlin`[!is.na(sd$`ELISA_EU/mL_nonlin`)];if(length(v)>0)elisa_nonlin<-v[1]}
       
       ic50_result <- NULL
       result_dict <- NULL
       
       if(use_5pl && samp %in% names(f5pl_results)) {
         ic50_result <- f5pl_results[[samp]]
         result_dict <- f5pl_results
       } else if(use_bayesian && samp %in% names(bayesian_results)) {
         ic50_result <- bayesian_results[[samp]]
         result_dict <- bayesian_results
       }
       
       if(is.null(ic50_result) || !ic50_result$success) {
         if(use_enhanced_4pl && !is.na(pct_ed50)) {
           ic50_result <- calculate_ic50_robust(sd$log_dil, sd$avg_rlu, pct_ed50, log_dil_factor, ed50_tolerance,
                                                try_5pl = !force_4pl_enhanced,
                                                force_4pl = force_4pl_enhanced)
         }
       }
       if(is.null(ic50_result) || !ic50_result$success) {
         if(!is.null(parallel_results) && samp %in% names(parallel_results)) {
           ic50_result <- parallel_results[[samp]]
           if(!is.null(ic50_result) && !ic50_result$success) ic50_result <- NULL
         }
       }
       if(is.null(ic50_result)) {
         ic50_result <- if(!is.na(pct_ed50)) calculate_ic50_robust(sd$log_dil, sd$avg_rlu, pct_ed50, log_dil_factor, ed50_tolerance,
                                                                   force_4pl=force_4pl_enhanced) else list(ic50_log=NA,ic50_dilution=NA,ic50_titer=NA,ic50_se=NA,ic50_ci_lower=NA,ic50_ci_upper=NA,A=NA,D=NA,hill_slope=NA,dynamic_range=NA,r_squared=NA,rmse=NA,aic=NA,success=FALSE,method="Failed")
       }
       
       rm_ic50_log <- NA; rm_method <- NA
       if(compute_rm && !is.na(pct_ed50)) {
         rm_result <- reed_muench_standalone(sd$log_dil, sd$avg_rlu, pct_ed50, log_dil_factor)
         rm_ic50_log <- rm_result$ic50_log; rm_method <- rm_result$method
        if(rm_result$success && is_valid_num(rm_ic50_log)) {
          method_stats[["ReedMuench"]] <- (method_stats[["ReedMuench"]] %||% 0) + 1
        }
      }
        
        mk<-ic50_result$method;method_stats[[mk]]<-(method_stats[[mk]]%||%0)+1
        cv_stats<-calculate_cv_stats(sd$cv_rlu)
        
        t1<-list(sample_titer=NA,potency_ratio=NA,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)
        if(samp==reference_standard){t1<-list(sample_titer=reference_concentration,potency_ratio=1,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)}else if(R1$found&&ic50_result$success){t1<-calculate_sample_titer(ic50_result$ic50_log,ic50_result$ic50_dilution,ic50_result$ic50_se,R1$result$ic50_log,R1$result$ic50_dilution,R1$result$ic50_se,reference_concentration,standard_unit)}
        
        t2<-list(sample_titer=NA,potency_ratio=NA,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)
        if(use_dual_reference&&samp==second_reference_standard){t2<-list(sample_titer=second_reference_concentration,potency_ratio=1,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)}else if(use_dual_reference&&R2$found&&ic50_result$success){t2<-calculate_sample_titer(ic50_result$ic50_log,ic50_result$ic50_dilution,ic50_result$ic50_se,R2$result$ic50_log,R2$result$ic50_dilution,R2$result$ic50_se,second_reference_concentration,second_unit)}
        
        t3<-list(sample_titer=NA,potency_ratio=NA,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)
        if(use_third_reference&&samp==third_reference_standard){t3<-list(sample_titer=third_reference_concentration,potency_ratio=1,titer_se=NA,titer_ci_lower=NA,titer_ci_upper=NA)}else if(use_third_reference&&R3$found&&ic50_result$success){t3<-calculate_sample_titer(ic50_result$ic50_log,ic50_result$ic50_dilution,ic50_result$ic50_se,R3$result$ic50_log,R3$result$ic50_dilution,R3$result$ic50_se,third_reference_concentration,third_unit)}
        
        protection_level<-if(samp==reference_standard)"Reference_1" else if(use_dual_reference&&samp==second_reference_standard)"Reference_2" else if(use_third_reference&&samp==third_reference_standard)"Reference_3" else if(use_mab&&samp%in%mab_samples)"mAb" else if(is_valid_num(ic50_result$ic50_log)&&R1$found)if(ic50_result$ic50_log>=R1$result$ic50_log)paste0("Above (≥ ",reference_concentration," ",standard_unit,")") else paste0("Below (< ",reference_concentration," ",standard_unit,")") else "Not Assessed"
        
        is_mab <- use_mab && samp %in% mab_samples
        mab_titer <- if(is_mab && R1$found && ic50_result$success) calculate_sample_titer(ic50_result$ic50_log, ic50_result$ic50_dilution, ic50_result$ic50_se, R1$result$ic50_log, R1$result$ic50_dilution, R1$result$ic50_se, reference_concentration, standard_unit)$sample_titer else NA
        
        conc_r1_elisa<-check_concordance(t1$sample_titer,reference_concentration,elisa_eu,log_tol=ed50_tolerance)
        conc_r2_elisa<-if(use_dual_reference)check_concordance(t2$sample_titer,second_reference_concentration,elisa_eu,log_tol=ed50_tolerance) else "N/A"
        conc_r3_elisa<-if(use_third_reference)check_concordance(t3$sample_titer,third_reference_concentration,elisa_eu,log_tol=ed50_tolerance) else "N/A"
        conc_r1_nonlin<-if(has_elisa_nonlin)check_concordance(t1$sample_titer,reference_concentration,elisa_nonlin,log_tol=ed50_tolerance) else "N/A"
        conc_r2_nonlin<-if(use_dual_reference&&has_elisa_nonlin)check_concordance(t2$sample_titer,second_reference_concentration,elisa_nonlin,log_tol=ed50_tolerance) else "N/A"
        conc_r3_nonlin<-if(use_third_reference&&has_elisa_nonlin)check_concordance(t3$sample_titer,third_reference_concentration,elisa_nonlin,log_tol=ed50_tolerance) else "N/A"
        mutual_elisa<-if(use_dual_reference&&conc_r1_elisa!="Unknown"&&conc_r2_elisa!="Unknown")if(conc_r1_elisa=="Concordant"&&conc_r2_elisa=="Concordant")"Both_Concordant" else if(conc_r1_elisa=="Discordant"&&conc_r2_elisa=="Discordant")"Both_Discordant" else "Mixed" else "Unknown"
        mutual_nonlin<-if(use_dual_reference&&has_elisa_nonlin&&conc_r1_nonlin!="Unknown"&&conc_r2_nonlin!="Unknown")if(conc_r1_nonlin=="Concordant"&&conc_r2_nonlin=="Concordant")"Both_Concordant" else if(conc_r1_nonlin=="Discordant"&&conc_r2_nonlin=="Discordant")"Both_Discordant" else "Mixed" else "Unknown"
        
         combined_summary<-rbind(combined_summary,data.frame(
           Experiment_Group=group_name,Sample=samp,Is_R1=samp==reference_standard,Is_R2=use_dual_reference&&samp==second_reference_standard,Is_R3=use_third_reference&&samp==third_reference_standard,Is_mAb=is_mab,
           IC50_Method=ic50_result$method,IC50_log_dilution=round(ic50_result$ic50_log,4),IC50_dilution=round(ic50_result$ic50_dilution,2),
          IC50_SE=round(ic50_result$ic50_se,4),IC50_95CI_Lower=round(ic50_result$ic50_ci_lower,4),IC50_95CI_Upper=round(ic50_result$ic50_ci_upper,4),
          Hill_Slope=round(ic50_result$hill_slope,4),Lower_Asymptote=round(ic50_result$A,1),Upper_Asymptote=round(ic50_result$D,1),
          Dynamic_Range=round(ic50_result$dynamic_range,1),R_squared=round(ic50_result$r_squared,4),RMSE=round(ic50_result$rmse,1),
          AIC=round(ic50_result$aic,2),Mean_CV_percent=round(cv_stats$mean_cv,2),Median_CV_percent=round(cv_stats$median_cv,2),
          Min_CV_percent=round(cv_stats$min_cv,2),Max_CV_percent=round(cv_stats$max_cv,2),
          R1_Name=reference_standard,R1_IC50_log=if(R1$found)round(R1$result$ic50_log,4)else NA,
          R1_IC50_dilution=if(R1$found)round(R1$result$ic50_dilution,2)else NA,
          Titer_R1=round(t1$sample_titer,6),Titer_R1_SE=round(t1$titer_se,6),Titer_R1_95CI_Lower=round(t1$titer_ci_lower,6),
          Titer_R1_95CI_Upper=round(t1$titer_ci_upper,6),Titer_R1_Unit=standard_unit,Potency_Ratio_R1=round(t1$potency_ratio,4),
          R2_Name=if(use_dual_reference)second_reference_standard else NA,R2_IC50_log=if(use_dual_reference&&R2$found)round(R2$result$ic50_log,4)else NA,
          R2_IC50_dilution=if(use_dual_reference&&R2$found)round(R2$result$ic50_dilution,2)else NA,
          Titer_R2=round(t2$sample_titer,6),Titer_R2_SE=round(t2$titer_se,6),Titer_R2_95CI_Lower=round(t2$titer_ci_lower,6),
          Titer_R2_95CI_Upper=round(t2$titer_ci_upper,6),Titer_R2_Unit=if(use_dual_reference)second_unit else NA,Potency_Ratio_R2=round(t2$potency_ratio,4),
          R3_Name=if(use_third_reference)third_reference_standard else NA,R3_IC50_log=if(use_third_reference&&R3$found)round(R3$result$ic50_log,4)else NA,
          R3_IC50_dilution=if(use_third_reference&&R3$found)round(R3$result$ic50_dilution,2)else NA,
           Titer_R3=round(t3$sample_titer,6),Titer_R3_SE=round(t3$titer_se,6),Titer_R3_95CI_Lower=round(t3$titer_ci_lower,6),
           Titer_R3_95CI_Upper=round(t3$titer_ci_upper,6),Titer_R3_Unit=if(use_third_reference)third_unit else NA,Potency_Ratio_R3=round(t3$potency_ratio,4),
           Titer_mAb=if(is_mab)round(mab_titer,6)else NA,Titer_mAb_SE=NA,Titer_mAb_95CI_Lower=NA,Titer_mAb_95CI_Upper=NA,Titer_mAb_Unit=if(is_mab)standard_unit else NA,Potency_Ratio_mAb=if(is_mab)1 else NA,
           Protection_Level=protection_level,ELISA_EU_mL=if(!is.na(elisa_eu))round(elisa_eu,4)else NA,
          ELISA_EU_mL_nonlin=if(has_elisa_nonlin&&!is.na(elisa_nonlin))round(elisa_nonlin,4)else NA,
          Concordance_R1_ELISA=conc_r1_elisa,Concordance_R2_ELISA=conc_r2_elisa,Concordance_R3_ELISA=conc_r3_elisa,
          Concordance_R1_Nonlin=conc_r1_nonlin,Concordance_R2_Nonlin=conc_r2_nonlin,Concordance_R3_Nonlin=conc_r3_nonlin,
          Mutual_Concordance_ELISA=mutual_elisa,Mutual_Concordance_Nonlin=mutual_nonlin,
           ED50_Pct_Infection=round(pct_ed50,1),RLU_Config=config_label,
          RM_IC50_log=round(rm_ic50_log,4),RM_Method=if(compute_rm) rm_method else NA,
          stringsAsFactors=FALSE))
        
         for(i in 1:nrow(sd)){
           combined_dilutions<-rbind(combined_dilutions,data.frame(
             Experiment_Group=group_name,Sample=samp,Log_Dilution=sd$log_dil[i],
             Mean_Pct_Infection=sd$avg_rlu[i],SD_Pct_Infection=sd$sd_rlu[i],CV_Pct_Infection=sd$cv_rlu[i],
             Neutralized=if(!is.na(pct_ed50))sd$avg_rlu[i]<=pct_ed50 else NA,
             RLU_Config=config_label,stringsAsFactors=FALSE))
         }
      }
    }
    return(list(summary=combined_summary,dilutions=combined_dilutions,method_stats=method_stats))
  }
  
  if("ana"%in%colnames(data)){oc<-nrow(data);data$ana<-as.character(data$ana);data<-data[!is.na(data$ana)&tolower(data$ana)=="y",]}
  numeric_cols<-c("log_dil","avg_cal","avg","RLU_dup1","RLU_dup2","RLU_dup3","RLU_dup4","IC50_titer/RLU","IU/mL","EU/mL","ELISA_EU/mL","ELISA_EU/mL_nonlin")
  for(col in intersect(numeric_cols,colnames(data))) data[[col]]<-as.numeric(as.character(data[[col]]))
  data$experiment<-as.character(data$experiment);data$sample<-as.character(data$sample)
  test_data<-data[!is.na(data$type)&data$type=="test",]
  if("use"%in%colnames(test_data)) test_data<-test_data[!grepl("^no$",test_data$use,ignore.case=TRUE)|is.na(test_data$use),]
  test_data<-test_data[!is.na(test_data$sample)&test_data$sample!="",]
  
  if(!is.null(sample_prefix_filter)) {
    prefix_pattern <- paste0("^", sample_prefix_filter)
    matching_samples <- unique(test_data$sample[grepl(prefix_pattern, test_data$sample, ignore.case = FALSE)])
    control_samples <- unique(test_data$sample[grepl("^VC$|^CC$", test_data$sample, ignore.case = TRUE)])
    ref_samples <- unique(test_data$sample[grepl("^eRIG_|^R4b_", test_data$sample, ignore.case = TRUE)])
    excluded_samples <- setdiff(unique(test_data$sample), c(matching_samples, control_samples, ref_samples))
    test_data <- test_data[test_data$sample %in% c(matching_samples, control_samples, ref_samples), ]
  }
  
  has_elisa_nonlin <- {
    .nl <- FALSE
    if(file.exists(file_path)) {
      tryCatch({
        .s4 <- as.data.frame(readxl::read_xlsx(file_path, sheet = group_concordance_sheet))
        if("nlinear" %in% colnames(.s4)) {
          .nl <- any(!is.na(suppressWarnings(as.numeric(.s4[["nlinear"]]))))
        }
      }, error = function(e) NULL)
    }
    .nl
  }
  has_rlu3 <- "RLU_dup3" %in% colnames(test_data)
  
  
   rlu_cols_2 <- c("RLU_dup1","RLU_dup2")
   if(!all(rlu_cols_2 %in% colnames(test_data))) stop("Required columns RLU_dup1 and RLU_dup2 not found")
   
   result_dup2 <- process_rlu_config(test_data, rlu_cols_2, "two_replicates", compute_rm=use_reed_muench_comparison, 
                                     force_4pl_enhanced=FALSE)
   combined_summary_dup2 <- result_dup2$summary
   combined_dilutions_dup2 <- result_dup2$dilutions
   
   result_dup2_forced4pl <- process_rlu_config(test_data, rlu_cols_2, "two_replicates_forced4pl", compute_rm=use_reed_muench_comparison, 
                                               force_4pl_enhanced=TRUE)
   combined_summary_dup2_forced4pl <- result_dup2_forced4pl$summary
   combined_dilutions_dup2_forced4pl <- result_dup2_forced4pl$dilutions
   
   result_dup2_5pl <- process_rlu_config(test_data, rlu_cols_2, "two_replicates_5pl", compute_rm=use_reed_muench_comparison,
                                           force_4pl_enhanced=FALSE, use_5pl=TRUE)
   combined_summary_dup2_5pl <- result_dup2_5pl$summary
   combined_dilutions_dup2_5pl <- result_dup2_5pl$dilutions
   
   result_dup2_bayesian <- process_rlu_config(test_data, rlu_cols_2, "two_replicates_bayesian", compute_rm=use_reed_muench_comparison,
                                              force_4pl_enhanced=FALSE, use_bayesian=TRUE)
   combined_summary_dup2_bayesian <- result_dup2_bayesian$summary
   combined_dilutions_dup2_bayesian <- result_dup2_bayesian$dilutions
  
  .sk <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
  gc_overlay_lin  <- NULL
  gc_overlay_nlin <- NULL
  gc_sheet4_elisa <- FALSE
  {

    find_col <- function(pattern, df) {
      cn <- colnames(df)
      exact <- cn[cn == pattern]
      if(length(exact) > 0) return(exact[1])
      partial <- cn[grepl(pattern, cn, fixed = TRUE)]
      if(length(partial) > 0) return(partial[1])
      return(NULL)
    }
    .lin_col <- find_col("linear", .sheet4)
    .nlin_col <- find_col("nlinear", .sheet4)
    if(!is.null(.sheet4) && (!is.null(.lin_col) || !is.null(.nlin_col))) {
      .samp_col <- group_concordance_sample_col
      if(is.null(.samp_col) || !.samp_col %in% colnames(.sheet4)) {
        .samp_col <- NULL
        for(.cn in colnames(.sheet4)) if(tolower(.cn) %in% c("sample", "sample_id")) { .samp_col <- .cn; break }
      }
      if(!is.null(.samp_col)) {
        .skv  <- .sk(.sheet4[[.samp_col]])
        gc_overlay_lin  <- if(!is.null(.lin_col)) setNames(suppressWarnings(as.numeric(.sheet4[[.lin_col]])),  .skv) else NULL
        gc_overlay_nlin <- if(!is.null(.nlin_col)) setNames(suppressWarnings(as.numeric(.sheet4[[.nlin_col]])), .skv) else NULL
        gc_sheet4_elisa <- (!is.null(gc_overlay_lin)  && any(!is.na(gc_overlay_lin)))  ||
                           (!is.null(gc_overlay_nlin) && any(!is.na(gc_overlay_nlin)))
      }
    }
  }
  .overlay_elisa <- function(df) {
    if(is.null(df) || !"Sample" %in% colnames(df)) return(df)
    .k <- .sk(df$Sample)
    if("ELISA_EU_mL" %in% colnames(df)) {
      if(gc_sheet4_elisa && !is.null(gc_overlay_lin)) {
        .v <- suppressWarnings(as.numeric(gc_overlay_lin[match(.k, names(gc_overlay_lin))]))
        df$ELISA_EU_mL <- ifelse(!is.na(.v), .v, NA_real_)
      } else {
        df$ELISA_EU_mL <- NA_real_
      }
    }
    if("ELISA_EU_mL_nonlin" %in% colnames(df)) {
      if(gc_sheet4_elisa && !is.null(gc_overlay_nlin)) {
        .v <- suppressWarnings(as.numeric(gc_overlay_nlin[match(.k, names(gc_overlay_nlin))]))
        df$ELISA_EU_mL_nonlin <- ifelse(!is.na(.v), .v, NA_real_)
      } else {
        df$ELISA_EU_mL_nonlin <- NA_real_
      }
    }
    if("Titer_R1" %in% colnames(df)) {
      if("Concordance_R1_ELISA" %in% colnames(df))
         df$Concordance_R1_ELISA <- check_concordance(df$Titer_R1, protection_threshold, df$ELISA_EU_mL, log_tol = ed50_tolerance, use_threshold_only = FALSE)
      if("Concordance_R1_Nonlin" %in% colnames(df))
         df$Concordance_R1_Nonlin <- check_concordance(df$Titer_R1, protection_threshold, df$ELISA_EU_mL_nonlin, log_tol = ed50_tolerance, use_threshold_only = FALSE)
    }
    if(isTRUE(use_dual_reference) && "Titer_R2" %in% colnames(df)) {
      if("Concordance_R2_ELISA" %in% colnames(df))
         df$Concordance_R2_ELISA <- check_concordance(df$Titer_R2, protection_threshold, df$ELISA_EU_mL, log_tol = ed50_tolerance, use_threshold_only = FALSE)
      if("Concordance_R2_Nonlin" %in% colnames(df))
         df$Concordance_R2_Nonlin <- check_concordance(df$Titer_R2, protection_threshold, df$ELISA_EU_mL_nonlin, log_tol = ed50_tolerance, use_threshold_only = FALSE)
    }
    if(isTRUE(use_third_reference) && "Titer_R3" %in% colnames(df)) {
      if("Concordance_R3_ELISA" %in% colnames(df))
         df$Concordance_R3_ELISA <- check_concordance(df$Titer_R3, protection_threshold, df$ELISA_EU_mL, log_tol = ed50_tolerance, use_threshold_only = FALSE)
      if("Concordance_R3_Nonlin" %in% colnames(df))
         df$Concordance_R3_Nonlin <- check_concordance(df$Titer_R3, protection_threshold, df$ELISA_EU_mL_nonlin, log_tol = ed50_tolerance, use_threshold_only = FALSE)
    }
    df
  }
  combined_summary_dup2           <- .overlay_elisa(combined_summary_dup2)
  combined_summary_dup2_forced4pl <- .overlay_elisa(combined_summary_dup2_forced4pl)
  combined_summary_dup2_5pl       <- .overlay_elisa(combined_summary_dup2_5pl)
  combined_summary_dup2_bayesian  <- .overlay_elisa(combined_summary_dup2_bayesian)

   if(has_rlu3) {
     rlu_cols_3 <- c("RLU_dup1","RLU_dup2","RLU_dup3")
     result_dup3 <- process_rlu_config(test_data, rlu_cols_3, "three_replicates", compute_rm=use_reed_muench_comparison,
                                       force_4pl_enhanced=FALSE)
     combined_summary_dup3 <- result_dup3$summary
     combined_dilutions_dup3 <- result_dup3$dilutions
     
     result_dup3_forced4pl <- process_rlu_config(test_data, rlu_cols_3, "three_replicates_forced4pl", compute_rm=use_reed_muench_comparison,
                                                 force_4pl_enhanced=TRUE)
     combined_summary_dup3_forced4pl <- result_dup3_forced4pl$summary
     combined_dilutions_dup3_forced4pl <- result_dup3_forced4pl$dilutions
     
     result_dup3_5pl <- process_rlu_config(test_data, rlu_cols_3, "three_replicates_5pl", compute_rm=use_reed_muench_comparison,
                                             force_4pl_enhanced=FALSE, use_5pl=TRUE)
     combined_summary_dup3_5pl <- result_dup3_5pl$summary
     combined_dilutions_dup3_5pl <- result_dup3_5pl$dilutions
     
     result_dup3_bayesian <- process_rlu_config(test_data, rlu_cols_3, "three_replicates_bayesian", compute_rm=use_reed_muench_comparison,
                                                force_4pl_enhanced=FALSE, use_bayesian=TRUE)
     combined_summary_dup3_bayesian <- result_dup3_bayesian$summary
     combined_dilutions_dup3_bayesian <- result_dup3_bayesian$dilutions
 
     combined_summary_dup3           <- .overlay_elisa(combined_summary_dup3)
     combined_summary_dup3_forced4pl <- .overlay_elisa(combined_summary_dup3_forced4pl)
     combined_summary_dup3_5pl     <- .overlay_elisa(combined_summary_dup3_5pl)
      combined_summary_dup3_bayesian  <- .overlay_elisa(combined_summary_dup3_bayesian)

    comparison_df <- merge(
      combined_summary_dup2[,c("Sample","Experiment_Group","IC50_log_dilution","IC50_Method","Titer_R1","Titer_R2",
                               "Concordance_R1_ELISA","Concordance_R2_ELISA","Protection_Level","R_squared","ELISA_EU_mL")],
      combined_summary_dup3[,c("Sample","IC50_log_dilution","IC50_Method","Titer_R1","Titer_R2",
                               "Concordance_R1_ELISA","Concordance_R2_ELISA","Protection_Level","R_squared")],
      by="Sample", suffixes=c("_2rep","_3rep"))
    
    comparison_df$IC50_Difference <- round(comparison_df$IC50_log_dilution_2rep - comparison_df$IC50_log_dilution_3rep, 4)
    comparison_df$Titer_Difference_R1 <- round(comparison_df$Titer_R1_2rep - comparison_df$Titer_R1_3rep, 8)
    comparison_df$Concordance_Agreement_R1 <- comparison_df$Concordance_R1_ELISA_2rep == comparison_df$Concordance_R1_ELISA_3rep
    comparison_df$Protection_Agreement <- comparison_df$Protection_Level_2rep == comparison_df$Protection_Level_3rep
    comparison_df$Method_Agreement <- comparison_df$IC50_Method_2rep == comparison_df$IC50_Method_3rep
    
    titer_cor <- if(sum(complete.cases(comparison_df$Titer_R1_2rep, comparison_df$Titer_R1_3rep)) > 1) {
      cor(comparison_df$Titer_R1_2rep, comparison_df$Titer_R1_3rep, use="complete.obs")
    } else { NA }
    ic50_cor <- if(sum(complete.cases(comparison_df$IC50_log_dilution_2rep, comparison_df$IC50_log_dilution_3rep)) > 1) {
      cor(comparison_df$IC50_log_dilution_2rep, comparison_df$IC50_log_dilution_3rep, use="complete.obs")
    } else { NA }
    
    .cmp_df_agree <- comparison_df
    if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(.cmp_df_agree)) {
      .cmp_df_agree <- .cmp_df_agree[grepl(paste0("^", sample_prefix_filter), as.character(.cmp_df_agree$Sample), ignore.case = FALSE), ]
    }
    n_concord_agree <- sum(.cmp_df_agree$Concordance_Agreement_R1, na.rm=TRUE)
    n_prot_agree <- sum(.cmp_df_agree$Protection_Agreement, na.rm=TRUE)
    n_total <- nrow(.cmp_df_agree)
    rm(.cmp_df_agree)
    
    
    
    comparison_temp <- merge(
      combined_summary_dup2[,c("Sample","IC50_log_dilution","Titer_R1","Concordance_R1_ELISA","Protection_Level","R_squared")],
      combined_summary_dup2_forced4pl[,c("Sample","IC50_log_dilution","Titer_R1","Concordance_R1_ELISA","Protection_Level","R_squared")],
      by="Sample", suffixes=c("_ACM", "_Forced4PL"), all=TRUE
    )
    
    comparison_ACM_forced <- merge(
      comparison_temp,
      combined_summary_dup2_5pl[,c("Sample","IC50_log_dilution","Titer_R1","Concordance_R1_ELISA","Protection_Level","R_squared")],
      by="Sample", all=TRUE
    )
    
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="IC50_log_dilution"] <- "IC50_log_dilution_5PL"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Titer_R1"] <- "Titer_R1_5PL"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Concordance_R1_ELISA"] <- "Concordance_R1_ELISA_5PL"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Protection_Level"] <- "Protection_Level_5PL"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="R_squared"] <- "R_squared_5PL"
    
    comparison_ACM_forced <- merge(
      comparison_ACM_forced,
      combined_summary_dup2_bayesian[,c("Sample","IC50_log_dilution","Titer_R1","Concordance_R1_ELISA","Protection_Level","R_squared")],
      by="Sample", all=TRUE
    )
    
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="IC50_log_dilution"] <- "IC50_log_dilution_Bayesian"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Titer_R1"] <- "Titer_R1_Bayesian"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Concordance_R1_ELISA"] <- "Concordance_R1_ELISA_Bayesian"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="Protection_Level"] <- "Protection_Level_Bayesian"
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="R_squared"] <- "R_squared_Bayesian"
    
    comparison_ACM_forced <- merge(
      comparison_ACM_forced,
      combined_summary_dup2[,c("Sample","RM_IC50_log","R1_IC50_dilution","R2_IC50_dilution","ELISA_EU_mL")],
      by="Sample", all=TRUE
    )
    colnames(comparison_ACM_forced)[colnames(comparison_ACM_forced)=="RM_IC50_log"] <- "RM_IC50_log"
    comparison_ACM_forced$RM_Titer_R1 <- ifelse(!is.na(comparison_ACM_forced$RM_IC50_log) & !is.na(comparison_ACM_forced$R1_IC50_dilution) & comparison_ACM_forced$R1_IC50_dilution > 0, 
                                                    (10^comparison_ACM_forced$RM_IC50_log / comparison_ACM_forced$R1_IC50_dilution) * reference_concentration, NA)
    if(!is.null(comparison_ACM_forced$ELISA_EU_mL) && !is.null(comparison_ACM_forced$RM_Titer_R1)) {
      comparison_ACM_forced$Concordance_R1_ELISA_RM <- ifelse(
        !is.na(comparison_ACM_forced$RM_Titer_R1) & !is.na(comparison_ACM_forced$ELISA_EU_mL),
        check_concordance(comparison_ACM_forced$RM_Titer_R1, reference_concentration,
                          comparison_ACM_forced$ELISA_EU_mL, log_tol=ed50_tolerance), NA)
    }
    comparison_ACM_forced$RM_Titer_R2 <- ifelse(use_dual_reference & !is.na(comparison_ACM_forced$RM_IC50_log) & !is.na(comparison_ACM_forced$R2_IC50_dilution) & comparison_ACM_forced$R2_IC50_dilution > 0,
                                                    (10^comparison_ACM_forced$RM_IC50_log / comparison_ACM_forced$R2_IC50_dilution) * second_reference_concentration, NA)
    
    comparison_ACM_forced <- comparison_ACM_forced[
      !is.na(comparison_ACM_forced$IC50_log_dilution_ACM) | 
        !is.na(comparison_ACM_forced$IC50_log_dilution_Forced4PL) | 
        !is.na(comparison_ACM_forced$IC50_log_dilution_5PL) |
        !is.na(comparison_ACM_forced$IC50_log_dilution_Bayesian) |
        !is.na(comparison_ACM_forced$RM_IC50_log), ]
    
    if(nrow(comparison_ACM_forced) > 0) {
      comparison_ACM_forced$IC50_Diff_ACM_vs_Forced <- NA
      comparison_ACM_forced$IC50_Diff_ACM_vs_5PL <- NA
      comparison_ACM_forced$IC50_Diff_Forced_vs_5PL <- NA
      comparison_ACM_forced$IC50_Diff_ACM_vs_Bayesian <- NA
      comparison_ACM_forced$IC50_Diff_Forced_vs_Bayesian <- NA
      comparison_ACM_forced$IC50_Diff_5PL_vs_Bayesian <- NA
      comparison_ACM_forced$IC50_Diff_ACM_vs_RM <- NA
      comparison_ACM_forced$IC50_Diff_Forced_vs_RM <- NA
      comparison_ACM_forced$IC50_Diff_5PL_vs_RM <- NA
      comparison_ACM_forced$IC50_Diff_Bayesian_vs_RM <- NA

      
      comparison_ACM_forced$IC50_Diff_ACM_vs_Forced[valid_ACM_forced] <- 
        round(comparison_ACM_forced$IC50_log_dilution_ACM[valid_ACM_forced] - 
                comparison_ACM_forced$IC50_log_dilution_Forced4PL[valid_ACM_forced], 4)
      
      comparison_ACM_forced$IC50_Diff_ACM_vs_5PL[valid_ACM_5pl] <- 
        round(comparison_ACM_forced$IC50_log_dilution_ACM[valid_ACM_5pl] - 
                comparison_ACM_forced$IC50_log_dilution_5PL[valid_ACM_5pl], 4)
      
      comparison_ACM_forced$IC50_Diff_Forced_vs_5PL[valid_forced_5pl] <- 
        round(comparison_ACM_forced$IC50_log_dilution_Forced4PL[valid_forced_5pl] - 
                comparison_ACM_forced$IC50_log_dilution_5PL[valid_forced_5pl], 4)
      
      comparison_ACM_forced$IC50_Diff_ACM_vs_Bayesian[valid_ACM_bayesian] <- 
        round(comparison_ACM_forced$IC50_log_dilution_ACM[valid_ACM_bayesian] - 
                comparison_ACM_forced$IC50_log_dilution_Bayesian[valid_ACM_bayesian], 4)
      
      comparison_ACM_forced$IC50_Diff_Forced_vs_Bayesian[valid_forced_bayesian] <- 
        round(comparison_ACM_forced$IC50_log_dilution_Forced4PL[valid_forced_bayesian] - 
                comparison_ACM_forced$IC50_log_dilution_Bayesian[valid_forced_bayesian], 4)
      
      comparison_ACM_forced$IC50_Diff_5PL_vs_Bayesian[valid_5pl_bayesian] <- 
        round(comparison_ACM_forced$IC50_log_dilution_5PL[valid_5pl_bayesian] - 
                comparison_ACM_forced$IC50_log_dilution_Bayesian[valid_5pl_bayesian], 4)
      
      comparison_ACM_forced$IC50_Diff_ACM_vs_RM[valid_ACM_rm] <- 
        round(comparison_ACM_forced$IC50_log_dilution_ACM[valid_ACM_rm] - 
                comparison_ACM_forced$RM_IC50_log[valid_ACM_rm], 4)
      
      comparison_ACM_forced$IC50_Diff_Forced_vs_RM[valid_forced_rm] <- 
        round(comparison_ACM_forced$IC50_log_dilution_Forced4PL[valid_forced_rm] - 
                comparison_ACM_forced$RM_IC50_log[valid_forced_rm], 4)
      
      comparison_ACM_forced$IC50_Diff_5PL_vs_RM[valid_5pl_rm] <- 
        round(comparison_ACM_forced$IC50_log_dilution_5PL[valid_5pl_rm] - 
                comparison_ACM_forced$RM_IC50_log[valid_5pl_rm], 4)
      
      comparison_ACM_forced$IC50_Diff_Bayesian_vs_RM[valid_bayesian_rm] <- 
        round(comparison_ACM_forced$IC50_log_dilution_Bayesian[valid_bayesian_rm] - 
                comparison_ACM_forced$RM_IC50_log[valid_bayesian_rm], 4)
      
      valid_conc <- !is.na(comparison_ACM_forced$Concordance_R1_ELISA_ACM) & 
        comparison_ACM_forced$Concordance_R1_ELISA_ACM %in% c("Concordant", "Discordant")
      
      ACM_conc <- sum(comparison_ACM_forced$Concordance_R1_ELISA_ACM == "Concordant", na.rm=TRUE)
      forced_conc <- sum(comparison_ACM_forced$Concordance_R1_ELISA_Forced4PL == "Concordant", na.rm=TRUE)
      f5pl_conc <- sum(comparison_ACM_forced$Concordance_R1_ELISA_5PL == "Concordant", na.rm=TRUE)
      bayesian_conc <- sum(comparison_ACM_forced$Concordance_R1_ELISA_Bayesian == "Concordant", na.rm=TRUE)
      rm_conc <- if("Concordance_R1_ELISA_RM" %in% colnames(comparison_ACM_forced))
        sum(comparison_ACM_forced$Concordance_R1_ELISA_RM == "Concordant", na.rm=TRUE) else NA
      total_valid <- sum(valid_conc, na.rm=TRUE)
      total_valid_rm <- if("Concordance_R1_ELISA_RM" %in% colnames(comparison_ACM_forced))
        sum(!is.na(comparison_ACM_forced$Concordance_R1_ELISA_RM) &
              comparison_ACM_forced$Concordance_R1_ELISA_RM %in% c("Concordant","Discordant"), na.rm=TRUE) else total_valid
      
      if(total_valid > 0) {
         concordance_rates <- data.frame(
           Method = c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench"),
           Concordance_Count = c(ACM_conc, forced_conc, f5pl_conc, bayesian_conc, rm_conc),
           Total = c(total_valid, total_valid, total_valid, total_valid, total_valid_rm),
           Rate = c(ACM_conc/total_valid*100, forced_conc/total_valid*100, f5pl_conc/total_valid*100, bayesian_conc/total_valid*100,
                    ifelse(total_valid_rm > 0, rm_conc/total_valid_rm*100, NA)),
           stringsAsFactors = FALSE
         )
      }
    }
    
    combined_summary_all <- rbind(combined_summary_dup2, combined_summary_dup3, 
                                  combined_summary_dup2_forced4pl, combined_summary_dup3_forced4pl,
                                  combined_summary_dup2_5pl, combined_summary_dup3_5pl,
                                  combined_summary_dup2_bayesian, combined_summary_dup3_bayesian)
    combined_dilutions_all <- rbind(combined_dilutions_dup2, combined_dilutions_dup3,
                                    combined_dilutions_dup2_forced4pl, combined_dilutions_dup3_forced4pl,
                                    combined_dilutions_dup2_5pl, combined_dilutions_dup3_5pl,
                                    combined_dilutions_dup2_bayesian, combined_dilutions_dup3_bayesian)
    
  } else {
    if(use_rlu3 && !has_rlu3) 
    combined_summary_all <- combined_summary_dup2
    combined_dilutions_all <- combined_dilutions_dup2
    comparison_df <- NULL
    combined_summary_dup3 <- NULL
    combined_summary_dup3_forced4pl <- NULL
    combined_summary_dup3_5pl <- NULL
    comparison_ACM_forced <- NULL
    concordance_rates <- NULL
  }
  
  combined_summary <- combined_summary_dup2
  combined_dilutions <- combined_dilutions_dup2
  method_stats <- result_dup2$method_stats
  
  make_concordance_table <- function(data, conc_col) {
    vals <- data[[conc_col]]; vals <- vals[!is.na(vals) & vals!="Unknown" & vals!="N/A"]
    if(length(vals)==0) return(NULL)
    ccount <- sum(vals=="Concordant"); dcount <- sum(vals=="Discordant"); total <- ccount+dcount
    if(total==0) return(NULL)
    data.frame(Concordant=ccount, Discordant=dcount, Total=total, Agreement_pct=round(ccount/total*100,1), stringsAsFactors=FALSE)
  }
  
  conc_tables <- list()
  conc_tables[["R1_vs_ELISA"]] <- make_concordance_table(combined_summary, "Concordance_R1_ELISA")
  if(use_dual_reference) conc_tables[[paste0(R2_name, "_vs_ELISA")]] <- make_concordance_table(combined_summary, "Concordance_R2_ELISA")
  if(use_third_reference) conc_tables[[paste0(R3_name, "_vs_ELISA")]] <- make_concordance_table(combined_summary, "Concordance_R3_ELISA")
  if(has_elisa_nonlin) {
    conc_tables[[paste0(R1_name, "_vs_Nonlin")]] <- make_concordance_table(combined_summary, "Concordance_R1_Nonlin")
    if(use_dual_reference) conc_tables[[paste0(R2_name, "_vs_Nonlin")]] <- make_concordance_table(combined_summary, "Concordance_R2_Nonlin")
    if(use_third_reference) conc_tables[[paste0(R3_name, "_vs_Nonlin")]] <- make_concordance_table(combined_summary, "Concordance_R3_Nonlin")
  }
  concordance_summary <- do.call(rbind, lapply(names(conc_tables), function(nm) { if(is.null(conc_tables[[nm]])) return(NULL); cbind(Comparison=nm, conc_tables[[nm]]) }))
  
  tolerance_levels <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
  concordance_multi_tol <- data.frame()
  
  replicate_configs <- list(
    "2 replicates" = list(
      mixed = combined_summary_dup2,
      forced = combined_summary_dup2_forced4pl,
      sk = combined_summary_dup2_5pl,
      bayesian = combined_summary_dup2_bayesian
    )
  )
  if(exists("combined_summary_dup3") && !is.null(combined_summary_dup3)) {
    replicate_configs[["3 replicates"]] <- list(
      mixed = combined_summary_dup3,
      forced = combined_summary_dup3_forced4pl,
      sk = combined_summary_dup3_5pl,
      bayesian = combined_summary_dup3_bayesian
    )
  }
  
  for(rep_name in names(replicate_configs)) {
    rep_config <- replicate_configs[[rep_name]]
    for(tol in tolerance_levels) {
      for(method_name in c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")) {
         for(ref_name in c("R1", "R2", "R3")) {
            if(ref_name == "R2" && !use_dual_reference) next
            if(ref_name == "R3" && !use_third_reference) next
           for(elisa_type in c("linear", "nonlin")) {
             if(elisa_type == "nonlin" && !has_elisa_nonlin) next
             
             if(method_name == "ACM") {
               data_df <- rep_config$mixed
             } else if(method_name == "F4PL") {
               data_df <- rep_config$forced
              } else if(method_name == "F5PL") {
                data_df <- rep_config$sk
             } else if(method_name == "Bayes-4PL") {
               data_df <- rep_config$bayesian
             } else {
               data_df <- rep_config$mixed
             }
            
            titer_col <- if(ref_name == "R1") "Titer_R1" else if(ref_name == "R2") "Titer_R2" else "Titer_R3"
            elisa_col <- if(elisa_type == "linear") "ELISA_EU_mL" else "ELISA_EU_mL_nonlin"
            conc_col <- paste0("Concordance_", ref_name, ifelse(elisa_type == "linear", "_ELISA", "_Nonlin"))
            ref_conc <- if(ref_name == "R1") reference_concentration else if(ref_name == "R2") second_reference_concentration else third_reference_concentration
            
            if(method_name == "Reed-Muench") {
              rm_col <- "RM_IC50_log"
              ref_dil_col <- paste0(ref_name, "_IC50_dilution")
              ref_conc <- if(ref_name == "R1") reference_concentration else if(ref_name == "R2") second_reference_concentration else third_reference_concentration
              if(rm_col %in% colnames(data_df) && ref_dil_col %in% colnames(data_df)) {
                data_df$RM_Titer_temp <- ifelse(!is.na(data_df[[rm_col]]) & !is.na(data_df[[ref_dil_col]]) & data_df[[ref_dil_col]] > 0,
                                                (10^data_df[[rm_col]] / data_df[[ref_dil_col]]) * ref_conc, NA)
                titer_col <- "RM_Titer_temp"
              } else {
                next
              }
            }
            
            if(all(c(titer_col, elisa_col) %in% colnames(data_df))) {
               d <- data_df[!is.na(data_df[[titer_col]]) & !is.na(data_df[[elisa_col]]) & !data_df$Is_R1 & !data_df$Is_R2 & !data_df$Is_R3 & !data_df$Is_mAb, ]
              if(nrow(d) > 0) {
                log_titer <- log10(d[[titer_col]])
                log_elisa <- log10(d[[elisa_col]])
                tol_lower <- 10^(log_titer - tol)
                tol_upper <- 10^(log_titer + tol)
                concordant <- d[[elisa_col]] >= tol_lower & d[[elisa_col]] <= tol_upper
                concordant[is.na(log_titer) | is.na(log_elisa)] <- NA
                
                n_conc <- sum(concordant, na.rm = TRUE)
                n_total <- sum(!is.na(concordant))
                rate <- if(n_total > 0) round(n_conc / n_total * 100, 1) else NA
                
                pct_diff <- 100 * (d[[titer_col]] / d[[elisa_col]] - 1)
                mean_pct_diff <- if(n_total > 0) round(mean(pct_diff, na.rm = TRUE), 2) else NA
                
                concordance_multi_tol <- rbind(concordance_multi_tol, data.frame(
                  Replicates = rep_name,
                  Method = method_name,
                  Reference = ref_name,
                  ELISA_Type = ifelse(elisa_type == "linear", "Linear", "Non-linear"),
                  Tolerance = tol,
                  Concordant = n_conc,
                  Discordant = n_total - n_conc,
                  Total = n_total,
                  Rate = rate,
                  Mean_Titer_Difference_pct = mean_pct_diff,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    }
  }
  
  sensitivity_specificity <- data.frame()
  
  for(method_name in c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")) {
     for(ref_name in c("R1", "R2", "R3")) {
       if(ref_name == "R2" && !use_dual_reference) next
       if(ref_name == "R3" && !use_third_reference) next
      for(elisa_type in c("linear", "nonlin")) {
        if(elisa_type == "nonlin" && !has_elisa_nonlin) next
        
        if(method_name == "ACM") {
          data_df <- combined_summary_dup2
        } else if(method_name == "F4PL") {
          data_df <- combined_summary_dup2_forced4pl
        } else if(method_name == "F5PL") {
          data_df <- combined_summary_dup2_5pl
        } else if(method_name == "Bayes-4PL") {
          data_df <- combined_summary_dup2_bayesian
        } else {
          data_df <- combined_summary_dup2
        }
        
         titer_col <- if(ref_name == "R1") "Titer_R1" else if(ref_name == "R2") "Titer_R2" else "Titer_R3"
         elisa_col <- if(elisa_type == "linear") "ELISA_EU_mL" else "ELISA_EU_mL_nonlin"
         ref_conc <- if(ref_name == "R1") reference_concentration else if(ref_name == "R2") second_reference_concentration else third_reference_concentration
        
        if(method_name == "Reed-Muench") {
          rm_col <- "RM_IC50_log"
          ref_dil_col <- paste0(ref_name, "_IC50_dilution")
          if(rm_col %in% colnames(data_df) && ref_dil_col %in% colnames(data_df)) {
            data_df$RM_Titer_temp <- ifelse(!is.na(data_df[[rm_col]]) & !is.na(data_df[[ref_dil_col]]) & data_df[[ref_dil_col]] > 0,
                                            (10^data_df[[rm_col]] / data_df[[ref_dil_col]]) * ref_conc, NA)
            titer_col <- "RM_Titer_temp"
          } else {
            next
          }
        }
        
        if(all(c(titer_col, elisa_col) %in% colnames(data_df))) {
           d <- data_df[!is.na(data_df[[titer_col]]) & !is.na(data_df[[elisa_col]]) & !data_df$Is_R1 & !data_df$Is_R2 & !data_df$Is_R3 & !data_df$Is_mAb, ]
          if(nrow(d) > 0) {
            elisa_pos <- d[[elisa_col]] >= ref_conc
            titer_pos <- d[[titer_col]] >= ref_conc
            
            tp <- sum(elisa_pos & titer_pos, na.rm = TRUE)
            fp <- sum(!elisa_pos & titer_pos, na.rm = TRUE)
            tn <- sum(!elisa_pos & !titer_pos, na.rm = TRUE)
            fn <- sum(elisa_pos & !titer_pos, na.rm = TRUE)
            
            sensitivity <- if((tp + fn) > 0) round(tp / (tp + fn) * 100, 1) else NA
            specificity <- if((tn + fp) > 0) round(tn / (tn + fp) * 100, 1) else NA
            ppv <- if((tp + fp) > 0) round(tp / (tp + fp) * 100, 1) else NA
            npv <- if((tn + fn) > 0) round(tn / (tn + fn) * 100, 1) else NA
            accuracy <- if((tp + tn + fp + fn) > 0) round((tp + tn) / (tp + tn + fp + fn) * 100, 1) else NA
            
            sensitivity_specificity <- rbind(sensitivity_specificity, data.frame(
              Method = method_name,
              Reference = ref_name,
              ELISA_Type = ifelse(elisa_type == "linear", "Linear", "Non-linear"),
              Threshold = ref_conc,
              TP = tp, FP = fp, TN = tn, FN = fn,
              Sensitivity = sensitivity,
              Specificity = specificity,
              PPV = ppv,
              NPV = npv,
              Accuracy = accuracy,
              Total = tp + fp + tn + fn,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  write.csv(concordance_multi_tol, file.path(output_dir, "concordance_multi_tolerance.csv"), row.names = FALSE)
  write.csv(sensitivity_specificity, file.path(output_dir, "sensitivity_specificity.csv"), row.names = FALSE)
  
  
   excel_sheets <- list(
     Sample_Results_ACM = combined_summary_dup2,
     Sample_Results_Forced4PL = combined_summary_dup2_forced4pl,
     Sample_Results_5PL = combined_summary_dup2_5pl,
     Dilution_Data_two_replicates = combined_dilutions_dup2,
     Concordance_Summary = concordance_summary,
     Concordance_Multi_Tolerance = concordance_multi_tol,
     Sensitivity_Specificity = sensitivity_specificity,
     Method_Summary = data.frame(Method=names(method_stats),Count=as.numeric(unlist(method_stats)),stringsAsFactors=FALSE),
     ELISA_Scale = data.frame(ELISA_Scale=elisa_scale,Titer_Scale=titer_scale,stringsAsFactors=FALSE))
   excel_sheets[["Sample_Results_Bayes-4PL"]] <- combined_summary_dup2_bayesian
  
  if(use_reed_muench_comparison) {
     rm_comparison <- combined_summary[!is.na(combined_summary$RM_IC50_log), 
                                       c("Sample","Experiment_Group","IC50_log_dilution","IC50_Method","RM_IC50_log","RM_Method","R_squared","ED50_Pct_Infection","Titer_R1")]
    rm_comparison$IC50_Difference <- round(rm_comparison$IC50_log_dilution - rm_comparison$RM_IC50_log, 4)
    rm_comparison$Titer_Difference <- round(rm_comparison$Titer_R1 - (1 / (10^rm_comparison$RM_IC50_log)), 6)
    excel_sheets$RM_Comparison <- rm_comparison
  }
  
  if(exists("comparison_ACM_forced") && !is.null(comparison_ACM_forced) && nrow(comparison_ACM_forced)>0) {
    excel_sheets$Method_Comparison_All_Five <- comparison_ACM_forced
    if(!is.null(concordance_rates) && nrow(concordance_rates)>0) {
      excel_sheets$Concordance_Rates_Comparison <- concordance_rates
      
      compute_cell_r2 <- function(titer_v, elisa_v) {
        if(length(titer_v) < 3 || length(elisa_v) < 3) return(NA)
        if(any(is.na(titer_v)) || any(is.na(elisa_v))) return(NA)
        if(sd(titer_v) == 0 || sd(elisa_v) == 0) return(NA)
        cor(titer_v, elisa_v, method = "pearson")^2
      }
      
      method_cell_r2 <- list()
      if(exists("comparison_ACM_forced") && is.data.frame(comparison_ACM_forced)) {
        comp <- comparison_ACM_forced
        elisa_lin <- suppressWarnings(as.numeric(comp$ELISA_EU_mL))
        elisa_nonlin <- suppressWarnings(as.numeric(comp$ELISA_EU_mL_nonlin))
        
        method_titer_cols <- list(
          ACM = list(R1 = "Titer_R1_ACM", R2 = "Titer_R2_ACM"),
          Forced4PL = list(R1 = "Titer_R1_Forced4PL", R2 = "Titer_R2_Forced4PL"),
          Forced5PL = list(R1 = "Titer_R1_5PL", R2 = "Titer_R2_5PL"),
          Bayes4PL = list(R1 = "Titer_R1_Bayesian", R2 = "Titer_R2_Bayesian"),
          ReedMuench = list(R1 = "RM_Titer_R1", R2 = "RM_Titer_R2")
        )
        
        for(m_name in names(method_titer_cols)) {
          cols <- method_titer_cols[[m_name]]
          r2_values <- c()
          
          for(ref in c("R1", "R2")) {
            titer_col <- cols[[ref]]
            if(!titer_col %in% colnames(comp)) next
            titer_v <- suppressWarnings(as.numeric(comp[[titer_col]]))
            
            if(ref == "R1") {
              if(!all(is.na(elisa_lin))) r2_values <- c(r2_values, compute_cell_r2(titer_v, elisa_lin))
              if(!all(is.na(elisa_nonlin))) r2_values <- c(r2_values, compute_cell_r2(titer_v, elisa_nonlin))
            } else {
              if(!all(is.na(elisa_lin))) r2_values <- c(r2_values, compute_cell_r2(titer_v, elisa_lin))
              if(!all(is.na(elisa_nonlin))) r2_values <- c(r2_values, compute_cell_r2(titer_v, elisa_nonlin))
            }
          }
          
          method_cell_r2[[m_name]] <- mean(r2_values, na.rm = TRUE)
        }
      }
      
      method_names <- c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")
      method_keys <- c("ACM", "Forced4PL", "Forced5PL", "Bayes-4PL", "ReedMuench")
      
      comprehensive_r2 <- sapply(method_keys, function(k) {
        if(k %in% names(method_cell_r2)) method_cell_r2[[k]] else NA
      })
      
      ranking_df <- data.frame(
        Method = method_names,
        Concordance_Rate = concordance_rates$Rate[match(method_names, concordance_rates$Method)],
        Comprehensive_R2 = comprehensive_r2,
        stringsAsFactors = FALSE
      )
      
      ranking_df$Concordance_Score <- ranking_df$Concordance_Rate / 100
      ranking_df$R2_Score <- ranking_df$Comprehensive_R2
      
      ranking_df$Concordance_Rank <- rank(-ranking_df$Concordance_Score, ties.method = "min")
      ranking_df$R2_Rank <- rank(-ranking_df$R2_Score, ties.method = "min")
      ranking_df$Overall_Rank <- rank(-ranking_df$R2_Score, ties.method = "min")
      ranking_df$Best_Method <- ifelse(ranking_df$Overall_Rank == min(ranking_df$Overall_Rank, na.rm = TRUE), "BEST", "")
      
      ranking_df <- ranking_df[order(ranking_df$Overall_Rank), c("Method", "Concordance_Rate", "Comprehensive_R2", "Concordance_Rank", "R2_Rank", "Overall_Rank", "Best_Method")]
      
      excel_sheets$Method_Ranking_Summary <- ranking_df
    }
  }
  
  if(has_rlu3 && !is.null(combined_summary_dup3)) {
    excel_sheets$Sample_Results_ACM_3rep <- combined_summary_dup3
    excel_sheets[["Sample_Results_Bayes-4PL_3rep"]] <- combined_summary_dup3_bayesian
    excel_sheets$Sample_Results_Forced4PL_3rep <- combined_summary_dup3_forced4pl
    excel_sheets$Sample_Results_5PL_3rep <- combined_summary_dup3_5pl
    excel_sheets[["Sample_Results_Bayes-4PL_3rep"]] <- combined_summary_dup3_bayesian
    excel_sheets$Dilution_Data_three_replicates <- combined_dilutions_dup3
    excel_sheets$RLU_Config_Comparison <- comparison_df
    if(!is.null(comparison_df)) {
    }
  }
  
  excel_sheets <- Filter(function(x) is.data.frame(x) && ncol(x) > 0, excel_sheets)
  if(length(excel_sheets) > 0) {
    tryCatch({
      write_xlsx(excel_sheets, file.path(output_dir,"neutralization_results.xlsx"))
    }, error = function(e) {
    })
  }
  
  if(generate_publication_graphs && nrow(combined_summary)>0) {
    tryCatch({
      
      combined_summary <- ref_filter(combined_summary, sample_prefix_filter = sample_prefix_filter)
      combined_summary_dup2 <- ref_filter(combined_summary_dup2, sample_prefix_filter = sample_prefix_filter)
      combined_summary_dup2_forced4pl <- ref_filter(combined_summary_dup2_forced4pl, sample_prefix_filter = sample_prefix_filter)
      combined_summary_dup2_5pl <- ref_filter(combined_summary_dup2_5pl, sample_prefix_filter = sample_prefix_filter)
      combined_summary_dup2_bayesian <- ref_filter(combined_summary_dup2_bayesian, sample_prefix_filter = sample_prefix_filter)
      if(exists("combined_summary_dup3") && !is.null(combined_summary_dup3)) combined_summary_dup3 <- ref_filter(combined_summary_dup3, sample_prefix_filter = sample_prefix_filter)
      if(exists("combined_summary_dup3_forced4pl") && !is.null(combined_summary_dup3_forced4pl)) combined_summary_dup3_forced4pl <- ref_filter(combined_summary_dup3_forced4pl, sample_prefix_filter = sample_prefix_filter)
      if(exists("combined_summary_dup3_5pl") && !is.null(combined_summary_dup3_5pl)) combined_summary_dup3_5pl <- ref_filter(combined_summary_dup3_5pl, sample_prefix_filter = sample_prefix_filter)
      if(exists("combined_summary_dup3_bayesian") && !is.null(combined_summary_dup3_bayesian)) combined_summary_dup3_bayesian <- ref_filter(combined_summary_dup3_bayesian, sample_prefix_filter = sample_prefix_filter)
      if(exists("comparison_df") && !is.null(comparison_df)) comparison_df <- ref_filter(comparison_df, sample_prefix_filter = sample_prefix_filter)
      if(exists("comparison_ACM_forced") && !is.null(comparison_ACM_forced)) comparison_ACM_forced <- ref_filter(comparison_ACM_forced, sample_prefix_filter = sample_prefix_filter)
      
      if(nrow(combined_summary) > 0) {
        method_counts <- table(combined_summary$IC50_Method)
        method_counts <- method_counts[names(method_counts) != "Failed"]
        method_stats <- as.list(method_counts)
        names(method_stats) <- names(method_counts)
      }
      method_stats_forced <- list()
      if(nrow(combined_summary_dup2_forced4pl) > 0) {
        method_counts_forced <- table(combined_summary_dup2_forced4pl$IC50_Method)
        method_counts_forced <- method_counts_forced[names(method_counts_forced) != "Failed"]
        method_stats_forced <- as.list(method_counts_forced)
        names(method_stats_forced) <- names(method_counts_forced)
      }
      method_stats_f5pl <- list()
      if(nrow(combined_summary_dup2_5pl) > 0) {
        method_counts_f5pl <- table(combined_summary_dup2_5pl$IC50_Method)
        method_counts_f5pl <- method_counts_f5pl[names(method_counts_f5pl) != "Failed"]
        method_stats_f5pl <- as.list(method_counts_f5pl)
        names(method_stats_f5pl) <- names(method_counts_f5pl)
      }
      method_stats_bayesian <- list()
      if(nrow(combined_summary_dup2_bayesian) > 0) {
        method_counts_bayesian <- table(combined_summary_dup2_bayesian$IC50_Method)
        method_counts_bayesian <- method_counts_bayesian[names(method_counts_bayesian) != "Failed"]
        method_stats_bayesian <- as.list(method_counts_bayesian)
        names(method_stats_bayesian) <- names(method_counts_bayesian)
      }
      if(use_reed_muench_comparison && nrow(combined_summary_dup2) > 0) {
        method_counts_rm <- table(combined_summary_dup2$IC50_Method)
        method_counts_rm <- method_counts_rm[names(method_counts_rm) != "Failed"]
        method_stats_rm <- as.list(method_counts_rm)
        names(method_stats_rm) <- names(method_counts_rm)
      }
      
      if(exists("comparison_ACM_forced") && !is.null(comparison_ACM_forced) && nrow(comparison_ACM_forced) > 0) {
        valid_conc_plot <- !is.na(comparison_ACM_forced$Concordance_R1_ELISA_ACM) & 
          comparison_ACM_forced$Concordance_R1_ELISA_ACM %in% c("Concordant", "Discordant")
        ACM_conc_plot <- sum(comparison_ACM_forced$Concordance_R1_ELISA_ACM == "Concordant", na.rm=TRUE)
        forced_conc_plot <- sum(comparison_ACM_forced$Concordance_R1_ELISA_Forced4PL == "Concordant", na.rm=TRUE)
        f5pl_conc_plot <- sum(comparison_ACM_forced$Concordance_R1_ELISA_5PL == "Concordant", na.rm=TRUE)
        bayesian_conc_plot <- sum(comparison_ACM_forced$Concordance_R1_ELISA_Bayesian == "Concordant", na.rm=TRUE)
        rm_conc_plot <- if("Concordance_R1_ELISA_RM" %in% colnames(comparison_ACM_forced))
          sum(comparison_ACM_forced$Concordance_R1_ELISA_RM == "Concordant", na.rm=TRUE) else NA
        total_valid_plot <- sum(valid_conc_plot, na.rm=TRUE)
        total_valid_rm_plot <- if("Concordance_R1_ELISA_RM" %in% colnames(comparison_ACM_forced))
          sum(!is.na(comparison_ACM_forced$Concordance_R1_ELISA_RM) &
                comparison_ACM_forced$Concordance_R1_ELISA_RM %in% c("Concordant","Discordant"), na.rm=TRUE) else total_valid_plot
        if(total_valid_plot > 0) {
          concordance_rates <- data.frame(
            Method = c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench"),
            Concordance_Count = c(ACM_conc_plot, forced_conc_plot, f5pl_conc_plot, bayesian_conc_plot, rm_conc_plot),
            Total = c(total_valid_plot, total_valid_plot, total_valid_plot, total_valid_plot, total_valid_rm_plot),
            Rate = c(ACM_conc_plot/total_valid_plot*100, forced_conc_plot/total_valid_plot*100, f5pl_conc_plot/total_valid_plot*100,
                     bayesian_conc_plot/total_valid_plot*100, ifelse(total_valid_rm_plot > 0, rm_conc_plot/total_valid_rm_plot*100, NA))
          )
        } else {
          concordance_rates <- NULL
        }
      } else {
        concordance_rates <- NULL
      }
      
      method_colors <- c("Parallel_4PL"="#1B9E77","4PL_Good_Fit"="#228B22","4PL_Moderate_Fit"="#4169E1",
                         "Extrapolated_Below_Regression"="#FFA500","Extrapolated_Below_Conservative"="#DC143C",
                         "Extrapolated_Above_Regression"="#9370DB","Extrapolated_Above_Conservative"="#DC143C",
                         "Fallback_4PL"="#808080","Failed"="#808080",
                         "ACM_4PL_Good_Fit"="#2E86AB","ACM_5PL_Good_Fit"="#F18F01",
                         "ACM_4PL_Constrained_Fit"="#1A7F3A","ACM_5PL_Constrained_Fit"="#C67B0B",
                         "ACM_4PL5PL_Weighted"="#1B9E77","ACM_5PL_Corrected"="#D95F02",
                         "ACM_4PL_Good_Fit_Corrected"="#2E86AB",
                         "5PL_Good_Fit"="#6A3D9A","5PL_5PL_Good_Fit"="#FF7F00",
                         "5PL_Constrained_Fit"="#CAB2D6","5PL_5PL_Constrained_Fit"="#FDBF6F",
                         "Bayes_4PL_Good_Fit"="#E41A1C","Bayes_4PL_Constrained"="#377EB8",
                         "RM_Transition_Fallback"="#FF6B6B","ED50_Interpolation_Fallback"="#4ECDC4",
                         "Conservative_Fallback"="#45B7D1","4PL_ED50_Interpolation"="#4ECDC4",
                         "5PL_ED50_Fallback"="#45B7D1","Parallel_ED50_Interpolation"="#4ECDC4",
                         "Parallel_Fallback"="#808080")
      conc_colors <- c("Concordant"="#228B22","Discordant"="#DC143C","Unknown"="#808080")
      config_colors <- c("two_replicates"="#2E86AB","three_replicates"="#F18F01")
      method_compare_colors <- c("ACM"="#000080", "F4PL"="#F18F01", "F5PL"="#6A3D9A", "Bayes-4PL"="#E41A1C", "Reed-Muench"="#88CCEE")
      
      add_regression <- function() {
        geom_smooth(method="lm", se=TRUE, color="black", fill="grey70", alpha=0.3, linewidth=0.8)
      }
      
      fig2_x_axis_text <- if(fig2_show_x_labels) {
        element_text(angle=0, hjust=0.5, size=axis_text_size)
      } else {
        element_blank()
      }
      
      all_method_data <- list()
      method_main_names <- c("ACM", "F4PL", "F5PL", "Bayes-4PL")
      method_stats_list <- list(method_stats, if(exists("method_stats_forced")) method_stats_forced else list(), if(exists("method_stats_f5pl")) method_stats_f5pl else list(), if(exists("method_stats_bayesian")) method_stats_bayesian else list())
      
      for(i in 1:length(method_main_names)) {
        main_name <- method_main_names[i]
        stats <- method_stats_list[[i]]
        if(length(stats) > 0) {
          for(sub_name in names(stats)) {
            count <- stats[[sub_name]]
            if(!is.na(count) && count > 0) {
              all_method_data <- c(all_method_data, list(data.frame(
                Main_Method = main_name,
                Sub_Method = sub_name,
                Count = count
              )))
            }
          }
        }
      }
      
      if(length(all_method_data) > 0) {
        method_dist_df <- do.call(rbind, all_method_data)
        method_dist_df$Sub_Method <- factor(method_dist_df$Sub_Method, 
                                            levels = unique(method_dist_df$Sub_Method[order(method_dist_df$Count, decreasing = TRUE)]))
        
        p1 <- ggplot(method_dist_df, aes(x = Main_Method, y = Count, fill = Sub_Method)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
          scale_fill_manual(values = method_colors, name = "Sub-method") +
          geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
          labs(x = "", y = "Number of samples") +
          theme_publication_custom(base_font_size) +
          theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
                legend.position = "bottom", legend.text = element_text(size = 8))
        ggsave(file.path(plot_dir, "Fig1_Method_Distribution.png"), p1, width = plot_width_one_panel * 1.2, height = plot_height_one_panel, dpi = dpi)
      }
      
      ci_panels <- list()
      
      ci_data <- combined_summary[!is.na(combined_summary$IC50_95CI_Lower) & !combined_summary$Is_R1 & !combined_summary$Is_R2 & !combined_summary$Is_R3 & !combined_summary$Is_mAb,]
      if(nrow(ci_data)>0){
        ci_data$Sample_Num <- factor(1:nrow(ci_data), levels = 1:nrow(ci_data))
        p2 <- ggplot(ci_data, aes(x=Sample_Num, y=IC50_log_dilution, color=IC50_Method)) +
          geom_point(size=point_size) +
          geom_errorbar(aes(ymin=IC50_95CI_Lower, ymax=IC50_95CI_Upper), width=errorbar_width) +
          scale_color_manual(values=method_colors) +
          labs(title="ACM", x="", y="IC50 log10(dilution)") +
          theme_publication_custom(base_font_size) +
          theme(axis.text.x=fig2_x_axis_text,
                plot.margin = margin(t=20, r=20, b=40, l=20))
        ci_panels[[1]] <- p2
      }
      
      ci_data_forced <- combined_summary_dup2_forced4pl[!is.na(combined_summary_dup2_forced4pl$IC50_95CI_Lower) & !combined_summary_dup2_forced4pl$Is_R1 & !combined_summary_dup2_forced4pl$Is_R2 & !combined_summary_dup2_forced4pl$Is_R3 & !combined_summary_dup2_forced4pl$Is_mAb,]
      if(nrow(ci_data_forced)>0){
        ci_data_forced$Sample_Num <- factor(1:nrow(ci_data_forced), levels = 1:nrow(ci_data_forced))
        forced_sub_colors <- c("ACM_4PL_Good_Fit"="#2E86AB","ACM_5PL_Good_Fit"="#F18F01",
                               "ACM_4PL_Constrained_Fit"="#1A7F3A","ACM_5PL_Constrained_Fit"="#C67B0B",
                               "ACM_4PL5PL_Weighted"="#1B9E77","ACM_5PL_Corrected"="#D95F02",
                               "ACM_4PL_Good_Fit_Corrected"="#2E86AB",
                               "4PL_ED50_Interpolation"="#4ECDC4","5PL_ED50_Fallback"="#45B7D1",
                               "ED50_Fallback"="#4ECDC4","Fallback_4PL"="#808080","Failed"="#808080")
        p2b <- ggplot(ci_data_forced, aes(x=Sample_Num, y=IC50_log_dilution, color=IC50_Method)) +
          geom_point(size=point_size) +
          geom_errorbar(aes(ymin=IC50_95CI_Lower, ymax=IC50_95CI_Upper), width=errorbar_width) +
          scale_color_manual(values=forced_sub_colors, name="F4PL Sub-method") +
          labs(title="F4PL", x="Sample (numbered)", y="IC50 log10(dilution)") +
          theme_publication_custom(base_font_size) +
          theme(axis.text.x=fig2_x_axis_text,
                plot.margin = margin(t=20, r=20, b=40, l=20),
                legend.position="bottom") +
          coord_cartesian(ylim = c(1, 5))
        ci_panels[[2]] <- p2b
      }
      
      ci_data_f5pl <- combined_summary_dup2_5pl[!is.na(combined_summary_dup2_5pl$IC50_95CI_Lower) & !combined_summary_dup2_5pl$Is_R1 & !combined_summary_dup2_5pl$Is_R2 & !combined_summary_dup2_5pl$Is_R3 & !combined_summary_dup2_5pl$Is_mAb,]
      if(nrow(ci_data_f5pl)>0){
        ci_data_f5pl$Sample_Num <- factor(1:nrow(ci_data_f5pl), levels = 1:nrow(ci_data_f5pl))
        f5pl_sub_colors <- c("5PL_Good_Fit"="#6A3D9A","5PL_5PL_Good_Fit"="#FF7F00",
                              "5PL_Constrained_Fit"="#CAB2D6","5PL_5PL_Constrained_Fit"="#FDBF6F",
                              "5PL_ED50_Fallback"="#4ECDC4","5PL_ED50_Interpolation"="#45B7D1",
                              "Fallback_4PL"="#808080","Failed"="#808080")
        p2c <- ggplot(ci_data_f5pl, aes(x=Sample_Num, y=IC50_log_dilution, color=IC50_Method)) +
          geom_point(size=point_size) +
          geom_errorbar(aes(ymin=IC50_95CI_Lower, ymax=IC50_95CI_Upper), width=errorbar_width) +
          scale_color_manual(values=f5pl_sub_colors, name="5PL Sub-method") +
          labs(title="F5PL", x="Sample (numbered)", y="IC50 log10(dilution)") +
          theme_publication_custom(base_font_size) +
          theme(axis.text.x=fig2_x_axis_text,
                plot.margin = margin(t=20, r=20, b=40, l=20),
                legend.position="bottom")
        ci_panels[[3]] <- p2c
      }
      
      ci_data_bayesian <- combined_summary_dup2_bayesian[!is.na(combined_summary_dup2_bayesian$IC50_95CI_Lower) & !combined_summary_dup2_bayesian$Is_R1 & !combined_summary_dup2_bayesian$Is_R2 & !combined_summary_dup2_bayesian$Is_R3 & !combined_summary_dup2_bayesian$Is_mAb,]
      if(nrow(ci_data_bayesian)>0){
        ci_data_bayesian$Sample_Num <- factor(1:nrow(ci_data_bayesian), levels = 1:nrow(ci_data_bayesian))
        bayesian_sub_colors <- c("Bayes_4PL_Good_Fit"="#E41A1C","Bayes_4PL_Constrained"="#377EB8",
                                 "Bayes_ED50_Fallback"="#4ECDC4","Bayes_4PL_ED50_Interpolation"="#45B7D1",
                                 "Fallback_4PL"="#808080","Failed"="#808080")
        p2d <- ggplot(ci_data_bayesian, aes(x=Sample_Num, y=IC50_log_dilution, color=IC50_Method)) +
          geom_point(size=point_size) +
          geom_errorbar(aes(ymin=IC50_95CI_Lower, ymax=IC50_95CI_Upper), width=errorbar_width) +
          scale_color_manual(values=bayesian_sub_colors, name="Bayes-4PL Sub-method") +
          labs(title="Bayes-4PL", x="Sample (numbered)", y="IC50 log10(dilution)") +
          theme_publication_custom(base_font_size) +
          theme(axis.text.x=fig2_x_axis_text,
                plot.margin = margin(t=20, r=20, b=40, l=20),
                legend.position="bottom")
        ci_panels[[4]] <- p2d
      }
      
      if(length(ci_panels) > 0) {
        if(length(ci_panels) == 1) {
          ggsave(file.path(plot_dir, "Fig2_IC50_with_CI.png"), ci_panels[[1]], width=fig2_width, height=fig2_height/3, dpi=dpi)
        } else if(length(ci_panels) == 2) {
          ggsave(file.path(plot_dir, "Fig2_IC50_with_CI.png"), arrangeGrob(ci_panels[[1]], ci_panels[[2]], nrow=2), width=fig2_width, height=fig2_height*2/3, dpi=dpi)
        } else if(length(ci_panels) == 3) {
          ggsave(file.path(plot_dir, "Fig2_IC50_with_CI.png"), arrangeGrob(ci_panels[[1]], ci_panels[[2]], ci_panels[[3]], nrow=3), width=fig2_width, height=fig2_height, dpi=dpi)
        } else {
          ggsave(file.path(plot_dir, "Fig2_IC50_with_CI.png"), arrangeGrob(ci_panels[[1]], ci_panels[[2]], ci_panels[[3]], ci_panels[[4]], nrow=4), width=fig2_width, height=fig2_height*4/3, dpi=dpi)
        }
      }
      
      ic50_rm_panels <- list()
      if(use_reed_muench_comparison) {
        rm_comp <- combined_summary[!is.na(combined_summary$RM_IC50_log) & !combined_summary$Is_R1 & !combined_summary$Is_R2 & !combined_summary$Is_R3 & !combined_summary$Is_mAb, ]
        if(nrow(rm_comp) > 0) {
          p5 <- ggplot(rm_comp, aes(x=RM_IC50_log, y=IC50_log_dilution, color=IC50_Method)) +
            geom_point(size=point_size, alpha=0.8) +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
            scale_color_manual(values=method_colors) +
            labs(title="Figure 5: 4PL vs Reed-Muench IC50 (Reference Only)", 
                 x="Reed-Muench IC50_log", y="4PL IC50_log") +
            theme_publication_custom(base_font_size)
          ic50_rm_panels[["Fig5"]] <- p5
        }
      }
      
      if(has_rlu3 && !is.null(comparison_df) && nrow(comparison_df)>0) {
        p6a <- ggplot(comparison_df, aes(x=IC50_log_dilution_2rep, y=IC50_log_dilution_3rep)) +
          geom_point(size=point_size, alpha=0.8, color=config_colors["two_replicates"]) +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          add_regression() +
          labs(title="Figure 6a: IC50 – Two vs Three Replicates", 
               x="IC50_log (two replicates)", y="IC50_log (three replicates)") +
          theme_publication_custom(base_font_size)
        ggsave(file.path(plot_dir,"Fig6a_RLU_IC50_Correlation.png"), p6a, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        
        p6b <- ggplot(comparison_df, aes(x=Titer_R1_2rep, y=Titer_R1_3rep)) +
          geom_point(size=point_size, alpha=0.8, color=config_colors["three_replicates"]) +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          add_regression() +
          labs(title="Figure 6b: Titer – Two vs Three Replicates", 
               x=paste0("Titer (two replicates, ", standard_unit, ")"), 
               y=paste0("Titer (three replicates, ", standard_unit, ")")) +
          theme_publication_custom(base_font_size)
        ggsave(file.path(plot_dir,"Fig6b_RLU_Titer_Correlation.png"), p6b, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        
        if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(comparison_df)) {
          agree_subset <- comparison_df[grepl(paste0("^", sample_prefix_filter), comparison_df$Sample, ignore.case = FALSE), ]
        } else {
          agree_subset <- comparison_df
        }
        agree_df <- data.frame(
          Category=c("Concordance\nAgreement","Protection\nAgreement","Method\nAgreement"),
          Count=c(sum(agree_subset$Concordance_Agreement_R1, na.rm=TRUE),
                  sum(agree_subset$Protection_Agreement, na.rm=TRUE),
                  sum(agree_subset$Method_Agreement, na.rm=TRUE)),
          Total=nrow(agree_subset))
        agree_df$Pct <- round(agree_df$Count/agree_df$Total*100, 1)
        p6c <- ggplot(agree_df, aes(x=Category, y=Pct)) +
          geom_bar(stat="identity", fill="#2E86AB", alpha=0.85) +
          geom_text(aes(label=paste0(Count,"/",Total,"\n(",Pct,"%)")), vjust=-0.3, size=3.5) +
          ylim(0,105) +
          labs(title="Figure 6c: Agreement between two and three replicates", 
               x="", y="Agreement (%)") +
          theme_publication_custom(base_font_size)
        ggsave(file.path(plot_dir,"Fig6c_RLU_Agreement.png"), p6c, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        
        if(has_rlu3 && !is.null(comparison_df) && nrow(comparison_df)>0 && exists("combined_summary_dup3") && !is.null(combined_summary_dup3)) {
          p6d_data <- merge(combined_summary_dup2[,c("Sample","IC50_log_dilution","IC50_Method")], 
                            combined_summary_dup3[,c("Sample","IC50_log_dilution","IC50_Method")], 
                            by="Sample", suffixes=c("_2rep","_3rep"))
          p6d_data <- p6d_data[!is.na(p6d_data$IC50_log_dilution_2rep) & !is.na(p6d_data$IC50_log_dilution_3rep), ]
          if(nrow(p6d_data) > 1) {
            p6d <- ggplot(p6d_data, aes(x=IC50_log_dilution_2rep, y=IC50_log_dilution_3rep, color=IC50_Method_2rep)) +
              geom_point(size=point_size, alpha=0.8) +
              geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
              geom_smooth(method="lm", se=TRUE, color="black", fill="grey70", alpha=0.3, linewidth=0.8) +
              scale_color_manual(values=method_colors, name="Method (2 rep)") +
              labs(title="Figure 6d: IC50 per Method – Two vs Three Replicates", 
                   4)),
                   x="IC50_log (two replicates)", y="IC50_log (three replicates)") +
              theme_publication_custom(base_font_size)
            ggsave(file.path(plot_dir,"Fig6d_IC50_Correlation_Per_Method.png"), p6d, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
          
          p6e_data <- merge(combined_summary_dup2[,c("Sample","Titer_R1","IC50_Method")], 
                            combined_summary_dup3[,c("Sample","Titer_R1","IC50_Method")], 
                            by="Sample", suffixes=c("_2rep","_3rep"))
          p6e_data <- p6e_data[!is.na(p6e_data$Titer_R1_2rep) & !is.na(p6e_data$Titer_R1_3rep), ]
          if(nrow(p6e_data) > 1) {
            p6e <- ggplot(p6e_data, aes(x=Titer_R1_2rep, y=Titer_R1_3rep, color=IC50_Method_2rep)) +
              geom_point(size=point_size, alpha=0.8) +
              geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
              geom_smooth(method="lm", se=TRUE, color="black", fill="grey70", alpha=0.3, linewidth=0.8) +
              scale_color_manual(values=method_colors, name="Method (2 rep)") +
              labs(title="Figure 6e: Titer per Method – Two vs Three Replicates", 
                   4)),
                   x=paste0("Titer (two replicates, ", standard_unit, ")"), 
                   y=paste0("Titer (three replicates, ", standard_unit, ")")) +
              theme_publication_custom(base_font_size)
            ggsave(file.path(plot_dir,"Fig6e_Titer_Correlation_Per_Method.png"), p6e, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
          
          if(nrow(p6e_data) > 1) {
            p6e_data$mean_titer <- rowMeans(p6e_data[, c("Titer_R1_2rep", "Titer_R1_3rep")], na.rm=TRUE)
            p6e_data$diff_titer <- p6e_data$Titer_R1_2rep - p6e_data$Titer_R1_3rep
            p6f <- ggplot(p6e_data, aes(x=mean_titer, y=diff_titer)) +
              geom_point(size=point_size, alpha=0.8, color="#2E86AB") +
              geom_hline(yintercept=mean(p6e_data$diff_titer, na.rm=TRUE), linetype="dashed", color="blue") +
              geom_hline(yintercept=mean(p6e_data$diff_titer, na.rm=TRUE) + 1.96*sd(p6e_data$diff_titer, na.rm=TRUE), linetype="dotted", color="red") +
              geom_hline(yintercept=mean(p6e_data$diff_titer, na.rm=TRUE) - 1.96*sd(p6e_data$diff_titer, na.rm=TRUE), linetype="dotted", color="red") +
              labs(
                   4), 
                                   ", LoA: [", round(mean(p6e_data$diff_titer, na.rm=TRUE) - 1.96*sd(p6e_data$diff_titer, na.rm=TRUE), 4), ", ",
                                   round(mean(p6e_data$diff_titer, na.rm=TRUE) + 1.96*sd(p6e_data$diff_titer, na.rm=TRUE), 4), "]"),
                   x=paste0("Mean Titer ", R1_name, " (IU/mL)"), y="Difference (2 rep - 3 rep)") +
              theme_publication_custom(base_font_size)
            ggsave(file.path(plot_dir,"Fig6f_BlandAltman_Titer_2vs3.png"), p6f, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
        }
      }
      
      create_elisa_titer_plot <- function(data, elisa_col, titer_col, ref_conc, ref_name, elisa_label, titer_label, fig_num, method_label="ACM", threshold=NULL, conc_col_override=NULL, data_3rep=NULL) {
        if(is.null(threshold)) threshold <- ref_conc
        
        make_panel <- function(df, rep_suffix) {
          if(is.null(df) || nrow(df)==0) return(NULL)
          plot_data <- df[((!is.na(df[[elisa_col]])) | (!is.null(et_elisa_linear_map) & !is.na(et_elisa_linear_map[.norm_key(df$Sample)])) | (!is.null(et_elisa_nlinear_map) & !is.na(et_elisa_nlinear_map[.norm_key(df$Sample)]))) & !is.na(df[[titer_col]]) & !df$Is_R1 & !df$Is_R2 & !df$Is_R3 & !df$Is_mAb, ]
          if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(plot_data)) {
            prefix_pattern <- paste0("^", sample_prefix_filter)
            plot_data <- plot_data[grepl(prefix_pattern, plot_data$Sample, ignore.case = FALSE), ]
          }
          if(nrow(plot_data)==0) return(NULL)
          
          .elisa_type <- if(elisa_col == "ELISA_EU_mL_nonlin") "nonlin" else "linear"
          if(!is.null(et_elisa_linear_map) || !is.null(et_elisa_nlinear_map)) {
            .sub <- mapply(function(sid, fit) .et_elisa_value(sid, .elisa_type, fit),
                           plot_data$Sample, plot_data[[elisa_col]])
            plot_data[[elisa_col]] <- as.numeric(.sub)
          }
          
          elisa_trans <- transform_elisa(plot_data[[elisa_col]], elisa_scale)
          titer_trans <- transform_titer(plot_data[[titer_col]], titer_scale)
          
          conc_col <- if(!is.null(conc_col_override)) {
            conc_col_override
          } else if(elisa_col == "ELISA_EU_mL") {
            if(ref_name == R1_name) "Concordance_R1_ELISA" else "Concordance_R2_ELISA"
          } else if(elisa_col == "ELISA_EU_mL_nonlin") {
            if(ref_name == R1_name) "Concordance_R1_Nonlin" else "Concordance_R2_Nonlin"
          } else { NULL }
          if(is.null(conc_col) || !conc_col %in% colnames(plot_data)) return(NULL)
          
          plot_data$Concordance <- plot_data[[conc_col]]
          if(exists(".elisa_type")) {
            plot_data$.titer_v <- suppressWarnings(as.numeric(plot_data[[titer_col]]))
            plot_data$.elisa_v <- suppressWarnings(as.numeric(plot_data[[elisa_col]]))
            plot_data$Concordance <- ifelse(!is.na(plot_data$.titer_v) & !is.na(plot_data$.elisa_v),
                                            check_concordance(plot_data$.titer_v, ref_conc, plot_data$.elisa_v,
                                                              log_tol = ed50_tolerance),
                                            plot_data$Concordance)
          }
          plot_data <- plot_data[plot_data$Concordance != "Unknown" & plot_data$Concordance != "N/A", ]
          if(nrow(plot_data)==0) return(NULL)
          
          elisa_thresh <- transform_threshold(threshold, elisa_scale)
          titer_thresh <- transform_titer_threshold(ref_conc, titer_scale)
          
          lm_fit <- lm(titer_trans$values ~ elisa_trans$values, data=plot_data)
          r_sq <- round(summary(lm_fit)$r.squared, 3)
          p_val <- if(!is.null(summary(lm_fit)$coefficients[2,4])) {
            round(summary(lm_fit)$coefficients[2,4], 4)
          } else { NA }
          
          p <- ggplot(plot_data, aes(x=elisa_trans$values, y=titer_trans$values, color=Concordance)) +
            geom_point(size=point_size, alpha=0.7) +
            geom_vline(xintercept=elisa_thresh, linetype="dashed", color="grey40", size=0.5) +
            geom_hline(yintercept=titer_thresh, linetype="dashed", color="grey40", size=0.5) +
            scale_color_manual(values=conc_colors) +
            add_regression() +
            labs(x=elisa_trans$label, y=titer_trans$label,
                 
                  "/", nrow(plot_data),
                                  " (", round(sum(plot_data$Concordance=="Concordant")/nrow(plot_data)*100,1), "%),  R² = ", r_sq,
                                  if(!is.na(p_val)) paste0(", p = ", p_val) else "")) +
             theme_publication_custom(base_font_size) +
            theme(legend.position="bottom")
          
          log_tol <- ed50_tolerance
          x_r <- range(plot_data[[elisa_col]], na.rm = TRUE)
          if(length(x_r) == 2 && all(x_r > 0) && !any(is.infinite(x_r))) {
            x_seq <- seq(x_r[1], x_r[2], length.out = 100)
            band_df <- data.frame(x = x_seq, y_upper = log10(x_seq) + log_tol, y_lower = log10(x_seq) - log_tol)
            y_r <- range(c(elisa_trans$values, titer_trans$values), na.rm = TRUE)
            if(show_concordance_rule) {
              p <- p + 
                geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
                geom_ribbon(data = band_df, aes(x = x, ymin = y_lower, ymax = y_upper), 
                            fill = "grey20", alpha = 0.15, inherit.aes = FALSE) +
                annotate("text", x = x_r[2], y = log10(x_r[2]) + log_tol * 0.5, 
                         label = paste0("+", log_tol, " log tolerance"), 
                         hjust = 1, vjust = -0.5, size = 3, color = "grey40") +
                annotate("text", x = x_r[2], y = log10(x_r[2]) - log_tol * 0.5, 
                         label = paste0("-", log_tol, " log tolerance"), 
                         hjust = 1, vjust = 1.5, size = 3, color = "grey40") +
                annotate("text", x = x_r[1], y = y_r[2], 
                         label = paste0("Concordance rule: |log10(Titer) - log10(ELISA)| < ", log_tol, 
                                        " OR both on same side of threshold (", threshold, " ", standard_unit, ")"),
                         hjust = 0, vjust = 1, size = 3, color = "black", fontface = "bold")
            }
          }
          return(p)
        }
        
        p_left <- make_panel(data, " (2 replicates)")
        if(is.null(p_left)) return(NULL)
        
        if(!is.null(data_3rep)) {
          p_right <- make_panel(data_3rep, " (3 replicates)")
          if(is.null(p_right)) return(p_left)
          arrangeGrob(p_left, p_right, ncol=2)
        } else {
          p_left
        }
      }
      
      create_delta_plot <- function(data, elisa_col, titer_col, ref_conc, ref_name, elisa_label, titer_label, fig_num, method_label="ACM", threshold=NULL, conc_col_override=NULL, data_3rep=NULL) {
        if(is.null(threshold)) threshold <- ref_conc
        
        make_panel <- function(df, rep_suffix) {
          if(is.null(df) || nrow(df)==0) return(NULL)
          plot_data <- df[((!is.na(df[[elisa_col]])) | (!is.null(et_elisa_linear_map) & !is.na(et_elisa_linear_map[.norm_key(df$Sample)])) | (!is.null(et_elisa_nlinear_map) & !is.na(et_elisa_nlinear_map[.norm_key(df$Sample)]))) & !is.na(df[[titer_col]]) & !df$Is_R1 & !df$Is_R2 & !df$Is_R3 & !df$Is_mAb, ]
          if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(plot_data)) {
            prefix_pattern <- paste0("^", sample_prefix_filter)
            plot_data <- plot_data[grepl(prefix_pattern, plot_data$Sample, ignore.case = FALSE), ]
          }
          if(nrow(plot_data)==0) return(NULL)
          
          .elisa_type <- if(elisa_col == "ELISA_EU_mL_nonlin") "nonlin" else "linear"
          if(!is.null(et_elisa_linear_map) || !is.null(et_elisa_nlinear_map)) {
            .sub <- mapply(function(sid, fit) .et_elisa_value(sid, .elisa_type, fit),
                           plot_data$Sample, plot_data[[elisa_col]])
            plot_data[[elisa_col]] <- as.numeric(.sub)
          }
          
          plot_data$log_elisa <- log10(plot_data[[elisa_col]])
          plot_data$log_titer <- log10(plot_data[[titer_col]])
          plot_data$delta <- plot_data$log_titer - plot_data$log_elisa
          if(diff_pct_mode) plot_data$delta <- (10^plot_data$delta - 1) * 100
          
          conc_col <- if(!is.null(conc_col_override)) {
            conc_col_override
          } else if(elisa_col == "ELISA_EU_mL") {
            if(ref_name == R1_name) "Concordance_R1_ELISA" else "Concordance_R2_ELISA"
          } else if(elisa_col == "ELISA_EU_mL_nonlin") {
            if(ref_name == R1_name) "Concordance_R1_Nonlin" else "Concordance_R2_Nonlin"
          } else { NULL }
          if(is.null(conc_col) || !conc_col %in% colnames(plot_data)) return(NULL)
          
          plot_data$Concordance <- plot_data[[conc_col]]
          if(exists(".elisa_type")) {
            plot_data$.titer_v <- suppressWarnings(as.numeric(plot_data[[titer_col]]))
            plot_data$.elisa_v <- suppressWarnings(as.numeric(plot_data[[elisa_col]]))
            plot_data$Concordance <- ifelse(!is.na(plot_data$.titer_v) & !is.na(plot_data$.elisa_v),
                                            check_concordance(plot_data$.titer_v, ref_conc, plot_data$.elisa_v,
                                                              log_tol = ed50_tolerance),
                                            plot_data$Concordance)
          }
          plot_data <- plot_data[plot_data$Concordance != "Unknown" & plot_data$Concordance != "N/A", ]
          if(nrow(plot_data)==0) return(NULL)
          
          mean_diff <- mean(plot_data$delta, na.rm = TRUE)
          sd_diff <- sd(plot_data$delta, na.rm = TRUE)
          if(!is.finite(mean_diff) || !is.finite(sd_diff)) return(NULL)
          loa_upper <- mean_diff + 1.96 * sd_diff
          loa_lower <- mean_diff - 1.96 * sd_diff
          
          p <- ggplot(plot_data, aes(x = log_elisa, y = delta, color = Concordance)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
            geom_hline(yintercept = conc_tol_upper, linetype = "dotted", color = "grey40", linewidth = 0.6) +
            geom_hline(yintercept = conc_tol_lower, linetype = "dotted", color = "grey40", linewidth = 0.6) +
            geom_hline(yintercept = mean_diff, linetype = "solid", color = "blue", linewidth = 0.8) +
            geom_hline(yintercept = loa_upper, linetype = "dashed", color = "red", linewidth = 0.6) +
            geom_hline(yintercept = loa_lower, linetype = "dashed", color = "red", linewidth = 0.6) +
            geom_point(size = point_size, alpha = 0.7) +
            scale_color_manual(values = conc_colors) +
            labs(x = paste0("log10(", elisa_label, ")"), 
                 y = if(diff_pct_mode) "Percent Titer Difference vs ELISA (%)" else paste0("log10(Titer) - log10(", elisa_label, ")"),
                 
                 " | SD: ", round(sd_diff, 3), 
                                   " | LoA: [", round(loa_lower, 3), ", ", round(loa_upper, 3), "]")) +
            theme_publication_custom(base_font_size) +
            theme(legend.position = "bottom")
          return(p)
        }
        
        p_left <- make_panel(data, " (2 rep)")
        if(is.null(p_left)) return(NULL)
        
        if(!is.null(data_3rep)) {
          p_right <- make_panel(data_3rep, " (3 rep)")
          if(is.null(p_right)) return(p_left)
          arrangeGrob(p_left, p_right, ncol=2)
        } else {
          p_left
        }
      }
      
      et_all_methods <- c("ACM", "F4PL", "F5PL", "Bayes-4PL")
      et_rm_ok <- isTRUE(use_reed_muench_comparison) &&
        exists("combined_summary_dup2") && is.data.frame(combined_summary_dup2) &&
        any(!is.na(suppressWarnings(as.numeric(combined_summary_dup2$RM_IC50_log))))
      if(et_rm_ok) et_all_methods <- c(et_all_methods, "Reed-Muench")
      
       et_get_frame <- function(method, rep_tag) {
         is_rm <- identical(method, "Reed-Muench")
         if(rep_tag == "3rep") {
           df <- switch(method,
                        "ACM" = if(exists("combined_summary_dup3")) combined_summary_dup3 else NULL,
                        "F4PL"    = if(exists("combined_summary_dup3_forced4pl")) combined_summary_dup3_forced4pl else NULL,
                        "F5PL"     = if(exists("combined_summary_dup3_5pl")) combined_summary_dup3_5pl else NULL,
                        "Bayes-4PL"      = if(exists("combined_summary_dup3_bayesian")) combined_summary_dup3_bayesian else NULL,
                        "Reed-Muench"   = if(exists("combined_summary_dup3")) combined_summary_dup3 else NULL)
         } else {
           df <- switch(method,
                        "ACM" = combined_summary_dup2,
                        "F4PL"    = combined_summary_dup2_forced4pl,
                        "F5PL"     = combined_summary_dup2_5pl,
                        "Bayes-4PL"      = combined_summary_dup2_bayesian,
                        "Reed-Muench"   = combined_summary_dup2)
         }
        if(is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
        df <- df
        if(is_rm) {
          req <- c("RM_IC50_log", "R1_IC50_dilution", "R2_IC50_dilution")
          if(!all(req %in% colnames(df))) return(NULL)
          rm_log   <- suppressWarnings(as.numeric(df$RM_IC50_log))
          R1_dil <- suppressWarnings(as.numeric(df$R1_IC50_dilution))
          R2_dil <- suppressWarnings(as.numeric(df$R2_IC50_dilution))
          df$RM_Titer_R1 <- ifelse(!is.na(rm_log) & !is.na(R1_dil) & R1_dil > 0,
                                     (10^rm_log / R1_dil) * reference_concentration, NA_real_)
          df$RM_Titer_R2 <- ifelse(!is.na(rm_log) & !is.na(R2_dil) & R2_dil > 0,
                                     (10^rm_log / R2_dil) * second_reference_concentration, NA_real_)
          if("ELISA_EU_mL" %in% colnames(df))
             df$Concordance_R1_ELISA <- check_concordance(df$RM_Titer_R1, protection_threshold, df$ELISA_EU_mL, log_tol = ed50_tolerance, use_threshold_only = FALSE)
          if("ELISA_EU_mL_nonlin" %in% colnames(df))
             df$Concordance_R1_Nonlin <- check_concordance(df$RM_Titer_R1, protection_threshold, df$ELISA_EU_mL_nonlin, log_tol = ed50_tolerance, use_threshold_only = FALSE)
          if(isTRUE(use_dual_reference)) {
            if("ELISA_EU_mL" %in% colnames(df))
               df$Concordance_R2_ELISA <- check_concordance(df$RM_Titer_R2, protection_threshold, df$ELISA_EU_mL, log_tol = ed50_tolerance, use_threshold_only = FALSE)
            if("ELISA_EU_mL_nonlin" %in% colnames(df))
               df$Concordance_R2_Nonlin <- check_concordance(df$RM_Titer_R2, protection_threshold, df$ELISA_EU_mL_nonlin, log_tol = ed50_tolerance, use_threshold_only = FALSE)
          }
        }
        df
      }
      
      et_save_if <- function(p, fname_base) {
        if(is.null(p)) return(FALSE)
        dims <- get_elisa_titer_dims(p)
        ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_", fname_base, ".png")),
               p, width = dims$width, height = dims$height, dpi = dpi)
        fig_counter <<- fig_counter + 1
        TRUE
      }
      
       .norm_key <- function(v) gsub("[^A-Za-z0-9._-]+", "", tolower(as.character(v)))
       find_col <- function(pattern, df) {
         cn <- colnames(df)
         exact <- cn[cn == pattern]
         if(length(exact) > 0) return(exact[1])
         partial <- cn[grepl(pattern, cn, fixed = TRUE)]
         if(length(partial) > 0) return(partial[1])
         return(NULL)
       }
       et_elisa_linear_map  <- NULL
       et_elisa_nlinear_map <- NULL
       if(file.exists(file_path)) {
         tryCatch({
           .s4 <- as.data.frame(readxl::read_xlsx(file_path, sheet = group_concordance_sheet))
           .s4_key <- if(any(tolower(colnames(.s4)) %in% c("sample", "sample_id"))) {
             .s4[[colnames(.s4)[tolower(colnames(.s4)) %in% c("sample", "sample_id")][1]]]
           } else NULL
            if(!is.null(.s4_key)) {
              .sk <- .norm_key(.s4_key)
              .lin_col <- find_col("linear", .s4)
              .nlin_col <- find_col("nlinear", .s4)
              if(!is.null(.lin_col)) {
               .x <- setNames(suppressWarnings(as.numeric(.s4[[.lin_col]])), .sk)
               .x <- .x[!is.na(names(.x)) & names(.x) != ""]
               if(length(.x) > 0) et_elisa_linear_map <- .x
             }
             if(!is.null(.nlin_col)) {
               .x <- setNames(suppressWarnings(as.numeric(.s4[[.nlin_col]])), .sk)
               .x <- .x[!is.na(names(.x)) & names(.x) != ""]
               if(length(.x) > 0) et_elisa_nlinear_map <- .x
             }
                       "; non-linear:", length(et_elisa_nlinear_map))
           }
        }, error = function(e) )
      }
      .et_elisa_value <- function(sample_id, elisa_type, fitted_val) {
        .m <- if(identical(elisa_type, "nonlin")) et_elisa_nlinear_map else et_elisa_linear_map
        if(!is.null(.m) && length(.m) > 0) {
          .v <- .m[[.norm_key(sample_id)]]
          if(!is.na(.v)) return(.v)
        }
        NA_real_
      }
      
      fig_counter <- 7
      for(et_m in et_all_methods) {
        et_2rep <- et_get_frame(et_m, "2rep")
        et_3rep <- et_get_frame(et_m, "3rep")
        if(is.null(et_2rep) || nrow(et_2rep) == 0) {
          next
        }
        et_safe <- gsub("[^A-Za-z0-9._-]", "_", et_m)
        
        et_save_if(create_elisa_titer_plot(et_2rep, "ELISA_EU_mL", "Titer_R1",
                                           reference_concentration, R1_name, "ELISA", paste0("Titer (", standard_unit, ")"),
                                           fig_counter, et_m, data_3rep = et_3rep),
                   paste0("ELISA_vs_Titer_", R1_file, "_", et_safe))
        if(delta_plot) et_save_if(create_delta_plot(et_2rep, "ELISA_EU_mL", "Titer_R1",
                                     reference_concentration, R1_name, "ELISA", paste0("Titer (", standard_unit, ")"),
                                     fig_counter, et_m, data_3rep = et_3rep),
                   paste0("Delta_ELISA_vs_Titer_", R1_file, "_", et_safe))
        
        if(has_elisa_nonlin) {
          et_save_if(create_elisa_titer_plot(et_2rep, "ELISA_EU_mL_nonlin", "Titer_R1",
                                             reference_concentration, R1_name, "ELISA (non-linear)", paste0("Titer (", standard_unit, ")"),
                                             fig_counter, et_m, data_3rep = et_3rep),
                     paste0("ELISA_nonlin_vs_Titer_", R1_file, "_", et_safe))
          if(delta_plot) et_save_if(create_delta_plot(et_2rep, "ELISA_EU_mL_nonlin", "Titer_R1",
                                       reference_concentration, R1_name, "ELISA (non-linear)", paste0("Titer (", standard_unit, ")"),
                                       fig_counter, et_m, data_3rep = et_3rep),
                     paste0("Delta_ELISA_nonlin_vs_Titer_", R1_file, "_", et_safe))
        }
        
        if(use_dual_reference) {
          et_titer_r2 <- if(identical(et_m, "Reed-Muench")) "RM_Titer_R2" else "Titer_R2"
          et_save_if(create_elisa_titer_plot(et_2rep, "ELISA_EU_mL", et_titer_r2,
                                             second_reference_concentration, R2_name, "ELISA", paste0("Titer (", second_unit, ")"),
                                             fig_counter, et_m, data_3rep = et_3rep),
                     paste0("ELISA_vs_Titer_", R2_file, "_", et_safe))
          if(delta_plot) et_save_if(create_delta_plot(et_2rep, "ELISA_EU_mL", et_titer_r2,
                                       second_reference_concentration, R2_name, "ELISA", paste0("Titer (", second_unit, ")"),
                                       fig_counter, et_m, data_3rep = et_3rep),
                     paste0("Delta_ELISA_vs_Titer_", R2_file, "_", et_safe))
          if(has_elisa_nonlin) {
            et_save_if(create_elisa_titer_plot(et_2rep, "ELISA_EU_mL_nonlin", et_titer_r2,
                                               second_reference_concentration, R2_name, "ELISA (non-linear)", paste0("Titer (", second_unit, ")"),
                                               fig_counter, et_m, data_3rep = et_3rep),
                       paste0("ELISA_nonlin_vs_Titer_", R2_file, "_", et_safe))
            if(delta_plot) et_save_if(create_delta_plot(et_2rep, "ELISA_EU_mL_nonlin", et_titer_r2,
                                         second_reference_concentration, R2_name, "ELISA (non-linear)", paste0("Titer (", second_unit, ")"),
                                         fig_counter, et_m, data_3rep = et_3rep),
                       paste0("Delta_ELISA_nonlin_vs_Titer_", R2_file, "_", et_safe))
          }
        }
      }
      
      if(exists("comparison_ACM_forced") && !is.null(comparison_ACM_forced) && nrow(comparison_ACM_forced)>0) {
        
        valid_plot <- !is.na(comparison_ACM_forced$IC50_log_dilution_ACM) & !is.na(comparison_ACM_forced$IC50_log_dilution_5PL)
        if(sum(valid_plot) > 1) {
          p11a <- ggplot(comparison_ACM_forced[valid_plot,], aes(x=IC50_log_dilution_5PL, y=IC50_log_dilution_ACM)) +
            geom_point(size=point_size, alpha=0.8, color="#D55E00") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
             labs(title="Figure 11a: IC50 - ACM vs 5PL",
                  4)),
                  x="IC50_log (F5PL)", y="IC50_log (ACM)") +
            theme_publication_custom(base_font_size)
          ggsave(file.path(plot_dir, "Fig11a_IC50_Comparison_ACM_vs_5PL.png"), p11a, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
        
        valid_plot2 <- !is.na(comparison_ACM_forced$IC50_log_dilution_Forced4PL) & !is.na(comparison_ACM_forced$IC50_log_dilution_5PL)
        if(sum(valid_plot2) > 1) {
          p11b <- ggplot(comparison_ACM_forced[valid_plot2,], aes(x=IC50_log_dilution_5PL, y=IC50_log_dilution_Forced4PL)) +
            geom_point(size=point_size, alpha=0.8, color="#F18F01") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
             labs(title="Figure 11b: IC50 - F4PL vs 5PL",
                  4)),
                  x="IC50_log (F5PL)", y="IC50_log (F4PL)") +
            theme_publication_custom(base_font_size)
          ggsave(file.path(plot_dir, "Fig11b_IC50_Comparison_Forced_vs_5PL.png"), p11b, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
        
        
        if("IC50_Diff_ACM_vs_5PL" %in% colnames(comparison_ACM_forced)) {
          diff_data <- data.frame(
            IC50_Diff = comparison_ACM_forced$IC50_Diff_ACM_vs_5PL
          )
          diff_data <- diff_data[!is.na(diff_data$IC50_Diff), ]
          
          if(is.data.frame(diff_data) && nrow(diff_data) > 1) {
            p11d <- ggplot(diff_data, aes(x=IC50_Diff)) +
              geom_histogram(bins=20, fill="#D55E00", alpha=0.7, color="black", size=0.3) +
              geom_vline(xintercept=0, linetype="dashed", color="red") +
              labs(title="Figure 11d: IC50 Difference Distribution",
                   4), 
                                   ", SD: ", round(sd(diff_data$IC50_Diff), 4)),
                    x="IC50_log (ACM - 5PL)", y="Count") +
              theme_publication_custom(base_font_size)
            ggsave(file.path(plot_dir, "Fig11d_IC50_Difference_ACM_vs_5PL.png"), p11d, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
        }
        
        if("IC50_Diff_Forced_vs_5PL" %in% colnames(comparison_ACM_forced)) {
          diff_data2 <- data.frame(
            IC50_Diff = comparison_ACM_forced$IC50_Diff_Forced_vs_5PL
          )
          diff_data2 <- diff_data2[!is.na(diff_data2$IC50_Diff), ]
          
          if(is.data.frame(diff_data2) && nrow(diff_data2) > 1) {
            p11e <- ggplot(diff_data2, aes(x=IC50_Diff)) +
              geom_histogram(bins=20, fill="#F18F01", alpha=0.7, color="black", size=0.3) +
              geom_vline(xintercept=0, linetype="dashed", color="red") +
              labs(title="Figure 11e: IC50 Difference Distribution",
                   4), 
                                   ", SD: ", round(sd(diff_data2$IC50_Diff), 4)),
                    x="IC50_log (F4PL - 5PL)", y="Count") +
              theme_publication_custom(base_font_size)
            ggsave(file.path(plot_dir, "Fig11e_IC50_Difference_Forced_vs_5PL.png"), p11e, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
        }
        
      }
      
      plot_12a_data <- comparison_ACM_forced[!is.na(comparison_ACM_forced$Titer_R1_ACM) & !is.na(comparison_ACM_forced$Titer_R1_5PL), ]
      if(nrow(plot_12a_data) > 1) {
        p12a <- ggplot(plot_12a_data, aes(x=Titer_R1_5PL, y=Titer_R1_ACM)) +
          geom_point(size=point_size, alpha=0.8, color="#2E86AB") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          geom_smooth(method="lm", se=TRUE, color="black", fill="grey70", alpha=0.3, linewidth=0.8) +
          labs(3)),
               x=paste0("Titer ", R1_name, " F5PL (IU/mL)"), y=paste0("Titer ", R1_name, " ACM (IU/mL)")) +
          theme_publication_custom(base_font_size)
        ggsave(file.path(plot_dir, paste0("Fig12a_Titer_", R1_file, "_ACM_vs_5PL.png")), p12a, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
      }
      
      if(use_dual_reference) {
        plot_12b_data <- comparison_ACM_forced[!is.na(comparison_ACM_forced$Titer_R2_ACM) & !is.na(comparison_ACM_forced$Titer_R2_5PL), ]
        if(nrow(plot_12b_data) > 1) {
          p12b <- ggplot(plot_12b_data, aes(x=Titer_R2_5PL, y=Titer_R2_ACM)) +
            geom_point(size=point_size, alpha=0.8, color="#F18F01") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            geom_smooth(method="lm", se=TRUE, color="black", fill="grey70", alpha=0.3, linewidth=0.8) +
            labs(3)),
                 x=paste0("Titer ", R2_name, " F5PL (EU/mL)"), y=paste0("Titer ", R2_name, " ACM (EU/mL)")) +
            theme_publication_custom(base_font_size)
          ggsave(file.path(plot_dir, paste0("Fig12b_Titer_", R2_file, "_ACM_vs_5PL.png")), p12b, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
      }
      
      if(nrow(plot_12a_data) > 1) {
        plot_12a_data$mean_titer <- rowMeans(plot_12a_data[, c("Titer_R1_ACM", "Titer_R1_5PL")], na.rm=TRUE)
        plot_12a_data$diff_titer <- plot_12a_data$Titer_R1_ACM - plot_12a_data$Titer_R1_5PL
        p13 <- ggplot(plot_12a_data, aes(x=mean_titer, y=diff_titer)) +
          geom_point(size=point_size, alpha=0.8, color="#2E86AB") +
          geom_hline(yintercept=mean(plot_12a_data$diff_titer, na.rm=TRUE), linetype="dashed", color="blue") +
          geom_hline(yintercept=mean(plot_12a_data$diff_titer, na.rm=TRUE) + 1.96*sd(plot_12a_data$diff_titer, na.rm=TRUE), linetype="dotted", color="red") +
          geom_hline(yintercept=mean(plot_12a_data$diff_titer, na.rm=TRUE) - 1.96*sd(plot_12a_data$diff_titer, na.rm=TRUE), linetype="dotted", color="red") +
          labs(
               4), 
                               ", LoA: [", round(mean(plot_12a_data$diff_titer, na.rm=TRUE) - 1.96*sd(plot_12a_data$diff_titer, na.rm=TRUE), 4), ", ",
                               round(mean(plot_12a_data$diff_titer, na.rm=TRUE) + 1.96*sd(plot_12a_data$diff_titer, na.rm=TRUE), 4), "]"),
               x=paste0("Mean Titer ", R1_name, " (IU/mL)"), y="Difference (ACM - 5PL)") +
          theme_publication_custom(base_font_size)
        ggsave(file.path(plot_dir, paste0("Fig13_BlandAltman_ACM_vs_5PL_", R1_file, ".png")), p13, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
      }
      
      .method_concordance <- function(tbl, method, ref, elisatype) {
        if(is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(rep(NA_character_, 0))
        ref_conc <- if(ref == "R1") reference_concentration else second_reference_concentration
        elisa_col<- if(elisatype == "ELISA") "ELISA_EU_mL" else "ELISA_EU_mL_nonlin"
        titer_col<- if(ref == "R1") "Titer_R1" else "Titer_R2"
        dil_col  <- if(ref == "R1") "R1_IC50_dilution" else "R2_IC50_dilution"
        if(method == "Reed-Muench") {
          if(!use_reed_muench_comparison) return(rep(NA_character_, nrow(tbl)))
          if(!all(c("RM_IC50_log", dil_col, elisa_col) %in% colnames(tbl))) return(rep(NA_character_, nrow(tbl)))
          rm_log  <- suppressWarnings(as.numeric(tbl$RM_IC50_log))
          ref_dil <- suppressWarnings(as.numeric(tbl[[dil_col]]))
          elisa_v <- suppressWarnings(as.numeric(tbl[[elisa_col]]))
          titer_v <- ifelse(!is.na(rm_log) & !is.na(ref_dil) & ref_dil > 0,
                            (10^rm_log / ref_dil) * ref_conc, NA_real_)
        } else {
          if(!all(c(titer_col, elisa_col) %in% colnames(tbl))) return(rep(NA_character_, nrow(tbl)))
          titer_v <- suppressWarnings(as.numeric(tbl[[titer_col]]))
          elisa_v <- suppressWarnings(as.numeric(tbl[[elisa_col]]))
        }
        out <- ifelse(!is.na(titer_v) & !is.na(elisa_v),
                      check_concordance(titer_v, ref_conc, elisa_v, log_tol = ed50_tolerance),
                      NA_character_)
        if(all(c("Is_R1", "Is_R2", "Is_R3", "Is_mAb") %in% colnames(tbl))) {
          keep <- !(tbl$Is_R1 %in% TRUE) & !(tbl$Is_R2 %in% TRUE) & !(tbl$Is_R3 %in% TRUE) & !(tbl$Is_mAb %in% TRUE)
          out[!keep] <- NA
        }
        out
      }

      build_method_concordance_rates <- function(ref, elisatype, mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl) {
        col <- paste0("Concordance_", ref, "_", elisatype)
        elisa_col <- if(elisatype == "ELISA") "ELISA_EU_mL" else "ELISA_EU_mL_nonlin"
        
        has_elisa <- function(tbl) !is.null(tbl) && is.data.frame(tbl) && elisa_col %in% colnames(tbl)
        
        d_mixed    <- if(has_elisa(mixed_tbl))    .method_concordance(mixed_tbl, "ACM", ref, elisatype) else rep(NA_character_, if(!is.null(mixed_tbl) && is.data.frame(mixed_tbl)) nrow(mixed_tbl) else 0)
        d_forced   <- if(has_elisa(forced_tbl))   .method_concordance(forced_tbl, "F4PL",   ref, elisatype) else rep(NA_character_, if(!is.null(forced_tbl) && is.data.frame(forced_tbl)) nrow(forced_tbl) else 0)
        d_5pl      <- if(has_elisa(f5pl_tbl))      .method_concordance(f5pl_tbl, "F5PL",   ref, elisatype) else rep(NA_character_, if(!is.null(f5pl_tbl) && is.data.frame(f5pl_tbl)) nrow(f5pl_tbl) else 0)
        d_bayesian <- if(has_elisa(bayesian_tbl)) .method_concordance(bayesian_tbl, "Bayes-4PL", ref, elisatype) else rep(NA_character_, if(!is.null(bayesian_tbl) && is.data.frame(bayesian_tbl)) nrow(bayesian_tbl) else 0)
        d_rm       <- if(has_elisa(mixed_tbl))    .method_concordance(mixed_tbl, "Reed-Muench", ref, elisatype) else rep(NA_character_, if(!is.null(mixed_tbl) && is.data.frame(mixed_tbl)) nrow(mixed_tbl) else 0)

        .valid_status <- function(v) !is.na(v) & v %in% c("Concordant", "Discordant")
        methods <- c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")
        cols <- list(d_mixed, d_forced, d_5pl, d_bayesian, d_rm)
        conc <- sapply(cols, function(c) sum(c == "Concordant", na.rm = TRUE))
        tot  <- sapply(cols, function(c) sum(.valid_status(c)))
        rate <- ifelse(tot > 0, conc / tot * 100, NA)
        data.frame(Method = methods, Concordance_Count = conc, Total = tot, Rate = rate, stringsAsFactors = FALSE)
      }
      
      make_conc_rate_panel <- function(rate_df, title, subtitle) {
        if(is.null(rate_df) || nrow(rate_df) == 0) return(NULL)
        all_methods <- c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")
        rate_df <- merge(data.frame(Method = all_methods, stringsAsFactors = FALSE), rate_df, by = "Method", all.x = TRUE)
        rate_df$Rate[is.na(rate_df$Rate)] <- 0
        rate_df$Concordance_Count[is.na(rate_df$Concordance_Count)] <- 0
        rate_df$Total[is.na(rate_df$Total)] <- 0
        rate_df$Method <- factor(rate_df$Method, levels = all_methods)
        ggplot(rate_df, aes(x = Method, y = Rate, fill = Method)) +
          geom_bar(stat = "identity", alpha = 0.85, width = bar_width) +
          geom_text(aes(label = paste0(round(Rate, 1), "%\n(", Concordance_Count, "/", Total, ")")), vjust = -0.3, size = 4) +
          scale_fill_manual(values = method_compare_colors) +
          labs( x = "", y = "Concordance with ELISA (%)") +
          theme_publication_custom(base_font_size) +
          theme(legend.position = "none") +
          ylim(0, 105)
      }
      
      build_concordance_fourpanel <- function(mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl, file_suffix, rep_label) {
        conc_panels <- list()
        p_a <- make_conc_rate_panel(build_method_concordance_rates("R1", "ELISA", mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl),
                                    paste0("Concordance Rate by Method (", R1_name, ", Linear ELISA) - ", rep_label),
                                    paste0("Reference: ", R1_name, "; ELISA type: Linear; ", rep_label))
        if(!is.null(p_a)) conc_panels[["R1_Linear"]] <- p_a
        
        if(has_elisa_nonlin) {
          p_b <- make_conc_rate_panel(build_method_concordance_rates("R1", "Nonlin", mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl),
                                      paste0("Concordance Rate by Method (", R1_name, ", Non-linear ELISA) - ", rep_label),
                                      paste0("Reference: ", R1_name, "; ELISA type: Non-linear; ", rep_label))
          if(!is.null(p_b)) conc_panels[["R1_Nonlin"]] <- p_b
        }
        if(use_dual_reference) {
          p_c <- make_conc_rate_panel(build_method_concordance_rates("R2", "ELISA", mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl),
                                      paste0("Concordance Rate by Method (", R2_name, ", Linear ELISA) - ", rep_label),
                                      paste0("Reference: ", R2_name, "; ELISA type: Linear; ", rep_label))
          if(!is.null(p_c)) conc_panels[["R2_Linear"]] <- p_c
          if(has_elisa_nonlin) {
            p_d <- make_conc_rate_panel(build_method_concordance_rates("R2", "Nonlin", mixed_tbl, forced_tbl, f5pl_tbl, bayesian_tbl),
                                        paste0("Concordance Rate by Method (", R2_name, ", Non-linear ELISA) - ", rep_label),
                                        paste0("Reference: ", R2_name, "; ELISA type: Non-linear; ", rep_label))
            if(!is.null(p_d)) conc_panels[["R2_Nonlin"]] <- p_d
          }
        }
        
        if(length(conc_panels) > 0) {
          np <- length(conc_panels)
          out_file <- file.path(plot_dir, paste0("Fig_Concordance_FourPanel", file_suffix, ".png"))
          if(np == 1) {
            ggsave(out_file, conc_panels[[1]], width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
          } else if(np == 2) {
            ggsave(out_file, arrangeGrob(conc_panels[[1]], conc_panels[[2]], ncol = 2),
                   width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
          } else if(np == 3) {
            ggsave(out_file, arrangeGrob(conc_panels[[1]], conc_panels[[2]], conc_panels[[3]], ncol = 3),
                   width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
          } else {
            ggsave(out_file, arrangeGrob(conc_panels[[1]], conc_panels[[2]], conc_panels[[3]], conc_panels[[4]],
                                         nrow = 2, ncol = 2),
                   width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
          }
        }
      }
      
      build_concordance_fourpanel(combined_summary_dup2, combined_summary_dup2_forced4pl, combined_summary_dup2_5pl, combined_summary_dup2_bayesian,
                                  "", "2 replicates")
      
      if(has_rlu3 &&
         exists("combined_summary_dup3") && !is.null(combined_summary_dup3) &&
         exists("combined_summary_dup3_forced4pl") && !is.null(combined_summary_dup3_forced4pl) &&
         exists("combined_summary_dup3_5pl") && !is.null(combined_summary_dup3_5pl) &&
         exists("combined_summary_dup3_bayesian") && !is.null(combined_summary_dup3_bayesian)) {
        build_concordance_fourpanel(combined_summary_dup3, combined_summary_dup3_forced4pl, combined_summary_dup3_5pl, combined_summary_dup3_bayesian,
                                    "_3rep", "3 replicates")
      }
      
      if(nrow(combined_summary_dup2) > 0 && "Experiment_Group" %in% colnames(combined_summary_dup2)) {
        group_list <- unique(combined_summary_dup2$Experiment_Group)
        group_list <- group_list[!is.na(group_list) & group_list != ""]
        for(grp in group_list) {
          g_mixed <- combined_summary_dup2[combined_summary_dup2$Experiment_Group == grp, , drop = FALSE]
          g_forced <- combined_summary_dup2_forced4pl[combined_summary_dup2_forced4pl$Experiment_Group == grp, , drop = FALSE]
           g_sk <- combined_summary_dup2_5pl[combined_summary_dup2_5pl$Experiment_Group == grp, , drop = FALSE]
           g_bayesian <- combined_summary_dup2_bayesian[combined_summary_dup2_bayesian$Experiment_Group == grp, , drop = FALSE]
           if(nrow(g_mixed) >= 2) {
             build_concordance_fourpanel(g_mixed, g_forced, g_sk, g_bayesian,
                                        paste0("_", gsub("[^A-Za-z0-9]", "_", grp)),
                                        paste0(grp, " (2 rep)"))
          }
        }
      }
      
      if(has_rlu3 &&
         exists("combined_summary_dup3") && !is.null(combined_summary_dup3) &&
         exists("combined_summary_dup3_forced4pl") && !is.null(combined_summary_dup3_forced4pl) &&
         exists("combined_summary_dup3_5pl") && !is.null(combined_summary_dup3_5pl) &&
         exists("combined_summary_dup3_bayesian") && !is.null(combined_summary_dup3_bayesian)) {
        if(nrow(combined_summary_dup3) > 0 && "Experiment_Group" %in% colnames(combined_summary_dup3)) {
          group_list <- unique(combined_summary_dup3$Experiment_Group)
          group_list <- group_list[!is.na(group_list) & group_list != ""]
          for(grp in group_list) {
            g_mixed <- combined_summary_dup3[combined_summary_dup3$Experiment_Group == grp, , drop = FALSE]
            g_forced <- combined_summary_dup3_forced4pl[combined_summary_dup3_forced4pl$Experiment_Group == grp, , drop = FALSE]
            g_sk <- combined_summary_dup3_5pl[combined_summary_dup3_5pl$Experiment_Group == grp, , drop = FALSE]
            g_bayesian <- combined_summary_dup3_bayesian[combined_summary_dup3_bayesian$Experiment_Group == grp, , drop = FALSE]
            if(nrow(g_mixed) >= 2) {
              build_concordance_fourpanel(g_mixed, g_forced, g_sk, g_bayesian,
                                          paste0("_3rep_", gsub("[^A-Za-z0-9]", "_", grp)),
                                          paste0(grp, " (3 rep)"))
            }
          }
        }
      }
      
      rm_plot_data <- combined_summary_dup2
      rm_plot_data$Method <- "Reed-Muench"
      mixed_plot <- combined_summary_dup2
      mixed_plot$Method <- "ACM"
      forced_plot <- combined_summary_dup2_forced4pl
      forced_plot$Method <- "F4PL"
      best_plot <- combined_summary_dup2_5pl
      best_plot$Method <- "F5PL"
      bayesian_plot <- combined_summary_dup2_bayesian
      bayesian_plot$Method <- "Bayes-4PL"
      if(has_elisa_nonlin) {
        rm_plot_data <- combined_summary_dup2
        rm_plot_data$Method <- "Reed-Muench"
        
        mixed_plot <- combined_summary_dup2
        mixed_plot$Method <- "ACM"
        forced_plot <- combined_summary_dup2_forced4pl
        forced_plot$Method <- "F4PL"
        best_plot <- combined_summary_dup2_5pl
        best_plot$Method <- "F5PL"
        bayesian_plot <- combined_summary_dup2_bayesian
        bayesian_plot$Method <- "Bayes-4PL"
        
        multi_plot_data <- rbind(mixed_plot, forced_plot, best_plot, bayesian_plot, rm_plot_data)
        multi_plot_data$RM_Titer_R1 <- ifelse(!is.na(multi_plot_data$RM_IC50_log) & !is.na(multi_plot_data$R1_IC50_dilution) & multi_plot_data$R1_IC50_dilution > 0,
                                                (10^multi_plot_data$RM_IC50_log / multi_plot_data$R1_IC50_dilution) * reference_concentration, NA)
        multi_plot_data <- multi_plot_data[!is.na(multi_plot_data$ELISA_EU_mL) & !is.na(multi_plot_data$Titer_R1) & !multi_plot_data$Is_R1 & !multi_plot_data$Is_R2 & !multi_plot_data$Is_R3 & !multi_plot_data$Is_mAb, ]
        if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data)) {
          prefix_pattern <- paste0("^", sample_prefix_filter)
          multi_plot_data <- multi_plot_data[grepl(prefix_pattern, multi_plot_data$Sample, ignore.case = FALSE), ]
        }
        
        make_multi_method_plot <- function(df, suffix, concordance_filter = NULL, method_compare_alpha = 0.8) {
          if(nrow(df) == 0) return(NULL)
          plot_df <- df
          if(!is.null(concordance_filter)) {
            plot_df <- plot_df[plot_df$Concordance_R1_ELISA == concordance_filter, ]
          }
          if(nrow(plot_df) == 0) return(NULL)
          
           method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
           method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                             "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
           method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                              "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
           method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
          
          plot_df$Method <- factor(plot_df$Method, levels=method_order)
          plot_df <- plot_df[order(plot_df$Method), ]
          plot_df$Method_alpha <- method_alphas[as.character(plot_df$Method)]
          
          p <- ggplot(plot_df, aes(x=log10(ELISA_EU_mL), y=log10(Titer_R1), color=Method, shape=Method, size=Method, alpha=Method_alpha)) +
            geom_point() +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
            scale_color_manual(values=method_compare_colors) +
            scale_shape_manual(values=method_shapes) +
            scale_size_manual(values=method_sizes) +
            scale_alpha_identity() +
            guides(size="none", alpha="none") +
            labs( "All", concordance_filter)),
                 x="log10(ELISA EU/mL)", y=paste0("log10(Titer ", R1_name, " (", standard_unit, "))")) +
            theme_publication_custom(base_font_size) +
            theme(legend.position="bottom")
          return(p)
        }
        
        p15_a <- make_multi_method_plot(multi_plot_data, "A (2-rep Concordant)", "Concordant")
        p15_b <- make_multi_method_plot(multi_plot_data, "B (2-rep Discordant)", "Discordant")
        
        p15_c <- NULL
        p15_d <- NULL
        if(!is.null(combined_summary_dup3)) {
          rm_plot_data_3rep <- combined_summary_dup3
          rm_plot_data_3rep$Method <- "Reed-Muench"
          mixed_plot_3rep <- combined_summary_dup3
          mixed_plot_3rep$Method <- "ACM"
          forced_plot_3rep <- combined_summary_dup3_forced4pl
          forced_plot_3rep$Method <- "F4PL"
          best_plot_3rep <- combined_summary_dup3_5pl
          best_plot_3rep$Method <- "F5PL"
          bayesian_plot_3rep <- combined_summary_dup3_bayesian
          bayesian_plot_3rep$Method <- "Bayes-4PL"
          multi_plot_data_3rep <- rbind(mixed_plot_3rep, forced_plot_3rep, best_plot_3rep, bayesian_plot_3rep, rm_plot_data_3rep)
          multi_plot_data_3rep$RM_Titer_R1 <- ifelse(!is.na(multi_plot_data_3rep$RM_IC50_log) & !is.na(multi_plot_data_3rep$R1_IC50_dilution) & multi_plot_data_3rep$R1_IC50_dilution > 0,
                                                       (10^multi_plot_data_3rep$RM_IC50_log / multi_plot_data_3rep$R1_IC50_dilution) * reference_concentration, NA)
          multi_plot_data_3rep <- multi_plot_data_3rep[!is.na(multi_plot_data_3rep$ELISA_EU_mL) & !is.na(multi_plot_data_3rep$Titer_R1) & !multi_plot_data_3rep$Is_R1, ]
          if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_3rep)) {
            prefix_pattern <- paste0("^", sample_prefix_filter)
            multi_plot_data_3rep <- multi_plot_data_3rep[grepl(prefix_pattern, multi_plot_data_3rep$Sample, ignore.case = FALSE), ]
          }
          p15_c <- make_multi_method_plot(multi_plot_data_3rep, "C (3-rep Concordant)", "Concordant")
          p15_d <- make_multi_method_plot(multi_plot_data_3rep, "D (3-rep Discordant)", "Discordant")
        }
        
        all_panels <- list(p15_a, p15_b, p15_c, p15_d)
        all_panels <- all_panels[!sapply(all_panels, is.null)]
        
        if(length(all_panels) > 0) {
          if(length(all_panels) == 1) {
            ggsave(file.path(plot_dir, paste0("Fig15_ELISA_vs_Titer_All_Methods_", R1_file, ".png")), all_panels[[1]], width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          } else if(length(all_panels) == 2) {
            ggsave(file.path(plot_dir, paste0("Fig15_ELISA_vs_Titer_All_Methods_", R1_file, ".png")), arrangeGrob(all_panels[[1]], all_panels[[2]], ncol=2), width=plot_width_two_panel, height=plot_height_two_panel, dpi=dpi)
          } else {
            ggsave(file.path(plot_dir, paste0("Fig15_ELISA_vs_Titer_All_Methods_", R1_file, ".png")), arrangeGrob(all_panels[[1]], all_panels[[2]], all_panels[[3]], all_panels[[4]], ncol=2), width=plot_width_four_panel, height=plot_height_four_panel, dpi=dpi)
          }
        }
      }
      
      multi_plot_data_R2 <- rbind(mixed_plot, forced_plot, best_plot, bayesian_plot, rm_plot_data)
      multi_plot_data_R2$RM_Titer_R2 <- ifelse(!is.na(multi_plot_data_R2$RM_IC50_log) & !is.na(multi_plot_data_R2$R2_IC50_dilution) & multi_plot_data_R2$R2_IC50_dilution > 0,
                                                   (10^multi_plot_data_R2$RM_IC50_log / multi_plot_data_R2$R2_IC50_dilution) * second_reference_concentration, NA)
      multi_plot_data_R2 <- multi_plot_data_R2[!is.na(multi_plot_data_R2$ELISA_EU_mL) & !is.na(multi_plot_data_R2$Titer_R2) & !multi_plot_data_R2$Is_R1 & !multi_plot_data_R2$Is_R2 & !multi_plot_data_R2$Is_R3 & !multi_plot_data_R2$Is_mAb, ]
      if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_R2)) {
        prefix_pattern <- paste0("^", sample_prefix_filter)
        multi_plot_data_R2 <- multi_plot_data_R2[grepl(prefix_pattern, multi_plot_data_R2$Sample, ignore.case = FALSE), ]
      }
      if(use_dual_reference && nrow(multi_plot_data_R2) > 0) {
        method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
        method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                          "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
        method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                           "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
        method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
        multi_plot_data_R2$Method <- factor(multi_plot_data_R2$Method, levels=method_order)
        multi_plot_data_R2 <- multi_plot_data_R2[order(multi_plot_data_R2$Method), ]
        multi_plot_data_R2$Method_alpha <- method_alphas[as.character(multi_plot_data_R2$Method)]
        
        p15b_left <- ggplot(multi_plot_data_R2, aes(x=log10(ELISA_EU_mL), y=log10(Titer_R2), color=Method, shape=Method, size=Method, alpha=Method_alpha)) +
          geom_point() +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
          scale_color_manual(values=method_compare_colors) +
          scale_shape_manual(values=method_shapes) +
          scale_size_manual(values=method_sizes) +
          scale_alpha_identity() +
          guides(size="none", alpha="none") +
          labs(
               x="log10(ELISA EU/mL)", y=paste0("log10(Titer ", R2_name, " (", second_unit, "))")) +
          theme_publication_custom(base_font_size) +
          theme(legend.position="bottom")
        p15b_right <- NULL
        if(!is.null(combined_summary_dup3)) {
          rm_plot_data_3rep_R2 <- combined_summary_dup3
          rm_plot_data_3rep_R2$Method <- "Reed-Muench"
          mixed_plot_3rep_R2 <- combined_summary_dup3
          mixed_plot_3rep_R2$Method <- "ACM"
          forced_plot_3rep_R2 <- combined_summary_dup3_forced4pl
          forced_plot_3rep_R2$Method <- "F4PL"
          best_plot_3rep_R2 <- combined_summary_dup3_5pl
          best_plot_3rep_R2$Method <- "F5PL"
          bayesian_plot_3rep_R2 <- combined_summary_dup3_bayesian
          bayesian_plot_3rep_R2$Method <- "Bayes-4PL"
          multi_plot_data_R2_3rep <- rbind(mixed_plot_3rep_R2, forced_plot_3rep_R2, best_plot_3rep_R2, bayesian_plot_3rep_R2, rm_plot_data_3rep_R2)
          multi_plot_data_R2_3rep$RM_Titer_R2 <- ifelse(!is.na(multi_plot_data_R2_3rep$RM_IC50_log) & !is.na(multi_plot_data_R2_3rep$R2_IC50_dilution) & multi_plot_data_R2_3rep$R2_IC50_dilution > 0,
                                                            (10^multi_plot_data_R2_3rep$RM_IC50_log / multi_plot_data_R2_3rep$R2_IC50_dilution) * second_reference_concentration, NA)
           multi_plot_data_R2_3rep <- multi_plot_data_R2_3rep[!is.na(multi_plot_data_R2_3rep$ELISA_EU_mL) & !is.na(multi_plot_data_R2_3rep$Titer_R2) & !multi_plot_data_R2_3rep$Is_R1 & !multi_plot_data_R2_3rep$Is_R2 & !multi_plot_data_R2_3rep$Is_R3 & !multi_plot_data_R2_3rep$Is_mAb, ]
          if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_R2_3rep)) {
            prefix_pattern <- paste0("^", sample_prefix_filter)
            multi_plot_data_R2_3rep <- multi_plot_data_R2_3rep[grepl(prefix_pattern, multi_plot_data_R2_3rep$Sample, ignore.case = FALSE), ]
          }
          if(nrow(multi_plot_data_R2_3rep) > 0) {
            method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
            method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                              "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
            method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                               "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
            method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
            multi_plot_data_R2_3rep$Method <- factor(multi_plot_data_R2_3rep$Method, levels=method_order)
            multi_plot_data_R2_3rep <- multi_plot_data_R2_3rep[order(multi_plot_data_R2_3rep$Method), ]
            multi_plot_data_R2_3rep$Method_alpha <- method_alphas[as.character(multi_plot_data_R2_3rep$Method)]
            
            p15b_right <- ggplot(multi_plot_data_R2_3rep, aes(x=log10(ELISA_EU_mL), y=log10(Titer_R2), color=Method, shape=Method, size=Method, alpha=Method_alpha)) +
              geom_point() +
              geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
              geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
              scale_color_manual(values=method_compare_colors) +
              scale_shape_manual(values=method_shapes) +
              scale_size_manual(values=method_sizes) +
              scale_alpha_identity() +
              guides(size="none", alpha="none") +
              labs(
                   x="log10(ELISA EU/mL)", y=paste0("log10(Titer ", R2_name, " (", second_unit, "))")) +
              theme_publication_custom(base_font_size) +
              theme(legend.position="bottom")
          }
        }
        if(!is.null(p15b_right)) {
          ggsave(file.path(plot_dir, paste0("Fig15b_ELISA_vs_Titer_All_Methods_", R2_file, ".png")), arrangeGrob(p15b_left, p15b_right, ncol=2), width=plot_width_two_panel, height=plot_height_two_panel, dpi=dpi)
        } else {
          ggsave(file.path(plot_dir, paste0("Fig15b_ELISA_vs_Titer_All_Methods_", R2_file, ".png")), p15b_left, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
      }
      
      multi_plot_data_nonlin <- rbind(mixed_plot, forced_plot, best_plot, bayesian_plot, rm_plot_data)
      multi_plot_data_nonlin$RM_Titer_R1 <- ifelse(!is.na(multi_plot_data_nonlin$RM_IC50_log) & !is.na(multi_plot_data_nonlin$R1_IC50_dilution) & multi_plot_data_nonlin$R1_IC50_dilution > 0,
                                                     (10^multi_plot_data_nonlin$RM_IC50_log / multi_plot_data_nonlin$R1_IC50_dilution) * reference_concentration, NA)
      multi_plot_data_nonlin <- multi_plot_data_nonlin[!is.na(multi_plot_data_nonlin$ELISA_EU_mL_nonlin) & !is.na(multi_plot_data_nonlin$Titer_R1) & !multi_plot_data_nonlin$Is_R1, ]
      if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_nonlin)) {
        prefix_pattern <- paste0("^", sample_prefix_filter)
        multi_plot_data_nonlin <- multi_plot_data_nonlin[grepl(prefix_pattern, multi_plot_data_nonlin$Sample, ignore.case = FALSE), ]
      }
      if(nrow(multi_plot_data_nonlin) > 0) {
        method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
        method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                          "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
        method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                           "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
        method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
        multi_plot_data_nonlin$Method <- factor(multi_plot_data_nonlin$Method, levels=method_order)
        multi_plot_data_nonlin$Method_alpha <- method_alphas[as.character(multi_plot_data_nonlin$Method)]
        
        multi_plot_data_nonlin <- multi_plot_data_nonlin[order(multi_plot_data_nonlin$Method), ]
        
        p15c_left <- ggplot(multi_plot_data_nonlin, aes(x=log10(ELISA_EU_mL_nonlin), y=log10(Titer_R1), color=Method, shape=Method, size=Method), alpha=Method_alpha) +
          geom_point() +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
          scale_color_manual(values=method_compare_colors) +
          scale_shape_manual(values=method_shapes) +
          scale_size_manual(values=method_sizes) +
          guides(size="none") +
          scale_alpha_identity() +
          guides(size="none", alpha="none") +
          labs(
               x="log10(ELISA EU/mL non-linear)", y=paste0("log10(Titer ", R1_name, " (", standard_unit, "))")) +
          theme_publication_custom(base_font_size) +
          theme(legend.position="bottom")
        p15c_right <- NULL
        if(!is.null(combined_summary_dup3)) {
          rm_plot_data_3rep_nonlin <- combined_summary_dup3
          rm_plot_data_3rep_nonlin$Method <- "Reed-Muench"
          mixed_plot_3rep_nonlin <- combined_summary_dup3
          mixed_plot_3rep_nonlin$Method <- "ACM"
          forced_plot_3rep_nonlin <- combined_summary_dup3_forced4pl
          forced_plot_3rep_nonlin$Method <- "F4PL"
          best_plot_3rep_nonlin <- combined_summary_dup3_5pl
          best_plot_3rep_nonlin$Method <- "F5PL"
          bayesian_plot_3rep_nonlin <- combined_summary_dup3_bayesian
          bayesian_plot_3rep_nonlin$Method <- "Bayes-4PL"
          multi_plot_data_nonlin_3rep <- rbind(mixed_plot_3rep_nonlin, forced_plot_3rep_nonlin, best_plot_3rep_nonlin, bayesian_plot_3rep_nonlin, rm_plot_data_3rep_nonlin)
          multi_plot_data_nonlin_3rep$RM_Titer_R1 <- ifelse(!is.na(multi_plot_data_nonlin_3rep$RM_IC50_log) & !is.na(multi_plot_data_nonlin_3rep$R1_IC50_dilution) & multi_plot_data_nonlin_3rep$R1_IC50_dilution > 0,
                                                              (10^multi_plot_data_nonlin_3rep$RM_IC50_log / multi_plot_data_nonlin_3rep$R1_IC50_dilution) * reference_concentration, NA)
          multi_plot_data_nonlin_3rep <- multi_plot_data_nonlin_3rep[!is.na(multi_plot_data_nonlin_3rep$ELISA_EU_mL_nonlin) & !is.na(multi_plot_data_nonlin_3rep$Titer_R1) & !multi_plot_data_nonlin_3rep$Is_R1, ]
          if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_nonlin_3rep)) {
            prefix_pattern <- paste0("^", sample_prefix_filter)
            multi_plot_data_nonlin_3rep <- multi_plot_data_nonlin_3rep[grepl(prefix_pattern, multi_plot_data_nonlin_3rep$Sample, ignore.case = FALSE), ]
          }
          if(nrow(multi_plot_data_nonlin_3rep) > 0) {
            method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
            method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                              "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
            method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                               "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
            method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
            multi_plot_data_nonlin_3rep$Method <- factor(multi_plot_data_nonlin_3rep$Method, levels=method_order)
            multi_plot_data_nonlin_3rep$Method_alpha <- method_alphas[as.character(multi_plot_data_nonlin_3rep$Method)]
            
            multi_plot_data_nonlin_3rep <- multi_plot_data_nonlin_3rep[order(multi_plot_data_nonlin_3rep$Method), ]
            
            p15c_right <- ggplot(multi_plot_data_nonlin_3rep, aes(x=log10(ELISA_EU_mL_nonlin), y=log10(Titer_R1), color=Method, shape=Method, size=Method), alpha=Method_alpha) +
              geom_point() +
              geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
              geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
              scale_color_manual(values=method_compare_colors) +
              scale_shape_manual(values=method_shapes) +
              scale_size_manual(values=method_sizes) +
              guides(size="none") +
              scale_alpha_identity() +
              guides(size="none", alpha="none") +
              labs(
                   x="log10(ELISA EU/mL non-linear)", y=paste0("log10(Titer ", R1_name, " (", standard_unit, "))")) +
              theme_publication_custom(base_font_size) +
              theme(legend.position="bottom")
          }
        }
        if(!is.null(p15c_right)) {
          ggsave(file.path(plot_dir, paste0("Fig15c_ELISA_Nonlin_vs_Titer_All_Methods_", R1_file, ".png")), arrangeGrob(p15c_left, p15c_right, ncol=2), width=plot_width_two_panel, height=plot_height_two_panel, dpi=dpi)
        } else {
          ggsave(file.path(plot_dir, paste0("Fig15c_ELISA_Nonlin_vs_Titer_All_Methods_", R1_file, ".png")), p15c_left, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
      }
      
      if(use_dual_reference) {
        multi_plot_data_R2_nonlin <- rbind(mixed_plot, forced_plot, best_plot, bayesian_plot, rm_plot_data)
        multi_plot_data_R2_nonlin$RM_Titer_R2 <- ifelse(!is.na(multi_plot_data_R2_nonlin$RM_IC50_log) & !is.na(multi_plot_data_R2_nonlin$R2_IC50_dilution) & multi_plot_data_R2_nonlin$R2_IC50_dilution > 0,
                                                            (10^multi_plot_data_R2_nonlin$RM_IC50_log / multi_plot_data_R2_nonlin$R2_IC50_dilution) * second_reference_concentration, NA)
        multi_plot_data_R2_nonlin <- multi_plot_data_R2_nonlin[!is.na(multi_plot_data_R2_nonlin$ELISA_EU_mL_nonlin) & !is.na(multi_plot_data_R2_nonlin$Titer_R2) & !multi_plot_data_R2_nonlin$Is_R1 & !multi_plot_data_R2_nonlin$Is_R2 & !multi_plot_data_R2_nonlin$Is_R3 & !multi_plot_data_R2_nonlin$Is_mAb, ]
        if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_R2_nonlin)) {
          prefix_pattern <- paste0("^", sample_prefix_filter)
          multi_plot_data_R2_nonlin <- multi_plot_data_R2_nonlin[grepl(prefix_pattern, multi_plot_data_R2_nonlin$Sample, ignore.case = FALSE), ]
        }
        if(nrow(multi_plot_data_R2_nonlin) > 0) {
           method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
           method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                             "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
           method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                              "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
           method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
          multi_plot_data_R2_nonlin$Method <- factor(multi_plot_data_R2_nonlin$Method, levels=method_order)
          multi_plot_data_R2_nonlin$Method_alpha <- method_alphas[as.character(multi_plot_data_R2_nonlin$Method)]
          
          multi_plot_data_R2_nonlin <- multi_plot_data_R2_nonlin[order(multi_plot_data_R2_nonlin$Method), ]
          
          p15d_left <- ggplot(multi_plot_data_R2_nonlin, aes(x=log10(ELISA_EU_mL_nonlin), y=log10(Titer_R2), color=Method, shape=Method, size=Method), alpha=Method_alpha) +
            geom_point() +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
            scale_color_manual(values=method_compare_colors) +
            scale_shape_manual(values=method_shapes) +
            scale_size_manual(values=method_sizes) +
            guides(size="none") +
            scale_alpha_identity() +
            guides(size="none", alpha="none") +
            labs(
                 x="log10(ELISA EU/mL non-linear)", y=paste0("log10(Titer ", R2_name, " (", second_unit, "))")) +
            theme_publication_custom(base_font_size) +
            theme(legend.position="bottom")
          p15d_right <- NULL
          if(!is.null(combined_summary_dup3)) {
            rm_plot_data_3rep_R2_nonlin <- combined_summary_dup3
            rm_plot_data_3rep_R2_nonlin$Method <- "Reed-Muench"
            mixed_plot_3rep_R2_nonlin <- combined_summary_dup3
            mixed_plot_3rep_R2_nonlin$Method <- "ACM"
            forced_plot_3rep_R2_nonlin <- combined_summary_dup3_forced4pl
            forced_plot_3rep_R2_nonlin$Method <- "F4PL"
            best_plot_3rep_R2_nonlin <- combined_summary_dup3_5pl
            best_plot_3rep_R2_nonlin$Method <- "F5PL"
            bayesian_plot_3rep_R2_nonlin <- combined_summary_dup3_bayesian
            bayesian_plot_3rep_R2_nonlin$Method <- "Bayes-4PL"
            multi_plot_data_R2_nonlin_3rep <- rbind(mixed_plot_3rep_R2_nonlin, forced_plot_3rep_R2_nonlin, best_plot_3rep_R2_nonlin, bayesian_plot_3rep_R2_nonlin, rm_plot_data_3rep_R2_nonlin)
            multi_plot_data_R2_nonlin_3rep$RM_Titer_R2 <- ifelse(!is.na(multi_plot_data_R2_nonlin_3rep$RM_IC50_log) & !is.na(multi_plot_data_R2_nonlin_3rep$R2_IC50_dilution) & multi_plot_data_R2_nonlin_3rep$R2_IC50_dilution > 0,
                                                                     (10^multi_plot_data_R2_nonlin_3rep$RM_IC50_log / multi_plot_data_R2_nonlin_3rep$R2_IC50_dilution) * second_reference_concentration, NA)
             multi_plot_data_R2_nonlin_3rep <- multi_plot_data_R2_nonlin_3rep[!is.na(multi_plot_data_R2_nonlin_3rep$ELISA_EU_mL_nonlin) & !is.na(multi_plot_data_R2_nonlin_3rep$Titer_R2) & !multi_plot_data_R2_nonlin_3rep$Is_R1 & !multi_plot_data_R2_nonlin_3rep$Is_R2 & !multi_plot_data_R2_nonlin_3rep$Is_R3 & !multi_plot_data_R2_nonlin_3rep$Is_mAb, ]
            if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(multi_plot_data_R2_nonlin_3rep)) {
              prefix_pattern <- paste0("^", sample_prefix_filter)
              multi_plot_data_R2_nonlin_3rep <- multi_plot_data_R2_nonlin_3rep[grepl(prefix_pattern, multi_plot_data_R2_nonlin_3rep$Sample, ignore.case = FALSE), ]
            }
            if(nrow(multi_plot_data_R2_nonlin_3rep) > 0) {
              method_shapes <- c("ACM"=16, "F4PL"=17, "F5PL"=best_fit_shape, "Bayes-4PL"=15, "Reed-Muench"=18)
              method_sizes <- c("ACM"=mixed_symbol_size, "F4PL"=forced_symbol_size,
                                "F5PL"=bestfit_symbol_size, "Bayes-4PL"=bayesian_symbol_size, "Reed-Muench"=reedmuench_symbol_size)
              method_alphas <- c("ACM"=mixed_symbol_alpha, "F4PL"=forced_symbol_alpha,
                                 "F5PL"=bestfit_symbol_alpha, "Bayes-4PL"=bayesian_symbol_alpha, "Reed-Muench"=reedmuench_symbol_alpha)
              method_order <- c("Reed-Muench", "F5PL", "F4PL", "Bayes-4PL", "ACM")
              multi_plot_data_R2_nonlin_3rep$Method <- factor(multi_plot_data_R2_nonlin_3rep$Method, levels=method_order)
              multi_plot_data_R2_nonlin_3rep$Method_alpha <- method_alphas[as.character(multi_plot_data_R2_nonlin_3rep$Method)]
              
              multi_plot_data_R2_nonlin_3rep <- multi_plot_data_R2_nonlin_3rep[order(multi_plot_data_R2_nonlin_3rep$Method), ]
              
              p15d_right <- ggplot(multi_plot_data_R2_nonlin_3rep, aes(x=log10(ELISA_EU_mL_nonlin), y=log10(Titer_R2), color=Method, shape=Method, size=Method), alpha=Method_alpha) +
                geom_point() +
                geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
                geom_smooth(method="lm", se=TRUE, alpha=0.2, linewidth=0.8) +
                scale_color_manual(values=method_compare_colors) +
                scale_shape_manual(values=method_shapes) +
                scale_size_manual(values=method_sizes) +
                guides(size="none") +
                scale_alpha_identity() +
                guides(size="none", alpha="none") +
                labs(
                     x="log10(ELISA EU/mL non-linear)", y=paste0("log10(Titer ", R2_name, " (", second_unit, "))")) +
                theme_publication_custom(base_font_size) +
                theme(legend.position="bottom")
            }
          }
          if(!is.null(p15d_right)) {
            ggsave(file.path(plot_dir, paste0("Fig15d_ELISA_Nonlin_vs_Titer_All_Methods_", R2_file, ".png")), arrangeGrob(p15d_left, p15d_right, ncol=2), width=plot_width_two_panel, height=plot_height_two_panel, dpi=dpi)
          } else {
            ggsave(file.path(plot_dir, paste0("Fig15d_ELISA_Nonlin_vs_Titer_All_Methods_", R2_file, ".png")), p15d_left, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
          }
        }
      }
    }, error = function(e) {
    })
    
    
    
    
    plot_titer_delta_set <- function(diff_df, base_fig, ref_name, unit, ref_tag, elisa_tag, include_density) {
      diff_df$Method <- factor(diff_df$Method, levels = c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench"))
      
      p <- ggplot(diff_df, aes(x = Method, y = Delta, fill = Method)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
        geom_hline(yintercept = conc_tol_upper, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_hline(yintercept = conc_tol_lower, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
        scale_fill_manual(values = method_compare_colors) +
        labs(
             subtitle = diff_sub_main,
             x = "", y = diff_y_lab(unit)) +
        theme_publication_custom(base_font_size) +
        theme(legend.position = "none")
      ggsave(file.path(plot_dir, paste0("Fig", base_fig, "_Titer_Delta_Boxplot_", ref_tag, "_", elisa_tag, ".png")), p, width = plot_width_one_panel * 1.2, height = plot_height_one_panel, dpi = dpi)
      
      forest_data <- aggregate(Delta ~ Method, data = diff_df, FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = length(x)))
      forest_data <- do.call(data.frame, forest_data)
      bias_methods <- c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench")
      forest_data <- forest_data[forest_data$Method %in% bias_methods, ]
      forest_data$Method <- factor(forest_data$Method, levels = bias_methods)
      forest_data$CI_Lower <- forest_data$Delta.mean - 1.96 * forest_data$Delta.sd / sqrt(forest_data$Delta.n)
      forest_data$CI_Upper <- forest_data$Delta.mean + 1.96 * forest_data$Delta.sd / sqrt(forest_data$Delta.n)
      
      p <- ggplot(forest_data, aes(x = Method, y = Delta.mean, color = Method)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
        geom_hline(yintercept = conc_tol_upper, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_hline(yintercept = conc_tol_lower, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 1) +
        scale_color_manual(values = method_compare_colors) +
        labs(
             subtitle = "Mean difference with 95% CI; points show method bias relative to ELISA",
             x = "", y = diff_y_lab_mn(unit)) +
        theme_publication_custom(base_font_size) +
        theme(legend.position = "none")
      ggsave(file.path(plot_dir, paste0("Fig", base_fig, "b_Titer_Delta_Forest_", ref_tag, "_", elisa_tag, ".png")), p, width = plot_width_one_panel * 1.2, height = plot_height_one_panel, dpi = dpi)
      
      p <- ggplot(diff_df, aes(x = Method, y = Delta, fill = Method)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
        geom_hline(yintercept = conc_tol_upper, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_hline(yintercept = conc_tol_lower, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_violin(alpha = 0.6, scale = "width") +
        geom_boxplot(width = 0.1, alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
        scale_fill_manual(values = method_compare_colors) +
        labs(
             subtitle = diff_sub_dist,
             x = "", y = diff_y_lab(unit)) +
        theme_publication_custom(base_font_size) +
        theme(legend.position = "none")
      ggsave(file.path(plot_dir, paste0("Fig", base_fig, "c_Titer_Delta_Violin_", ref_tag, "_", elisa_tag, ".png")), p, width = plot_width_one_panel * 1.2, height = plot_height_one_panel, dpi = dpi)
      
      heatmap_data <- aggregate(cbind(Delta, Abs_Delta) ~ Method, data = diff_df, FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), median = median(x, na.rm = TRUE), iqr = IQR(x, na.rm = TRUE), max = max(x, na.rm = TRUE), min = min(x, na.rm = TRUE)))
      heatmap_data <- do.call(data.frame, heatmap_data)
      heatmap_wide <- reshape2::melt(heatmap_data, id.vars = "Method")
      heatmap_wide$Metric <- gsub("\\.mean$", " (Mean)", heatmap_wide$variable)
      heatmap_wide$Metric <- gsub("\\.sd$", " (SD)", heatmap_wide$Metric)
      heatmap_wide$Metric <- gsub("\\.median$", " (Median)", heatmap_wide$Metric)
      heatmap_wide$Metric <- gsub("\\.iqr$", " (IQR)", heatmap_wide$Metric)
      heatmap_wide$Metric <- gsub("\\.max$", " (Max)", heatmap_wide$Metric)
      heatmap_wide$Metric <- gsub("\\.min$", " (Min)", heatmap_wide$Metric)
      heatmap_wide$Metric <- gsub("Abs_Delta\\.", "", heatmap_wide$Metric)
      
      p <- ggplot(heatmap_wide, aes(x = Metric, y = Method, fill = value)) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_gradient2(low = "#DC143C", mid = "white", high = "#228B22", midpoint = 0, name = "Value") +
        geom_text(aes(label = round(value, 3)), size = 3, color = "black") +
        labs(
             x = "Statistic", y = "Method") +
        theme_publication_custom(base_font_size) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
      ggsave(file.path(plot_dir, paste0("Fig", base_fig, "d_Titer_Delta_Heatmap_", ref_tag, "_", elisa_tag, ".png")), p, width = plot_width_one_panel * 1.5, height = plot_height_one_panel, dpi = dpi)
      
      if(include_density && requireNamespace("ggridges", quietly = TRUE)) {
        p <- ggplot(diff_df, aes(x = Delta, y = Method, fill = Method)) +
          geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
          geom_vline(xintercept = conc_tol_upper, linetype = "dotted", color = "grey40", linewidth = 0.6) +
          geom_vline(xintercept = conc_tol_lower, linetype = "dotted", color = "grey40", linewidth = 0.6) +
          scale_fill_manual(values = method_compare_colors) +
          labs(
               subtitle = diff_sub_dens,
               x = diff_x_lab(unit),
               y = "Method") +
          theme_publication_custom(base_font_size) +
          theme(legend.position = "none")
        ggsave(file.path(plot_dir, paste0("Fig", base_fig, "e_Titer_Delta_Density_", ref_tag, "_", elisa_tag, ".png")), p, width = plot_width_one_panel * 1.2, height = plot_height_one_panel * 1.2, dpi = dpi)
      }
    }
    
    prepare_diff_data <- function(data_df, titer_col, elisa_col, ref_conc, method_name) {
      if(is.null(data_df) || !is.data.frame(data_df) || nrow(data_df) == 0) return(NULL)
      if(!all(c(titer_col, elisa_col) %in% colnames(data_df))) return(NULL)
      if(!all(c("Is_R1", "Is_R2", "Is_R3", "Is_mAb") %in% colnames(data_df))) return(NULL)
      d <- data_df[!is.na(data_df[[titer_col]]) & !is.na(data_df[[elisa_col]]) & !data_df$Is_R1 & !data_df$Is_R2 & !data_df$Is_R3 & !data_df$Is_mAb, ]
      if(nrow(d) == 0) return(NULL)
      if(!is.null(sample_prefix_filter) && "Sample" %in% colnames(d)) {
        prefix_pattern <- paste0("^", sample_prefix_filter)
        d <- d[grepl(prefix_pattern, d$Sample, ignore.case = FALSE), ]
      }
      if(nrow(d) == 0) return(NULL)
      d$Method <- method_name
      d$Delta <- log10(d[[titer_col]]) - log10(d[[elisa_col]])
      if(diff_pct_mode) d$Delta <- (10^d$Delta - 1) * 100
      d <- d[is.finite(d$Delta), ]
      if(nrow(d) == 0) return(NULL)
      d$Mean_Delta <- mean(d$Delta, na.rm = TRUE)
      d$SD_Delta <- sd(d$Delta, na.rm = TRUE)
      if(!is.finite(d$Mean_Delta[1]) || !is.finite(d$SD_Delta[1]) || d$SD_Delta[1] <= 0) return(NULL)
      d$LoA_Upper <- d$Mean_Delta + 1.96 * d$SD_Delta
      d$LoA_Lower <- d$Mean_Delta - 1.96 * d$SD_Delta
      d$Abs_Delta <- abs(d$Delta)
      d$Max_Discordance <- max(d$Abs_Delta, na.rm = TRUE)
      return(d)
    }
    
    tryCatch({
      d_mixed_r1 <- prepare_diff_data(combined_summary_dup2, "Titer_R1", "ELISA_EU_mL", reference_concentration, "ACM")
      d_forced_r1 <- prepare_diff_data(combined_summary_dup2_forced4pl, "Titer_R1", "ELISA_EU_mL", reference_concentration, "F4PL")
      d_best_r1 <- prepare_diff_data(combined_summary_dup2_5pl, "Titer_R1", "ELISA_EU_mL", reference_concentration, "F5PL")
      d_bayesian_r1 <- prepare_diff_data(combined_summary_dup2_bayesian, "Titer_R1", "ELISA_EU_mL", reference_concentration, "Bayes-4PL")
      d_rm_r1 <- NULL
      if(use_reed_muench_comparison && "RM_IC50_log" %in% colnames(combined_summary_dup2)) {
        temp_df <- combined_summary_dup2
        temp_df$RM_Titer_R1 <- ifelse(!is.na(temp_df$RM_IC50_log) & !is.na(temp_df$R1_IC50_dilution) & temp_df$R1_IC50_dilution > 0,
                                        (10^temp_df$RM_IC50_log / temp_df$R1_IC50_dilution) * reference_concentration, NA)
        d_rm_r1 <- prepare_diff_data(temp_df, "RM_Titer_R1", "ELISA_EU_mL", reference_concentration, "Reed-Muench")
      }
      
      all_diff_r1 <- rbind(d_mixed_r1, d_forced_r1, d_best_r1, d_bayesian_r1, d_rm_r1)
      if(!is.null(all_diff_r1) && nrow(all_diff_r1) > 0) {
        plot_titer_delta_set(all_diff_r1, 16, R1_name, standard_unit, "R1", "Linear", FALSE)
      }
      
      if(has_elisa_nonlin) {
        d_mixed_r1_nl <- prepare_diff_data(combined_summary_dup2, "Titer_R1", "ELISA_EU_mL_nonlin", reference_concentration, "ACM")
        d_forced_r1_nl <- prepare_diff_data(combined_summary_dup2_forced4pl, "Titer_R1", "ELISA_EU_mL_nonlin", reference_concentration, "F4PL")
        d_best_r1_nl <- prepare_diff_data(combined_summary_dup2_5pl, "Titer_R1", "ELISA_EU_mL_nonlin", reference_concentration, "F5PL")
        d_bayesian_r1_nl <- prepare_diff_data(combined_summary_dup2_bayesian, "Titer_R1", "ELISA_EU_mL_nonlin", reference_concentration, "Bayes-4PL")
        d_rm_r1_nl <- NULL
        if(use_reed_muench_comparison && "RM_IC50_log" %in% colnames(combined_summary_dup2)) {
          temp_df <- combined_summary_dup2
          temp_df$RM_Titer_R1 <- ifelse(!is.na(temp_df$RM_IC50_log) & !is.na(temp_df$R1_IC50_dilution) & temp_df$R1_IC50_dilution > 0,
                                          (10^temp_df$RM_IC50_log / temp_df$R1_IC50_dilution) * reference_concentration, NA)
          d_rm_r1_nl <- prepare_diff_data(temp_df, "RM_Titer_R1", "ELISA_EU_mL_nonlin", reference_concentration, "Reed-Muench")
        }
        
        all_diff_r1_nl <- rbind(d_mixed_r1_nl, d_forced_r1_nl, d_best_r1_nl, d_bayesian_r1_nl, d_rm_r1_nl)
        if(!is.null(all_diff_r1_nl) && nrow(all_diff_r1_nl) > 0) {
          plot_titer_delta_set(all_diff_r1_nl, 17, R1_name, standard_unit, "R1", "Nonlin", TRUE)
        }
      }
      
      if(use_dual_reference) {
        d_mixed_r2 <- prepare_diff_data(combined_summary_dup2, "Titer_R2", "ELISA_EU_mL", second_reference_concentration, "ACM")
        d_forced_r2 <- prepare_diff_data(combined_summary_dup2_forced4pl, "Titer_R2", "ELISA_EU_mL", second_reference_concentration, "F4PL")
        d_best_r2 <- prepare_diff_data(combined_summary_dup2_5pl, "Titer_R2", "ELISA_EU_mL", second_reference_concentration, "F5PL")
        d_bayesian_r2 <- prepare_diff_data(combined_summary_dup2_bayesian, "Titer_R2", "ELISA_EU_mL", second_reference_concentration, "Bayes-4PL")
        d_rm_r2 <- NULL
        if(use_reed_muench_comparison && "RM_IC50_log" %in% colnames(combined_summary_dup2)) {
          temp_df <- combined_summary_dup2
          temp_df$RM_Titer_R2 <- ifelse(!is.na(temp_df$RM_IC50_log) & !is.na(temp_df$R2_IC50_dilution) & temp_df$R2_IC50_dilution > 0,
                                          (10^temp_df$RM_IC50_log / temp_df$R2_IC50_dilution) * second_reference_concentration, NA)
          d_rm_r2 <- prepare_diff_data(temp_df, "RM_Titer_R2", "ELISA_EU_mL", second_reference_concentration, "Reed-Muench")
        }
        
        all_diff_r2 <- rbind(d_mixed_r2, d_forced_r2, d_best_r2, d_bayesian_r2, d_rm_r2)
        if(!is.null(all_diff_r2) && nrow(all_diff_r2) > 0) {
          plot_titer_delta_set(all_diff_r2, 18, R2_name, second_unit, "R2", "Linear", TRUE)
        }
        
        if(has_elisa_nonlin) {
          d_mixed_r2_nl <- prepare_diff_data(combined_summary_dup2, "Titer_R2", "ELISA_EU_mL_nonlin", second_reference_concentration, "ACM")
          d_forced_r2_nl <- prepare_diff_data(combined_summary_dup2_forced4pl, "Titer_R2", "ELISA_EU_mL_nonlin", second_reference_concentration, "F4PL")
          d_best_r2_nl <- prepare_diff_data(combined_summary_dup2_5pl, "Titer_R2", "ELISA_EU_mL_nonlin", second_reference_concentration, "F5PL")
          d_bayesian_r2_nl <- prepare_diff_data(combined_summary_dup2_bayesian, "Titer_R2", "ELISA_EU_mL_nonlin", second_reference_concentration, "Bayes-4PL")
          d_rm_r2_nl <- NULL
          if(use_reed_muench_comparison && "RM_IC50_log" %in% colnames(combined_summary_dup2)) {
            temp_df <- combined_summary_dup2
            temp_df$RM_Titer_R2 <- ifelse(!is.na(temp_df$RM_IC50_log) & !is.na(temp_df$R2_IC50_dilution) & temp_df$R2_IC50_dilution > 0,
                                            (10^temp_df$RM_IC50_log / temp_df$R2_IC50_dilution) * second_reference_concentration, NA)
            d_rm_r2_nl <- prepare_diff_data(temp_df, "RM_Titer_R2", "ELISA_EU_mL_nonlin", second_reference_concentration, "Reed-Muench")
          }
          
          all_diff_r2_nl <- rbind(d_mixed_r2_nl, d_forced_r2_nl, d_best_r2_nl, d_bayesian_r2_nl, d_rm_r2_nl)
          if(!is.null(all_diff_r2_nl) && nrow(all_diff_r2_nl) > 0) {
            plot_titer_delta_set(all_diff_r2_nl, 19, R2_name, second_unit, "R2", "Nonlin", TRUE)
          }
        }
      }
      
      if(nrow(sensitivity_specificity) > 0) {
        sens_data <- sensitivity_specificity
        sens_data$Method <- factor(sens_data$Method, levels = c("ACM", "F4PL", "F5PL", "Bayes-4PL", "Reed-Muench"))
        
        p20 <- ggplot(sens_data, aes(x = Method, y = Sensitivity, fill = Method)) +
          geom_bar(stat = "identity", alpha = 0.7, width = 0.6, position = position_dodge(width = 0.8)) +
          geom_bar(aes(y = Specificity), stat = "identity", alpha = 0.4, width = 0.6, position = position_dodge(width = 0.8), fill = "grey50") +
          geom_text(aes(label = paste0(Sensitivity, "%")), vjust = -0.5, size = 3.5, position = position_dodge(width = 0.8)) +
          scale_fill_manual(values = method_compare_colors) +
          labs(x = "", y = "Percentage (%)") +
          theme_publication_custom(base_font_size) +
          theme(legend.position = "none") +
          ylim(0, 105)
        ggsave(file.path(plot_dir, "Fig20_Sensitivity_Specificity_All_Methods.png"), p20, width = plot_width_one_panel * 1.2, height = plot_height_one_panel, dpi = dpi)
      }
      
    }, error = function(e) {
    })
    
    tryCatch({
      .sk <- function(v) gsub("[^A-Za-z0-9._-]+", "", tolower(as.character(v)))
      .grp_map <- NULL
      if(file.exists(file_path)) {
        .s4 <- tryCatch(as.data.frame(readxl::read_xlsx(file_path, sheet = group_concordance_sheet)), error = function(e) NULL)
        if(!is.null(.s4) && is.character(elisa_titer_delta_group) && elisa_titer_delta_group %in% colnames(.s4)) {
          .sc <- group_concordance_sample_col
          if(is.null(.sc) || !.sc %in% colnames(.s4)) {
            .sc <- NULL
            for(.cn in colnames(.s4)) if(tolower(.cn) %in% c("sample","sample_id")) { .sc <- .cn; break }
          }
          if(!is.null(.sc)) .grp_map <- setNames(trimws(as.character(.s4[[elisa_titer_delta_group]])), .sk(.s4[[.sc]]))
        }
      }
      .canon_order <- c("<0.125 EU/mL", "0.125 to <0.5 EU/mL", "0.5 to <2.0 EU/mL",
                        "2.0 to <4.0 EU/mL", ">=4.0 EU/mL", ">=4 EU/mL")
      grp_levels <- "Ungrouped"
      if(!is.null(.grp_map)) {
        .gvals <- trimws(as.character(.s4[[elisa_titer_delta_group]]))
        .gvals[tolower(.gvals) %in% c("na","nan","null","")] <- NA
        .raw <- unique(.gvals[!is.na(.gvals)])
        grp_levels <- unique(c(.canon_order[.canon_order %in% .raw], setdiff(.raw, .canon_order)))
      }
      .final_levels <- grp_levels
      if("Ungrouped" %in% .final_levels) .final_levels <- c(setdiff(.final_levels, "Ungrouped"), "Ungrouped")
      .okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E25","#CC79A7","#000000")
      .grp_cols <- if(length(.final_levels) <= length(.okabe)) .okabe[seq_along(.final_levels)] else colorRampPalette(.okabe)(length(.final_levels))
      .grp_cols_named <- setNames(.grp_cols, .final_levels)
      .dify <- function(titer, elisa, mode) {
        ok <- is.finite(titer) & is.finite(elisa) & !is.na(titer) & !is.na(elisa)
        if(mode == "linear")  d <- ifelse(ok, titer - elisa, NA_real_)
        else if(mode == "percent") d <- ifelse(ok & elisa > 0, 100 * (titer/elisa - 1), NA_real_)
        else d <- ifelse(ok & titer > 0 & elisa > 0, log10(titer/elisa), NA_real_)
        d
      }
      .ylab <- if(elisa_titer_delta_y == "linear") "Titer - ELISA (EU/mL)"
               else if(elisa_titer_delta_y == "percent") "Percent difference (Titer vs ELISA) (%)"
               else "Log10(Titer / ELISA)"
      .build <- function(data_df, titer_col, elisa_col, rep_tag) {
        if(is.null(data_df) || !"Sample" %in% colnames(data_df)) return(NULL)
        if(!all(c(titer_col, elisa_col, "Is_R1", "Is_R2", "Is_R3", "Is_mAb") %in% colnames(data_df))) return(NULL)
        d <- data_df[!data_df$Is_R1 & !data_df$Is_R2 & !data_df$Is_R3 & !data_df$Is_mAb, , drop = FALSE]
        d <- d[!is.na(d[[titer_col]]) & !is.na(d[[elisa_col]]), , drop = FALSE]
        if(!is.null(sample_prefix_filter)) d <- d[grepl(paste0("^", sample_prefix_filter), d$Sample, ignore.case = FALSE), , drop = FALSE]
        if(nrow(d) == 0) return(NULL)
        d$Diff <- .dify(d[[titer_col]], d[[elisa_col]], elisa_titer_delta_y)
        d <- d[is.finite(d$Diff), , drop = FALSE]
        if(nrow(d) == 0) return(NULL)
        d$Group <- .grp_map[match(.sk(d$Sample), names(.grp_map))]
        d$Group <- ifelse(is.na(d$Group) | d$Group == "" | tolower(d$Group) %in% c("na","nan","null"), "Ungrouped", d$Group)
        d$Replicates <- rep_tag
        d[, c("Sample","Group","Diff","Replicates")]
      }
      .mk_rm <- function(df, ref_conc, sec_conc, use_dual) {
        if(is.null(df)) return(NULL)
        df <- as.data.frame(df)
        if("RM_IC50_log" %in% colnames(df) && "R1_IC50_dilution" %in% colnames(df)) {
          df$RM_Titer_R1 <- ifelse(!is.na(df$RM_IC50_log) & !is.na(df$R1_IC50_dilution) & df$R1_IC50_dilution > 0,
                                     (10^df$RM_IC50_log / df$R1_IC50_dilution) * ref_conc, NA)
          if(use_dual && "R2_IC50_dilution" %in% colnames(df))
            df$RM_Titer_R2 <- ifelse(!is.na(df$RM_IC50_log) & !is.na(df$R2_IC50_dilution) & df$R2_IC50_dilution > 0,
                                       (10^df$RM_IC50_log / df$R2_IC50_dilution) * sec_conc, NA)
        }
        df
      }
      .method_bygroup <- list(
        list(data = combined_summary_dup2,
             data_3rep = combined_summary_dup3,
             name = "ACM", label = "ACM", R1 = "Titer_R1", R2 = "Titer_R2"),
        list(data = combined_summary_dup2_forced4pl,
             data_3rep = combined_summary_dup3_forced4pl,
             name = "F4PL", label = "Forced_4PL", R1 = "Titer_R1", R2 = "Titer_R2"),
        list(data = combined_summary_dup2_5pl,
             data_3rep = combined_summary_dup3_5pl,
             name = "F5PL", label = "5PL", R1 = "Titer_R1", R2 = "Titer_R2"),
        list(data = combined_summary_dup2_bayesian,
             data_3rep = combined_summary_dup3_bayesian,
             name = "Bayes-4PL", label = "Bayes-4PL", R1 = "Titer_R1", R2 = "Titer_R2"),
        list(data = .mk_rm(combined_summary_dup2, reference_concentration, second_reference_concentration, use_dual_reference),
             data_3rep = .mk_rm(combined_summary_dup3, reference_concentration, second_reference_concentration, use_dual_reference),
             name = "Reed-Muench", label = "ReedMuench", R1 = "RM_Titer_R1", R2 = "RM_Titer_R2")
      )

      .delta_sets <- list()
      for(.m in .method_bygroup) {
        for(.et in c("Linear", "Non-linear")) {
          .ecol <- if(.et == "Linear") "ELISA_EU_mL" else "ELISA_EU_mL_nonlin"
          .rep_cfgs <- list(
            list(tag = "2 replicates", suffix = "2rep", df = .m$data),
            list(tag = "3 replicates", suffix = "3rep", df = .m$data_3rep)
          )
          .panels <- list()
          for(.rc in .rep_cfgs) {
            .mdata <- .rc$df
            if(is.null(.mdata)) next
            for(.ref in c("R1", "R2")) {
              if(.ref == "R2" && !use_dual_reference) next
              .tcol <- if(.ref == "R1") .m$R1 else .m$R2
              if(!(.tcol %in% colnames(.mdata))) next
              .p <- .build(.mdata, .tcol, .ecol, .rc$tag)
              if(!is.null(.p)) {
                .p$Ref <- .ref
                .p$ELISA_Type <- .et
                .panels[[paste0(.ref, "_", .rc$suffix)]] <- .p
              }
            }
          }
          if(length(.panels) == 0) next
          dd <- do.call(rbind, .panels)
          if(is.null(dd) || nrow(dd) == 0) next
          dd$Group <- factor(dd$Group, levels = .final_levels)
          dd$Ref <- factor(dd$Ref, levels = c("R1", "R2"))
          dd$Replicates <- factor(dd$Replicates, levels = c("2 replicates", "3 replicates"))
          dd$ELISA_Type <- factor(dd$ELISA_Type, levels = c("Linear", "Non-linear"))
          .panel_levels <- c("R1 / 2 replicates", "R2 / 2 replicates",
                             "R1 / 3 replicates", "R2 / 3 replicates")
          dd$Panel <- factor(paste(dd$Ref, "/", dd$Replicates), levels = .panel_levels)
          .delta_sets[[length(.delta_sets) + 1]] <- list(m = .m, dd = dd, et = .et)
        }
      }
      .ymax_all <- if(length(.delta_sets) > 0) max(sapply(.delta_sets, function(s) max(s$dd$Diff, na.rm = TRUE))) else 1
      if(!is.finite(.ymax_all) || .ymax_all <= 0) .ymax_all <- 1
      .ymax_all <- .ymax_all * 1.05

      .render_delta_fig <- function(.m, dd, .ymax_all, et) {
        .panel_levels <- c("R1 / 2 replicates", "R2 / 2 replicates",
                           "R1 / 3 replicates", "R2 / 3 replicates")
        .panel_grobs <- list()
        for(.pl in .panel_levels) {
          .sub <- dd[dd$Panel == .pl, , drop = FALSE]
          if(nrow(.sub) == 0) next
          .ymin <- min(.sub$Diff, na.rm = TRUE)
          if(!is.finite(.ymin)) .ymin <- 0
          .ymin <- .ymin - 0.05 * (.ymax_all - .ymin)
          .pg <- ggplot(.sub, aes(x = Group, y = Diff, fill = Group)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
            { if(elisa_titer_delta_violin) geom_violin(alpha = 0.5, scale = "width", color = "black", linewidth = 0.3) } +
            geom_point(position = position_jitter(width = 0.1), size = 1.5, alpha = 0.7,
                       shape = 21, color = "black", stroke = 0.3) +
            scale_fill_manual(values = .grp_cols_named) +
            scale_x_discrete(labels = function(x) gsub(" EU/mL", "\nEU/mL", x)) +
            coord_cartesian(ylim = c(.ymin, .ymax_all)) +
            labs(
                 x = elisa_titer_delta_group, y = .ylab) +
            theme_publication_custom(base_font_size) +
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                  legend.position = "none",
                  plot.title = element_text(size = title_size, face = "bold", hjust = 0.5))
          .panel_grobs[[length(.panel_grobs) + 1]] <- .pg
        }
        .leg_tmp <- ggplot(dd, aes(x = Group, y = Diff, fill = Group)) +
          geom_point() + scale_fill_manual(values = .grp_cols_named) +
          guides(fill = guide_legend(title = elisa_titer_delta_group,
                                     override.aes = list(alpha = 1, shape = 21, size = 3))) +
          theme(legend.position = "bottom",
                legend.direction = "horizontal",
                legend.title = element_text(size = legend_title_size, face = "bold"),
                legend.text = element_text(size = legend_text_size))
        .leg_grob <- ggplotGrob(.leg_tmp)
        .leg_idx <- which(sapply(.leg_grob$grobs, function(g) grepl("guide-box", g$name, fixed = TRUE)))
        if(length(.leg_idx) > 0) .leg_grob <- .leg_grob$grobs[[.leg_idx[1]]] else .leg_grob <- NULL
        .np <- length(.panel_grobs)
        .R2_present <- any(as.character(dd$Ref) == "R2" & is.finite(dd$Diff))
        if(.R2_present) {
          .panels_arr <- do.call(gridExtra::arrangeGrob,
                                 c(list(grobs = .panel_grobs, nrow = 2, ncol = 2)))
        } else {
          .panels_arr <- do.call(gridExtra::arrangeGrob,
                                 c(list(grobs = .panel_grobs, nrow = 2, ncol = 1)))
        }
        .et_tag <- if(et == "Linear") "Linear" else "Nonlinear"
        .title_txt <- paste0("ELISA vs Titer delta by group (", .m$name, ", ", et, " ELISA, ", elisa_titer_delta_y, " scale)")
        .leg_list <- if(!is.null(.leg_grob)) list(.leg_grob) else list()
        .comb <- do.call(gridExtra::arrangeGrob, c(
          list(grobs = c(list(.panels_arr), .leg_list),
               top = grid::textGrob(.title_txt, gp = grid::gpar(fontsize = title_size, fontface = "bold"))),
          if(length(.leg_list) > 0) list(nrow = 2, heights = c(0.88, 0.12)) else list(nrow = 1, heights = c(1))
        ))
        ggsave(file.path(plot_dir, paste0("Fig_ELISA_Titer_Delta_byGroup_", .m$label, "_", .et_tag, ".png")), .comb,
               width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
        if(.m$label == "ACM" && et == "Linear")
          ggsave(file.path(plot_dir, paste0("Fig_ELISA_Titer_Delta_byGroup.png")), .comb,
                 width = plot_width_four_panel, height = plot_height_four_panel, dpi = dpi)
      }

      for(.s in .delta_sets) {
        .render_delta_fig(.s$m, .s$dd, .ymax_all, .s$et)
      }
    }, error = function(e) {
    })
    
    tryCatch({
      method_concordance_list <- list(
        list(data = combined_summary_dup2, data_3rep = combined_summary_dup3, label = "ACM", name = "ACM"),
        list(data = combined_summary_dup2_forced4pl, data_3rep = combined_summary_dup3_forced4pl, label = "Forced_4PL", name = "F4PL"),
        list(data = combined_summary_dup2_5pl, data_3rep = combined_summary_dup3_5pl, label = "5PL", name = "F5PL"),
        list(data = combined_summary_dup2_bayesian, data_3rep = combined_summary_dup3_bayesian, label = "Bayes-4PL", name = "Bayes-4PL"),
         list(data = combined_summary_dup2[,c("Sample","Experiment_Group","Is_R1","Is_R2","Is_R3","Is_mAb","ELISA_EU_mL","ELISA_EU_mL_nonlin","Titer_R1","Titer_R2","Titer_R3","Titer_mAb","Concordance_R1_ELISA","Concordance_R2_ELISA","Concordance_R3_ELISA","Concordance_R1_Nonlin","Concordance_R2_Nonlin","Concordance_R3_Nonlin","RM_IC50_log","R1_IC50_dilution","R2_IC50_dilution")], 
              data_3rep = if(exists("combined_summary_dup3") && !is.null(combined_summary_dup3)) combined_summary_dup3[,c("Sample","Experiment_Group","Is_R1","Is_R2","Is_R3","Is_mAb","ELISA_EU_mL","ELISA_EU_mL_nonlin","Titer_R1","Titer_R2","Titer_R3","Titer_mAb","Concordance_R1_ELISA","Concordance_R2_ELISA","Concordance_R3_ELISA","Concordance_R1_Nonlin","Concordance_R2_Nonlin","Concordance_R3_Nonlin","RM_IC50_log","R1_IC50_dilution","R2_IC50_dilution")] else NULL,
             label = "ReedMuench", name = "Reed-Muench")
      )
      
      for(meth in method_concordance_list) {
        meth_data <- meth$data
        meth_data_3rep <- meth$data_3rep
        meth_label <- meth$name
        
        titer_col <- if(meth$label == "ReedMuench") "RM_Titer" else "Titer_R1"
        r1_conc_override <- NULL
        r1_nl_conc_override <- NULL
        if(meth$label == "ReedMuench") {
          meth_data$RM_Titer_R1 <- ifelse(!is.na(meth_data$RM_IC50_log) & !is.na(meth_data$R1_IC50_dilution) & meth_data$R1_IC50_dilution > 0,
                                            (10^meth_data$RM_IC50_log / meth_data$R1_IC50_dilution) * reference_concentration, NA)
          if(any(!is.na(meth_data$RM_Titer_R1)) && any(!is.na(meth_data$ELISA_EU_mL))) {
            r1_conc_override <- "Concordance_R1_ELISA_RM"
            meth_data$Concordance_R1_ELISA_RM <- ifelse(!is.na(meth_data$RM_Titer_R1) & !is.na(meth_data$ELISA_EU_mL),
                                                          check_concordance(meth_data$RM_Titer_R1, reference_concentration, meth_data$ELISA_EU_mL, log_tol=ed50_tolerance), NA)
          }
          if(has_elisa_nonlin && any(!is.na(meth_data$RM_Titer_R1)) && any(!is.na(meth_data$ELISA_EU_mL_nonlin))) {
            r1_nl_conc_override <- "Concordance_R1_Nonlin_RM"
            meth_data$Concordance_R1_Nonlin_RM <- ifelse(!is.na(meth_data$RM_Titer_R1) & !is.na(meth_data$ELISA_EU_mL_nonlin),
                                                           check_concordance(meth_data$RM_Titer_R1, reference_concentration, meth_data$ELISA_EU_mL_nonlin, log_tol=ed50_tolerance), NA)
          }
          if(use_dual_reference) {
            meth_data$RM_Titer_R2 <- ifelse(!is.na(meth_data$RM_IC50_log) & !is.na(meth_data$R2_IC50_dilution) & meth_data$R2_IC50_dilution > 0,
                                              (10^meth_data$RM_IC50_log / meth_data$R2_IC50_dilution) * second_reference_concentration, NA)
          }
          titer_col <- "RM_Titer_R1"
          
          if(!is.null(meth_data_3rep)) {
            meth_data_3rep$RM_Titer_R1 <- ifelse(!is.na(meth_data_3rep$RM_IC50_log) & !is.na(meth_data_3rep$R1_IC50_dilution) & meth_data_3rep$R1_IC50_dilution > 0,
                                                   (10^meth_data_3rep$RM_IC50_log / meth_data_3rep$R1_IC50_dilution) * reference_concentration, NA)
            if(any(!is.na(meth_data_3rep$RM_Titer_R1)) && any(!is.na(meth_data_3rep$ELISA_EU_mL))) {
              meth_data_3rep$Concordance_R1_ELISA_RM <- ifelse(!is.na(meth_data_3rep$RM_Titer_R1) & !is.na(meth_data_3rep$ELISA_EU_mL),
                                                                 check_concordance(meth_data_3rep$RM_Titer_R1, reference_concentration, meth_data_3rep$ELISA_EU_mL, log_tol=ed50_tolerance), NA)
            }
            if(has_elisa_nonlin && any(!is.na(meth_data_3rep$RM_Titer_R1)) && any(!is.na(meth_data_3rep$ELISA_EU_mL_nonlin))) {
              meth_data_3rep$Concordance_R1_Nonlin_RM <- ifelse(!is.na(meth_data_3rep$RM_Titer_R1) & !is.na(meth_data_3rep$ELISA_EU_mL_nonlin),
                                                                  check_concordance(meth_data_3rep$RM_Titer_R1, reference_concentration, meth_data_3rep$ELISA_EU_mL_nonlin, log_tol=ed50_tolerance), NA)
            }
            if(use_dual_reference) {
              meth_data_3rep$RM_Titer_R2 <- ifelse(!is.na(meth_data_3rep$RM_IC50_log) & !is.na(meth_data_3rep$R2_IC50_dilution) & meth_data_3rep$R2_IC50_dilution > 0,
                                                     (10^meth_data_3rep$RM_IC50_log / meth_data_3rep$R2_IC50_dilution) * second_reference_concentration, NA)
            }
          }
        }
        
        p_meth_r1 <- create_elisa_titer_plot(meth_data, "ELISA_EU_mL", titer_col,
                                             reference_concentration, R1_name, "ELISA", 
                                             paste0("Titer (", standard_unit, ")"), fig_counter, meth_label,
                                             conc_col_override = r1_conc_override,
                                             data_3rep = meth_data_3rep)
        if(!is.null(p_meth_r1)) {
          dims <- get_elisa_titer_dims(p_meth_r1)
          ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_ELISA_vs_Titer_", R1_file, "_", meth$label, ".png")), p_meth_r1, width=dims$width, height=dims$height, dpi=dpi)
          fig_counter <- fig_counter + 1
        }
        p_meth_r1_delta <- if(delta_plot) create_delta_plot(meth_data, "ELISA_EU_mL", titer_col,
                                             reference_concentration, R1_name, "ELISA", 
                                             paste0("Titer (", standard_unit, ")"), fig_counter, meth_label,
                                             conc_col_override = r1_conc_override,
                                             data_3rep = meth_data_3rep)
        if(!is.null(p_meth_r1_delta)) {
          dims <- get_elisa_titer_dims(p_meth_r1_delta)
          ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_Delta_", R1_file, "_", meth$label, ".png")), p_meth_r1_delta, width=dims$width, height=dims$height, dpi=dpi)
          fig_counter <- fig_counter + 1
        }
        
        if(has_elisa_nonlin) {
          p_meth_r1_nl <- create_elisa_titer_plot(meth_data, "ELISA_EU_mL_nonlin", titer_col,
                                                  reference_concentration, R1_name, "ELISA (non-linear)", 
                                                  paste0("Titer (", standard_unit, ")"), fig_counter, meth_label,
                                                  conc_col_override = r1_nl_conc_override,
                                                  data_3rep = meth_data_3rep)
          if(!is.null(p_meth_r1_nl)) {
            dims <- get_elisa_titer_dims(p_meth_r1_nl)
            ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_ELISA_Nonlin_vs_Titer_", R1_file, "_", meth$label, ".png")), p_meth_r1_nl, width=dims$width, height=dims$height, dpi=dpi)
            fig_counter <- fig_counter + 1
          }
          p_meth_r1_nl_delta <- if(delta_plot) create_delta_plot(meth_data, "ELISA_EU_mL_nonlin", titer_col,
                                                  reference_concentration, R1_name, "ELISA (non-linear)", 
                                                  paste0("Titer (", standard_unit, ")"), fig_counter, meth_label,
                                                  conc_col_override = r1_nl_conc_override,
                                                  data_3rep = meth_data_3rep)
          if(!is.null(p_meth_r1_nl_delta)) {
            dims <- get_elisa_titer_dims(p_meth_r1_nl_delta)
            ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_Delta_Nonlin_", R1_file, "_", meth$label, ".png")), p_meth_r1_nl_delta, width=dims$width, height=dims$height, dpi=dpi)
            fig_counter <- fig_counter + 1
          }
        }
        
        if(use_dual_reference) {
          r2_titer_col <- if(meth$label == "ReedMuench") "RM_Titer_R2" else "Titer_R2"
          r2_conc_override <- NULL
          if(meth$label == "ReedMuench") {
            meth_data$Concordance_R2_ELISA_RM <- ifelse(!is.na(meth_data$RM_Titer_R2) & !is.na(meth_data$ELISA_EU_mL),
                                                          check_concordance(meth_data$RM_Titer_R2, second_reference_concentration, meth_data$ELISA_EU_mL, log_tol=ed50_tolerance), NA)
            r2_conc_override <- "Concordance_R2_ELISA_RM"
            if(!is.null(meth_data_3rep)) {
              meth_data_3rep$Concordance_R2_ELISA_RM <- ifelse(!is.na(meth_data_3rep$RM_Titer_R2) & !is.na(meth_data_3rep$ELISA_EU_mL),
                                                                 check_concordance(meth_data_3rep$RM_Titer_R2, second_reference_concentration, meth_data_3rep$ELISA_EU_mL, log_tol=ed50_tolerance), NA)
            }
          }
          p_meth_r2 <- create_elisa_titer_plot(meth_data, "ELISA_EU_mL", r2_titer_col,
                                               second_reference_concentration, R2_name, "ELISA", 
                                               paste0("Titer (", second_unit, ")"), fig_counter, meth_label,
                                               conc_col_override = r2_conc_override,
                                               data_3rep = meth_data_3rep)
          if(!is.null(p_meth_r2)) {
            dims <- get_elisa_titer_dims(p_meth_r2)
            ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_ELISA_vs_Titer_", R2_file, "_", meth$label, ".png")), p_meth_r2, width=dims$width, height=dims$height, dpi=dpi)
            fig_counter <- fig_counter + 1
          }
          p_meth_r2_delta <- if(delta_plot) create_delta_plot(meth_data, "ELISA_EU_mL", r2_titer_col,
                                               second_reference_concentration, R2_name, "ELISA", 
                                               paste0("Titer (", second_unit, ")"), fig_counter, meth_label,
                                               conc_col_override = r2_conc_override,
                                               data_3rep = meth_data_3rep)
          if(!is.null(p_meth_r2_delta)) {
            dims <- get_elisa_titer_dims(p_meth_r2_delta)
            ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_Delta_", R2_file, "_", meth$label, ".png")), p_meth_r2_delta, width=dims$width, height=dims$height, dpi=dpi)
            fig_counter <- fig_counter + 1
          }
          
          if(has_elisa_nonlin) {
            r2_nl_titer_col <- if(meth$label == "ReedMuench") "RM_Titer_R2" else "Titer_R2"
            r2_nl_conc_override <- NULL
            if(meth$label == "ReedMuench") {
              r2_nl_titer_col <- "RM_Titer_R2"
              meth_data$Concordance_R2_Nonlin_RM <- ifelse(!is.na(meth_data$RM_Titer_R2) & !is.na(meth_data$ELISA_EU_mL_nonlin),
                                                             check_concordance(meth_data$RM_Titer_R2, second_reference_concentration, meth_data$ELISA_EU_mL_nonlin, log_tol=ed50_tolerance), NA)
              r2_nl_conc_override <- "Concordance_R2_Nonlin_RM"
              if(!is.null(meth_data_3rep)) {
                meth_data_3rep$Concordance_R2_Nonlin_RM <- ifelse(!is.na(meth_data_3rep$RM_Titer_R2) & !is.na(meth_data_3rep$ELISA_EU_mL_nonlin),
                                                                    check_concordance(meth_data_3rep$RM_Titer_R2, second_reference_concentration, meth_data_3rep$ELISA_EU_mL_nonlin, log_tol=ed50_tolerance), NA)
              }
            }
            p_meth_r2_nl <- create_elisa_titer_plot(meth_data, "ELISA_EU_mL_nonlin", r2_nl_titer_col,
                                                    second_reference_concentration, R2_name, "ELISA (non-linear)", 
                                                    paste0("Titer (", second_unit, ")"), fig_counter, meth_label,
                                                    conc_col_override = r2_nl_conc_override,
                                                    data_3rep = meth_data_3rep)
            if(!is.null(p_meth_r2_nl)) {
              dims <- get_elisa_titer_dims(p_meth_r2_nl)
              ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_ELISA_Nonlin_vs_Titer_", R2_file, "_", meth$label, ".png")), p_meth_r2_nl, width=dims$width, height=dims$height, dpi=dpi)
              fig_counter <- fig_counter + 1
            }
            p_meth_r2_nl_delta <- if(delta_plot) create_delta_plot(meth_data, "ELISA_EU_mL_nonlin", r2_nl_titer_col,
                                                    second_reference_concentration, R2_name, "ELISA (non-linear)", 
                                                    paste0("Titer (", second_unit, ")"), fig_counter, meth_label,
                                                    conc_col_override = r2_nl_conc_override,
                                                    data_3rep = meth_data_3rep)
            if(!is.null(p_meth_r2_nl_delta)) {
              dims <- get_elisa_titer_dims(p_meth_r2_nl_delta)
              ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_Delta_Nonlin_", R2_file, "_", meth$label, ".png")), p_meth_r2_nl_delta, width=dims$width, height=dims$height, dpi=dpi)
              fig_counter <- fig_counter + 1
            }
          }
        }
      }
      
      if(use_reed_muench_comparison && !is.null(comparison_ACM_forced) && nrow(comparison_ACM_forced)>0) {
        
        valid_plot_rm1 <- !is.na(comparison_ACM_forced$IC50_log_dilution_ACM) & !is.na(comparison_ACM_forced$RM_IC50_log)
        if(sum(valid_plot_rm1, na.rm=TRUE) > 1) {
          p_rm1 <- ggplot(comparison_ACM_forced[valid_plot_rm1,], aes(x=RM_IC50_log, y=IC50_log_dilution_ACM)) +
            geom_point(size=point_size, alpha=0.8, color="#2E86AB") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
             labs(title="Figure 28: IC50 - ACM vs Reed-Muench",
                  4)),
                  x="IC50_log (Reed-Muench)", y="IC50_log (ACM)") +
            theme_publication_custom(base_font_size)
          ic50_rm_panels[["Fig28"]] <- p_rm1
        }
        
        valid_plot_rm2 <- !is.na(comparison_ACM_forced$IC50_log_dilution_Forced4PL) & !is.na(comparison_ACM_forced$RM_IC50_log)
        if(sum(valid_plot_rm2, na.rm=TRUE) > 1) {
          p_rm2 <- ggplot(comparison_ACM_forced[valid_plot_rm2,], aes(x=RM_IC50_log, y=IC50_log_dilution_Forced4PL)) +
            geom_point(size=point_size, alpha=0.8, color="#F18F01") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
            labs(title="Figure 29: IC50 - F4PL vs Reed-Muench",
                 4)),
                 x="IC50_log (Reed-Muench)", y="IC50_log (F4PL)") +
            theme_publication_custom(base_font_size)
          ic50_rm_panels[["Fig29"]] <- p_rm2
        }
        
        valid_plot_rm3 <- !is.na(comparison_ACM_forced$IC50_log_dilution_5PL) & !is.na(comparison_ACM_forced$RM_IC50_log)
        if(sum(valid_plot_rm3, na.rm=TRUE) > 1) {
          p_rm3 <- ggplot(comparison_ACM_forced[valid_plot_rm3,], aes(x=RM_IC50_log, y=IC50_log_dilution_5PL)) +
            geom_point(size=point_size, alpha=0.8, color="#6A3D9A") +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
            labs(title="Figure 30: IC50 - 5PL vs Reed-Muench",
                 4)),
                 x="IC50_log (Reed-Muench)", y="IC50_log (F5PL)") +
            theme_publication_custom(base_font_size)
          ic50_rm_panels[["Fig30"]] <- p_rm3
        }
        
        rm_data <- combined_summary[!is.na(combined_summary$RM_IC50_log) & !combined_summary$Is_R1 & !combined_summary$Is_R2 & !combined_summary$Is_R3 & !combined_summary$Is_mAb, ]
        rm_data$RM_Titer_R1 <- ifelse(!is.na(rm_data$RM_IC50_log) & !is.na(rm_data$R1_IC50_dilution) & rm_data$R1_IC50_dilution > 0,
                                        (10^rm_data$RM_IC50_log / rm_data$R1_IC50_dilution) * reference_concentration, NA)
        if(use_dual_reference) {
          rm_data$RM_Titer_R2 <- ifelse(!is.na(rm_data$RM_IC50_log) & !is.na(rm_data$R2_IC50_dilution) & rm_data$R2_IC50_dilution > 0,
                                          (10^rm_data$RM_IC50_log / rm_data$R2_IC50_dilution) * second_reference_concentration, NA)
        }
        if(nrow(rm_data) > 0) {
          p_rm_titer <- ggplot(rm_data, aes(x=RM_Titer_R1, y=Titer_R1, color=IC50_Method)) +
            geom_point(size=point_size, alpha=0.8) +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            add_regression() +
            scale_color_manual(values=method_colors) +
             labs(title="Reed-Muench vs 4PL Titer",
                 x="Reed-Muench Titer", y=paste0("4PL Titer (", standard_unit, ")")) +
            theme_publication_custom(base_font_size)
          ggsave(file.path(plot_dir, paste0("Fig", fig_counter, "_RM_vs_4PL_Titer.png")), p_rm_titer, width=plot_width_one_panel, height=plot_height_one_panel, dpi=dpi)
        }
      }
      
      if(length(ic50_rm_panels) >= 1) {
        panel_order <- c("Fig5", "Fig28", "Fig29", "Fig30")
        four_panels <- lapply(panel_order, function(nm) if(nm %in% names(ic50_rm_panels)) ic50_rm_panels[[nm]] else NULL)
        four_panels <- four_panels[!sapply(four_panels, is.null)]
        if(length(four_panels) == 2) {
          ggsave(file.path(plot_dir, "Fig_IC50_Comparison_FourPanel.png"),
                 arrangeGrob(four_panels[[1]], four_panels[[2]], ncol=2),
                 width=plot_width_four_panel, height=plot_height_four_panel*0.6, dpi=dpi)
        } else if(length(four_panels) == 3) {
          ggsave(file.path(plot_dir, "Fig_IC50_Comparison_FourPanel.png"),
                 arrangeGrob(four_panels[[1]], four_panels[[2]], four_panels[[3]], ncol=3),
                 width=plot_width_four_panel, height=plot_height_four_panel*0.5, dpi=dpi)
        } else {
          ggsave(file.path(plot_dir, "Fig_IC50_Comparison_FourPanel.png"),
                 arrangeGrob(four_panels[[1]], four_panels[[2]], four_panels[[3]], four_panels[[4]], nrow=2, ncol=2),
                 width=plot_width_four_panel, height=plot_height_four_panel, dpi=dpi)
        }
      }
    }, error = function(e) {
    })
  }
  
  
  
  
  result_list <- list(
    combined_summary = combined_summary,
    combined_summary_forced4pl = combined_summary_dup2_forced4pl,
    combined_summary_sk = combined_summary_dup2_5pl,
    combined_summary_bayesian = combined_summary_dup2_bayesian,
    combined_dilutions = combined_dilutions,
    concordance_summary = concordance_summary,
    method_stats = method_stats,
    method_stats_forced = if(exists("method_stats_forced")) method_stats_forced else list(),
    method_stats_f5pl = if(exists("method_stats_f5pl")) method_stats_f5pl else list(),
    method_stats_bayesian = if(exists("method_stats_bayesian")) method_stats_bayesian else list(),
    method_stats_rm = if(use_reed_muench_comparison) method_stats_rm else NULL
  )
  
  if(exists("comparison_ACM_forced") && !is.null(comparison_ACM_forced) && nrow(comparison_ACM_forced)>0) {
    result_list$method_comparison <- comparison_ACM_forced
    if(!is.null(concordance_rates) && nrow(concordance_rates)>0) {
      result_list$concordance_rates <- concordance_rates
    }
  }
  
  if(has_rlu3) {
    result_list$combined_summary_dup3 <- combined_summary_dup3
    result_list$combined_summary_dup3_forced4pl <- combined_summary_dup3_forced4pl
    result_list$combined_summary_dup3_5pl <- combined_summary_dup3_5pl
    result_list$comparison_df <- comparison_df
  }
  
   gc_stats <- NULL
   gc_per_sample <- NULL
   if(isTRUE(group_concordance_individual) || isTRUE(group_concordance_composite)) {
    tryCatch({
      
      gc_reg <- tolower(as.character(group_concordance_regression))[1]
      if(!gc_reg %in% c("per_group", "all_data", "both")) {
                  "' - falling back to 'per_group'")
        gc_reg <- "per_group"
      }
      gc_reg_group <- gc_reg %in% c("per_group", "both")
      gc_reg_all   <- gc_reg %in% c("all_data", "both")
      
      gc_out <- if(!is.null(group_concordance_outdir)) group_concordance_outdir else file.path(output_dir, "group_concordance")
      gc_plot_dir  <- file.path(gc_out, "plots")
      gc_bygrp_dir <- file.path(gc_plot_dir, "by_group")
      for(.d in c(gc_out, gc_plot_dir, gc_bygrp_dir)) if(!dir.exists(.d)) dir.create(.d, recursive = TRUE)
      
      gc_key <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
      
      gc_has_R1 <- all(c("Titer_R1", "Concordance_R1_ELISA") %in% colnames(combined_summary))
      gc_has_R2 <- all(c("Titer_R2", "Concordance_R2_ELISA") %in% colnames(combined_summary))
      gc_refs <- c("R1", "R2")[c(gc_has_R1, gc_has_R2)]
      if(length(gc_refs) == 0) {
      }
      
      gc_group_map <- NULL
      gc_elisa_linear_map <- NULL
      gc_elisa_nlinear_map <- NULL
      if(file.exists(file_path)) {
        gc_maps <- tryCatch({
          meta_df <- as.data.frame(readxl::read_xlsx(file_path, sheet = group_concordance_sheet))
          samp_col <- group_concordance_sample_col
          if(is.null(samp_col) || !samp_col %in% colnames(meta_df)) {
            samp_col <- NULL
            for(cn in colnames(meta_df)) if(tolower(cn) %in% c("sample", "sample_id")) { samp_col <- cn; break }
          }
           gwant <- gc_key(group_concordance_col)
           gcol <- NULL
           for(cn in colnames(meta_df)) if(gc_key(cn) == gwant) { gcol <- cn; break }
           find_col <- function(pattern, df) {
             cn <- colnames(df)
             exact <- cn[cn == pattern]
             if(length(exact) > 0) return(exact[1])
             partial <- cn[grepl(pattern, cn, fixed = TRUE)]
             if(length(partial) > 0) return(partial[1])
             return(NULL)
           }
           if(is.null(samp_col) || is.null(gcol)) {
                       group_concordance_sheet, "- grouping skipped")
             list(group = NULL, linear = NULL, nlinear = NULL)
           } else {
             sk <- gc_key(meta_df[[samp_col]])
             gv <- trimws(as.character(meta_df[[gcol]]))
             gv[gv == "" | tolower(gv) %in% c("na", "nan", "null")] <- NA
             mp <- setNames(gv, sk)
             mp <- mp[!is.na(names(mp)) & names(mp) != ""]
             mp <- mp[!duplicated(names(mp))]
             lin_col <- find_col("linear", meta_df)
             nlin_col <- find_col("nlinear", meta_df)
             lin_map <- if(!is.null(lin_col)) {
               x <- setNames(suppressWarnings(as.numeric(meta_df[[lin_col]])), sk)
               x[!is.na(names(x)) & names(x) != ""]
             } else NULL
             nlin_map <- if(!is.null(nlin_col)) {
               x <- setNames(suppressWarnings(as.numeric(meta_df[[nlin_col]])), sk)
               x[!is.na(names(x)) & names(x) != ""]
             } else NULL
                       group_concordance_sheet, "column", gcol,
                       "; linear ELISA titers:", length(lin_map), "; non-linear:", length(nlin_map))
             list(group = mp, linear = lin_map, nlinear = nlin_map)
           }
        }, error = function(e) {
          list(group = NULL, linear = NULL, nlinear = NULL)
        })
        gc_group_map <- gc_maps$group
        gc_elisa_linear_map <- gc_maps$linear
        gc_elisa_nlinear_map <- gc_maps$nlinear
      }
      
      if(length(gc_refs) == 0 || is.null(gc_group_map) || length(gc_group_map) == 0) {
        } else {
          gc_linear_dir  <- file.path(gc_out, "plots", "linear")
          gc_nonlin_dir  <- file.path(gc_out, "plots", "non-linear")
          gc_bygrp_dir   <- file.path(gc_out, "plots", "by_group")
          for(.d in c(gc_out, gc_linear_dir, gc_nonlin_dir, gc_bygrp_dir)) if(!dir.exists(.d)) dir.create(.d, recursive = TRUE)
          
          ref_params <- list(
            R1 = list(name = R1_name, conc = reference_concentration, unit = standard_unit, titer_col = "Titer_R1", dil_col = "R1_IC50_dilution"),
            R2 = list(name = R2_name, conc = second_reference_concentration, unit = second_unit, titer_col = "Titer_R2", dil_col = "R2_IC50_dilution")
          )
          
          all_ref_frames <- list()
          for(gc_ref_up in gc_refs) {
           
           gc_ref_name  <- if(gc_ref_up == "R1") R1_name else R2_name
           gc_ref_conc  <- if(gc_ref_up == "R1") reference_concentration else second_reference_concentration
           gc_ref_unit  <- if(gc_ref_up == "R1") standard_unit else second_unit
           gc_titer_col <- paste0("Titer_", gc_ref_up)
           gc_dil_col   <- paste0(gc_ref_up, "_IC50_dilution")
           gc_conc_colors <- c("Concordant" = "#228B22", "Discordant" = "#DC143C", "Unknown" = "#808080")
           
             gc_nonlin_ok <- ("ELISA_EU_mL_nonlin" %in% colnames(combined_summary) &&
                              any(!is.na(combined_summary$ELISA_EU_mL_nonlin)))
            gc_rm_ok <- all(c("RM_IC50_log", "R1_IC50_dilution", "R2_IC50_dilution") %in% colnames(combined_summary)) &&
              any(!is.na(combined_summary$RM_IC50_log))
             
               gc_methods <- c("ACM", "F4PL", "F5PL", "Bayes-4PL")
               if(gc_rm_ok) gc_methods <- c(gc_methods, "Reed-Muench")
            
            gc_combos <- list(
             list(ref = "R1", rep = "2rep", elisa = "Linear",     label = "R1 - 2 rep - linear ELISA"),
             list(ref = "R1", rep = "3rep", elisa = "Linear",     label = "R1 - 3 rep - linear ELISA"),
             list(ref = "R2", rep = "2rep", elisa = "Linear",     label = "R2 - 2 rep - linear ELISA"),
             list(ref = "R2", rep = "3rep", elisa = "Linear",     label = "R2 - 3 rep - linear ELISA")
           )
           
             gc_base_frame <- function(method, rep_tag) {
               if(rep_tag == "2rep") {
                 if(method == "F4PL")   return(combined_summary_dup2_forced4pl)
                 if(method == "F5PL") return(combined_summary_dup2_5pl)
                 if(method == "Bayes-4PL") return(combined_summary_dup2_bayesian)
                 return(combined_summary_dup2)
               }
               if(!exists("combined_summary_dup3") || is.null(combined_summary_dup3)) return(NULL)
               if(method == "F4PL")   return(if(exists("combined_summary_dup3_forced4pl")) combined_summary_dup3_forced4pl else NULL)
               if(method == "F5PL") return(if(exists("combined_summary_dup3_5pl")) combined_summary_dup3_5pl else NULL)
               if(method == "Bayes-4PL") return(if(exists("combined_summary_dup3_bayesian")) combined_summary_dup3_bayesian else NULL)
               combined_summary_dup3
             }
          
          gc_x_lab <- function(elisa_type) {
            nm <- if(identical(elisa_type, "Non-linear")) "ELISA non-linear EU/mL" else "ELISA EU/mL"
            if(elisa_scale == "log10") paste0("log10(", nm, ")") else nm
          }
          gc_y_lab <- function(ref_name = gc_ref_name, ref_unit = gc_ref_unit) {
            nm <- paste0("Titer ", ref_name, " (", ref_unit, ")")
            if(titer_scale == "log10") paste0("log10(", nm, ")") else nm
          }
          
           gc_frame <- function(method, rep_tag, elisa_type, ref = "R1", ref_conc = 0.5, titer_col = "Titer_R1", dil_col = "R1_IC50_dilution") {
             if(identical(elisa_type, "Non-linear") && !gc_nonlin_ok) return(NULL)
             base <- gc_base_frame(method, rep_tag)
             if(is.null(base) || !is.data.frame(base) || nrow(base) == 0) return(NULL)
             d <- base
             if(all(c("Is_R1", "Is_R2", "Is_R3", "Is_mAb") %in% colnames(d))) d <- d[!d$Is_R1 & !d$Is_R2 & !d$Is_R3 & !d$Is_mAb, , drop = FALSE]
             if(!is.null(sample_prefix_filter)) {
               d <- d[grepl(paste0("^", sample_prefix_filter), as.character(d$Sample), ignore.case = FALSE), , drop = FALSE]
             }
             if(nrow(d) == 0) return(NULL)
             
             elisa_col <- if(identical(elisa_type, "Non-linear")) "ELISA_EU_mL_nonlin" else "ELISA_EU_mL"
             if(!elisa_col %in% colnames(d)) return(NULL)
             elisa_v <- suppressWarnings(as.numeric(d[[elisa_col]]))
             
             conc_tol <- ed50_tolerance
             
              if(method == "Reed-Muench") {
                if(!all(c("RM_IC50_log", dil_col) %in% colnames(d))) return(NULL)
                rm_log  <- suppressWarnings(as.numeric(d$RM_IC50_log))
                ref_dil <- suppressWarnings(as.numeric(d[[dil_col]]))
                titer_v <- ifelse(!is.na(rm_log) & !is.na(ref_dil) & ref_dil > 0,
                                  (10^rm_log / ref_dil) * ref_conc, NA_real_)
                conc_v  <- ifelse(!is.na(titer_v) & !is.na(elisa_v),
                                  check_concordance(titer_v, ref_conc, elisa_v, log_tol = conc_tol),
                                  NA_character_)
                if(group_concordance_compute_absolute) {
                  conc_v_abs <- ifelse(!is.na(titer_v) & !is.na(elisa_v),
                                       check_concordance(titer_v, ref_conc, elisa_v, log_tol = 0),
                                       NA_character_)
                } else {
                  conc_v_abs <- rep(NA_character_, length(conc_v))
                }
              } else {
                if(!titer_col %in% colnames(d)) return(NULL)
                titer_v <- suppressWarnings(as.numeric(d[[titer_col]]))
                conc_v  <- ifelse(!is.na(titer_v) & !is.na(elisa_v),
                                  check_concordance(titer_v, ref_conc, elisa_v, log_tol = conc_tol),
                                  NA_character_)
                if(group_concordance_compute_absolute) {
                  conc_v_abs <- ifelse(!is.na(titer_v) & !is.na(elisa_v),
                                       check_concordance(titer_v, ref_conc, elisa_v, log_tol = 0),
                                       NA_character_)
                } else {
                  conc_v_abs <- rep(NA_character_, length(conc_v))
                }
              }
              grp <- gc_group_map[match(gc_key(d$Sample), names(gc_group_map))]
              out <- data.frame(Sample = as.character(d$Sample), Group = as.character(grp),
                                Concordance = as.character(conc_v), Concordance_Absolute = as.character(conc_v_abs),
                                ELISA = elisa_v, Titer = titer_v,
                                stringsAsFactors = FALSE)
              out <- out[!is.na(out$ELISA) & !is.na(out$Titer), , drop = FALSE]
             out$Group <- as.character(out$Group)
             out$Group[is.na(out$Group) | out$Group == ""] <- "Ungrouped"
             out <- out[out$Concordance %in% c("Concordant", "Discordant"), , drop = FALSE]
             if(isTRUE(group_concordance_concordant_only)) out <- out[out$Concordance == "Concordant", , drop = FALSE]
             if(nrow(out) == 0) return(NULL)
             out$x <- transform_elisa(out$ELISA, elisa_scale)$values
             out$y <- transform_titer(out$Titer, titer_scale)$values
             out <- out[is.finite(out$x) & is.finite(out$y), , drop = FALSE]
             if(nrow(out) == 0) return(NULL)
             out$Method     <- method
             out$Replicate  <- rep_tag
             out$ELISA_Type <- elisa_type
             out
           }
          
          gc_fit <- function(x, y) {
            keep <- is.finite(x) & is.finite(y)
            x <- x[keep]; y <- y[keep]
            res <- list(n = length(x), r = NA_real_, r2 = NA_real_, adj = NA_real_,
                        slope = NA_real_, intercept = NA_real_, p = NA_real_)
            if(res$n < 2 || length(unique(x)) < 2) return(res)
            res$r <- suppressWarnings(as.numeric(stats::cor(x, y)))
            fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
            if(is.null(fit)) return(res)
            sm <- summary(fit)
            res$r2  <- as.numeric(sm$r.squared)
            res$adj <- as.numeric(sm$adj.r.squared)
            cf <- stats::coef(fit)
            res$intercept <- as.numeric(cf[1])
            if(length(cf) >= 2) res$slope <- as.numeric(cf[2])
            if(!is.null(sm$coefficients) && nrow(sm$coefficients) >= 2 && ncol(sm$coefficients) >= 4) {
              res$p <- as.numeric(sm$coefficients[2, 4])
            }
            res
          }
          
           gc_stat_row <- function(df, method, elisa_type, rep_tag, group_label) {
             st <- gc_fit(df$x, df$y)
             cn <- sum(df$Concordance == "Concordant", na.rm = TRUE)
             ct <- sum(df$Concordance %in% c("Concordant", "Discordant"), na.rm = TRUE)
             cn_abs <- if("Concordance_Absolute" %in% colnames(df)) {
               sum(df$Concordance_Absolute == "Concordant", na.rm = TRUE)
             } else NA_real_
             ct_abs <- if("Concordance_Absolute" %in% colnames(df)) {
               sum(df$Concordance_Absolute %in% c("Concordant", "Discordant"), na.rm = TRUE)
             } else NA_real_
             data.frame(Method = method, ELISA = elisa_type, Replicate = rep_tag, Group = group_label,
                        N = st$n, Pearson_r = round(st$r, 4), R_squared = round(st$r2, 4),
                        Adj_R_squared = round(st$adj, 4), Slope = round(st$slope, 4),
                        Intercept = round(st$intercept, 4), P_value = signif(st$p, 4),
                        Concordant = cn, Discordant = ct - cn, Total = ct,
                        Concordance_Rate_pct = if(ct > 0) round(cn / ct * 100, 1) else NA_real_,
                        Concordant_Absolute = cn_abs, Discordant_Absolute = ct_abs - cn_abs, Total_Absolute = ct_abs,
                        Concordance_Rate_pct_Absolute = if(ct_abs > 0) round(cn_abs / ct_abs * 100, 1) else NA_real_,
                        Regression_Scope = gc_reg, Reference = gc_ref_name,
                        stringsAsFactors = FALSE)
           }
          
          gc_data <- list()
          for(m in gc_methods) gc_data[[m]] <- lapply(gc_combos, function(cb) {
            rp <- ref_params[[cb$ref]]
            gc_frame(m, cb$rep, cb$elisa, cb$ref, rp$conc, rp$titer_col, rp$dil_col)
          })
          all_ref_frames[[gc_ref_up]] <- gc_data
          all_ref_frames[[gc_ref_up]] <- gc_data
          
          .gc_raw <- unique(unlist(lapply(gc_data, function(lst)
            unlist(lapply(lst, function(fr) if(is.null(fr)) NULL else fr$Group)))))
          .gc_canon <- c("<0.125 EU/mL", "0.125 to <0.5 EU/mL", "0.5 to <2.0 EU/mL",
                         "2.0 to <4.0 EU/mL", ">=4.0 EU/mL", ">=4 EU/mL")
          gc_all_groups <- unique(c(.gc_canon[.gc_canon %in% .gc_raw], setdiff(.gc_raw, .gc_canon)))
          if("Ungrouped" %in% gc_all_groups) gc_all_groups <- c(setdiff(gc_all_groups, "Ungrouped"), "Ungrouped")
          
          if(length(gc_all_groups) == 0) {
          } else {
            gc_okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000",
                          "#88CCEE", "#CC6677", "#117733", "#332288", "#AA4499", "#E06C75", "#5E4FA2")
            gc_colors <- if(length(gc_all_groups) <= length(gc_okabe)) gc_okabe[seq_along(gc_all_groups)]
            else colorRampPalette(gc_okabe)(length(gc_all_groups))
            names(gc_colors) <- gc_all_groups
            if("Ungrouped" %in% gc_all_groups) gc_colors["Ungrouped"] <- "#808080"
            gc_shapes <- c("Concordant" = 16, "Discordant" = 17)
            gc_safe   <- function(s) gsub("[^A-Za-z0-9._-]+", "_", s)
            gc_blank  <- function(txt) grid::textGrob(txt, gp = grid::gpar(fontsize = max(9, base_font_size), col = "grey50"))
            gc_extract_legend <- function(p) {
              if(!inherits(p, "ggplot")) return(NULL)
              gb <- tryCatch(ggplot_build(p), error = function(e) NULL)
              if(is.null(gb)) return(NULL)
              ggt <- tryCatch(ggplot_gtable(gb), error = function(e) NULL)
              if(is.null(ggt)) return(NULL)
              idx <- which(sapply(ggt$grobs, function(g) g$name) == "guide-box")
              if(length(idx) > 0) ggt$grobs[[idx]] else NULL
            }
            
            gc_limits <- function(frames) {
              xs <- unlist(lapply(frames, function(fr) if(is.null(fr)) NULL else fr$x))
              ys <- unlist(lapply(frames, function(fr) if(is.null(fr)) NULL else fr$y))
              xs <- xs[is.finite(xs)]; ys <- ys[is.finite(ys)]
              if(length(xs) == 0 || length(ys) == 0) return(NULL)
              px <- diff(range(xs)) * 0.06; if(!is.finite(px) || px <= 0) px <- 0.2
              py <- diff(range(ys)) * 0.06; if(!is.finite(py) || py <= 0) py <- 0.2
              xlim <- c(min(xs) - px, max(xs) + px)
              ylim <- c(min(ys) - py, max(ys) + py)
              if(!is.null(group_concordance_xmin)) {
                xmin_t <- transform_elisa(group_concordance_xmin, elisa_scale)$values
                if(is.finite(xmin_t)) xlim[1] <- min(xlim[1], xmin_t)
              }
              if(!is.null(group_concordance_ymin)) {
                ymin_t <- transform_titer(group_concordance_ymin, titer_scale)$values
                if(is.finite(ymin_t)) ylim[1] <- min(ylim[1], ymin_t)
              }
              list(x = xlim, y = ylim)
            }
            
            gc_panel_allgroups <- function(fr, panel_title, lims, ref_name = gc_ref_name, ref_unit = gc_ref_unit) {
              if(is.null(fr) || nrow(fr) < 2) return(NULL)
              st <- gc_fit(fr$x, fr$y)
              cn <- sum(fr$Concordance == "Concordant", na.rm = TRUE)
              bits <- c(paste0("N = ", nrow(fr), " in ", length(unique(fr$Group)), " groups"),
                        paste0("concordant ", cn, "/", nrow(fr), " (", round(cn / nrow(fr) * 100, 1), "%)"))
              if(gc_reg_all) bits <- c(bits, paste0("all-data R2 = ", round(st$r2, 3), ", p = ", signif(st$p, 3)))
              if(gc_reg_group && !gc_reg_all) bits <- c(bits, "regression: one line per group")
              thx <- transform_elisa(group_concordance_threshold, elisa_scale)$values
              thy <- transform_titer(group_concordance_threshold, titer_scale)$values
              bits <- c(bits, paste0("dotted: ", group_concordance_threshold, " EU/mL ELISA (vert) / ",
                                     group_concordance_threshold, " U/mL titer (horiz)"))
              p <- ggplot(fr, aes(x = x, y = y)) +
                geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "red", linewidth = 0.4) +
                geom_vline(xintercept = thx, linetype = "dotted", colour = "#1f78b4", linewidth = 0.5) +
                geom_hline(yintercept = thy, linetype = "dotted", colour = "#33a02c", linewidth = 0.5) +
                geom_point(aes(colour = Group, shape = Concordance), size = point_size, alpha = 0.85)
              if(group_concordance_highlight_discordant) {
                p <- p + geom_point(data = subset(fr, Concordance == "Discordant"), colour = "red", size = point_size * 0.6, alpha = 0.9)
              }
              if(gc_reg_group) {
                ok_grp <- names(which(table(fr$Group) >= 3))
                fr_fit <- fr[fr$Group %in% ok_grp, , drop = FALSE]
                if(nrow(fr_fit) > 0) {
                  p <- p + geom_smooth(data = fr_fit, aes(colour = Group, group = Group),
                                       method = "lm", formula = y ~ x, se = FALSE,
                                       linewidth = 0.7, na.rm = TRUE)
                }
              }
              if(gc_reg_all && nrow(fr) >= 3) {
                p <- p + geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
                                     colour = "black", fill = "grey70", alpha = 0.25, linewidth = 0.9,
                                     linetype = if(gc_reg_group) "dashed" else "solid", na.rm = TRUE)
              }
              p <- p +
                scale_colour_manual(values = gc_colors, limits = gc_all_groups, breaks = gc_all_groups,
                                    drop = FALSE, name = "Titer category") +
                scale_shape_manual(values = gc_shapes, drop = FALSE, name = "Concordance") +
                 labs( x = gc_x_lab(fr$ELISA_Type[1]), y = gc_y_lab(ref_name, ref_unit)) +
                 theme_publication_custom(base_font_size) +
                 guides(colour = guide_legend(ncol = group_concordance_legend_cols, byrow = TRUE,
                                              override.aes = list(size = 2, alpha = 1)),
                        shape = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
                 theme(legend.position = "bottom", legend.box = "vertical",
                       legend.title = element_text(size = base_font_size - 1),
                       legend.text = element_text(size = base_font_size - 2),
                       legend.spacing.x = unit(0.08, "cm"),
                       legend.spacing.y = unit(0.02, "cm"),
                       legend.key.spacing.x = unit(0.04, "cm"),
                       legend.key.spacing.y = unit(0.01, "cm"),
                       legend.key.width = unit(0.5, "cm"),
                       legend.key.height = unit(0.18, "cm"),
                       legend.margin = margin(0, 0, 0, 0, "pt"),
                       legend.box.spacing = unit(0, "cm"),
                       legend.box.margin = margin(0, 0, 0, 0))
              if(!is.null(lims)) p <- p + coord_cartesian(xlim = lims$x, ylim = lims$y)
              p
            }
            
            gc_panel_group <- function(fr_all, group_name, panel_title, lims, ref_name = gc_ref_name, ref_unit = gc_ref_unit) {
              if(is.null(fr_all)) return(NULL)
              fr <- fr_all[!is.na(fr_all$Group) & fr_all$Group == group_name, , drop = FALSE]
              if(nrow(fr) == 0) return(NULL)
              st <- gc_fit(fr$x, fr$y)
              cn <- sum(fr$Concordance == "Concordant", na.rm = TRUE)
              bits <- c(paste0("N = ", nrow(fr)),
                        paste0("concordant ", cn, "/", nrow(fr), " (", round(cn / nrow(fr) * 100, 1), "%)"))
              if(gc_reg_group) bits <- c(bits, paste0("group R2 = ", round(st$r2, 3), ", p = ", signif(st$p, 3)))
              if(gc_reg_all) {
                sta <- gc_fit(fr_all$x, fr_all$y)
                bits <- c(bits, paste0("all-groups R2 = ", round(sta$r2, 3), ", p = ", signif(sta$p, 3)))
              }
              thx <- transform_elisa(group_concordance_threshold, elisa_scale)$values
              thy <- transform_titer(group_concordance_threshold, titer_scale)$values
              bits <- c(bits, paste0("dotted: ", group_concordance_threshold, " EU/mL ELISA (vert) / ",
                                     group_concordance_threshold, " U/mL titer (horiz)"))
              p <- ggplot(fr, aes(x = x, y = y)) +
                geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "red", linewidth = 0.4) +
                geom_vline(xintercept = thx, linetype = "dotted", colour = "#1f78b4", linewidth = 0.5) +
                geom_hline(yintercept = thy, linetype = "dotted", colour = "#33a02c", linewidth = 0.5) +
                geom_point(aes(colour = Concordance), size = point_size, alpha = 0.85)
              if(group_concordance_highlight_discordant) {
                p <- p + geom_point(data = subset(fr, Concordance == "Discordant"), colour = "red", size = point_size * 0.6, alpha = 0.9)
              }
              if(gc_reg_group && nrow(fr) >= 3) {
                p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black",
                                     fill = "grey70", alpha = 0.3, linewidth = 0.8, na.rm = TRUE)
              }
              if(gc_reg_all && nrow(fr_all) >= 3) {
                p <- p + geom_smooth(data = fr_all, aes(x = x, y = y), method = "lm", formula = y ~ x,
                                     se = FALSE, colour = "grey25", linewidth = 0.8,
                                     linetype = if(gc_reg_group) "dashed" else "solid",
                                     inherit.aes = FALSE, na.rm = TRUE)
              }
              p <- p +
                scale_colour_manual(values = gc_conc_colors, drop = FALSE, name = "Concordance") +
                 labs( x = gc_x_lab(fr$ELISA_Type[1]), y = gc_y_lab(ref_name, ref_unit)) +
                 theme_publication_custom(base_font_size) +
                 guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
                 theme(legend.position = "bottom",
                       legend.title = element_text(size = base_font_size - 1),
                       legend.text = element_text(size = base_font_size - 2),
                       legend.spacing.x = unit(0.08, "cm"),
                       legend.spacing.y = unit(0.02, "cm"),
                       legend.key.spacing.x = unit(0.04, "cm"),
                       legend.key.spacing.y = unit(0.01, "cm"),
                       legend.key.width = unit(0.5, "cm"),
                       legend.key.height = unit(0.18, "cm"),
                       legend.margin = margin(0, 0, 0, 0, "pt"),
                       legend.box.spacing = unit(0, "cm"),
                       legend.box.margin = margin(0, 0, 0, 0))
              if(!is.null(lims)) p <- p + coord_cartesian(xlim = lims$x, ylim = lims$y)
              p
            }
            
            gc_reg_note <- switch(gc_reg,
                                  per_group = "regression: individual line per group",
                                  all_data  = "regression: single line over all plotted datapoints",
                                  both      = "regression: per-group lines + pooled all-data line (dashed)")
            
             gc_stats_ref <- data.frame()
             gc_per_sample_ref <- data.frame()
              for(m in gc_methods) {
                for(i in seq_along(gc_combos)) {
                  cb <- gc_combos[[i]]
                  if(cb$ref != gc_ref_up) next
                  fr <- gc_data[[m]][[i]]
                  if(is.null(fr) || nrow(fr) == 0) next
                  for(g in setdiff(sort(unique(fr$Group)), "Ungrouped")) {
                    sub <- fr[fr$Group == g, , drop = FALSE]
                    if(nrow(sub) < 2) next
                    gc_stats_ref <- rbind(gc_stats_ref, gc_stat_row(sub, m, cb$elisa, cb$rep, g))
                  }
                  gc_stats_ref <- rbind(gc_stats_ref, gc_stat_row(fr, m, cb$elisa, cb$rep, "All groups"))
                  
                  if(nrow(fr) > 0) {
                    log_titer <- log10(fr$Titer)
                     tol_lower <- 10^(log_titer - ed50_tolerance)
                     tol_upper <- 10^(log_titer + ed50_tolerance)
                    interval_concordant <- fr$ELISA >= tol_lower & fr$ELISA <= tol_upper
                    interval_concordant[is.na(log_titer) | is.na(fr$ELISA)] <- NA
                    per_sample <- data.frame(
                      Sample = fr$Sample,
                      Group = fr$Group,
                      Method = m,
                      Replicate = cb$rep,
                      ELISA_Type = cb$elisa,
                      Reference = gc_ref_up,
                      Titer_Neut = fr$Titer,
                      ELISA = fr$ELISA,
                      Tolerance_Lower = tol_lower,
                      Tolerance_Upper = tol_upper,
                      Interval_Concordant = interval_concordant,
                      stringsAsFactors = FALSE
                    )
                    gc_per_sample_ref <- rbind(gc_per_sample_ref, per_sample)
                  }
                }
              }
             if(nrow(gc_stats_ref) > 0) {
               gc_stats <- rbind(gc_stats, gc_stats_ref)
             }
             if(nrow(gc_per_sample_ref) > 0) {
               gc_per_sample <- rbind(gc_per_sample, gc_per_sample_ref)
             }
           }
        }  # end for(gc_ref_up in gc_refs)
        
        if((isTRUE(group_concordance_composite) || isTRUE(group_concordance_individual)) && length(all_ref_frames) >= 1) {
          
          for(elisa_type in c("Linear", "Non-linear")) {
            if(elisa_type == "Non-linear" && !gc_nonlin_ok) next
            elisa_dir <- file.path(gc_out, "plots", tolower(elisa_type))
            
            elisa_idxs <- which(sapply(gc_combos, function(cb) cb$elisa == elisa_type))
            R1_idxs <- elisa_idxs[grep("^R1$", sapply(gc_combos[elisa_idxs], function(cb) cb$ref))]
            R2_idxs <- elisa_idxs[grep("^R2$", sapply(gc_combos[elisa_idxs], function(cb) cb$ref))]
            if(length(R1_idxs) < 2 || length(R2_idxs) < 2) next
            
            panel_idx <- c(R1_idxs[1], R2_idxs[1], R1_idxs[2], R2_idxs[2])
            panel_refs <- c("R1", "R2", "R1", "R2")
            panel_reps <- c("2rep", "2rep", "3rep", "3rep")
            panel_labels <- sapply(panel_idx, function(i) {
              cb <- gc_combos[[i]]
              paste0(cb$ref, " - ", cb$rep, " - ", tolower(elisa_type), " ELISA")
            })
            
            for(m in gc_methods) {
              frames <- lapply(panel_idx, function(idx) {
                r <- panel_refs[which(panel_idx == idx)][1]
                all_ref_frames[[r]][[m]][[idx]]
              })
              
              if(isTRUE(group_concordance_composite)) {
                lims <- gc_limits(frames)
                panels <- vector("list", 4)
                n_ok <- 0
                for(i in 1:4) {
                  fr <- frames[[i]]
                  rp <- ref_params[[panel_refs[i]]]
                  if(is.null(fr) || nrow(fr) < 2) {
                    panels[[i]] <- gc_blank(paste0(panel_labels[i], "\n(not available)"))
                  } else {
                    panels[[i]] <- gc_panel_allgroups(fr, panel_labels[i], lims, rp$name, rp$unit)
                    n_ok <- n_ok + 1
                  }
                }
                if(n_ok == 0) { ; next }
                
                panels_no_leg <- lapply(panels, function(p) {
                  if(inherits(p, "ggplot")) p + theme(legend.position = "none") else p
                })
                legend_grob <- NULL
                for(p in panels) {
                  if(inherits(p, "ggplot")) {
                    legend_grob <- gc_extract_legend(p)
                    if(!is.null(legend_grob)) break
                  }
                }
                if(is.null(legend_grob)) legend_grob <- gc_blank("")
                
                top_lab <- paste0("ELISA vs Titer concordance - all groups by Titer category (", m, "; ", elisa_type, ")",
                                  if(isTRUE(group_concordance_concordant_only)) " - concordant only" else "",
                                  "\n", gc_reg_note)
                arr <- gridExtra::arrangeGrob(
                  grobs = c(panels_no_leg, list(legend_grob)),
                  ncol = 2, nrow = 3,
                  layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
                  heights = unit(c(1, 1, 0.4), "null"),
                  top = grid::textGrob(top_lab, gp = grid::gpar(fontsize = title_size + 2, fontface = "bold"))
                )
                ggsave(file.path(elisa_dir, paste0("Concordance_AllGroups_", gc_safe(m), ".png")), arr,
                       width = min(48, group_concordance_width * 2),
                       height = min(48, group_concordance_height * 2 + 0.8),
                       dpi = group_concordance_dpi, limitsize = FALSE)
              }
              
              if(isTRUE(group_concordance_individual)) {
                lims <- gc_limits(frames)
                grp_present <- sort(unique(unlist(lapply(frames, function(fr) if(is.null(fr)) NULL else fr$Group))))
                grp_present <- setdiff(grp_present, "Ungrouped")
                if(length(grp_present) == 0) { ; next }
                panels <- list(); n_ok <- 0
                for(g in grp_present) {
                  for(i in 1:4) {
                    fr <- frames[[i]]
                    rp <- ref_params[[panel_refs[i]]]
                    ttl <- paste0(g, "\n", panel_labels[i])
                    if(is.null(fr) || nrow(fr) == 0) {
                      pp <- gc_blank(paste0(ttl, "\n(not available)"))
                    } else {
                      pp <- gc_panel_group(fr, g, ttl, lims, rp$name, rp$unit)
                    }
                    if(is.null(pp)) {
                      panels[[length(panels) + 1]] <- gc_blank(paste0(ttl, "\n(not available)"))
                    } else {
                      panels[[length(panels) + 1]] <- pp; n_ok <- n_ok + 1
                    }
                  }
                }
                if(n_ok == 0) next
                ncol_arr <- 4; nrow_arr <- length(grp_present)
                
                panels_no_leg <- lapply(panels, function(p) {
                  if(inherits(p, "ggplot")) p + theme(legend.position = "none") else p
                })
                legend_grob <- NULL
                for(p in panels) {
                  if(inherits(p, "ggplot")) {
                    legend_grob <- gc_extract_legend(p)
                    if(!is.null(legend_grob)) break
                  }
                }
                if(is.null(legend_grob)) legend_grob <- gc_blank("")
                
                top_lab <- paste0("ELISA vs Titer concordance by Titer category (", m, "; ", elisa_type, ")",
                                  if(isTRUE(group_concordance_concordant_only)) " - concordant only" else "",
                                  "\n", gc_reg_note)
                arr <- gridExtra::arrangeGrob(
                  grobs = c(panels_no_leg, list(legend_grob)),
                  ncol = ncol_arr, nrow = nrow_arr + 1,
                  layout_matrix = rbind(matrix(1:(nrow_arr * ncol_arr), nrow = nrow_arr, byrow = TRUE),
                                        rep(nrow_arr * ncol_arr + 1, ncol_arr)),
                  heights = unit(c(rep(1, nrow_arr), 0.35), "null"),
                  top = grid::textGrob(top_lab, gp = grid::gpar(fontsize = title_size + 2, fontface = "bold"))
                )
                ggsave(file.path(elisa_dir, paste0("Concordance_ByGroup_", gc_safe(m), ".png")), arr,
                       width = min(48, group_concordance_width * ncol_arr),
                       height = min(48, group_concordance_height * nrow_arr + 1),
                       dpi = group_concordance_dpi, limitsize = FALSE)
                          ncol_arr, "replicate/ELISA panels")
              }
            }
          }
        }
        
         if(!is.null(gc_stats) && is.data.frame(gc_stats) && nrow(gc_stats) > 0) {
           gc_stats <- gc_stats[order(gc_stats$Reference, gc_stats$Method, gc_stats$ELISA, gc_stats$Replicate, gc_stats$Group), ]
           gc_per_group <- gc_stats[gc_stats$Group != "All groups", , drop = FALSE]
           gc_by_group <- as.data.frame(
             gc_per_group %>%
               group_by(Reference, Method, Group) %>%
               summarise(Strata = dplyr::n(),
                         Total_N = sum(N, na.rm = TRUE),
                         Mean_R_squared = round(mean(R_squared, na.rm = TRUE), 4),
                         Min_P_value = if(all(is.na(P_value))) NA_real_ else signif(min(P_value, na.rm = TRUE), 4),
                         Mean_Concordance_Rate_pct = round(mean(Concordance_Rate_pct, na.rm = TRUE), 1),
                         .groups = "drop")
           )
           excel_sheets <- list(Stratified_By_Group = gc_stats, Summary_By_Group = gc_by_group)
           if(!is.null(gc_per_sample) && is.data.frame(gc_per_sample) && nrow(gc_per_sample) > 0) {
             gc_per_sample <- gc_per_sample[order(gc_per_sample$Reference, gc_per_sample$Method, gc_per_sample$Sample), ]
             excel_sheets$Per_Sample_Calculations <- gc_per_sample
           }
           if(requireNamespace("writexl", quietly = TRUE)) {
             writexl::write_xlsx(excel_sheets,
                                 file.path(gc_out, "grpconcrd_stat.xlsx"))
           }
         }
      }  # end gate else (length(gc_refs) == 0 ...)
    }, error = function(e) ))
  }
   if(!is.null(gc_stats) && is.data.frame(gc_stats) && nrow(gc_stats) > 0) {
     result_list$grpconcrd_stat <- gc_stats
   }
   if(!is.null(gc_per_sample) && is.data.frame(gc_per_sample) && nrow(gc_per_sample) > 0) {
     result_list$group_concordance_per_sample <- gc_per_sample
   }
   
   invisible(result_list)
}


















