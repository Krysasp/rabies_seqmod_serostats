library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(writexl)

source(file.path(Sys.getenv("CATCH_ROOT", "."), "scripts", "_core", "obfuscated_core.R"))
catch_source("scripts/utils/helpers.R")

analyze_predict_elisa_category <- function(file_path, sheet_name = "sheetNum", output_file = NULL) {
  .run_anova <- function(group_var, group_name) {
    results <- lapply(metrics, function(m) {
      df_sub <- df %>% filter(!is.na(.data[[group_var]]), !is.na(.data[[m]]))
      if(nrow(df_sub) < 2 || length(unique(df_sub[[group_var]])) < 2) {
        return(NULL)
      }
      fml <- as.formula(paste(m, '~', group_var))
      fit <- aov(fml, data = df_sub)
      tidy_fit <- tidy(fit) %>% filter(term != '(Intercept)') %>% mutate(Metric = m)
      return(tidy_fit)
    })
    bind_rows(results) %>%
      select(Metric, everything()) %>%
      mutate(Group = group_name)
  }
  
  .run_summary <- function(group_var, group_name) {
    df %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        N_Samples = n(),
        N_With_Reference_ELISA = sum(!is.na(Reference_Clean)),
        N_Correctly_Predicted = sum(Predicted_Clean == Reference_Clean, na.rm = TRUE),
        Probabilistic_Concordance_Rate_percent = round(mean(Predicted_Clean == Reference_Clean, na.rm = TRUE) * 100, 1),
        Mean_Confidence = round(mean(Pred_Coverage, na.rm = TRUE), 1),
        N_MultiCategory = sum(N_Probable > 1, na.rm = TRUE),
        Pct_MultiCategory = round(mean(N_Probable > 1, na.rm = TRUE) * 100, 1),
        N_Uncertain = sum(Pred_Coverage < 50, na.rm = TRUE),
        Pct_Uncertain = round(mean(Pred_Coverage < 50, na.rm = TRUE) * 100, 1),
        Mean_Titer_EU_mL = round(mean(Measured_Titer, na.rm = TRUE), 4),
        .groups = 'drop'
      )
  }
  
  .save_results <- function(anova_ref_df, anova_pred_df, anova_nprob_df, anova_all, summary_ref, summary_pred, summary_nprob) {
    if(!is.null(output_file)) {
      write_xlsx(list(
        ANOVA_Reference = anova_ref_df,
        ANOVA_Predicted = anova_pred_df,
        ANOVA_NProbable = anova_nprob_df,
        ANOVA_Combined = anova_all,
        Summary_Reference = summary_ref,
        Summary_Predicted = summary_pred,
        Summary_NProbable = summary_nprob
      ), path = output_file)
      message("Results saved to: ", output_file)
    }
  }
  
  df <- .load_and_prepare_data(file_path, sheet_name)
  metrics <- c('Pred_Coverage', 'Ref_Coverage', 
               'N_Probable', 'Measured_Titer', 'ELISA_Titer_val')
  
  anova_ref_df <- .run_anova('Reference_Clean', 'Reference Category')
  anova_pred_df <- .run_anova('Predicted_Clean', 'Predicted Category')
  anova_nprob_df <- .run_anova('N_Probable', 'N Probable Categories')
  anova_all <- bind_rows(anova_ref_df, anova_pred_df, anova_nprob_df)
  
  summary_ref <- .run_summary('Reference_Clean', 'Reference Category')
  summary_pred <- .run_summary('Predicted_Clean', 'Predicted Category')
  summary_nprob <- .run_summary('N_Probable', 'N Probable Categories')
  
  .save_results(anova_ref_df, anova_pred_df, anova_nprob_df, anova_all, summary_ref, summary_pred, summary_nprob)
  
  list(
    anova_combined = anova_all,
    summary_reference = summary_ref,
    summary_predicted = summary_pred,
    summary_nprobable = summary_nprob
  )
}

.load_and_prepare_data <- function(file_path, sheet_name) {
  df <- read_excel(file_path, sheet = sheet_name)
  
  df <- df %>%
    mutate(
      Tol_Lower = as.numeric(sub(' .*', '', `Tolerance range`)),
      Tol_Upper = as.numeric(sub('.*- ', '', `Tolerance range`)),
      is_concordant = `ELISA Titer` >= Tol_Lower & `ELISA Titer` <= Tol_Upper,
      Reference_Clean = ifelse(`Reference Category (U/mL)` == 'Discordant', NA, `Reference Category (U/mL)`),
      Predicted_Clean = `Predicted Category (U/mL)`,
      Pred_Coverage = `Predicted Interval Coverage (%)`,
      Ref_Coverage = `Reference Interval Coverage (%)`,
      N_Probable = `N Probable Categories`,
      Measured_Titer = `Measured Titer`,
      ELISA_Titer_val = `ELISA Titer`
    )
  
  df
}
