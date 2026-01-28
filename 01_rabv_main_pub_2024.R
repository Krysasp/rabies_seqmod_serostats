### SCRIPTS FOR SARAWAK RABV SEROLOGICAL ANALYSIS AND LARGE-SCALE RABV ANALYSIS ###
### Date created: 2022-10-15 ###
### Created by: JONATHAN CHAN ###
### This is the general packages to be loaded for RABV PHYLO WORKS
suppressPackageStartupMessages({
library(treeio)
library(ggtree)
library(ggplot2)
library(sf)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(ggrepel)
library(ggforce)
library(readr)
library(dplyr)
library(scales)
library(grid)
library(gridExtra)
library(viridis)
library(colorspace)
#library(rhierbaps)
library(ape)
library(repr)
library(padr)
library(reshape)
library(Seurat)
library(patchwork)
library(stringr)
library(mgcv)
library(RColorBrewer)
library(patchwork)
library(ggbeeswarm)
})
################################### 1. TESTING: HIERBAPS SNP GENERATION FROM ALIGNED FASTA #############################
snp.matrix <- load_fasta("01_src/swk_rabv_alignred_pub.fasta")
hb.results <- hierBAPS(snp.matrix, max.depth = 3, n.pops = 20, quiet = TRUE)
#write.csv(hb.results$partition.df, "output/rabv_hierbaps_out.csv",row.names = F)
head(hb.results$partition.df)
region_cols <- c("Asian"='#F2A922',"Cosmopolitan"='#267d56',"Arctic-related"='#1bbbd0',"Africa"='#f55a76',"Indian_subcontinent"='#c3699d')
################################### 2A. SEROLOGY STUDIES #############################
################ 2A. SEROLOGICAL STUDIES: AGE, VAX TYPES, DOG COUNT ANALYSIS ################
{
sero_df<- read.csv(file = "01_src/rabv_serostat.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
month_labels <- c("1" = "Month-1", "2" = "Month-2", "3" = "Month-3", "4" = "Month-4",
                  "5" = "Month-5", "6" = "Month-6", "7" = "Month-7", "8" = "Month-8",
                  "9" = "Month-9", "10" = "MOnth-10", "11" = "Month-11", "12" = "Month-12", "13" = "Month-13", 
                  "14" = "Month-14", "15" = "Month-15", "16" = "Month-16", "17" = "Month-17", "18" = "Month-18")
age_levels <- c("≥3-6m",">7m-12m", "≥12m-24m", "≥25m-<6y", "≥6y", "Unknown")
sero_df$age_G2 <- factor(sero_df$age_G2, levels = age_levels)
summary_data <- sero_df %>%
  group_by(age_G2, dose_type) %>%
  summarise(count = n(), unique_count = n_distinct(comb), .groups = "drop") %>%
  complete(
    age_G2 = age_levels,
    dose_type = c("primary", "booster"),
    fill = list(count = 0, unique_count = 0)
  )
summary_data$age_G2 <- factor(summary_data$age_G2, levels = age_levels)
# Plot grouped bar chart
page <- ggplot(summary_data, aes(x = age_G2, y = unique_count, fill = dose_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7, color = "black") +
  
  # Add count labels on top of bars
  geom_text(aes(label = paste0("n=", unique_count)),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 3.8) +
  scale_y_continuous(limits = c(0, 40), expand = expansion(mult = c(0, 0.05))) +
  # Labels and titles
  labs(
    title = "Age distribution of dogs at the \n time of first blood collection post-vaccination",
    x = "Age group",
    y = "Number of dogs",
    fill = "Vaccination status"
  ) +
  
  # Styling: axis titles bold and larger
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(size = 1, colour = "black")
  ) +
  
  # Ensure plot borders are thick and black
  theme(panel.grid = element_blank())
ggsave("01-rvsero_ibet_dogdemo_116dogs.png", page, width = 9, height = 5)
# Format 'month' as an ordered factor
month_group_levels <- c("1", "2-3", "4-5", "6-7", "8-9", "10-11", "12-13", "14-15", "16-17", "18")
month_group_all <- c("1", "2","3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
sero_df$month_group <- cut(
  sero_df$month,
  breaks = c(0, 1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 18.5),
  labels = month_group_levels,
  right = TRUE,
  include.lowest = TRUE
)
sero_df$month_group <- factor(sero_df$month_group, levels = month_group_levels)
}
################ 2B. SEROLOGICAL STUDIES: SCATTERPLOTS OF GROUPED DATA POINTS (CONTINUES WITH sero_df FROM 2A.) ################
### 2B1. MAIN ANALYSIS CODES: FUNCTIONS & DATA MANAGEMENT ###
suppressPackageStartupMessages({
# Create unique cell IDs for Seurat object
sero_df$cell_id <- make.unique(as.character(sero_df$no))
# Create dummy expression matrix (required for Seurat)
cell_ids <- sero_df$cell_id
dummy_counts <- matrix(rpois(100 * length(cell_ids), lambda = 10),
                       nrow = 308,
                       ncol = length(cell_ids),
                       dimnames = list(paste0("Gene", 1:308), cell_ids))
# Create Seurat object with metadata
seu <- CreateSeuratObject(counts = dummy_counts, meta.data = sero_df)
seu <- RenameCells(seu, new.names = sero_df$cell_id)
# Define log2 y-axis breaks and labels
log2_breaks <- c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
log2_labels <- as.character(log2_breaks)
log2_labels[1] <- "0.0"
log2_labels[6] <- "1.0"
log2_labels[7] <- "2.0"
log2_labels[8] <- "4.0"
log2_labels[9] <- "8.0"
log2_labels[10] <- "16.0"
log2_labels[11] <- "32.0"
# Function to compute geometric mean
geo_mean <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
geo_sd_factor <- function(x) {
  x <- x[x > 0]  # Remove zeros (log not defined)
  exp(sd(log(x), na.rm = TRUE))
}
brown_to_grey_palette <- colorRampPalette(c("#c69e3b", "#cdaa54",  "#d4b66c", "#dbc285", "#e2ce9d", "#f3ebd7"))  # dark brown to medium grey
month_colors <- setNames(
  brown_to_grey_palette(length(month_group_levels)),
  month_group_levels
)
month_colors_max <- setNames(
  brown_to_grey_palette(18),
  month_group_all
)
})
### 2B2. MAIN VIZUALIZE CODES: VIOLINPLOTS ###
{
plot_violin1m <- function(seu_obj, dose_type_filter = NULL, conc_vax_filter = NULL, set_type_filter = NULL, age_G2_filter = NULL, plot_title = NULL) {
  filter_expr <- rep(TRUE, ncol(seu_obj))  # Start with all TRUE
  if (!is.null(dose_type_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$dose_type == dose_type_filter)
  } if (!is.null(conc_vax_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$conc_vax == conc_vax_filter)
  } if (!is.null(set_type_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$set_type == set_type_filter)
  } if (!is.null(age_G2_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$age_G2 == age_G2_filter)
  }
  cells <- WhichCells(seu_obj, cells = colnames(seu_obj)[filter_expr])
  if (length(cells) == 0) {
    stop("No cells match the provided filter criteria.")}
  seu_sub <- subset(seu_obj, cells = cells)
  meta <- seu_sub@meta.data
  meta$month <- factor(meta$month, levels = 1:18)
  seu_sub@meta.data$month <- meta$month
  cat("Number of rows in seu_obj:", nrow(seu_obj), "\n")
  cat("Number of rows in meta:", nrow(meta), "\n")
  gm_df <- meta %>%
    group_by(month) %>%
    summarise(n = n(),
              geo_mean = geo_mean(adjusted),
              geo_sd_factor = geo_sd_factor(adjusted),
              .groups = "drop") %>%
    mutate(month_num = as.numeric(as.character(month)))
  #vline_positions <- (unique(gm_df$month_num) %>% sort()) + 0.5
  #vline_positions <- vline_positions[vline_positions < max(gm_df$month_num) + 0.5]
  #total_groups <- length(month)
  vline_positions <- seq(1.5, 18 - 0.5, by = 1)
VlnPlot(
    seu_sub,
    features = "adjusted",
    group.by = "month",
    pt.size = 1.9,
    alpha = 0.01,
    cols = rep("#F4F4F4", length(unique(seu_sub$month))))
p <-  ggplot() + geom_quasirandom(
  data = meta,
  aes(x = month, y = adjusted),
  width = 0.3,
  varwidth = TRUE,
  size = 2.2,
  alpha = 0.6,
  color = "#323232",
  method = "smiley",  # or "quasirandom", "pseudorandom"
  inherit.aes = FALSE) +
  ggtitle(NULL) + xlab('months since last vaccination') + ylab('anti-RABV-G IgG (EU/mL')+
    geom_segment(
      data = gm_df,
      aes(x = as.numeric(month) - 0.4, 
          xend = as.numeric(month) + 0.4,
          y = geo_mean,
          yend = geo_mean),
      color = "black",
      size = 1.2,
      inherit.aes = FALSE) +
    geom_vline(
      xintercept = vline_positions,
      linetype = "dashed",
      size = 0.5,
      color = "#191919") +
    scale_y_continuous(
      trans = "log2",
      breaks = log2_breaks,
      labels = log2_labels,
      limits = c(0.03125, 32.0)) +
    scale_x_discrete(drop = FALSE) +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.8, color = "black") +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(vjust = 0.5, size = 10, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none") +
    guides(color = guide_legend(nrow = 1))
  # Create summary table
  summary_df <- gm_df %>%
    mutate(
      # Convert Sample Size to integer and then to character (no decimals)
      `Sample Size (n)` = as.character(as.integer(round(n))),
      # Format Geometric Mean and Geometric SD Factor as strings with 3 decimals
      `Geometric Mean` = sprintf("%.3f", geo_mean),
      `Geometric SD Factor` = sprintf("%.3f", geo_sd_factor)
    ) %>%
    select(`Month` = month, `Sample Size (n)`, `Geometric Mean`, `Geometric SD Factor`) %>%
    pivot_longer(
      cols = -Month,
      names_to = "Months",
      values_to = "Value") %>%
    pivot_wider(
      names_from = Month,
      values_from = Value)
  my_theme <- ttheme_minimal(
    core = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")),
    rowhead = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")))
  table_grob <- tableGrob(summary_df, rows = NULL, theme = my_theme)
  for (i in seq_len(nrow(table_grob$layout))) {
    l <- table_grob$layout$l[i]
    r <- table_grob$layout$r[i]
    t <- table_grob$layout$t[i]
    b <- table_grob$layout$b[i]
    table_grob <- gtable::gtable_add_grob(
      table_grob,
      grobs = rectGrob(gp = gpar(col = "black", lwd = 1.0, fill = NA)),
      t = t, b = b, l = l, r = r
    )}
  combined_plot <- p / table_grob + plot_layout(heights = c(4, 1))
  if (!is.null(plot_title)) {
    combined_plot <- combined_plot +
      plot_annotation(
        title = plot_title,
        theme = theme(plot.title = element_text(size = 15, hjust = 0.5)) )}
  return(combined_plot)
  }
plot_violin2m <- function(seu_obj, dose_type_filter = NULL, conc_vax_filter = NULL, set_type_filter = NULL, age_G2_filter = NULL, plot_title = NULL) {
  # Filtering
  filter_expr <- rep(TRUE, ncol(seu_obj))
  if (!is.null(dose_type_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$dose_type == dose_type_filter)}
  if (!is.null(conc_vax_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$conc_vax == conc_vax_filter)}
  if (!is.null(set_type_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$set_type == set_type_filter)}
  if (!is.null(age_G2_filter)) {
    filter_expr <- filter_expr & (seu_obj@meta.data$age_G2 == age_G2_filter)}
  
  cells <- WhichCells(seu_obj, cells = colnames(seu_obj)[filter_expr])
  if (length(cells) == 0) {
    stop("No cells match the provided filter criteria.")}
  seu_sub <- subset(seu_obj, cells = cells)
  meta <- seu_sub@meta.data
  meta$month_group <- factor(meta$month_group, levels = month_group_levels)
  seu_sub@meta.data$month_group <- meta$month_group
  cat("Number of rows in seu_obj:", nrow(seu_obj), "\n")
  cat("Number of rows in meta:", nrow(meta), "\n")
  
  # Summary statistics
  gm_df <- meta %>%
    group_by(month_group) %>%
    summarise(
      n = n(),
      geo_mean = geo_mean(adjusted),
      geo_sd_factor = geo_sd_factor(adjusted),
      .groups = "drop"
    ) %>%
    mutate(month_group_num = as.numeric(factor(month_group, levels = month_group_levels)))
  #vline_positions <- (unique(gm_df$month_group_num) %>% sort()) + 0.5
  #vline_positions <- vline_positions[vline_positions < max(gm_df$month_group_num) + 0.5]
  total_groups <- length(month_group_levels)
  vline_positions <- seq(1.5, total_groups - 0.5, by = 1)
VlnPlot(
    seu_sub,
    features = "adjusted",
    group.by = "month_group",
    pt.size = 1.9,
    alpha = 0.01,
    cols = rep("#F4F4F4", length(unique(month_group_levels)))) 
  p <-  ggplot() + geom_quasirandom(
    data = meta,
    aes(x = month_group, y = adjusted),
    width = 0.3,
    varwidth = TRUE,
    size = 2.2,
    alpha = 0.6,
    color = "#323232",
    method = "smiley",  # or "quasirandom", "pseudorandom"
    inherit.aes = FALSE) +
    ggtitle(NULL) + xlab('months since last vaccination') + ylab('anti-RABV-G IgG (EU/mL')+
    geom_smooth(
      data = meta,
      aes(x = as.numeric(month_group), y = adjusted),
      method = "gam",
      formula = y ~ s(x),
      color = "black",
      fill = "grey70",
      size = 1,
      alpha = 0.3,
      se = TRUE,
      inherit.aes = FALSE) +
    geom_vline(
      xintercept = vline_positions,
      linetype = "dashed",
      size = 0.5,
      color = "#191919") +
    geom_segment(
      data = gm_df,
      aes(x = month_group_num - 0.4,
          xend = month_group_num + 0.4,
          y = geo_mean,
          yend = geo_mean),
      color = "black",
      size = 1.2,
      inherit.aes = FALSE) +
    scale_y_continuous(
      trans = "log2",
      breaks = log2_breaks,
      labels = log2_labels,
      limits = c(0.03125, 32)) +
    scale_x_discrete(drop = FALSE) +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.8, color = "black") +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(vjust = 0.5, size = 10, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "none") +
    guides(color = guide_legend(nrow = 1))
  summary_df <- gm_df %>%
    mutate(
      # Convert Sample Size to integer and then to character (no decimals)
      `Sample Size (n)` = as.character(as.integer(round(n))),
      # Format Geometric Mean and Geometric SD Factor as strings with 3 decimals
      `Geometric Mean` = sprintf("%.3f", geo_mean),
      `Geometric SD Factor` = sprintf("%.3f", geo_sd_factor)) %>%
    select(`Month` = month_group, `Sample Size (n)`, `Geometric Mean`, `Geometric SD Factor`) %>%
    pivot_longer(
      cols = -Month,
      names_to = "Months",
      values_to = "Value") %>%
    pivot_wider(
      names_from = Month,
      values_from = Value)
  # Table theme
  my_theme <- ttheme_minimal(
    core = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")),
    rowhead = list(
      fg_params = list(fontface = "bold", fontsize = 10, hjust = 0.5, x = 0.5),
      padding = unit(c(1, 2), "mm")))
  table_grob <- tableGrob(summary_df, rows = NULL, theme = my_theme)
  for (i in seq_len(nrow(table_grob$layout))) {
    l <- table_grob$layout$l[i]
    r <- table_grob$layout$r[i]
    t <- table_grob$layout$t[i]
    b <- table_grob$layout$b[i]
    table_grob <- gtable::gtable_add_grob(
      table_grob,
      grobs = rectGrob(gp = gpar(col = "black", lwd = 1.0, fill = NA)),
      t = t, b = b, l = l, r = r
    )}
  combined_plot <- p / table_grob + plot_layout(heights = c(4, 1))
  if (!is.null(plot_title)) {
    combined_plot <- combined_plot +
      plot_annotation(
        title = plot_title,
        theme = theme(plot.title = element_text(size = 15, hjust = 0.5)))}
  return(combined_plot)}
# Generate plots for both dose types
df_bydose_pr <- sero_df %>%  filter(dose_type == "primary")
df_bydose_bo <- sero_df %>%  filter(dose_type == "booster")}
### 2B3. CODES FOR GENERATION OF PUBLICATION-READY PLOTS ###
{
pa1 <- plot_violin1m(seu, dose_type_filter = "primary",
                plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected", 
                                    " (N = ", length(unique(df_bydose_pr$comb)), " dogs)")))
pa2 <- plot_violin1m(seu, dose_type_filter = "booster",
                plot_title = (paste("Anti-RABV-G IgG titers of vaccine-experienced dogs after RABV \n vaccination with ≥ 1 blood samples collected", 
                                    " (N = ", length(unique(df_bydose_bo$comb)), " dogs)")))
pa3 <- plot_violin2m(seu, dose_type_filter = "primary",
                plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected", 
                                    " (N = ", length(unique(df_bydose_pr$comb)), " dogs)")))
pa4 <- plot_violin2m(seu, dose_type_filter = "booster",
                plot_title = (paste("Anti-RABV-G IgG titers of vaccine-experienced dogs after RABV \n vaccination with ≥ 1 blood samples collected", 
                                    " (N = ", length(unique(df_bydose_bo$comb)), " dogs)")))

df_bydose_age1 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥3 - 6m")
pa5 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥3 - 6m",
                     plot_title = (paste("Anti-RABV-G IgG titers of ≥3 - 6m old vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_age1$comb)), " dogs)")))
df_bydose_age4 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == ">7m-12m")
pa6 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter=">7m-12m",
                     plot_title = (paste("Anti-RABV-G IgG titers of >7m-12m old vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_age4$comb)), " dogs)")))
df_bydose_age2 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥12m-24m")
pa7 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥12m-24m",
                     plot_title = (paste("Anti-RABV-G IgG titers of ≥12m-24m old vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_age2$comb)), " dogs)")))
df_bydose_age3 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥25m - <6y")
pa8 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥25m - <6y",
                     plot_title = (paste("Anti-RABV-G IgG titers of ≥25m - <6y old vaccine-naive dogs after RABV \n vaccination with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_age3$comb)), " dogs)")))
df_bydosetp_age1 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥3 - 6m")%>%filter(set_type == "multi")
pa9 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥3 - 6m", set_type_filter = "multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of ≥3 - 6m old vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosetp_age1$comb)), " dogs)")))
df_bydosetp_age4 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == ">7m-12m")%>%filter(set_type == "multi")
pa10 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter=">7m-12m", set_type_filter = "multi",
                      plot_title = (paste("Anti-RABV-G IgG titers of >7m-12m old vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                          " (N = ", length(unique(df_bydosetp_age4$comb)), " dogs)")))
df_bydosetp_age2 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥12m-24m")%>%filter(set_type == "multi")
pa11 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥12m-24m", set_type_filter = "multi",
                      plot_title = (paste("Anti-RABV-G IgG titers of ≥12m-24m old vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                          " (N = ", length(unique(df_bydosetp_age2$comb)), " dogs)")))
df_bydosetp_age3 <- sero_df %>% filter(dose_type == "primary") %>% filter(age_G2 == "≥25m - <6y")%>%filter(set_type == "multi")
pa12 <- plot_violin2m(seu, dose_type_filter = "primary", age_G2_filter="≥25m - <6y", set_type_filter = "multi",
                      plot_title = (paste("Anti-RABV-G IgG titers of ≥25m - <6y old vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                          " (N = ", length(unique(df_bydosetp_age3$comb)), " dogs)")))

df_bydosebytp_pr <- sero_df %>%  filter(dose_type == "primary") %>%  filter(set_type == "multi")
df_bydosebytp_bo <- sero_df %>%  filter(dose_type == "booster") %>%  filter(set_type == "multi")
pb1 <- plot_violin1m(seu, dose_type_filter = "primary", set_type_filter="multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosebytp_pr$comb)), " dogs)")))
pb2 <- plot_violin1m(seu, dose_type_filter = "booster", set_type_filter="multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-experienced dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosebytp_bo$comb)), " dogs)")))
pb3 <- plot_violin2m(seu, dose_type_filter = "primary", set_type_filter="multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosebytp_pr$comb)), " dogs)")))
pb4 <- plot_violin2m(seu, dose_type_filter = "booster", set_type_filter="multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-experienced dogs after RABV \n vaccination with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosebytp_bo$comb)), " dogs)")))

df_bydose_rab <- sero_df %>% filter(dose_type == "primary") %>% filter(conc_vax == "Rabisin")
pb5 <- plot_violin2m(seu, dose_type_filter = "primary", conc_vax_filter="Rabisin", 
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Rabisin dose with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_rab$comb)), " dogs)")))
df_bydose_rabb <- sero_df %>% filter(dose_type == "booster") %>% filter(conc_vax == "Rabisin")
pb6 <- plot_violin2m(seu, dose_type_filter = "booster", conc_vax_filter="Rabisin",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Rabisin dose with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_rabb$comb)), " dogs)")))
df_bydose_nob <- sero_df %>% filter(dose_type == "primary") %>% filter(conc_vax == "Nobivac")
pb7 <- plot_violin2m(seu, dose_type_filter = "primary", conc_vax_filter="Nobivac",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Nobivac dose with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_nob$comb)), " dogs)")))
df_bydose_nobb <- sero_df %>% filter(dose_type == "booster") %>% filter(conc_vax == "Nobivac")
pb8 <- plot_violin2m(seu, dose_type_filter = "booster", conc_vax_filter="Nobivac",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Nobivac dose with ≥ 1 blood samples collected",
                                         " (N = ", length(unique(df_bydose_nobb$comb)), " dogs)")))

df_bydosetp_rab <- sero_df %>% filter(dose_type == "primary") %>% filter(conc_vax == "Rabisin")%>%filter(set_type == "multi")
pb9 <- plot_violin2m(seu, dose_type_filter = "primary", conc_vax_filter="Rabisin",set_type_filter="multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Rabisin dose with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosetp_rab$comb)), " dogs)")))
df_bydosetp_rabb <- sero_df %>% filter(dose_type == "booster") %>% filter(conc_vax == "Rabisin")%>%  filter(set_type == "multi")
pb10 <- plot_violin2m(seu, dose_type_filter = "booster", conc_vax_filter="Rabisin", set_type_filter = "multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Rabisin dose with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosetp_rabb$comb)), " dogs)")))
df_bydosetp_nob <- sero_df %>% filter(dose_type == "primary") %>% filter(conc_vax == "Nobivac")%>%  filter(set_type == "multi")
pb11 <- plot_violin2m(seu, dose_type_filter = "primary", conc_vax_filter="Nobivac", set_type_filter = "multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Nobivac dose with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosetp_nob$comb)), " dogs)")))
df_bydosetp_nobb <- sero_df %>% filter(dose_type == "booster") %>% filter(conc_vax == "Nobivac")%>%  filter(set_type == "multi")
pb12 <- plot_violin2m(seu, dose_type_filter = "booster", conc_vax_filter="Nobivac", set_type_filter= "multi",
                     plot_title = (paste("Anti-RABV-G IgG titers of vaccine-naive dogs after received a \n Nobivac dose with ≥ 2 blood samples collected",
                                         " (N = ", length(unique(df_bydosetp_nobb$comb)), " dogs)")))
}
################ 2C. SEROLOGICAL STUDY: LONGITUDINAL ANTIBODY TITERS OVERTIME ################
### 2C1. MAIN SOURCE CODES: DEFINE NUMERICAL VALUES AND SOURCES ###
file_path<- read.csv(file = "01_src/rabv_serostat.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
### 2C2. MAIN ANALYSIS CODES: FUNCTIONS & DATA MANAGEMENT ###
library(tidyverse)
library(mgcv)
library(RColorBrewer)
library(patchwork)
library(grid)
library(gridExtra)
prepare_serology_data <- function(file_path, 
                                  filter_unassign = TRUE,
                                  month_range = c(1, 18),
                                  filter_group = NULL,
                                  filter_vaccine = NULL,
                                  filter_dose_type = NULL,
                                  filter_set_type = NULL,
                                  filter_tested_tp = NULL) {
  # Read and clean data
  raw_df <- read.csv(file = file_path, sep = ',', header = TRUE, 
                     stringsAsFactors = FALSE)
  
  sero_df <- raw_df %>%
    filter(!is.na(adjusted), !is.na(month), !is.na(comb), !is.na(group)) %>%
    mutate(
      month = as.numeric(month),
      adjusted_plot = ifelse(adjusted < 0.03125, 0.03125, adjusted)
    )
  
  # Optional filter for unassigned groups
  if (filter_unassign) {
    sero_df <- sero_df %>%
      filter(!grepl("unassign", group, ignore.case = TRUE))
  }
  
  # Apply additional filters if specified
  if (!is.null(filter_group)) {
    sero_df <- sero_df %>% filter(grepl(filter_group, group, ignore.case = TRUE))
  }
  
  if (!is.null(filter_vaccine)) {
    sero_df <- sero_df %>% filter(grepl(filter_vaccine, conc_vax, ignore.case = TRUE))
  }
  
  if (!is.null(filter_dose_type)) {
    sero_df <- sero_df %>% filter(grepl(filter_dose_type, dose_type, ignore.case = TRUE))
  }
  
  if (!is.null(filter_set_type)) {
    sero_df <- sero_df %>% filter(grepl(filter_set_type, set_type, ignore.case = TRUE))
  }
  
  if (!is.null(filter_tested_tp)) {
    sero_df <- sero_df %>% filter(!grepl(filter_tested_tp, tested_tp, ignore.case = TRUE))
  }
  
  # Filter by month range
  sero_df <- sero_df %>%
    filter(month >= month_range[1] & month <= month_range[2])
  
  return(sero_df)
}
##### 2C3. PLOTTING FUNCTIONS #####
# Common scale definitions
log2_breaks <- c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
log2_labels <- c("0.0", "0.0625", "0.125", "0.25", "0.5", "1.0", "2.0", "4.0", "8.0", "16.0", "32.0")
plot_serology <- function(data_group, 
                          label = NULL,
                          color_by = c("comb", "vaccine"),
                          show_title = FALSE,
                          show_gam = TRUE,
                          show_gridlines = TRUE,
                          x_limits = c(1, 18.5),
                          y_limits = c(0.0625, 32)) {
  
  color_by <- match.arg(color_by)
  group_name <- unique(data_group$group)
  unique_combs <- sort(unique(data_group$comb))
  n_comb <- length(unique_combs)
  
  # Prepare unique months and k_val
  unique_months <- length(unique(na.omit(data_group$month)))
  k_val <- max(3, min(10, unique_months - 1))
  
  # GAM model calculations (for title only)
  enough_data <- sum(!is.na(data_group$adjusted) & !is.na(data_group$month)) >= 3 && unique_months >= 2
  
  if (enough_data && show_gam) {
    gam_model <- tryCatch(
      mgcv::gam(adjusted ~ s(month, k = k_val), data = data_group, method = "REML", select = TRUE),
      error = function(e) NULL
    )
    
    if (!is.null(gam_model)) {
      gam_summary <- summary(gam_model)
      edf <- round(gam_summary$s.table[1, "edf"], 2)
      r_sq <- round(gam_summary$r.sq, 3)
      dev_expl <- round(gam_summary$dev.expl * 100, 1)
    } else {
      edf <- NA; r_sq <- NA; dev_expl <- NA
    }
  } else {
    gam_model <- NULL
    edf <- NA; r_sq <- NA; dev_expl <- NA
  }
  
  # Title (optional)
  title_text <- if (show_title) {
    paste0(
      "Serial anti-RABV-G IgG titers after primary dose, N=", n_comb, 
      " dog(s), \nGAM: s(month, k=", k_val, "), edf=", edf, 
      ", R²=", r_sq, ", DevExpl=", dev_expl, "%"
    )
  } else {
    NULL
  }
  
  # Prepare data
  data_group <- data_group %>% arrange(comb, month)
  
  # Color schemes
  if (color_by == "comb") {
    # Color palette for individual dogs
    full_palette <- brewer.pal(12, "Paired")
    no_yellow_indices <- setdiff(1:12, c(7, 9))
    filtered_palette <- full_palette[no_yellow_indices]
    palette_13 <- colorRampPalette(filtered_palette)(13)
    comb_colors <- setNames(palette_13[1:length(unique_combs)], unique_combs)
    color_var <- "comb"
    legend_title <- "Owner-dog identity"
  } else {
    # Color by vaccine type
    vaccine_colors <- c("Nobivac" = "#1c1c84", "Rabisin" = "red")
    color_var <- "conc_vax"
    legend_title <- "Vaccine type"
  }
  
  # Base theme selection
  if (show_gridlines) {
    base_theme <- theme_bw(base_size = 12)
  } else {
    base_theme <- theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
      )
  }
  
  # Base plot
  p <- ggplot(data_group, aes(x = month, y = adjusted_plot, 
                              group = comb, 
                              color = !!sym(color_var))) +
    geom_point(size = 2, alpha = 0.8) +
    geom_path(alpha = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.8, color = "black") +
    scale_y_continuous(
      trans = "log2",
      limits = y_limits,
      breaks = log2_breaks,
      labels = log2_labels
    ) +
    scale_x_continuous(
      breaks = 1:18,
      limits = x_limits,
      expand = expansion(mult = c(0.05, 0))
    ) +
    labs(
      x = paste0("months since last vaccination, ", group_name),
      y = "anti-RABG-G IgG titer (EU/mL)",
      title = title_text,
      color = legend_title
    ) +
    base_theme +
    theme(
      plot.title = if (show_title) element_text(size = 12, face = "bold", hjust = 0.5) else element_blank(),
      axis.text = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.5, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 1, l = 5)
    )
  
  # Add color scale
  if (color_by == "comb") {
    p <- p + scale_color_manual(values = comb_colors)
  } else {
    p <- p + scale_color_manual(values = vaccine_colors)
  }
  
  # Add GAM smooth if requested
  if (show_gam && enough_data) {
    p <- p + geom_smooth(
      aes(x = month, y = adjusted, group = 1),
      method = "gam",
      formula = y ~ s(x, k = k_val),
      method.args = list(select = TRUE),
      color = if (color_by == "comb") "#242424" else "#525252",
      size = 0.6,
      se = TRUE,
      inherit.aes = FALSE,
      data = data_group)
  }
  
  # Add panel label for publication
  if (!is.null(label)) {
    #p <- p + annotate(
    #  "text", x = -Inf, y = Inf, label = label,
    #  hjust = -0.2, vjust = 1.2,
    #  size = 6, fontface = "bold")
    p <- p + 
      coord_cartesian(clip = "off") +  # Allow labeling done outside of plot
      annotate(
        "text", 
        x = -Inf, 
        y = Inf, 
        label = paste0("  ", label, "  "),  # Add spaces for padding
        hjust = -0.2,      # Align left
        vjust = -0.5,      # Align to top
        size = 5.5,       # Slightly larger
        fontface = "bold",
        color = "black"
      ) +
      # Increase left and top margins to accommodate the label
      theme(plot.margin = margin(t = 30, r = 5, b = 5, l = 25, unit = "pt"))
  }
  
  return(p)
}

# Text box function
text_box_grob <- function(text, width = 120, lineheight = 1.0, fontsize = 1.0) {
  if (!is.null(width)) {
    wrapped_text <- strwrap(text, width = width)
    text_string <- paste(wrapped_text, collapse = "\n")
  } else {
    text_string <- text
  }
  
  grob <- grobTree(
    rectGrob(gp = gpar(fill = NA, col = "black", lwd = 0.7)),
    textGrob(
      label = text_string,
      x = 0.01, y = 1, just = c("left", "top"),
      gp = gpar(cex = fontsize, lineheight = lineheight)
    )
  )
  
  ggplot() +
    theme_void() +
    annotation_custom(grob) +
    theme(plot.margin = margin(0, 0, 0, 0))
}
##### 3. MULTI-PLOT GENERATION FUNCTION #####
generate_combined_plot <- function(data,
                                   group_column = "group",
                                   group_order = NULL,
                                   color_by = c("comb", "vaccine"),
                                   caption_text = NULL,
                                   ncol = 2,
                                   show_title = FALSE,
                                   show_gam = TRUE,
                                   show_gridlines = FALSE,
                                   save_path = NULL,
                                   width = 39,
                                   height = 39,
                                   units = "cm",
                                   caption_height = 0.2,
                                   caption_width = 120) {
  color_by <- match.arg(color_by)
  # Set group order
  if (!is.null(group_order)) {
    data[[group_column]] <- factor(data[[group_column]], levels = group_order)
  }
  # Split data by groups
  data_split <- data %>%
    split(.[[group_column]]) %>%
    keep(~ nrow(.x) > 0)  # Remove empty groups
  # Generate labels
  labels <- LETTERS[1:length(data_split)]
  # Create individual plots
  plots <- map2(data_split, labels, ~ plot_serology(.x, 
                                                    label = .y,
                                                    color_by = color_by,
                                                    show_title = show_title,
                                                    show_gam = show_gam,
                                                    show_gridlines = show_gridlines))
  # Combine plots
  combined_plot <- wrap_plots(plots, ncol = ncol)
  
  # Add caption if provided
  if (!is.null(caption_text)) {
    text_plot <- text_box_grob(caption_text, width = caption_width)
    final_plot <- combined_plot / text_plot +
      plot_layout(heights = c(1, caption_height))
  } else {
    final_plot <- combined_plot
  }
  
  # Save if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, 
           plot = final_plot, 
           width = width, 
           height = height, 
           units = units,
           limitsize = FALSE)
    message("Plot saved to: ", save_path)
  }
  return(final_plot)
}
##### 2C4. SPECIFIC PLOT CONFIGURATIONS #####
# Main data preparation
sero_df <- prepare_serology_data("01_src/rabv_serostat.csv")
# 4.1 Overall plot (all groups)
desired_order_overall <- c("below", "above", "decay", "anam1", "anam2", "anam3", 
                           "undesc", "B_up", "B_anam", "B_down")
overall_caption <- "Anti-RABV-G IgG trajectory after a dose of RABV vaccination in dogs that provided ≥ 3 sera, N=47dogs
A: Consistently below seroprotective level [n=13, 28%]
B: Consistently above seroprotective level [n=6, 13%]
C: Progressive waning to < seroprotective titer [n=13, 28%] 
D: Anamnestic response occurred within first 6 months of follow-up [n=5, 11%]
E: Anamnestic response occurred after first 6 months of follow-up [n=3, 6%]
F: Anamnestic response with high pre-exposure titer after over six months of follow-up [n=3, 6%]
G: Miscellaneous: Responses cannot be described as others [n=4, 8.5%]

Anti-RABV-G IgG trajectory after vaccine-experience dogs received RABV vaccination and provided ≥ 3 sera, N=9 dogs
H: Vaccine-experienced dogs: Above 0.5EU/mL [n=5, 56%]
I: Vaccine-experienced dogs: Anamnestic response [n=2, 22%]
J: Vaccine-experienced dogs: Below 0.5EU/mL[n=2, 22%]"
overall_plot <- generate_combined_plot(
  data = sero_df,
  group_column = "group",
  group_order = desired_order_overall,
  color_by = "comb",
  caption_text = overall_caption,
  show_title = FALSE,
  show_gam = TRUE,
  show_gridlines = TRUE,
  save_path = "00-overallAb_trajectories.png",
  width = 39,
  height = 39,
  caption_height = 0.23
)
# 4.2 Rabisin-only plots
rab_df <- prepare_serology_data("01_src/rabv_serostat.csv",
                                filter_vaccine = "Rab")
rab_order <- c("below", "above", "decay", "anam1", "anam2", "anam3", "undesc")
rab_caption <- "Longitudinal antibody trajectories after a dose of RABV vaccination in dogs that provided ≥ 3 sera, N=20 dogs 
A: Consistently below seroprotective level [n=6, 32% ]
B. Consistently above seroprotective level [n=3, 16%]
C: Progressive waning to < seroprotective titer [n=4, 20%]
D: Anamnestic response occurred within first 6 months of follow-up [n=3, 15%]
E: Anamnestic response occurred after first 6 months of follow-up [n=1, 5%]
F: Anamnestic response occurred after first 6 months of follow-up [n=2, 10%]
G: Miscellaneous: Response cannot be described as others [n=1, 5%]"
rab_plot <- generate_combined_plot(
  data = rab_df,
  group_column = "group",
  group_order = rab_order,
  color_by = "comb",
  caption_text = rab_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-rabisin_trajectories.png",
  width = 39,
  height = 39,
  caption_height = 0.17
)
# 4.3 Nobivac-only plots
nob_df <- prepare_serology_data("01_src/rabv_serostat.csv",
                                filter_vaccine = "Nob")
nob_order <- c("below", "above", "decay", "anam1", "anam2", "anam3", "undesc")
nob_caption <- "Longitudinal antibody trajectories after a dose of RABV vaccination in dogs that provided ≥ 3 sera, N=27 dogs 
A: Consistently below seroprotective level [n=7, 26% ]
B. Consistently above seroprotective level [n=3, 11%]
C: Progressive waning to < seroprotective titer [n=9, 33%]
D: Anamnestic response occurred within first 6 months of follow-up [n=2, 7.4%]
E: Anamnestic response occurred after first 6 months of follow-up [n=2, 7.4%]
F: Anamnestic response occurred after first 6 months of follow-up [n=1, 3%]
G: Miscellaneous: Responses cannot be described as others [n=3, 11%]"

nob_plot <- generate_combined_plot(
  data = nob_df,
  group_column = "group",
  group_order = nob_order,
  color_by = "comb",
  caption_text = nob_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-nobivac_trajectories.png",
  width = 39,
  height = 39,
  caption_height = 0.17
)
# 4.4 Vaccine-experienced dogs (Nobivac)
nobb_df <- prepare_serology_data("01_src/rabv_serostat.csv",
                                 filter_vaccine = "Nob")
nobb_order <- c("B_up", "B_anam", "B_down")
nobb_caption <- "Trajectory of anti-RABV-G IgG in vaccine-experienced dogs after a dose of RABV vaccination: distinct patterns, N=8 
A: Vaccine-experienced dog: Above/progressive decline but ≥0.5EU/mL [n=4, 50%]
B: Vaccine-experienced dog: Consistently below 0.5EU/m [n=2, 25%]; 
C: Vaccine-experienced dog: Anamnestic response at follow-up [n=2, 25%]"

nobb_plot <- generate_combined_plot(
  data = nobb_df,
  group_column = "group",
  group_order = nobb_order,
  color_by = "comb",
  caption_text = nobb_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-nobivac_boosted.png",
  width = 39,
  height = 39,
  caption_height = 0.17
)
# 4.5 Vaccine-experienced dogs (Rabisin, boosted)
raw_df <- read.csv("01_src/rabv_serostat.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
rabb_df <- raw_df %>%
  filter(!is.na(adjusted), !is.na(month), !is.na(comb), !is.na(group)) %>%
  filter(!grepl("Nob", conc_vax)) %>%
  filter(grepl("booster", dose_type)) %>%
  filter(grepl("multi", set_type)) %>%
  mutate(
    month = as.numeric(month),
    adjusted_plot = ifelse(adjusted < 0.03125, 0.03125, adjusted)
  ) %>%
  filter(month >= 1 & month <= 18)
rabb_order <- c("B_up", "unassign")
rabb_caption <- "Trajectory of anti-RABV-G IgG in vaccine-experienced dogs after a dose of RABV vaccination: distinct patterns, N=5
A: Vaccine-experienced dog: Above/progressive decline but ≥0.5EU/mL, [n=1 dog]
B: Vaccine-experienced dog: Unassigned response, only 2 sera [n=4 dogs]"
rabb_plot <- generate_combined_plot(
  data = rabb_df,
  group_column = "group",
  group_order = rabb_order,
  color_by = "comb",
  caption_text = rabb_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-rabisin_boosted.png",
  width = 39,
  height = 39,
  caption_height = 0.21
)
# 4.6 Vaccine-colored plots
desired_order_vax <- c("below", "above", "decay", "anam1", "anam2", "anam3", 
                       "undesc", "B_up", "B_anam", "B_down")
overall_vax_caption <- "Longitudinal antibody trajectories after a dose of \n 
RABV vaccination in dogs that provided ≥3 sera, N=47 dogs
A: Consistently below seroprotective level [n=13, 28% ]
B. Consistently above seroprotective level [n=6, 13%]
C: Progressive waning to < seroprotective titer [n=13, 28%]
D: Anamnestic response occurred within first 6 months of follow-up [n=5, 11%]
E: Anamnestic response occurred after first 6 months of follow-up [n=3, 6%]
F: Anamnestic response occurred after first 6 months of follow-up [n=3, 6%]
G: Miscellaneous: Responses cannot be described as others [n=4, 9%]

Anti-RABV-G IgG trajectory after vaccine-experience dogs received RABV vaccination and provided ≥ 3 sera, N=9 dogs
H: Vaccine-experienced dogs: Above 0.5EU/mL [n=5, 56%]
I: Vaccine-experienced dogs: Anamnestic response [n=2, 22%]; 
J: Vaccine-experienced dogs: Below 0.5EU/mL[n=2, 22%]"
overall_plot_vax <- generate_combined_plot(
  data = sero_df,
  group_column = "group",
  group_order = desired_order_vax,
  color_by = "vaccine",
  caption_text = overall_vax_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-overall_vaccine_colored.png",
  width = 39,
  height = 39,
  caption_height = 0.25
)

# 4.7 Month-1 filtered plots (missing from original)
raw_df <- read.csv("01_src/rabv_serostat.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
M1_df <- raw_df %>%
  group_by(comb) %>%
  filter(any(month == 1)) %>%
  filter(grepl("multi", set_type)) %>%
  filter(!grepl("2tp", tested_tp)) %>%
  filter(grepl("primary", dose_type)) %>%
  ungroup() %>%
  filter(!grepl("unassign", group)) %>%
  mutate(
    month = as.numeric(month),
    adjusted_plot = ifelse(adjusted < 0.03125, 0.03125, adjusted)
  ) %>%
  filter(month >= 1 & month <= 18)

M1_order <- c("1wane", "2fluct", "3anam", "4below")
M1_caption <- "Longitudinal antibody trajectories in vaccine-naïve dogs with anti-RABV-G IgG 
≥0.5 EU/mL at Month-1 post-vax: Had a total of ≥3 sera, N=15 dogs
A: Had ≥0.5EU/mL @Month-1 subsequently waned to BDL or 0.125EU/mL [n=8, 53%]
B: Had ≥0.5EU/mL @Month-1 & fluctuated around 0.5-1.0 EU/mL [n=4, 27%]
C: Had anamnestic response at follow-up [n=3, 20%]
D: Others: Consistently below seroprotective level [n=2, 100%]"

overall_plot_1M <- generate_combined_plot(
  data = M1_df,
  group_column = "group1M",
  group_order = M1_order,
  color_by = "comb",
  caption_text = M1_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-month1_trajectories.png",
  width = 38,
  height = 29,
  caption_height = 0.26
)

# 4.8 Month-1 vaccine-colored plots
M1_vax_caption <- "Longitudinal antibody trajectories in vaccine-naïve dogs with anti-RABV-G IgG
≥0.5 EU/mL at Month-1 post-vax: Had a total of ≥3 sera, N=15 dogs
A: Had ≥0.5EU/mL @Month-1 subsequently waned to BDL or 0.125EU/mL [n=8, 53%]
B: Had ≥0.5EU/mL @Month-1 & fluctuated around 0.5-1.0 EU/mL [n=4, 27% ]
C: Had anamnestic response at follow-up [n=3, 20% ]
D: Others: Consistently below seroprotective level [n=2, 100% ]"

overall_plot_1Mvax <- generate_combined_plot(
  data = M1_df,
  group_column = "group1M",
  group_order = M1_order,
  color_by = "vaccine",
  caption_text = M1_vax_caption,
  show_title = FALSE,
  show_gam = TRUE,
  save_path = "00-month1_vaccine_colored.png",
  width = 38,
  height = 29,
  caption_height = 0.26
)

##### 5. BATCH PROCESSING EXAMPLE #####
# Example: Generate all plots without GAM and without gridlines
all_configurations <- list(
  list(name = "overall_no_gam", data = sero_df, order = desired_order_overall, 
       color = "comb", caption = overall_caption, gam = FALSE, gridlines = FALSE),
  list(name = "rabisin_no_gam", data = rab_df, order = rab_order,
       color = "comb", caption = rab_caption, gam = FALSE, gridlines = FALSE),
  list(name = "nobivac_no_gam", data = nob_df, order = nob_order,
       color = "comb", caption = nob_caption, gam = FALSE, gridlines = FALSE)
)

for (config in all_configurations) {
  generate_combined_plot(
    data = config$data,
    group_column = "group",
    group_order = config$order,
    color_by = config$color,
    caption_text = config$caption,
    show_title = FALSE,
    show_gam = config$gam,
    show_gridlines = config$gridlines,
    save_path = paste0("00-", config$name, ".png"),
    width = 39,
    height = 45,
    caption_height = 0.2
  )
}

##### 6. CUSTOM PLOT GENERATION (for specific needs) #####
# Example: Create a custom plot with specific options
custom_plot <- function() {
  # Prepare specific dataset
  custom_data <- sero_df %>%
    filter(group %in% c("below", "above", "decay")) %>%
    mutate(group = factor(group, levels = c("below", "above", "decay")))
  
  # Generate plot without GAM, without gridlines, with custom size
  plot <- generate_combined_plot(
    data = custom_data,
    group_column = "group",
    color_by = "vaccine",
    caption_text = "Custom subset: Below, Above, and Decay groups only",
    show_title = FALSE,
    show_gam = FALSE,
    show_gridlines = FALSE,
    save_path = "00-custom_subset.png",
    width = 30,
    height = 20,
    units = "cm"
  )
  
  return(plot)
}

# Generate custom plot
custom_plot_result <- custom_plot()



################ 2D. SEROLOGICAL STUDY: STACKED BAR STATISTICS AROUND SEROPROTECTIVE LEVELS ################
### 2D1. MAIN SOURCE CODES: DEFINE NUMERICAL VALUES AND SOURCES ###
{
library(scales)
library(ggnewscale)
library(cowplot)
raw_df<- read.csv(file = "01_src/rabv_serostat.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
geo_mean <- function(x) {
  exp(mean(log(x[x > 0]), na.rm = TRUE))
}
}
### 2D2. MAIN ANALYSIS CODES: FUNCTIONS & DATA MANAGEMENT ###
suppressPackageStartupMessages({
plot_titer_stacked_bar <- function(data, filter_by = list(), sample_sz = 1, agegroup = NULL, vaxtype = NULL) {
  gm_ratio <- 25  # Transformation ratio: 1 / 0.04, for 0–100 to 0–4.0 on secondary axis
  required_cols <- c("month", "adjusted", "comb")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  if (length(filter_by) > 0) {
    for (var in names(filter_by)) {
      if (!var %in% names(data)) stop(paste("Column", var, "not found in dataset."))
      data <- dplyr::filter(data, .data[[var]] == filter_by[[var]])
    }
  }
  if (nrow(data) == 0) stop("No data left after filtering.")
  # --- Count total dogs ---
  n_dogs <- dplyr::n_distinct(data$comb)
  # --- Prepare data ---
  data_prepped <- data %>%
    mutate(
      month_group = factor(case_when(
        month == 1 ~ "1",
        month %in% 2:3 ~ "2-3",
        month %in% 4:5 ~ "4-5",
        month %in% 6:7 ~ "6-7",
        month %in% 8:9 ~ "8-9",
        month %in% 10:11 ~ "10-11",
        month %in% 12:13 ~ "12-13",
        month %in% 14:15 ~ "14-15",
        month %in% 16:17 ~ "16-17",
        month == 18 ~ "18"
      ), levels = c("1", "2-3", "4-5", "6-7", "8-9", "10-11", "12-13", "14-15", "16-17", "18")),
      adj_category = factor(case_when(
        adjusted >= 0.5 ~ "Had ≥ 0.5EU/mL",
        adjusted < 0.5 & adjusted > 0.063 ~ "Had < 0.5EU/mL but > BDL",
        adjusted == 0.063 ~ "IgG level is BDL",
        TRUE ~ NA_character_
      ), levels = c("IgG level is BDL", "Had < 0.5EU/mL but > BDL", "Had ≥ 0.5EU/mL"))
    ) %>%
    filter(!is.na(adj_category), !is.na(month_group))
  if (nrow(data_prepped) == 0) stop("No valid data after processing.")
  # --- Geometric mean function ---
  geo_mean <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
  # --- Bar data ---
  bar_data <- data_prepped %>%
    group_by(month_group, adj_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(month_group) %>%
    mutate(percent = count / sum(count) * 100)
  # --- Geometric mean data ---
  geo_data <- data_prepped %>%
    group_by(month_group, adj_category) %>%
    summarise(gm = geo_mean(adjusted), .groups = "drop") %>%
    filter(adj_category %in% c("Had ≥ 0.5EU/mL", "Had < 0.5EU/mL but > BDL"))
  # --- Dog/sample counts for labeling ---
  dog_counts <- data %>%
    mutate(month_group = factor(case_when(
      month == 1 ~ "1",
      month %in% 2:3 ~ "2-3",
      month %in% 4:5 ~ "4-5",
      month %in% 6:7 ~ "6-7",
      month %in% 8:9 ~ "8-9",
      month %in% 10:11 ~ "10-11",
      month %in% 12:13 ~ "12-13",
      month %in% 14:15 ~ "14-15",
      month %in% 16:17 ~ "16-17",
      month == 18 ~ "18"
    ), levels = levels(data_prepped$month_group))) %>%
    group_by(month_group) %>%
    summarise(D = n_distinct(comb), .groups = "drop")
  sd_labels <- data_prepped %>%
    group_by(month_group) %>%
    summarise(S = sum(!is.na(adjusted)), .groups = "drop") %>%
    left_join(dog_counts, by = "month_group") %>%
    mutate(label = paste0("S: ", S, "/D: ", D)) %>%
    left_join(
      bar_data %>% group_by(month_group) %>% summarise(y_pos = sum(percent), .groups = "drop"),
      by = "month_group"
    ) %>%
    mutate(y_pos = y_pos + 5)
  
  # --- Colors ---
  fill_colors <- c(
    "IgG level is BDL" = "#4c4c4c",
    "Had < 0.5EU/mL but > BDL" = "#FFCC33",
    "Had ≥ 0.5EU/mL" = "#4FA861"
  )
  line_colors <- c(
    "Had ≥ 0.5EU/mL" = "#e41a1c",
    "Had < 0.5EU/mL but > BDL" = "#377eb8"
  )
  plot_title <- paste0(
    "Distribution of IgG titers in ", agegroup, " vaccine-naive dogs with ≥", sample_sz,
    " blood \n samples collected after a dose of ", vaxtype, " RABV vaccination, N = ", n_dogs, " dogs"
  )
}})
### 2D3. CODES FOR GENERATION OF PUBLICATION-READY PLOTS ###
suppressPackageStartupMessages({
  # --- Final plot ---
  p <- ggplot() +
    geom_bar(
      data = bar_data,
      aes(x = month_group, y = percent, fill = adj_category),
      stat = "identity", position = "stack", width = 0.8
    ) +
    geom_text(
      data = sd_labels,
      aes(x = month_group, y = y_pos, label = label),
      size = 2.6, fontface = "bold", color = "black"
    ) +
    scale_fill_manual(
      values = fill_colors,
      name = "IgG Ab Response Range\n(S: sera datapoints, D: dogs)"
    ) +
    scale_y_continuous(
      name = "Sample proportion (%)",
      limits = c(0, 100),
      sec.axis = sec_axis(~ . * 0.04, name = "Geometric Mean IgG Titer (EU/mL)", breaks = seq(0, 4, 0.5))
    ) +
    ggnewscale::new_scale_color() +
    geom_line(
      data = geo_data %>% filter(adj_category == "Had ≥ 0.5EU/mL"),
      aes(x = month_group, y = gm * gm_ratio, group = 1, color = "Had ≥ 0.5EU/mL"),
      size = 1.2
    ) +
    geom_point(
      data = geo_data %>% filter(adj_category == "Had ≥ 0.5EU/mL"),
      aes(x = month_group, y = gm * gm_ratio, color = "Had ≥ 0.5EU/mL"),
      size = 2.5
    ) +
    geom_text(
      data = geo_data %>% filter(adj_category == "Had ≥ 0.5EU/mL"),
      aes(x = month_group, y = gm * gm_ratio - 3, label = round(gm, 2)),
      color = "#6c1010", size = 3.5, fontface = "bold"
    ) +
    geom_line(
      data = geo_data %>% filter(adj_category == "Had < 0.5EU/mL but > BDL"),
      aes(x = month_group, y = gm * gm_ratio, group = 1, color = "Had < 0.5EU/mL but > BDL"),
      linetype = "dashed", size = 1.2
    ) +
    geom_point(
      data = geo_data %>% filter(adj_category == "Had < 0.5EU/mL but > BDL"),
      aes(x = month_group, y = gm * gm_ratio, color = "Had < 0.5EU/mL but > BDL"),
      size = 2.5
    ) +
    geom_text(
      data = geo_data %>% filter(adj_category == "Had < 0.5EU/mL but > BDL"),
      aes(x = month_group, y = gm * gm_ratio - 3, label = round(gm, 2)),
      color = "#003747", size = 3.5, fontface = "bold"
    ) +
    scale_color_manual(
      values = line_colors,
      name = "Geometric Mean Lines"
    ) +
    labs(
      x = "Months since last vaccination",
      title = plot_title
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(face = "bold", size = 9),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5)
    )
  return(p)
text_box_grob <- function(text, width = 80) {
  if (!is.null(width)) {
    wrapped_text <- strwrap(text, width = width)
    text_string <- paste(wrapped_text, collapse = "\n")
  } else {
    text_string <- text  # Use as-is if width not specified
  }
  # Create grob with tighter line spacing
  grob <- grobTree(
    rectGrob(gp = gpar(fill = NA, col = "black", lwd = 0.7)),
    textGrob(
      label = text_string,
      x = 0.01, y = 1, just = c("left", "top"),
      gp = gpar(cex = 1.0, lineheight = 1)  # Reduced line spacing
    )
  )
  # Wrap grob into a ggplot object
  ggplot() +
    theme_void() +
    annotation_custom(grob) +
    theme(plot.margin = margin(0, 0, 0, 0),legend.position = "none")
}
})
### 2D4. SAVING PLOTS WITH GGSAVE IN PNG: WRAP, DEFINE PATH, AND FILENAME  ###
suppressPackageStartupMessages({
ps_101 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary"), sample_sz = 1)
ps_101_comb <- ps_101 / text_box_grob("Interpret: Gradual early decline reflects natural waning of vaccine-induced immunity, subsequent titer rise and fall may be explained by cohort effects as a result of including other dogs sampled in a less longitudinal structure.") +
  plot_layout(heights = c(1, 0.15)) 
ps_73 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", set_type = "multi"), sample_sz = 2) #Have to change plot title text
ps_73_comb <- ps_73 / text_box_grob("Interpret: Gradual early decline consistent with waning of humoral immunity. Subsequent titer rise could be cohort effects or natural re-exposure of antigen.") +
  plot_layout(heights = c(1, 0.15)) 
ps_age01 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", age_G2 = "≥3 - 6m"), agegroup = '≥3 - 6m old ')
ps_age01_comb <- ps_age01 / text_box_grob("Interpret: Dog remain protected during the time of sampling, data did not suggest pattern of a booster response or natural re-exposure to antigen.") +
  plot_layout(heights = c(1, 0.15)) 
ps_age02 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", age_G2 = "≥12m-24m"), agegroup = '≥12m-24m old ')
ps_age02_comb <- ps_age02 / text_box_grob("Interpret: Small-sample bias in this age cohort after 9-months post-vaccination.") +
  plot_layout(heights = c(1, 0.08)) 
ps_age03 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", age_G2 = "≥25m - <6y"), agegroup = '≥25m - <6y old ')
ps_age03_comb <- ps_age03 / text_box_grob("Interpret: Strong primary immune response with expected waning of humoral immunity. The next rise in titers may due to natural re-exposure to antigen.") +
  plot_layout(heights = c(1, 0.12)) 
ps_age04 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", age_G2 = ">7m-12m"), agegroup = '>7m-12m old ')
ps_age04_comb <- ps_age04 / text_box_grob("Interpret: Small-sample bias.") +
  plot_layout(heights = c(1, 0.12)) 
ps_age05 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", age_G2 = "≥6y "), agegroup = '>6 years old ')
ps_age05_comb <- ps_age05 / text_box_grob("Interpret: Small-sample bias.") +
  plot_layout(heights = c(1, 0.12)) 
ps_rab06 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", conc_vax = "Rabisin"), vaxtype =' Rabisin ')
ps_rab06_comb <- ps_rab06 / text_box_grob("Interpret: Strong primary immune response with a consisent trend.") +
  plot_layout(heights = c(1, 0.12)) 
ps_nob07 <- plot_titer_stacked_bar(data = raw_df, filter_by = list(dose_type = "primary", conc_vax = "Nobivac"), vaxtype =' Nobivac ')
ps_nob07_comb <- ps_nob07 / text_box_grob("Interpret: Initial peak of primary immune response. Seroconversion failure was higher in the first 6 months than Rabisin.") +
  plot_layout(heights = c(1, 0.12)) 
ps_116 <- plot_titer_stacked_bar(data = raw_df, sample_sz = 1, agegroup = 'vaccine-experience & ')
ps_116_comb <- ps_116 / text_box_grob("Interpret: Immunity above protective threshold, GM is 3/4 times higher, no significant waning. Irrespective of vaccine-naive or -experienced dogs.") +
  plot_layout(heights = c(1, 0.12)) 
plots_stacked <- ggarrange(ps_101_comb, ps_73_comb, ps_age01_comb, ps_age02_comb, ps_age03_comb, ps_age04_comb, ps_age05_comb, ps_rab06_comb, ps_nob07_comb, ps_116_comb, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), nrow = 5, ncol = 2, common.legend = FALSE)
legend <- get_legend(ps_101 + theme(legend.position = "bottom"))
final_plot <- plot_grid(
  plots_stacked,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.10)
)
})
ggsave(final_plot, filename = '00-rvsero_ibet_resplvls_116dogsv1.png', width = 36, height = 63, units = "cm",limitsize = FALSE, device = "png")

################ 2E. SEROLOGICAL STUDY: GEOMMAPPING VACCINATION ROLLOUT BY IBET ################
### 2E1. MAIN SOURCE CODES: DEFINE NUMERICAL VALUES AND SOURCES  ###
suppressPackageStartupMessages({
  library(sf)
library(dplyr)
library(leaflet)
library(lubridate)
library(osrm)
library(readr)
library(RColorBrewer)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggrepel)
library(viridis)
library(webshot2)
library(htmlwidgets)
library(leaflet.extras)

# 1. Load CSV
data<- read.csv(file = "01_src/rabv_serostat.csv")
data <- data[, c("dst_site", "WKT", "cdate", "comb", "groupW", "dose_type")]
data <- data[order(data$cdate), ]

# --- Toggle label visibility for dst_site ---
show_labels <- FALSE  # Set to TRUE to display destination site labels
# --- Load and preprocess data ---
swk_divisions_sf <- st_read("swk_maps/swk_divmap_2022.shp") %>%
  st_transform(crs = 4326) %>%
  filter(!st_is_empty(.))
if (!"ADM2_EN" %in% names(swk_divisions_sf)) {
  stop("Column 'ADM2_EN' not found in shapefile.")}
})
### 2E2. MAIN ANALYSIS CODES: FUNCTIONS & DATA MANAGEMENT ###
{
# --- Prepare and clean the route data ---
data2 <- data %>%
  select(dst_site, WKT, cdate, comb, groupW) %>%
  mutate(cdate = as.Date(cdate),
         month_label = format(cdate, "%Y-%b")  # Keep lowercase here for sorting
  ) %>%
  arrange(cdate) %>%
  distinct()
data2 <- data %>%
  mutate(
    cdate = as.Date(cdate),
    month_label = format(cdate, "%Y-%b")  # Keep lowercase here for sorting
  ) %>%
  arrange(cdate) %>%
  distinct()
# 2. Convert to sf object
routes_sf <- st_as_sf(data2, wkt = "WKT", crs = 4326) %>%
  arrange(cdate)
# 3. Extract and sort unique month labels based on actual date order
month_levels <- routes_sf %>%
  distinct(month_label, cdate) %>%
  arrange(cdate) %>%
  pull(month_label)
# Optional: Capitalize the labels for display
routes_sf$month_label <- toupper(routes_sf$month_label)
month_levels <- unique(routes_sf$month_label)
palette_colors <- viridis(length(month_levels) + 2, option = "inferno")[1:length(month_levels)]
month_color_df <- data.frame(
  month_label = month_levels,
  color = palette_colors,
  stringsAsFactors = FALSE
)
routes_sf <- routes_sf %>%
  left_join(month_color_df, by = "month_label")
# --- Extract coordinates + distance ---
get_coords <- function(geom) {
  coords <- st_coordinates(geom)[, 1:2]
  list(
    start = coords[1, ],
    end   = coords[nrow(coords), ],
    mid   = coords[round(nrow(coords)/2), ],
    coords = coords
  )}
coords_list <- lapply(st_geometry(routes_sf), get_coords)
routes_sf <- routes_sf %>%
  mutate(
    start_lon = sapply(coords_list, function(x) x$start[1]),
    start_lat = sapply(coords_list, function(x) x$start[2]),
    end_lon   = sapply(coords_list, function(x) x$end[1]),
    end_lat   = sapply(coords_list, function(x) x$end[2]),
    mid_lon   = sapply(coords_list, function(x) x$mid[1]),
    mid_lat   = sapply(coords_list, function(x) x$mid[2]),
    coords    = coords_list
  ) %>%
  rowwise() %>%
  mutate(
    distance_km = round(distHaversine(c(start_lon, start_lat), c(end_lon, end_lat)) / 1000)
  ) %>%
  ungroup()
# --- Offset overlapping polylines to reduce clutter ---
offset_polyline <- function(coords, index, offset_m = 300) {
  coords <- as.matrix(coords)
  if (ncol(coords) > 2) coords <- coords[, 1:2]
  b <- bearing(coords[1, ], coords[nrow(coords), ])
  shift_b <- ifelse(index %% 2 == 0, b + 95, b + 85)
  shifted_coords <- t(apply(coords, 1, function(pt) destPoint(pt, b = shift_b, d = offset_m * index)))
  st_linestring(shifted_coords)
}
offset_geoms <- mapply(
  offset_polyline,
  coords = lapply(routes_sf$coords, function(x) x$coords),
  index = seq_along(routes_sf$coords),
  SIMPLIFY = FALSE
)
routes_sf$geometry <- st_sfc(offset_geoms, crs = 4326)
# --- Summarize by destination for deduplicated distance labels ---
set.seed(42)
jitter_amt <- 0.1
site_stats <- routes_sf %>%
  st_drop_geometry() %>%
  group_by(dst_site, end_lon, end_lat) %>%
  summarise(
    freq = n(),
    distance_km = first(distance_km),
    .groups = "drop")
site_stats2 <- routes_sf %>%
  st_drop_geometry() %>%
  group_by(end_lon, end_lat) %>%
  summarise(
    freq = n(),
    unique_comb = n_distinct(comb),
    primary_count = n_distinct(comb[dose_type == "primary"]),
    boost_count = n_distinct(comb[dose_type == "booster"]),
    anam_count = n_distinct(comb[groupW == "anam"]),
    dose_type = as.character(dose_type),
    .groups = "drop"
  ) %>%
  mutate(
    jittered = duplicated(paste(end_lon, end_lat)) | duplicated(paste(end_lon, end_lat), fromLast = TRUE),
    jitter_lon = ifelse(jittered, end_lon + runif(n(), -jitter_amt, jitter_amt), end_lon),
    jitter_lat = ifelse(jittered, end_lat + runif(n(), -jitter_amt, jitter_amt), end_lat),
    jitter_square_lon = jitter_lon + 0.04,  # fixed offset for square
    primary_count = as.integer(primary_count),
    boost_count = as.integer(boost_count),# ensure numeric
    anam_count = as.integer(anam_count),        # ensure numeric
    label_text = as.character(primary_count)# ensure character
  )

# --- Apply jitter only if overlapping coordinates ---
coords <- site_stats2 %>% select(end_lon, end_lat)
dup_coords <- duplicated(coords) | duplicated(coords, fromLast = TRUE)
site_stats2 <- site_stats2 %>%
  mutate(
    jitter_circle_lon = ifelse(dup_coords, end_lon - jitter_amt, end_lon),
    jitter_square_lon = ifelse(dup_coords, end_lon + jitter_amt, end_lon)
  )
}
### 2E3. MAIN ANALYSIS CODES: MAPPING ###
{
basemap_tiles <- "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}{r}.png"
attrib <- '&copy; <a href="https://www.openstreetmap.org/">OSM</a>'
m <- leaflet() %>%
  addTiles(urlTemplate = basemap_tiles, attribution = attrib) %>%
  setView(lng = 114.0, lat = 1.5, zoom = 6)
# --- Add division boundaries with dashed lines ---
m1 <- m %>%
  addPolylines(
    data = swk_divisions_sf,
    color = "black",
    weight = 2.5,
    opacity = 0.8,
    dashArray = "4,4"
  )
# --- Add division names using internal points (avoids overlapping lines) ---
label_points <- st_point_on_surface(swk_divisions_sf)
label_coords <- st_coordinates(label_points)
div_df <- data.frame(
  label = swk_divisions_sf$ADM2_EN,
  lon = label_coords[, 1],
  lat = label_coords[, 2],
  stringsAsFactors = FALSE
)
m1 <- m1 %>%
  addLabelOnlyMarkers(
    lng = div_df$lon,
    lat = div_df$lat,
    label = div_df$label,
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "center",
      textOnly = TRUE,
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "12px",
        "background" = "rgba(255,255,255,0.6)",
        "padding" = "2px 4px",
        "border" = "1px solid #00000066"
      )
    )
  )
# --- Add route polylines ---
m2 <- m1 %>%
  addPolylines(
    data = routes_sf,
    color = ~color,
    weight = 3.5,
    opacity = 0.9,
    popup = ~paste0(
      "<b>", dst_site, "</b><br>",
      "Date: ", month_label, "<br>",
      "Distance: ", distance_km, " km"
    )
  )
# --- Add destination markers sized by frequency ---
m2 <- m2 %>%
  addCircleMarkers(
    data = site_stats,
    lng = ~end_lon,
    lat = ~end_lat,
    radius = ~2.2 + sqrt(freq) * 1.5,
    color = "#141414",
    weight = 1,
    fillColor = "#3b3b3b",
    fillOpacity = 1,
    stroke = TRUE
  )
set.seed(123)
jitter_amt <- 0.02
site_stats2 <- site_stats2 %>%
  group_by(end_lon, end_lat) %>%
  mutate(
    jittered = n() > 1,
    jitter_lon = ifelse(jittered, end_lon + runif(1, -jitter_amt, jitter_amt), end_lon),
    jitter_lat = ifelse(jittered, end_lat + runif(1, -jitter_amt, jitter_amt), end_lat),
    jitter_square_lon = jitter_lon + 0.02,
    jitter_square_lat = jitter_lat + 0.02# offset square to the right
  ) %>%
  ungroup()

m4 <- m1 %>%
  addCircleMarkers(
    data = site_stats2 %>% filter(primary_count > 0),
    lng = ~jitter_lon,
    lat = ~jitter_lat,
    radius = ~6 + sqrt(freq) * 2,
    color = "#000000",
    weight = 1,
    # change fillColor based on dose_type
    fillColor = "#e0e0e0",
    fillOpacity = 1,
    stroke = TRUE,
    label = ~as.character(primary_count),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "center",
      textOnly = TRUE,
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "11px",
        "text-align" = "center"
      )
    )
  )
m4 <- m4 %>%
  addCircleMarkers(
    data = site_stats2 %>% filter(boost_count > 0),
    lng = ~jitter_lon+0.01,
    lat = ~jitter_lat+0.01,
    radius = ~6 + sqrt(freq) * 2,
    color = "#000000",
    weight = 1,
    # change fillColor based on dose_type
    fillColor = "#33a02c",
    fillOpacity = 1,
    stroke = TRUE,
    label = ~as.character(boost_count),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "center",
      textOnly = TRUE,
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "11px",
        "text-align" = "center"
      )
    )
  )
# Square marker icon
square_icon <- makeIcon(
  iconUrl = "https://img.icons8.com/?size=100&id=4knsel5keE42&format=png&color=000000",
  iconWidth = 20,
  iconHeight = 20
)
# Add square markers (for 'anam' counts)
m4 <- m4 %>%
  addMarkers(
    data = site_stats2 %>% filter(anam_count > 0),
    lng = ~jitter_square_lon,
    lat = ~jitter_square_lat,
    icon = square_icon,
    label = ~as.character(anam_count),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "right",
      textOnly = TRUE,
      style = list(
        "color" = "darkgreen",
        "font-weight" = "bold",
        "font-size" = "12px",
        "text-align" = "left"
      )
    )
  )
# --- Add deduplicated distance labels with jitter ---
site_stats <- site_stats %>%
  mutate(
    jitter_lon = end_lon + runif(n(), -jitter_amt, jitter_amt),
    jitter_lat = end_lat + runif(n(), -jitter_amt, jitter_amt)
  )
m2 <- m2 %>%
  addLabelOnlyMarkers(
    data = site_stats,
    lng = ~jitter_lon,
    lat = ~jitter_lat,
    label = ~paste0(distance_km, " km"),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "auto",
      textOnly = TRUE,
      style = list(
        "color" = "darkblue",
        "font-style" = "italic",
        "font-size" = "9px"
      )
    )
  )
# --- Add Markas UKPS label at the start point ---
m2 <- m2 %>%
  addLabelOnlyMarkers(
    lng = routes_sf$start_lon[1],
    lat = routes_sf$start_lat[1],
    label = "Markas UKPS",
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "right",
      textOnly = TRUE,
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "14px",
        "background" = "white",
        "padding" = "4px 8px",
        "border" = "2px solid black"
      )
    )
  )
m5 <- m4 %>%
  addLabelOnlyMarkers(
    lng = routes_sf$start_lon[1],
    lat = routes_sf$start_lat[1],
    label = "Markas UKPS",
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "right",
      textOnly = TRUE,
      style = list(
        "color" = "black",
        "font-weight" = "bold",
        "font-size" = "14px",
        "background" = "white",
        "padding" = "4px 8px",
        "border" = "2px solid black"
      )
    )
  )
}
### 2E4. MAIN PLOTTING CODES: MAPPING ###
{
# --- Create individual dots per 'comb' (m6 map) ---
set.seed(42)  # For reproducible jitter
max_jitter_radius <- 0.05  # Max radius (degrees) for enclosing circles

# 1. Get number of combs per site
site_comb_counts <- routes_sf %>%
  st_drop_geometry() %>%
  group_by(end_lon, end_lat) %>%
  filter(dose_type == 'primary') %>%  ###remove this if regardless of booster or primary
  summarise(n_comb = n_distinct(comb), .groups = "drop")

# 2. Calculate radius to enclose all dots (empirical formula or fixed max)
# We scale the radius up slightly based on the number of dots
site_comb_counts <- site_comb_counts %>%
  mutate(circle_radius = pmin(0.005 + sqrt(n_comb) * 0.007, max_jitter_radius))

# 3. Join radius back to routes_sf
routes_with_radius <- routes_sf %>%
  st_drop_geometry() %>%
  left_join(site_comb_counts, by = c("end_lon", "end_lat")) %>%
  mutate(site_id = paste(end_lon, end_lat, sep = "_"))

individual_points <- routes_with_radius %>%
  group_by(site_id) %>%
  mutate(
    # Random angle between 0 and 2*pi
    angle = runif(n(), 0, 2 * pi),
    # Random radius scaled down to 70% of enclosing radius
    radius_use = runif(n(), 0.2, 0.8) * circle_radius[1],
    # Calculate jittered positions inside circle
    jitter_lon = end_lon + radius_use * cos(angle),
    jitter_lat = end_lat + radius_use * sin(angle),
    dot_color = ifelse(groupW == "anam", "#e41a1c", "#1f78b4")) %>% 
  filter(dose_type == 'primary') %>%  ###remove this if regardless of booster or primary
  ungroup()

# --- Create map m6 with enclosing circles and individual dots ---
m6 <- m1 %>%
  # Add enclosing circles at each destination site
  addCircles(
    data = site_comb_counts,
    lng = ~end_lon,
    lat = ~end_lat,
    radius = ~circle_radius * 111000,  # convert degrees to meters (~111km per degree)
    fillColor = "#f0f0f0",
    fillOpacity = 0.2,
    color = "#555555",
    weight = 1
  ) %>%
  # Add individual jittered dots per 'comb'
  addCircleMarkers(
    data = individual_points,
    lng = ~jitter_lon,
    lat = ~jitter_lat,
    radius = 3,
    color = "#000000",
    fillColor = "#1f78b4",
    fillOpacity = 0.8,
    stroke = TRUE,
    weight = 0.5,
    popup = ~paste0("comb: ", comb, "<br>dst_site: ", dst_site)
  )

# Optional: Add site labels
if (show_labels) {
  m6 <- m6 %>%
    addLabelOnlyMarkers(
      data = site_comb_counts,
      lng = ~end_lon,
      lat = ~end_lat,
      label = ~paste0("n=", n_comb),
      labelOptions = labelOptions(
        noHide = TRUE,
        direction = 'top',
        textOnly = TRUE,
        style = list(
          "color" = "black",
          "font-weight" = "bold",
          "font-size" = "10px"
        )))
}
# Assign colors based on groupW
individual_points <- individual_points %>%
  mutate(
    dot_color = ifelse(groupW == "anam", "#e41a1c", "#1f78b4")  # red for anam
  )
m7 <- m1 %>%
  # Add enclosing circles for each site
  addCircles(
    data = site_comb_counts,
    lng = ~end_lon,
    lat = ~end_lat,
    radius = ~circle_radius * 111000,
    fillColor = "#f0f0f0",
    fillOpacity = 0.2,
    color = "#555555",
    weight = 1
  ) %>%
  # Add individual dots, with color based on groupW
  addCircleMarkers(
    data = individual_points,
    lng = ~jitter_lon,
    lat = ~jitter_lat,
    radius = 3,
    color = "#000000",
    fillColor = ~dot_color,
    fillOpacity = 0.8,
    stroke = TRUE,
    weight = 0.5,
    popup = ~paste0("comb: ", comb, "<br>dst_site: ", dst_site, "<br>groupW: ", groupW)
  )

# Optional: Add labels again if enabled
if (show_labels) {
  m7 <- m7 %>%
    addLabelOnlyMarkers(
      data = site_comb_counts,
      lng = ~end_lon,
      lat = ~end_lat,
      label = ~paste0("n=", n_comb),
      labelOptions = labelOptions(
        noHide = TRUE,
        direction = 'top',
        textOnly = TRUE,
        style = list(
          "color" = "black",
          "font-weight" = "bold",
          "font-size" = "10px"
        )
      )
    )
}

# --- Add legend based on filtered months used ---
m2 <- m2 %>%
  addLegend(
    "topright",
    colors = month_color_df$color,
    labels = month_color_df$month_label,
    title = "Trip Month",
    opacity = 1
  )
# ----------------- Arrow-pointed dog unique count on map -------------------------------- #
arrow_length_m <- 10000  # adjust for visibility
}
### 2E5. MAIN SAVE CODES: SAVING IN HTML WITH SETVIEW ###
{
m3 <- m2 %>%
  setView(lng = 110.18732995643, lat = 1.422862652800, zoom = 9.5)
m4 <- m4 %>%
  setView(lng = 110.18732995643, lat = 1.422862652800, zoom = 9.5)
m6 <- m6 %>%
  setView(lng = 110.18732995643, lat = 1.422862652800, zoom = 9.5)
m7 <- m7 %>%
  setView(lng = 110.18732995643, lat = 1.422862652800, zoom = 9.5)
m1 <- m1 %>%
  setView(lng = 112.9268672476329, lat = 2.531578652732783, zoom = 8)
m0 <- m1 %>%
  setView(lng = 108.655368573478, lat = 3.69083071746549, zoom = 6)
saveWidget(m3, "temp_leaflet1_map.html", selfcontained = TRUE)
saveWidget(m4, "temp_leafletsamp_map.html", selfcontained = TRUE)
saveWidget(m6, "temp_leafletsamp_map2.html", selfcontained = TRUE)
saveWidget(m7, "temp_leafletsamp_map3.html", selfcontained = TRUE)
saveWidget(m1, "temp_divisions_map.html", selfcontained = TRUE)
saveWidget(m0, "temp_states_map.html", selfcontained = TRUE)
# Take screenshot and save as high-resolution PNG
webshot("temp_divisions_map.html", 
        file = "a-rvsero_ibet_swkpmap.png", 
        vwidth = 1060,    # width in pixels
        vheight = 800,
        zoom = 8)
webshot("temp_leaflet1_map.html", 
        file = "b-rvsero_ibet_vaxtripmap.png", 
        vwidth = 1060,    # width in pixels
        vheight = 800,
        zoom = 9.0)
webshot("temp_leafletsamp_map2.html", 
        file = "b-rvsero_ibet_sampmap2.png", 
        vwidth = 1060,    # width in pixels
        vheight = 800,
        zoom = 9.0)
webshot("temp_states_map.html", 
        file = "b-rvsero_ibet_statemap.png", 
        vwidth = 1060,    # width in pixels
        vheight = 800,
        zoom = 6)
}

################################### 2F. SEROLOGICAL STUDY: Longitudinal Timeline with Death Events ##############################
### 2F1. MAIN SOURCE CODES: DEFINE NUMERICAL VALUES AND SOURCES  ###
suppressPackageStartupMessages({
library(ggplot2)
library(dplyr)
library(ggimage)
library(forcats)
library(lubridate)
df <- read.csv(file = "01_src/rabv_serostat.csv")
# --- 1. Data Preparation ---
df <- df %>%
  mutate(
    vdate = as.Date(vdate),
    died_on = as.Date(died_on),
    sample_date = vdate %m+% months(month)
  )

# Filter to dogs with death date
df_dead <- df %>%
  filter(!is.na(died_on))

# Vaccine type lookup
comb_vax_lookup <- df_dead %>%
  select(comb, conc_vax) %>%
  distinct()
})
### 2F2. MAIN ANALYSIS CODES: FUNCTIONS & DATA MANAGEMENT ###
{
samples_df <- df_dead %>%
  select(comb, sample_date, adjusted) %>%
  distinct() %>%
  rename(date = sample_date)
# Build timeline events
events <- df_dead %>%
  group_by(comb) %>%
  summarise(
    vdate = min(vdate),
    died_on = unique(died_on),
    sample_dates = list(unique(sample_date)),
    .groups = "drop"
  ) %>%
  left_join(comb_vax_lookup, by = "comb")
# Flatten into long format
timeline_df <- events %>%
  rowwise() %>%
  mutate(
    events = list(data.frame(
      comb = comb,
      conc_vax = conc_vax,
      event = c("vaccinated", rep("sample", length(sample_dates)), "died"),
      date = c(vdate, sample_dates, died_on),
      stringsAsFactors = FALSE
    ))) %>%
  pull(events) %>%
  bind_rows()

# Join `adjusted` values to sample rows only
timeline_df <- timeline_df %>%
  left_join(samples_df, by = c("comb", "date")) %>%
  mutate(
    icon = case_when(
      event == "vaccinated" ~ "https://img.icons8.com/?size=100&id=ORLbe5bmWXza&format=png&color=000000",
      event == "sample"     ~ "https://img.icons8.com/?size=100&id=81267&format=png&color=000000",
      event == "died"       ~ "https://freemiumicons.com/wp-content/uploads/2022/05/dead-dog-icon.png"
    ),
    event = factor(event, levels = c("vaccinated", "sample", "died"))
  )
}
### 2F3. MAIN PLOTTING AND SAVE CODES: PLOTTING ###
{
library(ggtext)
# --- Step 3: Plotting Function ---
plot_timeline_by_vaccine <- function(data, vaccine_label) {
  data_filtered <- data %>% filter(conc_vax == vaccine_label)
  
  # Calculate month differences
  event_diff_labels <- data_filtered %>%
    arrange(comb, date) %>%
    group_by(comb) %>%
    mutate(
      next_date = lead(date),
      label_months = round(interval(date, next_date) / months(1)),
      mid_date = date + (next_date - date) / 2
    ) %>%
    filter(!is.na(label_months))
  data_samples <- data_filtered %>% 
    filter(event == "sample") %>%
    mutate(
      rich_label = paste0("<span style='font-size:8pt;'>", adjusted,
                          "</span><span style='font-size:5pt;'> EU/mL</span>")
    )
  # Plot
  ggplot(data_filtered, aes(x = date, y = fct_rev(comb))) +
    geom_line(aes(group = comb), color = "gray60", linewidth = 0.6) +
    geom_image(aes(image = icon), size = 0.03) +
    
    # Add month gap labels
    geom_text(
      data = event_diff_labels,
      aes(x = mid_date, y = fct_rev(comb), label = paste0(label_months, " mo")),
      size = 3.2, vjust = -0.5
    ) +
    
    # Add adjusted value above each sample
    ggtext::geom_richtext(
      data = data_samples,
      aes(label = rich_label),
      fill = NA, label.color = NA,  # no background
      vjust = -0.15,
      size = 2.8
    ) +
    
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.04))) +
    labs(
      title = paste("Life course of a vaccine-naive dog received ", vaccine_label, 'before death'),
      x = "Date",
      y = "Dog-Owner Identity"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10) 
    )
}

# --- Step 4: Generate plots for each vaccine ---
plot_rabisin <- plot_timeline_by_vaccine(timeline_df, "Rabisin")
plot_nobivac <- plot_timeline_by_vaccine(timeline_df, "Nobivac")
library(ggpubr)
# --- Step 5: Combine plots using ggarange ---
combvax_plot <- ggarrange(
  plot_rabisin,
  plot_nobivac,
  nrow = 2,
  labels = c("A", "B")
)

ggsave(combvax_plot, filename = 'swk_IBETseroevents_2025.png',width = 31, height = 36, units = "cm",limitsize = FALSE, device = 'png')
}

################################### 3A. RABV PHYLOGENY: MAIN ANALYSIS STARTS HERE ##############################
### 3A1. MAIN SOURCE CODES: DEFINE CLADE, LOCATIONS, REGIONS, VALUES ###
suppressPackageStartupMessages({
library(readxl)
library(colorspace)
  subclade_colors <- c(
    ## ---- Africa ----
    Africa_2 = "#1B9E77",   # teal green
    Africa_3 = "#0F766E",   # deep teal
    ## ---- Arctic-related ----
    Arctic_related_AL1b = "#4575B4",  # cold blue
    Arctic_related_AL2  = "#313695",  # deep arctic blue
    ## ---- Asian (SEA) ----
    Asian_SEA1a = "#E6AB02",  # golden yellow
    Asian_SEA1b = "#D95F02",  # orange
    Asian_SEA1c = "#B35806",  # darker orange
    Asian_SEA1c_q = "#F1A340",# lighter variant (?)
    Asian_SEA2a = "#66A61E",  # green
    Asian_SEA2b = "#4D9221",  # dark green
    Asian_SEA3  = "#1B7837",  # forest green
    Asian_SEA4  = "#01665E",  # teal
    Asian_SEA5  = "#35978F",  # blue-green
    Asian_SEA5_q = "#80CDC1", # lighter variant (?)
    ## ---- Cosmopolitan (AF) ----
    Cosmopolitan_AF1a = "#E7298A",  # magenta
    Cosmopolitan_AF1b = "#CE1256",  # deep magenta
    Cosmopolitan_AF1c = "#91003F",  # dark wine
    Cosmopolitan_AF4  = "#DD1C77",  # pink-red
    ## ---- Cosmopolitan (CA) ----
    Cosmopolitan_CA1 = "#7570B3",   # purple-blue
    Cosmopolitan_CA2 = "#542788",   # deep purple
    ## ---- Cosmopolitan (Vac) ----
    Cosmopolitan_Vac = "#FB8072",   # coral
    ## ---- Indian Subcontinent ----
    Indian_Subcontinent = "#A6761D" # warm brown
  )
  sublin_colors <- c(
    ## ---- Africa ----
    Africa_2 = "#1B9E77",   # teal green
    Africa_3 = "#0F766E",   # deep teal
    
    ## ---- Arctic-related ----
    Bat = "#4575B4",  # cold blue
    Arctic_related  = "#313695",  # deep arctic blue
    
    ## ---- Asian (SEA) ----
    SEA1a = "#E6AB02",  # golden yellow
    SEA1b = "#D95F02",  # orange
    SEA1c = "#B35806",  # darker orange

    SEA2 = "#66A61E",  # green
    SEA3  = "#1B7837",  # forest green
    SEA4  = "#01665E",  # teal
    SEA5  = "#35978F",  # blue-green
    ## ---- Cosmopolitan (AF) ----
    MexSK_1 = "#E7298A",  # magenta
    RAC = "#CE1256",  # deep magenta
    SCSK = "#91003F",  # dark wine
    Cosmopolitan_wild_type = "#542788",   # deep purple
    Cosmopolitan_vaccine = "#FB8072",   # coral
    ## ---- Indian Subcontinent ----
    Indian_Subcontinent = "#A6761D" # warm brown
  )
  country_colors <- c(
    "Sarawak"            = "#5E3C99",  # deep violet
    "China"          = "#B22222",  # deep firebrick red
    "France"         = "#1F78B4",  # strong deep blue
    "Japan"          = "#A9447A",  # deep rose-purple
    "Hungary"        = "#2F7F3E",  # deep forest green
    "Taiwan"         = "#A0792A",  # deep golden-brown
    "USA"            = "#5A4A42",  # deep warm gray
    "VietNam"        = "#156F88",  # deep blue-teal
    "Finland"        = "#345E9E",  # cooler deep blue
    "Italy"          = "#2E8B57",  # deep sea green
    "Madagascar"     = "#8B5A2B",  # saddle brown (deep)
    "SouthKorea"     = "#6C3B8E",  # deep purple
    "Thailand"       = "#FF0000",  # deep jade green
    "Australia"      = "#E6AB02",  # deep copper
    "India"          = "#B05E3A",  # deep terracotta
    "Venezuela"      = "#E7298A",  # deep amethyst
    "Spain"          = "#7C8F2A",  # deep olive-lime
    "Denmark"        = "#8C4638",  # deep brick
    "Germany"        = "#3A6F62",  # deep teal-gray
    "UnitedKingdom"  = "#D95F02",  # deep slate-blue
    "Sweden"         = "#3C7F49",  # deep green
    "Russia"         = "#0047AB",  # deep cobalt blue
    "Philippines"    = "#96CA2D",  # rich mustard brown
    "Argentina"      = "#4F7A9A"   # deep ice blue
  )
  continent_colors <- c(
    "Asia"          = "#FF0000",  # deep red
    "Europe"        = "#0047AB",  # deep indigo-blue
    "Africa"        = "#27632A",  # deep green
    "Oceania"       = "#6C3B8E",  # deep royal purple
    "NorthAmerica" = "#8B5A2B",  # deep brown-copper
    "SouthAmerica" = "#9A7E10",  # deep golden-olive
    "Unknown"       = "#4E4E7E"   # deep gray-blue
  )
  linclade <- c(
    "SEA1a_SC1" = "#D95F02",
    "SEA1a_SC1a" = "#27632A", 
    "SEA1a_SC1b" = "#6C3B8E", 
    "SEA1a_SC1c" = "#B22222",
    "SEA1a_SC2"  = "#0047AB", 
    "Unknown"    = "#4E4E7E"
  )
  root_tip = "NC_001542|rabv|NA|ND"
  min_clade_size = 5
  mrsd_auto = NULL
  nodelabsz = 1.5
  nodeptsz = 0.3
  tipptsz = 1.4
  hjustsz = -0.5
  ggtlinesz = 0.1
  color_branches = TRUE
  color_ancestral_branches = TRUE
  root_branch_black = TRUE
  tip_style = "subgeno"
  near_tip_threshold = 0.02
  label_clades = TRUE
  geno_based = "linshort"
  max_label_depth = 0.6
  treescale_x = 0.2
  treescale_y = 100
  treescale_fsz = 3
  treescale_ofst = 150
  radial_scale = 120
  hl_mySarawak = TRUE
  hl_dx = TRUE
  hl_val = 'southwestern'
  hl_col = 'city'
  hl_fill = sublin_colors
  subsrate = 0.05
  subgeno_choice = NULL
  show_tips = TRUE
  show_only_hl_dx = FALSE
  show_node_labels = FALSE
  show_node_points = FALSE
  branch_color_mode = TRUE
  radial_mode = FALSE
  branch_stretch_factor = 1
  fruit_column = 'Serotypes'
  fruitoffset = 0.03
  plot_tt = TRUE
  radial_center_text = NULL
  radial_center_text_size = 6
  radial_center_text_font = "bold"
  fruit_ribbons = FALSE
  fruit_ribbon_width =0.8
  fruit_colors = sublin_colors
  fruit_text = FALSE
  fruit_text_around_ribbon = FALSE
  cladogram_mode = FALSE
  x_axis_break = NULL
  country_colors = continent_colors
lin_val <- c("SEA1a_SC1","SEA1a_SC2","SEA1a_SC1a","SEA1a_SC1b","SEA1a_SC1c")
lin_col <- c("SEA1a_SC1"= "#D95F02","SEA1a_SC2"= "#0047AB","SEA1a_SC1a"= "#27632A","SEA1a_SC1b"= "#91003F","SEA1a_SC1c"= "#FF0000")

city_val <- c("Kuching", "Kuching_Bau", "Kuching_Lundu", "Bintulu", "Kapit", "Limbang",
              "Miri", "Mukah", "Pontianak", "Samarahan", "Sarikei", "Serian", "Sibu", "Sri_Aman")
city_col <- c(
"Kuching" = "#D95F02", 
"Kuching_Bau" = "#D95F02", 
"Kuching_Lundu" = "#D95F02", 
"Bintulu" = "#91003F", 
"Kapit" = "#8B5A2B", 
"Limbang" = "#0047AB", 
"Miri" = "#3C7F49", 
"Mukah" = "#35978F", 
"Pontianak" = "#E7298A", 
"Samarahan" = "#27632A", 
"Sarikei" = "#6C3B8E", 
"Serian" = "#FF0000", 
"Sibu" = "#0047AB", 
"Sri_Aman" = "#5A4A42")
})
### 3A2. MAIN SOURCE CODES: DEFINE SOURCE, FILES, PARAMETERS ###
{
city_col <- lighten(city_col, amount = 0.2)  # lighten by 30%
rabvseq_read <- read_excel("01_src/rabv_glob_metadata.xlsx", sheet = 1)
meta_data = rabvseq_read
legend_labels <- function(capitalize = FALSE) {
  if (capitalize) {
    function(x) toupper(x)
  } else {
    waiver()
  }
}
}
### 3A3. MAIN SOURCE CODES: RABV CDS, G, N PHYLOGENIES ###
{
process_tree_pair <- function(
    iqtree_file,
    tt_file = NULL,
    meta_data,
    root_tip = "NC_001542|rabv|NA|ND",
    min_clade_size = 5,
    mrsd_auto = NULL,
    nodelabsz = 1.5,
    nodeptsz = 0.3,
    tipptsz = 1.4,
    hjustsz = -0.5,
    ggtlinesz = 0.1,
    return_basic_iq_tree = FALSE,
    plot_basic_iq_tree = FALSE,
    capitalize_legend = FALSE,
    color_branches = TRUE,
    color_ancestral_branches = TRUE,
    root_branch_black = TRUE,
    tip_style = "subgeno",
    near_tip_threshold = 0.02,
    label_clades = TRUE,
    geno_based = "lin",
    max_label_depth = 0.6,
    treescale_x = 0.2,
    treescale_y = 100,
    treescale_fsz = 3,
    treescale_ofst = 150,
    radial_scale = 120,
    hl_mySarawak = TRUE,
    hl_dx = TRUE,
    hl_val = city_val,
    hl_col = "sar_city",
    hl_fill = city_col,
    subsrate = 0.05,
    subgeno_choice = NULL,
    show_tips = TRUE,
    show_only_hl_dx = FALSE,
    show_node_labels = FALSE,
    show_node_points = FALSE,
    branch_color_mode = TRUE,
    radial_mode = FALSE,
    branch_stretch_factor = 1,
    fruit_column = 'lin',
    fruitoffset = 0.03,
    plot_tt = TRUE,
    radial_center_text = NULL,
    radial_center_text_size = 6,
    radial_center_text_font = "bold",
    fruit_ribbons = FALSE,
    fruit_ribbon_width =0.8,
    fruit_colors = subclade_colors,
    fruit_text = FALSE,
    fruit_text_around_ribbon = FALSE,
    cladogram_mode = FALSE
){
  suppressPackageStartupMessages({
    library(ape); library(phytools); library(readr); library(dplyr)
    library(ggtree); library(ggforce); library(ggtreeExtra); library(tidyverse); library(lubridate)
    library(phangorn);library(treeio)
    library(tidyverse); library(lubridate); library(colorspace)
    library(ggrepel); library(ggnewscale); library(viridis); library(patchwork)
  })

  # Ensure tips are always shown if we want to highlight
  if (show_only_hl_dx) show_tips <- TRUE
  
  read_tree <- function(file, type){
    if (type == "newick") read.tree(file) else read.nexus(file)
  }
  
  iq_phylo <- read_tree(iqtree_file, "newick")
  tt_phylo <- if (!is.null(tt_file)) read_tree(tt_file, "nexus") else NULL
  
  if (!is.null(tt_phylo)) {
    shared <- intersect(iq_phylo$tip.label, tt_phylo$tip.label)
    iq_phylo <- keep.tip(iq_phylo, shared)
    tt_phylo <- keep.tip(tt_phylo, shared)
  }
  
  root_tree_safe <- function(phy, root_tip, tree_name = "tree") {
    if (is.null(root_tip)) return(phy)
    if (root_tip %in% phy$tip.label) {
      phy <- root(phy, outgroup = root_tip, resolve.root = TRUE)
      phy <- reorder.phylo(phy, order = "cladewise")
    } else {
      warning(sprintf("Root tip '%s' not found in %s", root_tip, tree_name))
    }
    phy
  }
  
  iq_phylo <- root_tree_safe(iq_phylo, root_tip, "IQ-TREE")
  if (return_basic_iq_tree || plot_basic_iq_tree) {
    
    if (is.null(iq_phylo$node.label)) {
      iq_phylo$node.label <- rep(NA, iq_phylo$Nnode)
    }
    
    # convert node labels to numeric for filtering
    node_labels_num <- suppressWarnings(as.numeric(iq_phylo$node.label))
    # set labels <=75 or NA to NA
    node_labels_filtered <- ifelse(!is.na(node_labels_num) & node_labels_num > 75,
                                   node_labels_num, NA)
    
    # create basic ggtree object
    basic_tree_plot <- ggtree(iq_phylo) +
      #geom_tiplab(size = tipptsz) +
      geom_nodepoint(aes(subset = !isTip), size = nodeptsz) +
      geom_nodelab(size = nodelabsz)
    
    # if only plotting is requested, return both tree and plot
    if (plot_basic_iq_tree) {
      return(list(
        basic_iq_tree = iq_phylo,
        basic_iq_plot = basic_tree_plot
      ))
    } else {
      # only returning tree
      return(list(basic_iq_tree = iq_phylo))
    }
  }
  if (!is.null(tt_phylo)) tt_phylo <- root_tree_safe(tt_phylo, root_tip, "TT tree")
  
  get_clade_tips <- function(tr){
    m <- length(tr$tip.label)
    n <- tr$Nnode
    lapply((m+1):(m+n), function(i){
      sort(tr$tip.label[Descendants(tr, i, "tips")[[1]]])
    })
  }
  
  transfer_bootstrap <- function(src, tgt){
    src_clades <- get_clade_tips(src)
    tgt_clades <- get_clade_tips(tgt)
    match_idx <- match(sapply(tgt_clades, paste, collapse=","), sapply(src_clades, paste, collapse=","))
    tgt$node.label <- rep(NA, tgt$Nnode)
    for (i in seq_along(match_idx)){
      j <- match_idx[i]
      if (!is.na(j) && length(src$node.label) >= j && !is.null(src$node.label[j])){
        tgt$node.label[i] <- src$node.label[j]
      }
    }
    tgt
  }
  
  if (!is.null(tt_phylo)) tt_phylo <- transfer_bootstrap(iq_phylo, tt_phylo)
  
  iq_td <- as.treedata(iq_phylo)
  tt_td <- if (!is.null(tt_phylo)) as.treedata(tt_phylo) else NULL
  if (is.null(mrsd_auto)) mrsd_auto <- max(ymd(meta_data$col_date), na.rm = TRUE)
  
  if (is.null(subgeno_choice)) subgeno_colors <- c("na"="grey70") else {
    subgeno_colors <- subgeno_choice
    if (!("na" %in% names(subgeno_colors))) subgeno_colors <- c(subgeno_colors, na="grey70")
  }
  branch_colors <- sapply(subgeno_colors, function(z) lighten(z, 0.25))
  
  prepare_tbl <- function(td){
    tbl <- fortify(td)
    
    # Merge meta_data
    tip_info <- meta_data %>%
      filter(seq_name %in% tbl$label) %>%
      select(seq_name, all_of(geno_based), country, continent, !!hl_col)
    colnames(tip_info)[1] <- "label"
    tbl <- tbl %>% left_join(tip_info, by="label")
    
    # Clade size calculation
    tbl$clade_size <- 1L
    ntip <- sum(tbl$isTip)
    for (i in seq(nrow(tbl), ntip+1)){
      kids <- which(tbl$parent == tbl$node[i])
      tbl$clade_size[i] <- sum(tbl$clade_size[kids])
    }
    
    # Fill internal nodes with dominant geno_based
    for (i in which(!tbl$isTip)) {
      desc_tips <- Descendants(td@phylo, tbl$node[i], type = "tips")[[1]]
      vals <- tbl[[geno_based]][desc_tips]
      vals <- vals[!is.na(vals) & vals != "na"]
      tbl[[geno_based]][i] <- if(length(vals)>0) names(which.max(table(vals))) else names(subgeno_colors)[1]
    }
    
    # Propagate downwards
    root_node <- min(tbl$node)
    queue <- root_node
    while(length(queue) > 0){
      nd <- queue[1]; queue <- queue[-1]
      mygeno <- tbl[[geno_based]][tbl$node == nd]
      kids <- tbl$node[tbl$parent == nd]
      for (c in kids){
        if (is.na(tbl[[geno_based]][tbl$node == c])) tbl[[geno_based]][tbl$node == c] <- mygeno
        queue <- c(queue, c)
      }
    }
    
    # Factorize geno_based
    tbl[[geno_based]][is.na(tbl[[geno_based]])] <- "na"
    tbl[[geno_based]][!(tbl[[geno_based]] %in% names(subgeno_colors))] <- "na"
    tbl[[geno_based]] <- factor(tbl[[geno_based]], levels=names(subgeno_colors))
    
    # Safe dx_flag creation
    if (!hl_col %in% colnames(tbl)) tbl[[hl_col]] <- NA
    tbl <- tbl %>% mutate(
      country_flag = ifelse(country=="Sarawak","Sarawak","Other"),
      dx_flag      = ifelse(!is.na(.data[[hl_col]]) & .data[[hl_col]] %in% hl_val,
                            .data[[hl_col]], "Other")
    )
    
    tbl
  }
  tree_tbl_iq <- prepare_tbl(iq_td)
  tree_tbl_tt <- if (!is.null(tt_td)) prepare_tbl(tt_td) else NULL
  
  stretch_tree <- function(td, fact){
    if (fact == 1) return(td)
    phy <- td@phylo
    nd  <- node.depth.edgelength(phy)
    md  <- max(nd)
    parent_depth <- nd[phy$edge[,1]]
    phy$edge.length <- phy$edge.length * (1 + (fact-1)*(parent_depth/md))
    td@phylo <- phy
    td
  }
  
  rescale_tree_for_radial <- function(td, scale_factor){
    phy <- td@phylo
    phy$edge.length <- phy$edge.length * scale_factor
    td@phylo <- phy
    td
  }
  
  iq_td <- stretch_tree(iq_td, branch_stretch_factor)
  if (!is.null(tt_td)) tt_td <- stretch_tree(tt_td, branch_stretch_factor)
  if (radial_mode) {
    iq_td <- rescale_tree_for_radial(iq_td, radial_scale)
    if (!is.null(tt_td)) tt_td <- rescale_tree_for_radial(tt_td, radial_scale)
  }
  
  color_edges <- function(tbl){
    tbl$edge_color <- as.character(tbl[[geno_based]])
    tbl$edge_color[is.na(tbl$edge_color)] <- "na"
    if (!color_ancestral_branches) tbl$edge_color[!tbl$isTip] <- "black"
    if (root_branch_black) {
      root_edge_idx <- which(tbl$parent == 0)
      if (length(root_edge_idx) == 1) tbl$edge_color[root_edge_idx] <- "black"
    }
    tbl
  }
  
  tree_tbl_iq <- color_edges(tree_tbl_iq)
  tree_tbl_tt <- if (!is.null(tree_tbl_tt)) color_edges(tree_tbl_tt) else NULL
  
  make_tree_plot <- function(td, tbl){
    phy <- td@phylo
    nd  <- node.depth.edgelength(phy)
    tip_depths <- nd[1:length(phy$tip.label)]
    tip_pos <- setNames(tip_depths, phy$tip.label)
    
    # Initial tip positions
    tbl <- tbl %>% mutate(
      x_tip = ifelse(isTip, tip_pos[label], NA_real_),
      y_tip = y
    )
    
    # Base ggtree object
    if(cladogram_mode) p <- ggtree(td, branch.length = "none") %<+% tbl
    else p <- ggtree(td) %<+% tbl
    
    # Reassign tip positions for cladogram
    if (cladogram_mode && nrow(tbl) > 0){
      tip_coords <- p$data %>% filter(isTip)
      tbl <- tbl %>% left_join(
        tip_coords %>% select(label, x, y),
        by = "label",
        suffix = c("", "_clad")
      )
      tbl <- tbl %>% mutate(
        x_tip = ifelse(isTip, x_clad, x_tip),
        y_tip = ifelse(isTip, y_clad, y_tip)
      )
    }
    
    # Branch coloring
    if (color_branches || branch_color_mode) {
      p <- p +
        geom_tree(
          aes(color = edge_color),
          linewidth = ggtlinesz,
          lineend = "round"
        ) +
        scale_color_manual(
          values = branch_colors,
          labels = legend_labels(capitalize_legend),
          drop = FALSE
        )
    } else {
      p <- p + geom_tree(linewidth = ggtlinesz, color = "white")
    }
    # Add treescale for both radial and non-radial trees
    if (!is.null(tbl) && nrow(tbl) > 0) {
      
      if (radial_mode) {
        # radial tree: place scale near center, offset label for visibility
        p <- p + theme(axis.line.x = element_blank()) +
          geom_treescale(
            x = treescale_x,             # small offset from center
            y = treescale_y,                # center
            width = subsrate,         # 0.02 substitutions per site
            linesize = 0.05,
            fontsize = treescale_fsz,
            color = "black",
            offset = treescale_ofst,        # offset pushes label slightly above the bar
            label = " substitutions per site"
          )
      } else {
        # non-radial tree: bottom-left, offset label slightly
        p <- p + geom_treescale(
          x = min(tbl$x_tip, na.rm = TRUE),
          y = min(tbl$y, na.rm = TRUE),
          width = subsrate,
          linesize = 0.05,
          fontsize = 3,
          color = "black",
          offset = treescale_ofst,        # offset pushes label slightly above the bar
          label = " substitutions per site"
        )
      }
    }
    if (hl_dx) {
      # Prepare a palette: named vector linking each hl_val to a color
      if (length(hl_val) > 1) {
        if (is.null(names(hl_fill))) {
          # If hl_fill is unnamed, assign colors in order
          if (length(hl_fill) < length(hl_val)) {
            stop("hl_fill must have at least as many colors as hl_val values")
          }
          hl_palette <- setNames(hl_fill[seq_along(hl_val)], hl_val)
        } else {
          hl_palette <- hl_fill
        }
      } else {
        # single value
        hl_palette <- setNames(hl_fill, hl_val)
      }
      
      if (show_only_hl_dx) {
        tip_df <- tbl %>%
          filter(isTip & country_flag == "Sarawak" & .data[[hl_col]] %in% hl_val)
      } else {
        tip_df <- tbl %>% filter(isTip & country_flag == "Sarawak")
      }
      
      if (nrow(tip_df) > 0) {
        tip_df$tip_color <- ifelse(
          tip_df[[hl_col]] %in% hl_val,
          hl_palette[tip_df[[hl_col]]],
          "white"
        )
        
        p <- p + geom_point(
          data = tip_df,
          aes(x = x_tip, y = y_tip, fill = tip_color),
          shape = 21, color = "black",
          size = tipptsz, stroke = 0.3
        ) +
          scale_fill_identity()
      }
      
    } else if (show_tips) {
      tip_df <- tbl %>% filter(isTip)
      
      tip_fill <- case_when(
        tip_style == "clinical_dx" & tip_df[[hl_col]] %in% hl_val ~ hl_palette[tip_df[[hl_col]]],
        TRUE ~ as.character(tip_df[[geno_based]])
      )
      
      if (nrow(tip_df) > 0) {
        p <- p + geom_point(
          data = tip_df,
          aes(x = x_tip, y = y_tip, fill = tip_fill),
          shape = 21, size = tipptsz, color = "black"
        ) +
          scale_fill_manual(values = c(subgeno_colors, hl_palette),
                            labels = legend_labels(capitalize_legend))
      }
    }
    
    # Fruit ribbons / annotation
    if ((fruit_ribbons || fruit_text) && !is.null(fruit_column) && nrow(tbl)>0){
      if (is.null(fruit_colors)){
        fruit_colors <- structure(distinct(tbl, !!sym(fruit_column))[[fruit_column]],
                                  .Names=distinct(tbl, !!sym(fruit_column))[[fruit_column]])
      }
      
      if (fruit_ribbons){
        p <- p + geom_fruit(
          geom=geom_tile,
          mapping=aes_string(fill=fruit_column),
          offset=fruitoffset,
          width=fruit_ribbon_width,
          color=NA
        )+
          scale_fill_manual(
            values = fruit_colors,
            labels = legend_labels(capitalize_legend)
          )
      } else {
        p <- p + geom_fruit(
          geom=geom_point,
          mapping=aes_string(color=fruit_column),
          offset=fruitoffset,
          size=tipptsz*0.8
        ) + scale_color_manual(
          values = fruit_colors,
          labels = legend_labels(capitalize_legend)
        )
      }
      
      if (fruit_text & fruit_text_around_ribbon){
        fruit_mid <- tbl %>% group_by(.data[[fruit_column]]) %>%
          summarise(x_mid = mean(x_tip), y_mid = mean(y_tip))
        fruit_mid <- fruit_mid %>% rowwise() %>%
          mutate(
            angle = atan2(y_mid, x_mid) * 180/pi,
            radius = sqrt(x_mid^2 + y_mid^2),
            unit_x = x_mid/radius, unit_y = y_mid/radius,
            x_mid = x_mid + (fruit_ribbon_width+fruitoffset)*unit_x,
            y_mid = y_mid + (fruit_ribbon_width+fruitoffset)*unit_y,
            angle_tangent = angle + 90,
            angle_tangent = ifelse(angle_tangent>90, angle_tangent-180,
                                   ifelse(angle_tangent < -90, angle_tangent+180, angle_tangent)),
            hjust = ifelse(angle_tangent>90 | angle_tangent< -90, 1,0)
          )
        p <- p + geom_text(data=fruit_mid,
                           aes(x=x_mid, y=y_mid, label = if (capitalize_legend)
                             toupper(.data[[fruit_column]])
                             else
                               .data[[fruit_column]],
                               angle=angle_tangent, hjust=hjust),
                           size=2.5, fontface="bold", color="black", vjust=0.5)
      }
    }
    
    # Radial layout
    if (radial_mode) p <- p + layout_circular() + theme_void()
    if (!is.null(radial_center_text)) p <- p + annotate("text", x=0, y=0,
                                                        label=radial_center_text,
                                                        size=radial_center_text_size,
                                                        fontface=radial_center_text_font)
    p
  }
  
  p_iq <- make_tree_plot(iq_td, tree_tbl_iq)
  p_tt <- if (!is.null(tt_td) && plot_tt) make_tree_plot(tt_td, tree_tbl_tt)
  
  list(
    iq_plot = p_iq,
    tt_plot = p_tt,
    iq_tree = iq_phylo,
    tt_tree = tt_phylo,
    iq_table = tree_tbl_iq,
    tt_table = tree_tbl_tt
  )
}
  library(ggplot2)
  library(ggpubr)
res1 <- process_tree_pair(
  iqtree_file = "05_tree/rabv_384globloc_cds_align.fasta.treefile",
  tt_file = "05_tree/2026-01-02_rabv384.nexus", subgeno_choice = sublin_colors,root_tip ='NC_001542|rabv|NA|ND',show_only_hl_dx = FALSE,
  treescale_x = 0.5, treescale_y = -10, treescale_fsz= 1.8, treescale_ofst = 120,cladogram_mode = TRUE,
  meta_data = meta_data, geno_based = 'linshort',radial_mode = TRUE, show_tips = TRUE, hl_dx = TRUE, radial_scale = 10, subsrate = 2,
  nodelabsz = 1.6,nodeptsz =0.2,tipptsz = 1.0, ggtlinesz = 0.4, hjustsz=-0.2,max_label_depth = 3, branch_stretch_factor = 0.0,
  color_branches = TRUE, tip_style = "subgenotype", near_tip_threshold = 0.2, hl_val = lin_val,hl_col = "linclade",hl_fill = lin_col,
)

res2 <- process_tree_pair(
  iqtree_file = "05_tree/rabv_384globloc_G_align.fasta.treefile",
  tt_file = "05_tree/2026-01-02_rabv384.nexus", subgeno_choice = sublin_colors,root_tip ='NC_001542|rabv|NA|ND',show_only_hl_dx = FALSE,
  treescale_x = 0.5, treescale_y = -10, treescale_fsz= 1.8, treescale_ofst = 120,cladogram_mode = TRUE,
  meta_data = meta_data, geno_based = 'linshort',radial_mode = TRUE, show_tips = TRUE, hl_dx = TRUE, radial_scale = 10, subsrate = 2,
  nodelabsz = 1.6,nodeptsz =0.2,tipptsz = 1.0, ggtlinesz = 0.4, hjustsz=-0.2,max_label_depth = 3, branch_stretch_factor = 0.0,
  color_branches = TRUE, tip_style = "subgenotype", near_tip_threshold = 0.2, hl_val = lin_val,hl_col = "linclade",hl_fill = lin_col,
)

res3 <- process_tree_pair(
  iqtree_file = "05_tree/rabv_417globloc_N_align.fasta.treefile",
  tt_file = "05_tree/2026-01-02_rabv384.nexus", subgeno_choice = sublin_colors,root_tip ='NC_001542|rabv|NA|ND',show_only_hl_dx = FALSE,
  treescale_x = 0.5, treescale_y = -10, treescale_fsz= 1.8, treescale_ofst = 120,cladogram_mode = TRUE,
  meta_data = meta_data, geno_based = 'linshort',radial_mode = TRUE, show_tips = TRUE, hl_dx = TRUE, radial_scale = 10, subsrate = 2,
  nodelabsz = 1.6,nodeptsz =0.2,tipptsz = 1.0, ggtlinesz = 0.4, hjustsz=-0.2,max_label_depth = 3, branch_stretch_factor = 0.0,
  color_branches = TRUE, tip_style = "subgenotype", near_tip_threshold = 0.2, hl_val = lin_val,hl_col = "linclade",hl_fill = lin_col,
)
ggsave("output/rabv_384globloc_cds_plot1b.png", plot=res1$iq_plot, width=9.6, height=12, dpi=300)
ggsave("output/rabv_384globloc_G_plot1b.png", plot=res2$iq_plot, width=9.6, height=12, dpi=300)
ggsave("output/rabv_384globloc_N_plot1b.png", plot=res3$iq_plot, width=9.6, height=12, dpi=300)
phyloplot <- ggarrange(res1$iq_plot, res2$iq_plot, res3$iq_plot, 
          ncol = 3, nrow = 1,
          labels = c("A", "B", "C"),         # labels for plots
          widths = c(1, 1, 1),                # relative widths of plots
          common.legend = TRUE, legend = "bottom",
          hjust = -0.5)
ggsave("output/rabv_384globloc_phylocomb.png", plot=phyloplot, width=14, height=5, dpi=300)


}

### 3A3. MAIN SOURCE CODES: RABV MCC CDS, G, N PHYLOGENY ###
{
seq_col="seq_name"
root_tree=FALSE
root_tip=NULL
branch_lab="linclade"
tip_lab="sar_city"
branch_palette=linclade
tip_palette=city_col
tmrca_date=2010.126
show_node_dates=TRUE
node_label_mode="subclade"
extra_annot_nodes=NULL
show_hpd_ribbon=FALSE
show_state_pies=FALSE
state_col="state"
state_prob_col="state.prob"
radial_mode=FALSE
cladogram_mode=FALSE
ggtlinesz=0.3
tipptsz=2
odelabsz=3
fruit_column="linclade"
fruit_colors=linclade
fruit_offset=0.05
fruit_width=0.3
col_darken=0.1
hpd_digits=2
tr_linesz = 0.3
x_axis_interval=2
node_label_hjust=0.9
node_label_vjust=-2
node_label_size = 3.4
node_label_color = "branch"
node_label_repel = TRUE
basic_node_tree=FALSE
return_data=TRUE

process_mcc_tree <- function(
    beast_tree_file, meta_data=NULL, iqtree_file=NULL, seq_col="seq_name",
    root_tree=FALSE, root_tip=NULL,
    branch_lab="linclade", tip_lab="linclade",
    branch_palette=NULL, tip_palette=NULL,
    tmrca_date=2010.126,
    show_node_dates=TRUE,
    node_label_mode=c("all","internal","subclade"),
    extra_annot_nodes=NULL,
    show_tippoints = TRUE,
      show_hpd_ribbon=FALSE,
    show_state_pies=FALSE,
    state_col="state", state_prob_col="state.prob",
    radial_mode=FALSE, cladogram_mode=FALSE,
    ggtlinesz=0.3, tipptsz=2, nodelabsz=3,
    fruit_column=NULL, fruit_colors=NULL,
    fruit_offset=0.05, fruit_width=0.3,
    col_darken=0.1,
    hpd_digits=2,
    tr_linesz = 0.3,
    x_axis_interval=2,
    node_label_hjust=0.9,
    node_label_vjust=-2,
    node_label_size = 3.4,
    node_label_color = "branch",
    node_label_repel = TRUE,
    basic_node_tree=FALSE,
    return_data=TRUE
){
  suppressPackageStartupMessages({
    library(treeio); library(ggtree); library(ggtreeExtra)
    library(ape); library(phytools)
    library(dplyr); library(tibble)
    library(RColorBrewer); library(scales); library(colorspace)
    library(ggrepel)
  })
  node_label_mode <- match.arg(node_label_mode)
  auto_palette <- function(x,pal="Set2"){
    ux <- sort(unique(na.omit(x))); if(length(ux)==0) return(NULL)
    cols <- brewer.pal(max(3,min(8,length(ux))),pal); setNames(cols[seq_along(ux)],ux)
  }
  sanitize_phylo <- function(phy){
    phy <- reorder.phylo(phy,"postorder")
    if(is.null(phy$edge.length)) phy$edge.length <- rep(1e-6,nrow(phy$edge))
    phy$edge.length[is.na(phy$edge.length)|phy$edge.length<=0] <- 1e-6
    if(is.null(phy$tip.label)||any(is.na(phy$tip.label))) stop("Tree tips missing labels")
    if(is.null(phy$node.label)||length(phy$node.label)!=phy$Nnode) phy$node.label <- rep(NA_character_,phy$Nnode)
    phy
  }
  # ---------------- Read and root tree ----------------
  td <- read.beast(beast_tree_file); phy <- td@phylo
  if(root_tree && !is.null(root_tip) && root_tip %in% phy$tip.label)
    phy <- root(phy,root_tip,resolve.root=TRUE)
  phy <- sanitize_phylo(phy)
  
  # ---------------- Basic node-only tree ----------------
  basic_tree_plot <- NULL
  if(basic_node_tree){
    basic_tbl <- as_tibble(ggtree::fortify(phy))
    basic_tree_plot <- ggtree(phy, branch.length=if(cladogram_mode) "none" else "branch.length") +
      geom_text2(data=dplyr::filter(basic_tbl,!isTip), aes(label=node), size=nodelabsz) +
      theme_tree2() + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
    if(radial_mode) basic_tree_plot <- basic_tree_plot + layout_circular()
  }
  
  # ---------------- Annotations ----------------
  tbl_anno <- td@data %>% as.data.frame(check.names=FALSE) %>%
    {colnames(.)<-make.unique(colnames(.));.} %>% as_tibble()
  if(!"node"%in%colnames(tbl_anno)) tbl_anno <- rownames_to_column(tbl_anno,"node")
  tbl_anno$node <- as.integer(tbl_anno$node)
  
  parse_hpd <- function(x){
    x <- gsub("[c\\(\\)\"]","",as.character(x))
    nums <- as.numeric(unlist(strsplit(x,"[,\\s]+")))
    if(length(nums)<2) nums <- c(nums[1],nums[1])
    nums[1:2]
  }
  
  hpd_cols <- grep("HPD|range",colnames(tbl_anno),value=TRUE)
  for(col in hpd_cols){
    vals <- lapply(tbl_anno[[col]],parse_hpd)
    base <- gsub("\\.","_",col)
    tbl_anno[[paste0(base,"_low")]] <- sapply(vals,`[`,1)
    tbl_anno[[paste0(base,"_high")]] <- sapply(vals,`[`,2)
  }
  tbl_anno <- tbl_anno %>% select(where(~!is.list(.))) %>% distinct(node,.keep_all=TRUE)
  
  # ---------------- Convert heights to calendar years ----------------
  if("height_median"%in%colnames(tbl_anno)){
    max_h <- max(tbl_anno$height_median,na.rm=TRUE)
    tbl_anno <- tbl_anno %>% mutate(
      node_date = as.numeric(tmrca_date) + (max_h - height_median)
    )
    if(all(c("height_0_95_HPD_low","height_0_95_HPD_high")%in%colnames(tbl_anno))){
      tbl_anno <- tbl_anno %>% mutate(
        year_low  = as.numeric(tmrca_date) + (max_h - height_0_95_HPD_high),
        year_high = as.numeric(tmrca_date) + (max_h - height_0_95_HPD_low),
        node_date_label = paste0(round(node_date,0)," [",round(year_low,hpd_digits),"-",round(year_high,hpd_digits),"]")
      )
    } else tbl_anno <- tbl_anno %>% mutate(node_date_label = as.character(round(node_date,0)))
  }
  
  # ---------------- Optional IQ-TREE labels ----------------
  if(!is.null(iqtree_file)){
    iq_phy <- sanitize_phylo(read.tree(iqtree_file))
    shared <- intersect(iq_phy$tip.label,phy$tip.label)
    phy <- keep.tip(phy,shared); iq_phy <- keep.tip(iq_phy,shared)
    get_clades <- function(tr) lapply((length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode), function(i) sort(tr$tip.label[Descendants(tr,i,"tips")[[1]]]))
    tgt <- get_clades(phy); src <- get_clades(iq_phy)
    m <- match(sapply(tgt,paste,collapse=","),sapply(src,paste,collapse=","))
    phy$node.label <- rep(NA_character_,phy$Nnode)
    for(i in seq_along(m)) if(!is.na(m[i])) phy$node.label[i] <- iq_phy$node.label[m[i]]
  }
  
  # ---------------- Fortify main tree ----------------
  tbl <- as_tibble(ggtree::fortify(phy))
  tbl <- left_join(tbl, tbl_anno, by="node")
  if(!is.null(meta_data)&&"label"%in%colnames(tbl)&&seq_col%in%colnames(meta_data))
    tbl <- left_join(tbl, meta_data, by=c("label"=seq_col))
  
  # ---------------- Dominant branch assignment ----------------
  tbl$dominant_branch <- sapply(tbl$node,function(n){
    if(n<=length(phy$tip.label)) return(as.character(tbl[[branch_lab]][tbl$node==n]))
    tips <- Descendants(phy,n,"tips")[[1]]; vals <- tbl[[branch_lab]][tips]; vals <- vals[!is.na(vals)]
    if(length(vals)==0) return(NA_character_)
    tab <- table(vals); dom <- names(tab)[which.max(tab)]
    if(tab[dom]/sum(tab) > 0.5) dom else NA_character_
  })
  
  # ---------------- Palettes ----------------
  if(is.null(branch_palette)) branch_palette <- auto_palette(tbl$dominant_branch)
  if(is.null(tip_palette)) tip_palette <- auto_palette(tbl[[tip_lab]])
  if(col_darken!=0){branch_palette <- darken(branch_palette,col_darken); tip_palette <- darken(tip_palette,col_darken)}
  
  # ---------------- Main tree plot ----------------
  p <- ggtree(phy, branch.length=if(cladogram_mode)"none" else "branch.length", linewidth=tr_linesz) %<+% tbl +
    geom_tree(aes(color=dominant_branch), linewidth=ggtlinesz) +
    scale_color_manual(values=branch_palette, na.value="grey70")
  if(show_tippoints){ p <- p + 
    geom_tippoint(aes(fill=.data[[tip_lab]]), shape=21, size=tipptsz, color="black") +
    scale_fill_manual(values=tip_palette)
  }
  # ---------------- Node date labels ----------------
  if(show_node_dates && "node_date_label"%in%colnames(tbl)){
    lbl <- tbl %>% filter(!isTip)
    if(node_label_mode=="subclade") lbl <- lbl %>% filter(!duplicated(dominant_branch))
    if(!is.null(extra_annot_nodes)) lbl <- bind_rows(lbl, tbl %>% filter(node %in% extra_annot_nodes))
    lbl <- distinct(lbl,node,.keep_all=TRUE)
    lbl$label_color <- if(node_label_color=="branch") lbl$dominant_branch else node_label_color
    
    if(node_label_repel){
      p <- p + ggrepel::geom_text_repel(
        data = lbl,
        aes(x=x, y=y, label=node_date_label, color=label_color),
        size=node_label_size,
        hjust=node_label_hjust,
        vjust=node_label_vjust,
        fontface="bold",
        segment.color="black",
        segment.alpha=0.5,
        show.legend=FALSE
      )
    } else {
      p <- p + geom_text(
        data = lbl,
        aes(x=x, y=y, label=node_date_label, color=label_color),
        size=node_label_size,
        hjust=node_label_hjust,
        show.legend=FALSE
      )
    }
  }
  
  # ---------------- HPD ribbon ----------------
  if(show_hpd_ribbon && all(c("year_low","year_high")%in%colnames(tbl)))
    p <- p + geom_segment(data=tbl %>% filter(!isTip),
                          aes(x=year_low, xend=year_high, y=y, yend=y), linewidth=1, alpha=0.4)
  
  # ---------------- State pies ----------------
  if(show_state_pies && all(c(state_col,state_prob_col)%in%colnames(tbl)))
    p <- p + geom_point(data=tbl %>% filter(!isTip),
                        aes(fill=.data[[state_col]], size=.data[[state_prob_col]], x=x),
                        shape=21, alpha=0.9) + scale_size(range=c(1,5))
  
  # ---------------- Fruit annotations ----------------
  if(!is.null(fruit_column) && fruit_column %in% colnames(tbl)){
    
    if(is.null(fruit_colors)) fruit_colors <- auto_palette(tbl[[fruit_column]])
    p <- p + ggnewscale::new_scale_fill() +
      geom_fruit(
        geom = geom_tile,
        mapping = aes_string(fill = fruit_column),
        offset = fruit_offset,
        width = fruit_width
      ) +
      scale_fill_manual(values = fruit_colors)
  }
  
  # ---------------- X-axis scale by calendar year ----------------
  if(!radial_mode && !is.null(x_axis_interval) && "x"%in%colnames(tbl)){
    x_min <- min(tbl$x, na.rm=TRUE)
    x_max <- max(tbl$x, na.rm=TRUE)
    x_breaks <- seq(floor(x_min), ceiling(x_max), by=x_axis_interval)
    
    # Correct: tip x=0 = tmrca_date
    year_labels <- as.integer(as.numeric(tmrca_date) + x_breaks)
    
    p <- p + theme_tree2() +
      scale_x_continuous(
        name = "Year",
        breaks = x_breaks,
        labels = year_labels
      )
  }
  if(radial_mode) p <- p + layout_circular() + theme_void()
  if(return_data) return(list(plot=p, tree=phy, table=tbl, basic_tree=basic_tree_plot))
  return(p)
}

rabv_mcc <- process_mcc_tree(
  beast_tree_file = "05_tree/rabv_240loccds_BSkyL_GTRG_RCLN_300m15k.treeanno.txt",
  meta_data = readxl::read_xlsx("01_src/rabv_glob_metadata.xlsx", sheet = 1),
  tr_linesz = 0.4,
  tip_lab = "sar_city",
  tip_palette = city_col,
  branch_lab = "linclade",
  branch_palette = linclade,
  tmrca_date = "2010.126",
  col_darken=0.2,
  fruit_column = "linclade",
  fruit_colors = linclade,
  node_label_mode = "subclade",
  extra_annot_nodes = c(262, 395, 450),
  x_axis_interval = 2
)
ggsave("output/rabv_240loccds_BSkyL_GTRG_RCLN_treeplot.png", plot=rabv_mcc$p, width=9.5, height=11, dpi=300)
}

################################### 3B. RABV PHYLOGENY: ANALYSIS FOR SMALL, LOCAL DATASETS ##############################
### 3B1. MAIN SOURCE CODES: DEFINE PATHS TO TREE FILES, SOURCES  ###
suppressPackageStartupMessages({
read_df<- read.table(file = "04_fN/fN-gene_localredglobalred_meta.tsv", sep = '\t', header = TRUE)
metadata_df <- as.data.frame(read_df)
read_fG<- read.table(file = "05_fG/fG_localredglobalred_meta.tsv", sep = '\t', header = TRUE)
metadata_fG <- as.data.frame(read_fG)
tree <- read.newick("02_bootstraps/pN-gene_locgloberabvsub1_2022.fasta.treefile", node.label='support')
tree <- ape::root(tree, outgroup = "1970s/U22627/unknown/Egypt/Africa4")
treeN <- read.newick("02_bootstraps/fN-gene_locgloberabv_m7k2_ml.nwk", node.label='support')
treeN <- ape::root(treeN, outgroup = "1970s/U22627/unknown/Egypt/Africa4")
treeG <- read.newick("02_bootstraps/fcds_locgloberabv_2023.fasta.treefile", node.label='support')
treeG <- ape::root(treeG, outgroup = "1950/KF154998/Dog/Israel/1950")
fasta.file.name <- system.file("extdata", "01_src/swk_rabv_alignred_pub.fasta", package = "rhierbaps")
})
### 3B2. MAIN ANALYSIS CODES: RABV_COMPLETE CDS DISTREE  ###
{
options(max.print=1000)
p0 <- ggtree(treeG, size=0.75, show.legend = TRUE) %<+% metadata_fG + theme_tree2()+
  geom_tippoint(aes(color=Lineage), size=3.5)+
  scale_x_continuous(breaks=c(0.00,0.10,0.20,0.30), labels=c("0.00", "0.10", "0.20","0.30"))+
  expand_limits(x = c(0, 0.35), y=c(0,133))+
  geom_nodelab(geom='label', aes(label=as.numeric(support), subset=as.numeric(support) > 70),size = 2.5, label.padding = unit(0.1, "lines"),
               alpha=0.7, nudge_x = -0.005)+
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, hjust=-0.02, size = 2.8,fontface = "bold")+ 
  theme(axis.text = element_text(vjust=0.5,size=10, face = "bold", angle = 90), 
        legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.direction="vertical", 
        #legend.key = element_blank(),
        legend.position = "left", legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.1,'cm'),
        legend.margin = margin(2, -4, 14, 1, "cm"),
        #axis.ticks.length = unit(0.25, "cm")+
        plot.margin = unit(c(-10,0,0,0), "lines"))+  
        geom_strip('RV74/2-1','RV59-2-2', offset = 0.195, barsize=2.5, color='#ff854d', 
        label="RABV isolated \n from Sarawak", offset.text=-0.15,size = 8.5, face="bold")+
  guides(fill=guide_legend(title="Rabies Lineages"))+
  coord_cartesian(xlim = c(-0.02, 0.48), # This focuses the x-axis on the range of interest
                  clip = 'off')
panel0 <- p0 %<+% metadata_fG +
  geom_tippoint(
    aes(subset=(Category=='Asian')),fill='#F2A922',size=3.0,
    color='#704c07',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Cosmopolitan')),fill='#267d56',size=3.0,
                color='#1a563b',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Arctic-related')),fill='#1bbbd0',size=3.0,
                color='#157aa3',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Africa')),fill='#f55a76',size=3.0,
                color='#83081e',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Indian_subcontinent')),fill='#c3699d',size=3.0,
                color='#bc575f',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  coord_cartesian(xlim = c(0, 0.26), # This focuses the x-axis on the range of interest
                  clip = 'off')+
  scale_fill_manual(values=c("#F2A922","#267d56","#FF0000",
                             "#1bbbd0", "#f55a76","#c3699d"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)
scaleClade(p0, 353, .3) %>% collapse(353, 'min')
ggsave('output/swk_fG_globloc_2023.png',width = 23, height = 36, units = "cm",limitsize = FALSE, device = 'png')
}
### 3B3. MAIN ANALYSIS CODES: RABV COMPLETE NUCLEOPROTEIN GENE DISTREE  ###
{
p2 <- ggtree(treeN, size=0.35, show.legend = TRUE) %<+% metadata_df + theme_tree2()+
  geom_tippoint(aes(fill=Lineage), size=2.2,stroke=0.4,color='#704c07',shape=21)+
  scale_color_manual(values=region_cols)+
  scale_x_continuous(breaks = seq(from = 0, to = 0.26, by = 0.03))+
  expand_limits(x = c(0, 0.30), y=c(-3,176))+
  geom_nodelab(geom='label', aes(label=as.numeric(support)*100, subset=as.numeric(support) > 0.7),size = 1.2, label.padding = unit(0.1, "lines"),
               alpha=0.7, nudge_x = -0.005)+ 
  #geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, hjust=-0.02, size = 2.8,fontface = "bold")+ 
  theme(axis.text = element_text(vjust=0.5,size=10, face = "bold", angle = 90), 
        legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.direction="vertical", 
        legend.position = "left", legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.1,'cm'),
        legend.margin = margin(15, -4, 14, 1, "cm"),
        plot.margin = unit(c(-8,0,5,0), "lines"))+  
  geom_strip('2018/RV30/Human/Sarawak','2018/RV46/Human/Sarawak', offset = 0.180, barsize=2.5, color='#ff854d', 
             label="Sarawak \n samples", offset.text=-0.05,size = 2.5, fontface =2, fontface="bold")+
  guides(fill=guide_legend(title="Rabies Lineages"))+
  coord_cartesian(xlim = c(-0.02, 0.48), # This focuses the x-axis on the range of interest
                  clip = 'off')

p2.1 <- scaleClade(p2, 343, scale=0.2) %>% collapse(343, 'max', fill="darkorange") +
  geom_cladelab(node=343, label="Sarawak strains", align=TRUE, vjust=-0.5, geom='label',size=4, fill='darkorange') +
  geom_cladelab(node=336, label="Indonesia strains", align=TRUE, vjust=-0.5, geom='label',size=4, fill='lightblue')
  
p2.2 <- viewClade(p2, node = MRCA(p2, "RVD1159/dbs/MKM/2024", "1997/AB154220/Dog/Indonesia_WestJava")) %>% 
    scaleClade(343, scale=0.1) %>% collapse(343, 'max', fill="orange") +
    geom_cladelab(node=343, label="Sarawak strains", align=TRUE, hjust=0.2, size=7.5, vjust=-4, geom='label',face = "bold", fill='darkorange')+
  scale_color_manual(values=region_cols)+
  geom_tiplab(align=FALSE, linetype='dashed', linesize=.3, hjust=-0.02, 
              size = 2.8, fontface = "bold", show.legend = FALSE)+
  theme(legend.position = "none", plot.margin = unit(c(-5,-30,8,0), "lines"))

  #geom_strip('2017/99/Human/Sarawak','2018/RV28/Human/Sarawak', offset = 0.182, barsize=2.5, color='#ff854d', 
  #           label="Sarawak \n samples", offset.text=-0.056,size = 3.0, fontface =2) + 
  #geom_strip('2019/RV56/Human/Sarawak', '2018/RV19/Human/Sarawak', offset = 0.182, barsize=2.9, color='#ff854d')+ 

plotT1 <- ggarrange(p2.1, p2.2, labels = c("A", "B"), nrow = 1, ncol = 2, widths = c(0.4, 0.5))
ggsave(plotT1, file ='output/global_fN_swkrabv_2024.png',width = 19, height = 22, units = "cm",limitsize = FALSE, device = 'png')
}
### 3B4. MAIN ANALYSIS CODES: RABV PARTIAL NUCLEOPROTEIN GENE DISTREE  ###
{
p3 <- ggtree(tree, size=0.75, show.legend = TRUE) %<+% metadata_df + theme_tree2()+
  geom_tippoint(aes(color=Lineage), size=3.5)+
  scale_color_manual(values=region_cols)+
    scale_x_continuous(breaks=c(0.00,0.10,0.20,0.30), labels=c("0.00", "0.10", "0.20","0.30"))+
    expand_limits(x = c(0, 0.35), y=c(0,105))+
  geom_nodelab(geom='label', aes(label=as.numeric(support), subset=as.numeric(support) > 70),size = 2.5, label.padding = unit(0.1, "lines"),
               alpha=0.7, nudge_x = -0.005)+
  #geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, hjust=-0.02, size = 2.8,fontface = "bold")+ 
  scale_color_manual(values=c(2018/RV119/Human/Sarawak = "#e71e00"))+
  theme(axis.text = element_text(vjust=0.5,size=10, face = "bold", angle = 90), 
        legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.direction="vertical", 
        #legend.key = element_blank(),
        legend.position = "left", legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.1,'cm'),
        legend.margin = margin(2, -4, 14, 1, "cm"),
        #axis.ticks.length = unit(0.25, "cm")+
  plot.margin = unit(c(-10,-30,0,0), "lines"))+  
  guides(fill=guide_legend(title="Rabies Lineages"))+
  coord_cartesian(xlim = c(-0.02, 0.43), # This focuses the x-axis on the range of interest
                  clip = 'off')
p3.1 <- collapse(p3, 106, 'mixed')
p3.2 <- viewClade(p3, node = MRCA(p, "2019/RV56/Human/Sarawak", "1997/AB154220/Dog/Indonesia_WestJava"))+
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, hjust=-0.02, 
              size = 2.8, fontface = "bold", show.legend = FALSE)+
  geom_strip('2017/99/Human/Sarawak','2018/RV28/Human/Sarawak', offset = 0.182, barsize=2.5, color='#ff854d', 
             label="Sarawak \n samples", offset.text=-0.056,size = 3.0, fontface =2) + 
  geom_strip('2019/RV56/Human/Sarawak', '2018/RV19/Human/Sarawak', offset = 0.182, barsize=2.9, color='#ff854d')+ 
theme(legend.position = "none", plot.margin = unit(c(-5,-30,0,0), "lines"))

p3.1 + p3.2
ggsave('output/global_pN_swkrabv_2022.png',width = 21, height = 21, units = "cm",limitsize = FALSE, device = 'png')

panelA <- p3 %<+% metadata_df +
  geom_tippoint(
    aes(subset=(Category=='Asian')),fill='#F2A922',size=3.0,
                            color='#704c07',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Cosmopolitan')),fill='#267d56',size=3.0,
                color='#1a563b',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Arctic-related')),fill='#1bbbd0',size=3.0,
                color='#157aa3',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Africa')),fill='#f55a76',size=3.0,
                color='#83081e',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  geom_tippoint(aes(subset=(Category=='Indian_subcontinent')),fill='#c3699d',size=3.0,
                color='#bc575f',stroke=0.1,shape=21,x_offset = 0.01, y_offset = 0.01, sigma =3)+
  coord_cartesian(xlim = c(0, 0.26), # This focuses the x-axis on the range of interest
                  clip = 'off')+
scale_fill_manual(values=c("#F2A922","#267d56","#FF0000",
                           "#1bbbd0", "#f55a76","#c3699d"),
                  guide=guide_legend(keywidth = 0.5, 
                                     keyheight = 0.5, order=1,
                                     override.aes=list(starshape=15)),
                  na.translate=FALSE)
panelA
ggsave('output/global_swk_rabv_2022.png',width = 18, height = 23, units = "cm",limitsize = FALSE, device = 'png')
}
### 3B5. MAIN ANALYSIS CODES: RABV GENOME TIMETREE  ###
#### 3B5a. MAIN ANALYSIS CODES: RABV GENOME TIMETREE part1: Define file paths & Source  ####
{
  library(ggplot2)
  library("readxl")
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
  library(cowplot)
  library(ggthemes)
  library(viridis)
  library(ggalluvial)
  library("lubridate")
  library("scales")
  library(ggtree)
  library(ggbreak)
  library(svglite)
  library(tidytree)
  library(ape)
  library(ggfx)
  library(treeio)
  library(choroplethr)
  library(choroplethrMaps)
  swk_rabv_mcc <-read.beast("07_phylogeo/rabv_nhm_cRRW_HKY-G_RCLN_CoalesBayesSG_500m5k.trees.anno")
  swk_rabv_time <-read.nexus("2024-10-30_indo_ml_treetime/timetree_bootstrap.nexus")
  seq_meta <- read_excel("01_src/rabv_glob_metadata.csv", sheet = 1)
  tree_meta <- seq_meta[!(is.na(seq_meta$col_date) | seq_meta$col_date=="" | seq_meta$col_date=="na"), ]
  tree_meta <- tree_meta[(tree_meta$country=="Sarawak" | tree_meta$country=="Indonesia"), ]
  tree_meta <- tree_meta[,-(1),drop=FALSE] 
  custom1<-c("KUCHING"='#F2A922',"KCH-PADAWAN"='#D2B55B',"KCH-BAU"='#BEAE92', "SIBU"='#83083e',"SERIAN"='#3BBF82',"SIBURAN"='#ff8c00', "SAMARAHAN"='#1c7da3', 
             "SARIKEI"='#a70000',"SIMUNJAN"='#085c82',"MUKAH"='#f55270',"KCH-LUNDU"='#f45548',"LIMBANG"='#006666', "ENGKILILI"='#00b300',
             "JULAU"='#4d2f75',"MIRI"='#bb8fce',"SRI AMAN"='dodgerblue1',"SIMUNJAN"='hotpink2',"BINTULU"='mediumorchid3', "MERADONG"='palegreen3', "(REF) Kalimantan" ='#b8e3e9')
  dtclade <- data.frame(Swkrv01=c("RVD687/dbs/KCH/2020", "RVD402/dbs/KCH/2019"), Swkrv02=c("SVDL2096_22/dbs/LUN/2022", "RVD831/dbs/KCH/2022"), 
                        Swkrv03=c("RVD503/dbs/KCH/2019", "RVD172/dbs/SER/2018"), Swkrv04=c("RVD401/dbs/SER/2019", "RVD348/dbs/KCH/2019"), 
                        Swkrv05=c("SVDL2108_22/dbs/BTU/2022", "SVDL0103_22/dbs/SAR/2022"), Swkrv06=c("RVD686/dbs/KCH/2020", "RVC199/cbs/KCH/2018"), 
                        group=c('Swkrv01', 'Swkrv02', 'Swkrv03', 'Swkrv04', 'Swkrv05', 'Swkrv06'))
  
  dtclade <- data.frame(ta1=c("RVD687/dbs/KCH/2020", "RVD402/dbs/KCH/2019"), ta2=c("SVDL2096_22/dbs/LUN/2022", "RVD831/dbs/KCH/2022"), 
                        ta3=c("RVD503/dbs/KCH/2019", "RVD172/dbs/SER/2018"), ta4=c("RVD401/dbs/SER/2019", "RVD348/dbs/KCH/2019"),
                        group=c('Swkrv01', 'Swkrv02', 'Swkrv03', 'Swkrv03'))
  Swkrv07=c("RVD1038/dsal/KCH/2024", "RVD555")
}
#### 3B5b. MAIN ANALYSIS CODES: RABV GENOME TIMETREE part2: Plotting & Saving  ####
{
  ptt2 <- ggtree(swk_rabv_time, mrsd="2024-07-11", as.Date=TRUE, color='#6a777b', size=0.3) +
    geom_tiplab(size=0.8)+
    scale_x_break(c(as.Date('1920-01-01'), as.Date('2001-12-31')), space = 0.5)+
    #scale_x_break(c(as.Date('1992-01-01'), as.Date('2002-12-31')), space = 0.5)+
    scale_x_break(c(as.Date('2005-01-01'), as.Date('2015-12-31')), space = 0.5)
  
  ptt2 <- ggtree(swk_rabv_time, mrsd="2024-07-11", as.Date=TRUE, linetype ='twodash',color='#6a777b', size=0.65)%<+% tree_meta + theme_tree2()+
    scale_x_date(date_labels = "%Y",date_breaks = "1 year",limits = as.Date(c("1918-01-01", "2024-07-31"))) +
    scale_x_break(c(as.Date('1920-01-01'), as.Date('2001-12-31')), space = 0.5)+
    scale_x_break(c(as.Date('2005-01-01'), as.Date('2015-12-31')), space = 0.5)+
    expand_limits(y = c(-3, 226))+
    #geom_tiplab(=0.8)+
    geom_tippoint(aes(fill=location), size=0.6,stroke=0.4,color='#704c07',shape=21)+
    geom_nodepoint(geom='label', size = 0.9, shape=23, 
                   fill="#f45948", aes(label=as.numeric(label), subset=as.numeric(label) > 0.70), label.padding = unit(0.1, "lines"), alpha=0.7)+
    geom_nodelab(geom='label', aes(label=as.numeric(label) * 100, subset=as.numeric(label) > 0.70), size = 1.2, label.padding = unit(0.1, "lines"),
                 alpha=0.7, nudge_x = -350)+
    scale_fill_manual("location", values= custom1)+
    theme(axis.text = element_text(vjust=0.5,size=9, face = "bold", angle = 90), 
          legend.title = element_text(hjust = 0.5, vjust = -0.5, size=10, face = "bold"),
          legend.text = element_text(size = 8),
          legend.direction="vertical", 
          legend.position = "left", legend.key.size = unit(0.3, 'cm'),
          legend.spacing = unit(0.1,'cm'),
          legend.margin = margin(15, -5, 2, 0, "cm"),
          plot.margin = unit(c(5,0,0,0), "lines"))+ 
    guides(fill=guide_legend(title="Location", ncol=2))
  ptt2.1 <- ggtree(swk_rabv_time, mrsd="2024-09-11", as.Date=FALSE)%<+% tree_meta +
    geom_strip('RVD687/dbs/KCH/2020','RVD402/dbs/KCH/2019', offset = -5, barsize=1.5, color='#c1464c', label="swkC1", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('SVDL2096_22/dbs/LUN/2022', 'RVD831/dbs/KCH/2022', offset = -5, barsize=1.5, color='#333d64', label="swkC2", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('RVD503/dbs/KCH/2019','RVD172/dbs/SER/2018', offset = -5, barsize=1.5, color='#765ed3', label="swkC3", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('RVD401/dbs/SER/2019','RVD348/dbs/KCH/2019', offset = -5, barsize=1.5, color='#fa95ad', label="swkC4", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('SVDL2108_22/dbs/BTU/2022','SVDL0103_22/dbs/SAR/2022', offset = -5, barsize=1.5, color='#ff854d', label="swkC5", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('RVD686/dbs/KCH/2020','RVC199/cbs/KCH/2018', offset = -5, barsize=1.5, color='#006699', label="swkC6", offset.text=0.5,size = 2.5, fontface="bold")+
    geom_strip('RVD1038/dsal/KCH/2024','RVD555/dbs/PDWN/2019', offset = -5, barsize=1.5, color='#00cc99', label="swkC7", offset.text=0.5,size = 2.5, fontface="bold")+
    theme(plot.margin = unit(c(0,5,0,-35), "lines"))
  ggsave(ptt2, filename = 'output/swk_cds_indo_timetreeML_2024_tip.png',width = 26, height = 23, units = "cm",limitsize = FALSE, device = 'png')
  ggsave(ptt2.1, filename = 'output/swk_cds_indo_timetreeclade_2024.png',width = 26, height = 23, units = "cm",limitsize = FALSE, device = 'png')
  
  rabv_Re_gridded <- read.table("2024-10-30_indo_ml_treetime/skyline.tsv", sep="\t", header= TRUE)
  rabv_Re_gridded = subset(rabv_Re_gridded, rabv_Re_gridded$date  > 1975)  #remove Limbang entry to simplify analysis
  p_locseqre <- ggplot()+
    geom_ribbon(data = rabv_Re_gridded, mapping = aes(ymin=lower, ymax=upper, x =date), fill="#ff981f", alpha=0.5)+
    geom_line(data = rabv_Re_gridded, aes(x=date, y=N_e), size=0.6, col="#814500", alpha=0.8)+
    geom_hline(yintercept = 1, linetype = 'dotdash', color = '#494949') +
    theme_classic()+
    ylab("Effective reproductive number, Re")+xlab("Year")+coord_cartesian(clip = "off")+
    theme(legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
          legend.text = element_text(size = 10),
          axis.title.x = element_text(color="black", size=12, face="bold"),
          axis.text.x = element_text(color="black", size=8,vjust=0.5,angle = 90),
          axis.title.y = element_text(color="black", size=12, face="bold"),
          axis.text.y = element_text(color="black", size=8))
  ggsave(p_locseqre, filename = "output/swk_cds_indo_timetreeRE_2024.png", width = 12, height = 10, units = "cm",device = "png")
  
  
  php2 <- ggtree(swk_rabv_mcc, layout="fan",linetype ='twodash',color='#6a777b', size=0.65) %<+% hb.results$partition.df + theme_tree2()+
    geom_tippoint(aes(color = factor(`level 2`))) +
    #expand_limits(y = c(-3, 150))+
    geom_nodepoint(geom='label', aes(label=round(as.numeric(posterior), 2),subset=!isTip & as.numeric(posterior)> 0.90),size = 1.5, shape=23, 
                   fill="#f45948", label.padding = unit(0.1, "lines"), alpha=0.7, nudge_x = -0.005)+
    scale_fill_manual("location", values= custom1)+
    theme(axis.text = element_text(vjust=0.5,size=10, face = "bold", angle = 90), 
          legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
          legend.text = element_text(size = 8),
          legend.direction="vertical", 
          legend.position = "left", legend.key.size = unit(0.5, 'cm'),
          legend.spacing = unit(0.1,'cm'),
          legend.margin = margin(-5, -2, 2, 0, "cm"),
          plot.margin = unit(c(5,0,0,0), "lines"))+ 
    guides(fill=guide_legend(title=""))
  php2
  ggsave('output/swk_fGhp_swkrabv_2024.png',width = 21, height = 30, units = "cm",limitsize = FALSE, device = 'png')
  ############################### UNUSED_CUSTOMIZATION ############################### 
  panelB <- p2 %<+% swk_rvmeta + 
    geom_tippoint(aes(subset=(location=='KUCHING')), fill='#F2A922',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SIBU')), fill='#1c7da3',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SERIAN')), fill='#2b8c60',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SIBURAN')), fill='#3BBF82',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SAMARAHAN')), fill='#36561a',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='BAU')), fill='#681c3f',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SIMUNJAN')), fill='#085c82',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='MUKAH')), fill='#f55270',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='LUNDU')), fill='#83081e',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='ENGKILILI')), fill='#A62D65',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='JULAU')), fill='#4d2f75',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='MIRI')), fill='#68409f',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='BINTULU')), fill='#9a78c9',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SRI AMAN')), fill='dodgerblue1',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SIMUNJAN')), fill='hotpink2',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='SELANGAU')), fill='mediumorchid3',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    geom_tippoint(aes(subset=(location=='MERADONG')), fill='palegreen3',size=3.5,stroke=0.4,color='#704c07', shape=21)+
    
    panelB
  ggsave('output/local_swk_rabv_bayesSG_mcc.svg',width = 36, height = 25, units = "cm",limitsize = FALSE, device = 'svg')
  
}

################################### 4A. RABV EPIDEMIOLOGY: EPID CASES & TRENDS -TO BE MODIFIED IF MORE DATA OBTAINED ##############################
### 4A1. MAIN SOURCE CODES: DEFINE PATHS TO CSV, XLSX FILES, SOURCES  ###
{
library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggthemes)
library(wesanderson)
library(RColorBrewer)
library(ggalt)
library(slider)
library(lubridate)
library(HDInterval)
library(smoother)
library(padr)
library(readr)
library(reprex)
ymax = 110000
swkstat <- read_excel("source/swk_rabv_stat.xlsx")
swkstat <- rename(swkstat, c("cummu_Date"="date")) %>%
  mutate(cummu_Date = as_date(cummu_Date))%>%
  pad(start_val = ymd("2017-07-12"), end_val = ymd("2020-08-01"))
p_swk_rabv <- ggplot(swkstat) +
  geom_bar(aes(x = cummu_Date, y= case_bite),colour="#3b3b3b",width = 1*1*2,fill='#a95a00',position = position_nudge(x = 1*1),size=0.05)
p_swk_rabv
geom_bar(mapping = aes(x = cummu_Date, y = case_bite),size=0.8,color = "#e03974", alpha =0.8) +
  geom_bar(mapping = aes(y = vac_dog, x = cummu_Date),size=0.8,color = "#ef8000", alpha =0.8) +
  geom_bar(mapping = aes(y = vac_cat, x = cummu_Date),size=0.8,color = "#00adb4", alpha =0.8) +
  theme_classic()
p_swk_rabv
ggsave(p_swk_rabv, filename = "Epid_re_n_rt_2022.png", width = 38, height = 20, units = "cm",device = "png")
}

############################### 4B. GENERATE EPID TREND MAPS ############################### 
### 4B1. MAIN SOURCE CODES: DEFINE PATHS TO FILES, VALUES, DATA SORTING  ###
{
  options(scipen = 999)
  options(digits=6)
  swk_sf <- st_read("swk_maps/swk_divmap_2022.shp")
  metamap <- read_csv("01_src/swk_mapmeta.csv")
  # Map base
  
  swk_sf <- swk_sf %>% mutate(centroid = map(geometry, st_centroid), 
                              coords = map(centroid, st_coordinates), 
                              coords_x = map_dbl(coords, 1), coords_y = map_dbl(coords, 2))
  #swk_sf <- swk_sf[-c(6:15)]
  target <- c("Kuching","Samarahan","Serian","Sri Aman","Betong","Sarikei","Sibu","Mukah", "Bintulu","Kapit","Miri","Limbang")
  swk_sf <- swk_sf[match(target, swk_sf$ADM2_EN),]
  swk_sf2 <- do.call("rbind", replicate(25, swk_sf, simplify = FALSE))
  swk_sf2 <- tibble::rowid_to_column(swk_sf2, "ID")
  metamap <- metamap[-c(2)]
  swk_tab <- left_join(swk_sf2, metamap, by = "ID") %>%
    #select(OBJECTID,ADM2_EN,coords,case_prop,coords_x,coords_y,case_per100k,col_dateend,col_year,type,value,value_c,population,Region,geometry) %>%
    mutate(col_dateend = as.Date(col_dateend, format = "%m/%d/%Y"))
  swk_tabsamp<-subset(swk_tab,type=='total')
  swk_tabpos<-subset(swk_tab,type=='pos')
  swk_tabngs<-subset(swk_tab,type=='seq_ngs')
  
  samp_ranges <- c(
    "<10", "11-50","51-100","101-200","201-300", "300-500", ">500")
  swk_samprgn1 <- 
    swk_tabsamp %>%
    mutate(samp_range = factor(case_when(
      `value_c` < 10 ~ samp_ranges[1],
      between(`value_c`, 11, 50) ~ samp_ranges[2],
      between(`value_c`, 51, 100) ~ samp_ranges[3],
      between(`value_c`, 101, 200) ~ samp_ranges[4],
      between(`value_c`, 201, 300) ~ samp_ranges[5],
      between(`value_c`, 301, 500) ~ samp_ranges[6],
      `value_c` >= 501 ~ samp_ranges[7]), levels = samp_ranges))
  colors_map <- diverging_hcl(7, palette = "Blue-Red 3", rev = F)
}
### 4B2. MAIN PLOTTING CODES: PLOTTING & SAVING  ###
{
  plot1 <- swk_samprgn1 %>% 
    filter(col_dateend == as.Date('2018-12-31') | col_dateend == as.Date('2019-12-31') | col_dateend == as.Date('2020-12-31') | col_dateend == as.Date('2021-12-31') | col_dateend == as.Date('2022-12-31') | col_dateend == as.Date('2023-12-31') | col_dateend == as.Date('2024-07-11')) %>%
    ggplot() +  geom_sf(color = "black", size = 0.1) +
    aes(fill = samp_range) +
    #coord_sf(datum = NA)+
    scale_fill_manual(values = colors_map) +
    theme_void()+
    #geom_point(data = swk_tab1, aes(coords_x, coords_y, size = population, label = "Population size"), color = "#2c9264")+
    #geom_label_repel(data = swk_tab1, aes(coords_x, coords_y, label = seq_sum), size = 4, color = "#9e3566",fontface = "bold",nudge_y=0.9)+
    labs(fill = "Cummulative samples",title="Number and divisions where animal samples were collected")+
    theme(legend.position="right", strip.text = element_text(size = 10, face="bold"),legend.key.size = unit(0.3, 'cm'),
          legend.title= element_text(size = 9,face = "bold"),
          plot.title = element_text(hjust = 0, size = 12, face = "bold",margin=margin(0,0,8,0)),plot.margin = unit(c(-1, 0, -2, 0), "cm")) + facet_wrap("col_dateend", ncol = 4)
  
  pos_ranges <- c("<10", "11-30", "31-50", "51-100", "101-200", "201-300", ">300")
  swk_samprgn2 <- 
    swk_tabpos %>%
    mutate(pos_range = factor(case_when(
      `value_c` < 10 ~ pos_ranges[1],
      between(`value_c`, 11, 30) ~ pos_ranges[2],
      between(`value_c`, 31, 50) ~ pos_ranges[3],
      between(`value_c`, 51, 100) ~ pos_ranges[4],
      between(`value_c`, 101, 200) ~ pos_ranges[5],
      between(`value_c`, 201, 300) ~ pos_ranges[6],
      `value_c` >= 301 ~ pos_ranges[7]), levels = pos_ranges))
  colors_map <- diverging_hcl(7, palette = "Green-Brown", rev = F)
  plot2 <- swk_samprgn2 %>% 
    filter(col_dateend == as.Date('2018-12-31') | col_dateend == as.Date('2019-12-31') | col_dateend == as.Date('2020-12-31') | col_dateend == as.Date('2021-12-31') | col_dateend == as.Date('2022-12-31') | col_dateend == as.Date('2023-12-31') | col_dateend == as.Date('2024-07-11')) %>%
    ggplot() +  geom_sf(color = "black", size = 0.1) +
    aes(fill = pos_range) +
    #coord_sf(datum = NA)+
    scale_fill_manual(values = colors_map) +
    theme_void()+
    #geom_point(data = swk_tab1, aes(coords_x, coords_y, size = population, label = "Population size"), color = "#2c9264")+
    #geom_label_repel(data = swk_tab1, aes(coords_x, coords_y, label = seq_sum), size = 4, color = "#9e3566",fontface = "bold",nudge_y=0.9)+
    labs(fill = "Cummulative samples",title="Number of PCR positive samples by different divisions")+
    theme(legend.position="right", strip.text = element_text(size = 10, face="bold"),legend.key.size = unit(0.3, 'cm'),
          legend.title= element_text(size = 9,face = "bold"),
          plot.title = element_text(hjust = 0, size = 12, face = "bold",margin=margin(0,0,8,0)),plot.margin = unit(c(-1, 0, -2, 0), "cm")) + facet_wrap("col_dateend", ncol = 4)
  
  seq_ranges <- c("<5", "6-10","11-30","31-50","51-80","81-100", ">100")
  swk_samprgn3 <- 
    swk_tabngs %>%
    mutate(seq_range = factor(case_when(
      `value_c` < 5 ~ seq_ranges[1],
      between(`value_c`, 6, 10) ~ seq_ranges[2],
      between(`value_c`, 11, 30) ~ seq_ranges[3],
      between(`value_c`, 31, 50) ~ seq_ranges[4],
      between(`value_c`, 51, 80) ~ seq_ranges[5],
      between(`value_c`, 81, 100) ~ seq_ranges[6],
      `value_c` >= 101 ~ seq_ranges[7]), levels = seq_ranges))
  colors_map <- diverging_hcl(7, palette = "Green-Brown", rev = F)
  plot3 <- swk_samprgn3 %>% 
    filter(col_dateend == as.Date('2018-12-31') | col_dateend == as.Date('2019-12-31') | col_dateend == as.Date('2020-12-31') | col_dateend == as.Date('2021-12-31') | col_dateend == as.Date('2022-12-31') | col_dateend == as.Date('2023-12-31') | col_dateend == as.Date('2024-07-11')) %>%
    ggplot() +  geom_sf(color = "black", size = 0.1) +
    aes(fill = seq_range) +
    #coord_sf(datum = NA)+
    scale_fill_manual(values = colors_map) +
    theme_void()+
    #geom_point(data = swk_tab1, aes(coords_x, coords_y, size = population, label = "Population size"), color = "#2c9264")+
    #geom_label_repel(data = swk_tab1, aes(coords_x, coords_y, label = seq_sum), size = 4, color = "#9e3566",fontface = "bold",nudge_y=0.9)+
    labs(fill = "Cummulative sequences",title="Number of genome sequences by different divisions")+
    theme(legend.position="right", strip.text = element_text(size = 10, face="bold"),legend.key.size = unit(0.3, 'cm'),
          legend.title= element_text(size = 9,face = "bold"),
          plot.title = element_text(hjust = 0, size = 12, face = "bold",margin=margin(2,0,8,0)),plot.margin = unit(c(0, 0, -1, 0), "cm")) + facet_wrap("col_dateend", ncol = 4)
  arrange <- ggarrange(plot1, plot2, plot3,font.label=list(color="black",size=12),
                       ncol = 1, nrow = 3)
  ggsave("rabv_casemapplot_2024.png", arrange, width = 13, height = 14)
}
############################### 4C. GENERATE EPID TREND GRAPH + TRACING MAPS ############################### 
{
swk_samprgn <- swk_tab[swk_tab$type =="pos" & !(swk_tab$value=="0"), ]
swk_samprgnp <- as.data.frame(swk_samprgn) %>% group_by(col_year) %>%
  summarize(value_p = sum(value))
swk_samprgnt <- swk_tab[swk_tab$type =="total" & !(swk_tab$value=="0"), ]
swk_samprgnt <- as.data.frame(swk_samprgnt) %>% group_by(col_year) %>%
  summarize(value_t = sum(value))
swk_samprgn4n <- left_join(swk_samprgnt, swk_samprgnp, by="col_year") %>%
  mutate(col_year = as.numeric(col_year)) %>%
  group_by(col_year, value_t, value_p) %>%
  summarise(test_pct = value_p / value_t)
plot4 <- ggplot() +
  geom_bar(data=swk_samprgn, aes(fill=ADM2_EN, x = col_year, y = value), position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  scale_y_continuous(breaks = seq(100, 800, by = 50))+
  scale_x_continuous(breaks = seq(2018, 2024, by = 1))+
  geom_line(data=swk_samprgn4n, aes(x=col_year, y = test_pct *100),size=1.5, color='#FF6327') +
  scale_y_continuous(sec.axis=sec_axis( 
    ~.*1, name="Percentage of cases tested positive")) +
  #ggtitle("C. Total number of positive animals by RT-PCR") +
  labs(y = "Case number", x = "Year of collection") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        text=element_text(family="Arial"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11))+ guides(fill=guide_legend(title="Division"))
}
############################### 4D.  GENERATE CASE MAP, COVERAGE VS CT VALS (3 PANELS) ############################### 
{
writeLines(gsub("\t", ",", readLines("01_src/ngs_all_readstats.txt")), "01_src/ngs_all_readstats.csv")
epistat <- read_csv("01_src/ngs_all_readstats.csv")
seq_meta<- read_csv("01_src/rabv_glob_metadata.csv")
seq_meta <- seq_meta[!(is.na(seq_meta$col_date) | seq_meta$col_date=="" | seq_meta$col_date=="na"), ]
seq_meta <- seq_meta[seq_meta$country=="Sarawak", ]
swk_seqstats <- left_join(epistat, seq_meta, by = "sample_id") %>%
  mutate(col_date = as.Date(col_date, format = "%m/%d/%Y")) 
swk_seqstats <- swk_seqstats[c(2:4,6:10,16:18,30:34)]
fill = c("#4652B4", "#FF76B4")
plot5 <- ggplot(swk_seqstats, aes(x = ct_val, y = coverage.x, size = meandepth, fill = seq_method)) +
  geom_point(shape = 21) +
  #ggtitle("D. Genome coverage and mean read depth reflected by sample ct values") +
  labs(x = "Ct values", y = "Genome coverage") +
  scale_y_continuous(trans = log2_trans(), limit = c(93,100), breaks = c(93,98,99,99.5,100))+
  scale_x_continuous(breaks = seq(0, 40, by = 5))+
  scale_fill_manual(values = fill) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        text=element_text(family="Arial"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11))+ guides(fill=guide_legend(title="Primer scheme"))

custom1<-c("KUCHING"='#F2A922',"SIBU"='#83081e',"SERIAN"='#3BBF82',"SIBURAN"='#ff8c00', "SAMARAHAN"='#36561a', 
           "BAU"='#681c3f',"SIMUNJAN"='#085c82',"MUKAH"='#f55270',"LUNDU"='#1c7da3',"ENGKILILI"='#A62D65',
           "JULAU"='#4d2f75',"MIRI"='#68409f',"BINTULU"='#BEAE92',"SRI AMAN"='dodgerblue1',"SIMUNJAN"='hotpink2',"SELANGAU"='mediumorchid3', "MERADONG"='palegreen3')
library(tidyverse)
library(sf)
library(ggmap)
swk_samprg<- read.table("01_src/rabv_swk_samptstat.txt", sep="\t", header=TRUE)
swk_samprgmap <- swk_samprg[c(2:5,7,10:12,15:20,24:26)]
swk_samprgmap <- swk_samprgmap[!(is.na(swk_samprgmap$lat) & swk_samprgmap$lat=="" & swk_samprgmap$uniq==""), ]
swk_samprgmap <- swk_samprgmap[ swk_samprgmap$result=="pos" | swk_samprgmap$result=="gel_pos" | swk_samprgmap$result=="w_pos" , ]
swk_samprgmapseq <- swk_samprgmap[swk_samprgmap$seq_done=="y", ]
register_stadiamaps(key = "16011e6c-4b8e-4ace-99f4-8b4657994543")

map <- get_stadiamap( bbox = c(left = 109.8, bottom = 0.78, right = 113.5, top = 3.4), zoom = 8, maptype = "alidade_smooth")
plot6 <-  ggmap(map) + 
  geom_point(data = swk_samprgmap, aes(x = as.numeric(lon), y = as.numeric(lat), fill = div), size = 2.2, shape = 21, alpha=0.6) +
  theme_void() + 
  scale_fill_manual("div", values= custom1)+
  #ggtitle("A. Positive samples by location") +
  labs(y = "Case number", x = "Year of collection") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+ guides(fill=guide_legend(title="Division"))

plot7 <-  ggmap(map) + 
  geom_point(data = swk_samprgmapseq, aes(x = as.numeric(lon), y = as.numeric(lat), fill = div), size = 2.2, shape = 21, alpha=0.6) +
  theme_void() + 
  scale_fill_manual("div", values= custom1)+
  #ggtitle("B. Virus genomes sequenced from positive samples by location") +
  labs(y = "Case number", x = "Year of collection") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+ guides(fill=guide_legend(title="Division"))
mapzm <- get_stadiamap( bbox = c(left = 110.12, bottom = 1.36, right = 110.45, top = 1.596), zoom = 12, maptype = "alidade_smooth")
plot8 <-  ggmap(mapzm) + 
  geom_point(data = swk_samprgmapseq, aes(x = as.numeric(lon), y = as.numeric(lat), fill = div), size = 2.2, shape = 21, alpha=0.5) +
  theme_void() + 
  scale_fill_manual("div", values= custom1)+
  #ggtitle("C. Closer look at sample collected from  Kuching") +
  labs(y = "Case number", x = "Year of collection") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )+ guides(fill=guide_legend(title="Division"))
plot6_7_8 <- ggarrange(
  plot6, plot7, plot8, labels = c("A", "B", "C"), nrow = 1, ncol = 3,
  common.legend = TRUE, legend = "bottom"
)
plot4_5 <- ggarrange(
  plot4, plot5, labels = c("D", "E"), nrow = 1, ncol = 2,
  common.legend = FALSE, legend = "bottom"
)
layout <- c(
  area(1, 1, 4, 10),
  area(5, 1, 8, 10)
)
swk_epid_trend <- plot6_7_8 /plot4_5 + plot_layout(design = layout)
ggsave(swk_epid_trend, filename = 'rabv_tracemap_2024', width = 31, height = 19, units = "cm",limitsize = FALSE,device = "png")
}
############################### 4E.  SAMPLE COUNT, GENOME SEQ DEPTH (2 PANELS) ############################### 
{
  type_cols1 <- c("#33658A","#F26157", "#A54657","#582630","#F7EE7F","#4DAA57","#F9ECCC","#679289","#33658A",
                  "#F6AE2D","#86BBD8","#9b5fce", "#e6799a")
  rabv_demog <- read.table(file = "01_src/rabv_swk_loc_cds_nhm_demog.txt", sep = ',', header = TRUE)
  #group by each input feature
  yr_spec <- rabv_demog %>%
    group_by(year, species) %>%
    summarise(n = n()) %>%
    mutate(pct= prop.table(n) * 100)
  
  yr_div <- rabv_demog %>%
    group_by(year, division) %>%
    summarise(n = n()) %>%
    mutate(pct= prop.table(n) * 100)
  pyr_spec <- ggplot(yr_spec, aes(x = year, y = n , fill= species)) +
    geom_bar(position = position_stack(), stat = "identity", width = .7) + 
    geom_label(data=subset(yr_spec, pct>=30), aes(label=paste0(sprintf("%1.1f", pct),"%")),
               position=position_stack(vjust=0.75),show.legend = FALSE, size=3)+
    theme(axis.text.x = element_text(angle = 0, size = 12, colour = "black", face= "bold"), 
          axis.title = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 9, face = "bold", colour = "black"), 
          axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(x = "Collection year", y = "Sample count", fill = "Species") + 
    scale_fill_manual(values = type_cols1)
  pyr_div <- ggplot(yr_div, aes(x = year, y = n , fill= division)) +
    geom_bar(position = position_stack(), stat = "identity", width = .7) +
    geom_label(data=subset(yr_div, pct>=20), aes(label=paste0(sprintf("%1.1f", pct),"%")),
               position=position_stack(vjust=0.75),show.legend = FALSE, size=3)+
    theme(axis.text.x = element_text(angle = 0, size = 12, colour = "black", face= "bold"), 
          axis.title = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 9, face = "bold", colour = "black"), 
          axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(x = "Collection year", y = "Sample count", fill = "Location") +
    scale_fill_brewer(palette = "Spectral")
  
  plot9_10 <- ggarrange(pyr_div, pyr_spec, labels = c("A", "B"), nrow = 2, ncol = 1,
                        common.legend = FALSE, legend = "bottom")
  plot9_10_C <- ggarrange(plot9_10, seqdepth, labels = c("", "C"), nrow = 1, ncol = 2,
                          common.legend = FALSE, legend = "bottom")
  ggsave(plot9_10_C, filename = 'output/rabv_combostats_2024.png', width = 31, height = 19, units = "cm",limitsize = FALSE,device = "png")
}
############################### 2C. Selection PRESSURE ANALYSIS - HYPHY ###############################
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(tibble)
excel_path = "01_src/rabv_glob_metadata.xlsx"
sheet_name = NULL
make_virus_selection_plots <- function(
    excel_path,
    sheet_name = NULL,
    sheet_number = 4,
    highlight_all_pos = TRUE, # TRUE: highlight all pos from each method
    # FALSE: only highlight sites pos in all three methods
    show_nonpos_points = TRUE, # Control whether to show non-positive sites as dark gray
    gene_bar_width = 0.25,
    gene_bar_nudge_x = 0,
    gene_label_size = 3.4,
    point_shape = 21,
    point_size = 2,
    highlight_shape = 21,
    highlight_size = 2.7,
    output_dir = "plots"
) {
  
  # Load required libraries
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read data from Excel
  if (!is.null(sheet_name)) {
    df <- read_excel(excel_path, sheet = sheet_name)
  } else {
    df <- read_excel(excel_path, sheet = sheet_number)
  }
  
  # Ensure required columns exist
  required_cols <- c("virus", "gene", "site1_comb", "site2", "P[dN/dS>1]", 
                     "slac_sel", "P[A<B]", "fubar_sel", "fel_aphabeta_dif", "fel_sel")
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for dN and dS columns for dN/dS calculation
  has_dNdS <- all(c("dN", "dS") %in% colnames(df))
  
  # Rename columns for easier handling
  df <- df %>%
    rename(
      site = site1_comb,
      gene_site = site2,
      slac_pval = `P[dN/dS>1]`,
      fubar_pval = `P[A<B]`,
      fel_dif = fel_aphabeta_dif
    )
  
  # Calculate dN/dS ratio for each site if columns exist
  if (has_dNdS) {
    df <- df %>%
      mutate(dN_dS = ifelse(dS > 0, dN / dS, NA_real_))
  }
  
  # Standardize selection columns: convert all non-"pos" values to "nonpos"
  df <- df %>%
    mutate(
      slac_sel = ifelse(slac_sel == "pos", "pos", "nonpos"),
      fubar_sel = ifelse(fubar_sel == "pos", "pos", "nonpos"),
      fel_sel = ifelse(fel_sel == "pos", "pos", "nonpos")
    )
  
  # Filter out 'cds' gene if it exists
  if ("cds" %in% unique(df$gene)) {
    df <- df %>% filter(gene != "cds")
  }
  
  # Get individual genes and concatenate them
  # First, get the gene order based on their starting position in the CDS
  gene_order <- df %>%
    group_by(virus, gene) %>%
    summarise(
      gene_start_cds = min(site, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(virus) %>%
    arrange(gene_start_cds) %>%
    mutate(
      gene_order = seq_len(n())
    ) %>%
    ungroup()
  
  # Join the gene order back to the main data and concatenate genes
  processed_df <- df %>%
    left_join(gene_order, by = c("virus", "gene")) %>%
    group_by(virus) %>%
    arrange(gene_order, gene_site) %>%
    mutate(
      plot_x = row_number(),  # Simple sequential numbering for concatenated genes
      feature = gene
    ) %>%
    ungroup()
  
  # Calculate gene positions for annotation
  gene_positions <- processed_df %>%
    group_by(virus, gene) %>%
    summarise(
      concat_start = min(plot_x, na.rm = TRUE),
      concat_end = max(plot_x, na.rm = TRUE),
      concat_mid = (min(plot_x) + max(plot_x)) / 2,
      cds_start = min(site, na.rm = TRUE),
      cds_end = max(site, na.rm = TRUE),
      .groups = "drop"
    )
  # Add a column for sites that are positive in all three methods
  processed_df <- processed_df %>%
    mutate(
      all_three_pos = ifelse(
        slac_sel == "pos" & fubar_sel == "pos" & fel_sel == "pos",
        "Positive",
        NA
      )
    )
  
  # Generate positive sites table
  positive_sites_table <- processed_df %>%
    filter(all_three_pos == "Positive") %>%
    arrange(virus, gene, site) %>%
    select(virus, gene, site, gene_site) %>%
    rename(
      Serotype = virus,
      Gene = gene,
      `CDS Position` = site,
      `Gene Position` = gene_site
    ) %>%
    group_by(Serotype, Gene) %>%
    summarise(
      `CDS Positions` = paste(sort(unique(`CDS Position`)), collapse = ", "),
      `Gene Positions` = paste(sort(unique(`Gene Position`)), collapse = ", "),
      `Number of positive sites` = n(),
      .groups = "drop"
    ) %>%
    arrange(Serotype, Gene)
  
  # Calculate comprehensive gene dN/dS statistics if dN/dS data exists
  gene_dnds_statistics <- NULL
  if (has_dNdS) {
    gene_dnds_statistics <- processed_df %>%
      group_by(virus, gene) %>%
      summarise(
        n_sites = n(),
        mean_dN = mean(dN, na.rm = TRUE),
        mean_dS = mean(dS, na.rm = TRUE),
        mean_dNdS = ifelse(mean_dS > 0, mean_dN / mean_dS, NA_real_),
        median_dNdS = median(dN_dS, na.rm = TRUE),
        min_dNdS = min(dN_dS, na.rm = TRUE),
        max_dNdS = max(dN_dS, na.rm = TRUE),
        sd_dNdS = sd(dN_dS, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        across(c(mean_dN, mean_dS, mean_dNdS, median_dNdS, min_dNdS, max_dNdS, sd_dNdS),
               ~ round(., 3))
      ) %>%
      rename(
        Serotype = virus,
        Gene = gene
      ) %>%
      arrange(Serotype, Gene)
  }
  
  # Generate plots for each virus
  virus_list <- unique(processed_df$virus)
  all_plots <- list()
  
  for (virus_name in virus_list) {
    virus_data <- processed_df %>%
      filter(virus == virus_name)
    
    # Get gene positions for this virus
    gene_pos <- gene_positions %>%
      filter(virus == virus_name) %>%
      select(feature = gene, start = concat_start, end = concat_end, 
             original_start = cds_start, original_end = cds_end,
             mid = concat_mid) %>%
      arrange(start)
    
    # Create function to determine y-axis limits automatically
    calculate_y_limits <- function(data, y_var, method_name) {
      y_values <- as.numeric(data[[y_var]])
      y_values <- y_values[!is.na(y_values)]
      
      if (length(y_values) == 0) {
        return(c(0, 1, 0.2))
      }
      
      # Handle extreme values (like -99999 for FEL)
      if (method_name == "FEL") {
        y_values <- y_values[y_values > -1000]
      }
      
      y_min <- floor(min(y_values, na.rm = TRUE) * 2) / 2
      y_max <- ceiling(max(y_values, na.rm = TRUE) * 2) / 2
      
      # Ensure reasonable limits for each method
      if (method_name == "SLAC") {
        y_min <- max(y_min, 0)
        y_max <- min(y_max, 3)
        if (y_max - y_min < 0.5) y_max <- y_min + 1
      } else if (method_name == "FEL") {
        y_min <- max(y_min, -2)
        y_max <- min(y_max, 3)
        if (y_max - y_min < 0.5) y_max <- y_min + 1
      } else if (method_name == "FUBAR") {
        y_min <- max(y_min, 0)
        y_max <- min(y_max, 1)
        if (y_max - y_min < 0.1) y_max <- y_min + 0.2
      }
      
      # Calculate nice breaks
      range_diff <- y_max - y_min
      if (range_diff <= 1) {
        breaks_by <- 0.2
      } else if (range_diff <= 2) {
        breaks_by <- 0.5
      } else {
        breaks_by <- 1
      }
      
      return(c(y_min, y_max, breaks_by))
    }
    
    # Create individual plots for each method
    create_method_plot <- function(method_data, method_name, y_var, sel_var, 
                                   y_label, show_x = FALSE) {
      
      # Calculate y-axis limits automatically
      y_limits <- calculate_y_limits(method_data, y_var, method_name)
      
      # Prepare data
      if (highlight_all_pos) {
        highlight_data <- method_data %>%
          filter(get(sel_var) == "pos")
        nonpos_data <- method_data %>%
          filter(get(sel_var) == "nonpos")
      } else {
        highlight_data <- method_data %>%
          filter(all_three_pos == "Positive")
        nonpos_data <- method_data %>%
          filter(is.na(all_three_pos) | all_three_pos != "Positive")
      }
      # Start plot
      p <- ggplot()
      
      # Add vertical dotted lines to separate genes
      if (nrow(gene_pos) > 1) {
        for (i in 1:(nrow(gene_pos) - 1)) {
          gene_end <- gene_pos$end[i]
          p <- p +
            geom_vline(
              xintercept = gene_end + 0.5,
              linetype = "dotted",
              color = "gray50",
              linewidth = 0.3,
              alpha = 0.7
            )
        }
      }
      # Add non-positive points as dark gray dots if show_nonpos_points is TRUE
      if (show_nonpos_points && nrow(nonpos_data) > 0) {
        p <- p +
          geom_point(
            data = nonpos_data,
            aes(x = plot_x, y = as.numeric(get(y_var))),
            shape = point_shape,
            color = "#222222",
            stroke = 0.1,
            fill = "gray30",
            size = point_size,
            alpha = 0.7
          )
      }
      
      # Add positive points (red circles)
      if (nrow(highlight_data) > 0) {
        p <- p +
          geom_point(
            data = highlight_data,
            aes(x = plot_x, y = as.numeric(get(y_var))),
            shape = point_shape,
            color = "#222222",
            stroke = 0.1,
            fill = "darkred",
            size = point_size
          ) +
          geom_point(
            data = highlight_data,
            aes(x = plot_x, y = as.numeric(get(y_var))),
            shape = highlight_shape,
            fill = NA,
            color = "darkred",
            stroke = 0.3,
            size = highlight_size
          )
      }
      # Set axis limits
      y_min <- y_limits[1]
      y_max <- y_limits[2]
      breaks_by <- y_limits[3]
      p <- p +
        scale_y_continuous(
          limits = c(y_min, y_max),
          breaks = seq(y_min, y_max, by = breaks_by),
          expand = expansion(mult = 0.05)
        ) +
        scale_x_continuous(
          limits = c(min(method_data$plot_x), max(method_data$plot_x)),
          expand = expansion(mult = 0.02)
        ) +
        ylab(y_label) +
        theme_bw(base_size = 11) +
        theme(
          panel.border = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 5)),
          plot.margin = margin(5, 10, 5, 5),
          legend.position = "none"
        )
      
      if (show_x) {
        p <- p + 
          xlab("Codon position (genes concatenated)") +
          theme(
            axis.text.x = element_text(size = 9, color = "black"),
            axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 5))
          )
      } else {
        p <- p +
          xlab(NULL) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
      }
      return(p)
    }
    # Create plots
    slac_plot <- create_method_plot(
      virus_data,
      "SLAC",
      "slac_pval",
      "slac_sel",
      y_label = "p-value, dN/dS > 1"
    )
    
    fel_plot <- create_method_plot(
      virus_data,
      "FEL",
      "fel_dif",
      "fel_sel",
      y_label = "β minus α"
    )
    fubar_plot <- create_method_plot(
      virus_data,
      "FUBAR",
      "fubar_pval",
      "fubar_sel",
      y_label = "p-value, α < β",
      show_x = TRUE
    )
    # Create gene structure bar
    n_genes <- nrow(gene_pos)
    gene_colors <- hue_pal()(n_genes)
    
    # Calculate label positions with automatic nudge
    gene_pos <- gene_pos %>%
      mutate(
        # Calculate gene width
        gene_width = end - start,
        # Calculate distance to next gene label
        next_mid = lead(mid, default = max(mid) + 1000),
        prev_mid = lag(mid, default = -1000),
        # Calculate minimum distance to neighboring labels
        min_distance = pmin(mid - prev_mid, next_mid - mid),
        # Determine if nudge is needed (if labels would be too close)
        needs_nudge = gene_width < 100 | min_distance < 100,  # Adjust thresholds as needed
        # Calculate nudge direction and amount
        nudge_x = ifelse(needs_nudge,
                         ifelse(mid - prev_mid < next_mid - mid, 
                                -gene_width * 0.4,  # Nudge left if closer to left neighbor
                                gene_width * 0.4),   # Nudge right if closer to right neighbor
                         0),
        # Alternative: Simple staggered approach for very close genes
        stagger = ifelse(row_number() %% 2 == 0, 0.2, -0.2),  # Alternate up/down
        # Final x position for label
        label_x = ifelse(gene_width < 50, mid + nudge_x, mid),
        # Use staggered y-position for very narrow genes
        label_y = ifelse(gene_width < 50, 
                         gene_bar_width / 2 + (gene_bar_width * stagger),
                         gene_bar_width / 2)
      )
    gene_bar <- ggplot(gene_pos) +
      geom_rect(
        aes(
          xmin = start,
          xmax = end,
          ymin = 0,
          ymax = gene_bar_width,
          fill = feature
        ),
        color = "black",
        linewidth = 0.2
      ) +
      # Add gene labels with automatic positioning
      geom_text(
        aes(
          x = label_x,
          y = label_y,
          label = feature
        ),
        fontface = "bold",
        size = gene_label_size,
        vjust = 0.5,
        hjust = 0.5,
        color = "black",
        check_overlap = TRUE  # This prevents ggplot from drawing overlapping labels
      ) +
      # Alternative: Use geom_text_repel for automatic label repulsion (requires ggrepel)
      # geom_text_repel(
      #   aes(
      #     x = mid,
      #     y = gene_bar_width / 2,
      #     label = feature
      #   ),
      #   fontface = "bold",
      #   size = gene_label_size,
      #   color = "white",
      #   segment.color = NA,  # Remove connecting lines
      #   force = 0.5,
      #   max.overlaps = Inf
      # ) +
      scale_fill_manual(values = setNames(gene_colors, gene_pos$feature)) +
      scale_x_continuous(
        limits = c(min(virus_data$plot_x), max(virus_data$plot_x)),
        expand = expansion(mult = 0.02)
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)
      )
    
    # Combine plots
    panel_c <- ggarrange(
      fubar_plot,
      gene_bar,
      ncol = 1,
      heights = c(4, 1),
      align = "v"
    )
    
    combined_plot <- ggarrange(
      slac_plot,
      fel_plot,
      panel_c,
      labels = c("A", "B", "C"),
      font.label = list(size = 12, face = "bold"),
      nrow = 3,
      heights = c(1, 1, 1.2),
      align = "v"
    )
    
    # Add title
    combined_plot <- annotate_figure(
      combined_plot,
      top = text_grob(
        paste("Virus:", virus_name, "(", nrow(gene_pos), "genes concatenated)"),
        face = "bold",
        size = 14
      )
    )
    
    # Save plot
    filename <- paste0(
      output_dir, "/", 
      gsub("[^[:alnum:]]", "_", virus_name), 
      "_highlight", ifelse(highlight_all_pos, "All", "Consensus"),
      "_nonpos", ifelse(show_nonpos_points, "On", "Off"),
      ".png"
    )
    
    ggsave(
      filename,
      plot = combined_plot,
      width = 12,
      height = 10,
      dpi = 300,
      bg = "white"
    )
    
    # Store in list
    all_plots[[virus_name]] <- list(
      slac = slac_plot,
      fel = fel_plot,
      fubar = fubar_plot,
      gene_bar = gene_bar,
      combined = combined_plot,
      filename = filename,
      gene_positions = gene_pos
    )
    
    cat("Plot saved for virus:", virus_name, "->", filename, "\n")
  }
  
  # Save tables to CSV
  table_filename <- paste0(output_dir, "/positive_sites_table.csv")
  write.csv(positive_sites_table, table_filename, row.names = FALSE)
  cat("\nPositive sites table saved to:", table_filename, "\n")
  
  if (!is.null(gene_dnds_statistics)) {
    gene_stats_file <- paste0(output_dir, "/gene_dnds_statistics.csv")
    write.csv(gene_dnds_statistics, gene_stats_file, row.names = FALSE)
    cat("Gene dN/dS statistics saved to:", gene_stats_file, "\n")
  }
  
  # Return comprehensive results
  return(list(
    plots = all_plots,
    positive_sites_table = positive_sites_table,
    gene_dnds_statistics = gene_dnds_statistics,
    processed_data = processed_df,
    table_files = list(
      positive_sites = table_filename,
      gene_dnds = if (!is.null(gene_dnds_statistics)) gene_stats_file else NULL
    )
  ))
}

create_publication_table <- function(selection_results, 
                                     caption_pos = "Sites under positive selection identified by all three methods",
                                     caption_dnds = "Gene-level dN/dS statistics",
                                     label_pos = "tab:positive_sites",
                                     label_dnds = "tab:gene_dnds_stats",
                                     output_format = "html") {
  
  # Extract data from selection_results
  positive_sites_table <- selection_results$positive_sites_table
  gene_dnds_statistics <- selection_results$gene_dnds_statistics
  
  # Create a list to store all outputs
  result <- list(
    positive_sites_data = positive_sites_table,
    gene_dnds_data = gene_dnds_statistics,
    formatted_positive_sites = NULL,
    formatted_gene_dnds = NULL
  )
  
  # Create formatted tables if knitr and kableExtra are available
  if (requireNamespace("knitr", quietly = TRUE) && 
      requireNamespace("kableExtra", quietly = TRUE)) {
    
    # 1. Format positive sites table
    if (nrow(positive_sites_table) > 0) {
      result$formatted_positive_sites <- positive_sites_table %>%
        knitr::kable(
          format = output_format,
          caption = caption_pos,
          align = c("l", "l", "l", "l", "c"),
          booktabs = TRUE
        ) %>%
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE,
          position = "center"
        ) %>%
        kableExtra::column_spec(1, bold = TRUE) %>%
        kableExtra::column_spec(3, width = "3cm") %>%
        kableExtra::column_spec(4, width = "2.5cm") %>%
        kableExtra::row_spec(0, bold = TRUE, background = "#f0f0f0")
    }
    
    # 2. Format gene dN/dS statistics table
    if (!is.null(gene_dnds_statistics) && nrow(gene_dnds_statistics) > 0) {
      # Clean column names for display
      gene_dnds_display <- gene_dnds_statistics %>%
        rename(
          `Number of sites` = n_sites,
          `Mean dN` = mean_dN,
          `Mean dS` = mean_dS,
          `Mean dN/dS` = mean_dNdS,
          `Median dN/dS` = median_dNdS,
          `Minimum dN/dS` = min_dNdS,
          `Maximum dN/dS` = max_dNdS,
          `SD dN/dS` = sd_dNdS
        )
      
      result$formatted_gene_dnds <- gene_dnds_display %>%
        knitr::kable(
          format = output_format,
          caption = caption_dnds,
          align = c("l", "l", "c", "c", "c", "c", "c", "c", "c"),
          booktabs = TRUE
        ) %>%
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE,
          position = "center"
        ) %>%
        kableExtra::column_spec(1, bold = TRUE) %>%
        kableExtra::row_spec(0, bold = TRUE, background = "#f0f0f0")
    }
    
  } else {
    warning("knitr and/or kableExtra not installed. Returning data frames only.")
  }
  
  # Return all components
  return(result)
}

# Example usage:
results <- make_virus_selection_plots(
  excel_path = "01_src/rabv_glob_metadata.xlsx",
  sheet_number = 4,
  #show_neutral_points = TRUE,
  #use_cds = FALSE,  # Set to FALSE for individual genes
  highlight_all_pos = TRUE,  # Set to FALSE for only all-three-positive sites
  gene_bar_width = 0.6,
  gene_bar_nudge_x = 0,
  output_dir = "06_output/rabv_selection_plots"
)
tables <- create_publication_table(
  selection_results = results,
  caption_pos = "Positively selected sites across all three methods",
  caption_dnds = "Gene-level evolutionary statistics"
)
# Access the formatted tables
if (!is.null(tables$formatted_positive_sites)) {
  cat(as.character(tables$formatted_positive_sites))
}
if (!is.null(tables$formatted_gene_dnds)) {
  cat(as.character(tables$formatted_gene_dnds))
}
# Save tables to HTML files
if (!is.null(tables$formatted_positive_sites)) {
  writeLines(as.character(tables$formatted_positive_sites), 
             "results/positive_sites_table.html")
}
if (!is.null(tables$formatted_gene_dnds)) {
  writeLines(as.character(tables$formatted_gene_dnds), 
             "results/gene_dnds_statistics.html")
}
################################### 5. RABV PHYLOGENY: PHYLODYNAMICS & PHLOGEOGRAPHY ##############################
##### 5A. MAIN SOURCE CODES: PHYLOGEOGRAPHY WITH SERAPHIM #####
{

library(diagram)
require(seraphim)
library(rgdal)
library(lubridate)
library(doMC)
source("07_phylogeo/mccExtractions.r")
library("readxl")
library(dplyr)
library(maptools)
library(phyloch)
locmeta <- read.csv("07_phylogeo/rabv_swk_loc_cds_nhm_latlon.txt")
# Note: The following process describes how to generate the phylogeography maps like in Figure 4 for the cluster C.1.
#To generate the same for other cluster, just substitute the .trees, .tree, and .asc files
# 1. Extracting the spatio-temporal information contained in posterior trees
#treefiles<- "rabv_nhmprecov_cRRW_HKY-G_RCLN_CoalesBayesSG_500m5k.trees.txt"
treefiles<- "rabv_nhmprecov_cRRW_StrictC_RCLN_CoalesBayesSG_500m5k.trees.txt"
allTrees = scan(file=treefiles, what="", sep="\n", quiet=T)
localTreesDirectory = "07_phylogeo/swk_rabv_extr_postcovStrict"
burnIn = 0.1
randomSampling = FALSE
nberOfTreesToSample = 85
mostRecentSamplingDatum = decimal_date(ymd("2021-10-13"))
coordinateAttributeName = "location"
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)

# 2. Extracting the spatio-temporal information embedded in the MCC tree
treefile <- "07_phylogeo/rabv_nhm_cRRW_HKY-G_RCLN_CoalesBayesSG_500m5k.treeanno.txt"
#beast_tree <- read.beast(treefile)
#beast2 <- beast_tree[-(1:1), ]
#beast2 <- apply(beast2,2,as.character)
#write.csv(beast2, file="output/rabv_nhm_cRRW_HKY.csv")
mostRecentSamplingDatum = decimal_date(ymd("2024-10-13"))
mcc_tre = readAnnotatedNexus(treefile)
mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
mcc_tabdf <- as.data.frame(mcc_tab)
mcc_tabdf$startYear_round <- round(mcc_tabdf$startYear ,digit=0)
mcc_tabdf = subset(mcc_tabdf, mcc_tabdf$startYear_round > 2017)
mcc_tabdf = subset(mcc_tabdf, mcc_tabdf$endLat  < 4.5)  #remove Limbang entry to simplify analysis
# 3. Estimating the HPD region for each time slice
nberOfExtractionFiles = nberOfTreesToSample
prob = 0.95; precision = 0.75
startDatum = min(mcc_tabdf[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


# 4.1 spatial boundaries
my_spdf <- readOGR( 
  dsn= paste0("swk_maps/") , 
  layer="swk_divmap_2022",
  verbose=FALSE
)
capital <- my_spdf[my_spdf@data$ADM2_EN=='Kuching']

template_raster = raster("swk_maps/swk_divmap_raster.asc")
#borders = crop(getData("GADM", country="ZAF", level=1), extent(template_raster))
borders = crop(capital, extent(template_raster))

# 4.2 Defining the different colour scales
minYear = min(mcc_tabdf[,"startYear"]); maxYear = max(mcc_tabdf[,"endYear"])
endYears_indices = (((mcc_tabdf[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
##colors have to go to max time slice - AT LEAST
n_number_colours_needed<- max(round(endYears_indices))
n_repeats_discrete<- 10
c1<- rev(brewer.pal(9,"YlOrBr"))
#c1<- rev(brewer.pal(4,"PuRd"))
c2<- (brewer.pal(4,"Blues"))
colours<- rev(rep(c(c1,c2), each=n_repeats_discrete))
colour_scale<- colorRampPalette(colours)(n_number_colours_needed)

endYears_colours = colour_scale[round(endYears_indices)]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"55")
}

# 5. Generating the dispersal history plot
pdf('map_swk_rabv_postcovStrict_dispphist.pdf',width=12, height=7,bg="white")
ptsize<- 0.35
pitjit<- 0.01
par(mar=c(0,0,0,0), oma=c(1.2,2.4,1,0), mgp=c(0,0.2,0), lwd=0.1, bty="o")
plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
plot(borders, add=T, lwd=0.1, border="gray10")
#plot(borders, lwd=0.1, border="gray10")
for (i in length(polygons):0.1)
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tabdf)[1])
{
  # curvedarrow(cbind(mcc_tabdf[i,"startLon"],mcc_tabdf[i,"startLat"]), cbind(mcc_tabdf[i,"endLon"],mcc_tabdf[i,"endLat"]), arr.length=0,
  # arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tabdf[i,"startLon"],mcc_tabdf[i,"startLat"]), cbind(mcc_tabdf[i,"endLon"],mcc_tabdf[i,"endLat"]), arr.length=0,
              arr.width=0.00, lwd=1.6*1.1, lty=1, lcol="grey22", arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tabdf[i,"startLon"],mcc_tabdf[i,"startLat"]), cbind(mcc_tabdf[i,"endLon"],mcc_tabdf[i,"endLat"]), arr.length=0,
              arr.width=0.00, lwd=1.0*1.1, lty=1, lcol=endYears_colours[i], arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
}
for (i in dim(mcc_tabdf)[1]:1)
{
  xs<- mcc_tabdf[i,"startLon"]
  ys<- mcc_tabdf[i,"startLat"]
  xe<- jitter(mcc_tabdf[i,"endLon"],pitjit)
  ye<- jitter(mcc_tabdf[i,"endLat"],pitjit)
  if (i == 1)
  {
    points(xs, ys, pch=16, col=colour_scale[1], cex=ptsize)
    points(xs, ys, pch=1, col="gray10", cex=ptsize)
  }
  points(xe, ye, pch=16, col=endYears_colours[i], cex=ptsize)
  points(xe, ye, pch=1, col="gray10", cex=ptsize)
}

xrange<- c(xmin(template_raster), xmax(template_raster))
yrange<- c(ymin(template_raster), ymax(template_raster))
rect(xrange[1], yrange[1], xrange[2], yrange[2], xpd=T, lwd=0.2)
axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tabdf[,"startYear"]); rast[2] = max(mcc_tabdf[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.40,0.6,0.62),
     legend.args=list(text="", cex=0.1, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.08,0)))
pointLabel(x= locmeta$lat, y= locmeta$lon, labels = locmeta$location)
a<-dev.off()
}
##### 5B. MAIN SOURCE CODES: PHYLOGEOGRAPHY WITH SERAPHIM VERSION 2 #####
suppressPackageStartupMessages({
  library(seraphim); library(diagram); library(rgdal); library(maptools)
  library(raster); library(lubridate); library(dplyr); library(RColorBrewer); library(shape)
  library(phyloch); library(sf); library(ggplot2); library(gganimate); library(raster)
  library(remotes)
})
### Section 1: MCC Tree Processing Functions
treefiles = "07_phylogeo/rabv_nhm_cRRW_HKY-G_RCLN_CoalesBayesSG_500m5k.trees.txt"
mcc_treefiles = "07_phylogeo/rabv_nhm_cRRW_HKY-G_RCLN_CoalesBayesSG_500m5k.treeanno.txt"
localTreesDirectory = "07_phylogeo/trees/"
nberOfExtractionFiles = 85
prob = 0.95
precision = 0.75
convert_to_sf = FALSE
color_scheme = "plasma"
borders_fill_col = 'ADM2_REF'

shp_dsn = "swk_maps"
shp_layer = "swk_divmap_2022"
raster_file = "swk_maps/swk_divmap_raster.asc"
color_scheme = "plasma"
borders_fill_col = 'ADM2_REF'
borders_palette = "plasma"
borders_alpha = 0.8
arrow_length = 0.08
arrow_width = 0.06
curve = 0.3
n_points = 100
pitjit = 0.01
convert_to_sf = FALSE
borders_color_reverse = TRUE
intensity_threshold = 0.95
min_distance_km = 50
show_uncertainty_legend = TRUE
polygon_shape = "ellipse"
polygon_alpha = 0.4
uncertainty_metric = "std_dev"
kernel_bandwidth = 0.1
uncertainty_legend_position = "topright"
min_events_for_uncertainty = 2

curvedarrow_lwd = 1.1
arrowhead_length = 0.08
arrowhead_width = 0.06
point_cex = 0.35
annotate_year = FALSE
intensity_threshold = 0.95
solid_lwd_multiplier = 2.0
dashed_lwd_multiplier = 1.0
dashed_lty = 2
show_dashed_lines = FALSE
dashed_line_alpha = 0.6
use_curvedarrow = TRUE
label_lines = "none"
show_polygon_legend = TRUE
polygon_color = "steelblue3"

polygon_style = "sphere_events"
sphere_alpha = 0.4
sphere_scaling = 1.0
bandwidth_months = 2
distance_ribbon_alpha = 0.25
mostRecentSamplingDate = 2022.7
raster_file
startYearFilter = 2014
maxLatitude = 4.5
intensity_threshold = 0.95
output_pdf = NULL
plot_temporal = TRUE
save_distance_plot = NULL
### Section 1: Data Extraction
mccTreeExtraction <- function(mcc_tre, mostRecentSamplingDatum) {
  # Validate input
  if (is.null(mcc_tre) || is.null(mcc_tre$edge)) {
    warning("Invalid MCC tree structure")
    return(data.frame())
  }
  # Initialize matrix with standard columns
  n_edges <- nrow(mcc_tre$edge)
  mcc_tab <- matrix(nrow = n_edges, ncol = 14)
  colnames(mcc_tab) <- c("node1", "node2", "length", "startLon", "startLat", 
                         "endLon", "endLat", "startNodeL", "endNodeL", 
                         "startYear", "endYear", "location_code", "duration", "posterior_prob")
  # Fill basic edge information
  mcc_tab[, 1:2] <- mcc_tre$edge
  mcc_tab[, "length"] <- mcc_tre$edge.length
  
  # Initialize with default values
  mcc_tab[, "posterior_prob"] <- 1.0
  mcc_tab[, "location_code"] <- NA_character_
  mcc_tab[, "duration"] <- mcc_tab[, "length"]
  
  # Extract information from annotations
  for (i in seq_len(n_edges)) {
    # Get annotation for this edge
    if (i <= length(mcc_tre$annotations) && !is.null(mcc_tre$annotations[[i]])) {
      ann <- mcc_tre$annotations[[i]]
    } else {
      ann <- list()
    }
    
    # Extract coordinates safely
    if (!is.null(ann$location2) && !is.null(ann$location1)) {
      mcc_tab[i, "endLon"] <- as.numeric(ann$location2)
      mcc_tab[i, "endLat"] <- as.numeric(ann$location1)
    } else if (!is.null(ann$longitude) && !is.null(ann$latitude)) {
      mcc_tab[i, "endLon"] <- as.numeric(ann$longitude)
      mcc_tab[i, "endLat"] <- as.numeric(ann$latitude)
    } else {
      mcc_tab[i, c("endLon", "endLat")] <- NA_real_
    }
    
    # Extract node height
    if (!is.null(ann$height)) {
      mcc_tab[i, "endNodeL"] <- as.numeric(ann$height)
    } else if (!is.null(ann$node_height)) {
      mcc_tab[i, "endNodeL"] <- as.numeric(ann$node_height)
    } else {
      mcc_tab[i, "endNodeL"] <- NA_real_
    }
    
    # Extract posterior probability
    if (!is.null(ann$posterior)) {
      mcc_tab[i, "posterior_prob"] <- as.numeric(ann$posterior)
    } else if (!is.null(ann$probability)) {
      mcc_tab[i, "posterior_prob"] <- as.numeric(ann$probability)
    } else if (!is.null(ann$prob)) {
      mcc_tab[i, "posterior_prob"] <- as.numeric(ann$prob)
    } else if (!is.null(ann$posterior_probability)) {
      mcc_tab[i, "posterior_prob"] <- as.numeric(ann$posterior_probability)
    }
    
    # Create location code
    if (!is.na(mcc_tab[i, "endLon"]) && !is.na(mcc_tab[i, "endLat"])) {
      mcc_tab[i, "location_code"] <- paste(
        sprintf("%.4f", as.numeric(mcc_tab[i, "endLon"])),
        sprintf("%.4f", as.numeric(mcc_tab[i, "endLat"])),
        sep = "_"
      )
    }
  }
  
  # Trace back to parent nodes to get start positions
  for (i in seq_len(n_edges)) {
    parent_idx <- which(mcc_tab[, "node2"] == mcc_tab[i, "node1"])
    
    if (length(parent_idx) > 0) {
      # Copy from parent's end position
      parent_row <- parent_idx[1]
      mcc_tab[i, "startLon"] <- mcc_tab[parent_row, "endLon"]
      mcc_tab[i, "startLat"] <- mcc_tab[parent_row, "endLat"]
      mcc_tab[i, "startNodeL"] <- mcc_tab[parent_row, "endNodeL"]
      mcc_tab[i, "location_code"] <- mcc_tab[parent_row, "location_code"]
    } else {
      # Root node - get from root annotation
      ann <- mcc_tre$root.annotation
      if (!is.null(ann)) {
        if (!is.null(ann$location2) && !is.null(ann$location1)) {
          mcc_tab[i, "startLon"] <- as.numeric(ann$location2)
          mcc_tab[i, "startLat"] <- as.numeric(ann$location1)
        } else if (!is.null(ann$longitude) && !is.null(ann$latitude)) {
          mcc_tab[i, "startLon"] <- as.numeric(ann$longitude)
          mcc_tab[i, "startLat"] <- as.numeric(ann$latitude)
        } else {
          mcc_tab[i, c("startLon", "startLat")] <- NA_real_
        }
        
        if (!is.null(ann$height)) {
          mcc_tab[i, "startNodeL"] <- as.numeric(ann$height)
        } else if (!is.null(ann$node_height)) {
          mcc_tab[i, "startNodeL"] <- as.numeric(ann$node_height)
        } else {
          mcc_tab[i, "startNodeL"] <- NA_real_
        }
        
        # Create location code for root
        if (!is.na(mcc_tab[i, "startLon"]) && !is.na(mcc_tab[i, "startLat"])) {
          mcc_tab[i, "location_code"] <- paste(
            sprintf("%.4f", as.numeric(mcc_tab[i, "startLon"])),
            sprintf("%.4f", as.numeric(mcc_tab[i, "startLat"])),
            sep = "_"
          )
        }
      } else {
        mcc_tab[i, c("startLon", "startLat", "startNodeL")] <- NA_real_
      }
    }
    
    # Calculate years
    if (!is.na(mcc_tab[i, "startNodeL"])) {
      mcc_tab[i, "startYear"] <- mostRecentSamplingDatum - as.numeric(mcc_tab[i, "startNodeL"])
    } else {
      mcc_tab[i, "startYear"] <- NA_real_
    }
    
    if (!is.na(mcc_tab[i, "endNodeL"])) {
      mcc_tab[i, "endYear"] <- mostRecentSamplingDatum - as.numeric(mcc_tab[i, "endNodeL"])
    } else {
      mcc_tab[i, "endYear"] <- NA_real_
    }
  }
  
  # Convert to data frame
  mcc_tab <- as.data.frame(mcc_tab, stringsAsFactors = FALSE)
  
  # Convert columns to appropriate types
  numeric_cols <- c("node1", "node2", "length", "startLon", "startLat", 
                    "endLon", "endLat", "startNodeL", "endNodeL", 
                    "startYear", "endYear", "duration", "posterior_prob")
  
  for (col in numeric_cols) {
    if (col %in% names(mcc_tab)) {
      mcc_tab[[col]] <- as.numeric(as.character(mcc_tab[[col]]))
    }
  }
  
  # Remove rows with missing essential data
  mcc_tab <- mcc_tab[complete.cases(mcc_tab[, c("startLon", "startLat", "endLon", "endLat")]), ]
  
  if (nrow(mcc_tab) == 0) {
    return(data.frame())
  }
  
  # Order by startYear, then by endYear
  mcc_tab <- mcc_tab[order(mcc_tab$startYear, mcc_tab$endYear), ]
  rownames(mcc_tab) <- NULL
  
  return(mcc_tab)
}
extract_mcc_data <- function(mcc_treefiles, mostRecentSamplingDate, 
                             startYearFilter = 2010, maxLatitude = 4.5) {
  # Load required package
  if (!requireNamespace("seraphim", quietly = TRUE)) stop("Package 'seraphim' required")
  
  mcc_list <- list()
  for (i in seq_along(mcc_treefiles)) {
    tree_file <- mcc_treefiles[i]
    message(sprintf("  Processing MCC tree %d/%d: %s", i, length(mcc_treefiles), basename(tree_file)))
    
    tryCatch({
      tre <- seraphim::readAnnotatedNexus(tree_file)
      df <- mccTreeExtraction(tre, mostRecentSamplingDate)
      
      if (nrow(df) > 0) {
        df$tree_id <- paste0("Tree_", i)
        df$tree_file <- basename(tree_file)
        mcc_list[[i]] <- df
        message(sprintf("    Extracted %d dispersal events", nrow(df)))
      } else {
        warning(sprintf("    No data extracted from %s", basename(tree_file)))
      }
    }, error = function(e) {
      warning(sprintf("    Error processing %s: %s", basename(tree_file), e$message))
    })
  }
  mcc_list <- mcc_list[!sapply(mcc_list, is.null)]
  if (length(mcc_list) == 0) {
    stop("No valid data extracted from any MCC tree files.")
  }
  mcc_tabdf <- do.call(rbind, mcc_list)
  if (nrow(mcc_tabdf) == 0) {
    stop("No dispersal events extracted from MCC trees.")
  }
  # Filter data
  mcc_tabdf <- mcc_tabdf |>
    dplyr::mutate(startYear_round = round(startYear)) |>
    dplyr::filter(startYear_round > startYearFilter, 
                  endLat < maxLatitude,
                  !is.na(startLon), !is.na(startLat),
                  !is.na(endLon), !is.na(endLat))
  message(sprintf("Total: %d dispersal events from %d MCC trees", 
                  nrow(mcc_tabdf), length(mcc_list)))
  return(mcc_tabdf)
}
mcc_tabdf <- extract_mcc_data(mcc_treefiles, mostRecentSamplingDate)
### Section 2: Data Processing and Preparation
process_spatial_data <- function(mcc_tabdf, shp_dsn, shp_layer, raster_file,
                                 localTreesDirectory = "07_phylogeo/trees/",
                                 nberOfExtractionFiles = 85,
                                 prob = 0.95,
                                 precision = 0.75,
                                 convert_to_sf = FALSE,
                                 color_scheme = "plasma",
                                 borders_fill_col = 'ADM2_REF',
                                 borders_alpha = 0.8,
                                 borders_color_reverse = FALSE,
                                 intensity_threshold = 0.95,
                                 min_distance_km = 0,
                                 polygon_alpha = 0.33) {  # New parameter for polygon transparency
  
  # Load required packages
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' required")
  if (!requireNamespace("raster", quietly = TRUE)) stop("Package 'raster' required")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Package 'viridis' required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required")
  if (!requireNamespace("geosphere", quietly = TRUE)) stop("Package 'geosphere' required")
  if (!requireNamespace("rgdal", quietly = TRUE)) stop("Package 'rgdal' required")
  library("seraphim")
  library("ks")
  ## ---- Calculate dispersal distances ----
  if (min_distance_km > 0) {
    message(sprintf("Calculating dispersal distances with minimum threshold: %s km", min_distance_km))
    
    mcc_tabdf$dispersal_distance_km <- NA_real_
    
    for (i in seq_len(nrow(mcc_tabdf))) {
      if (!is.na(mcc_tabdf$startLon[i]) && !is.na(mcc_tabdf$startLat[i]) &&
          !is.na(mcc_tabdf$endLon[i]) && !is.na(mcc_tabdf$endLat[i])) {
        
        dist_m <- geosphere::distGeo(
          c(mcc_tabdf$startLon[i], mcc_tabdf$startLat[i]),
          c(mcc_tabdf$endLon[i], mcc_tabdf$endLat[i])
        )
        mcc_tabdf$dispersal_distance_km[i] <- dist_m / 1000
      }
    }
    # Create mask for events to show
    n_before <- nrow(mcc_tabdf)
    show_event <- !is.na(mcc_tabdf$dispersal_distance_km) & 
      mcc_tabdf$dispersal_distance_km >= min_distance_km
    mcc_tabdf$show_event <- show_event
    
    n_after <- sum(show_event, na.rm = TRUE)
    message(sprintf("Filtering dispersal events: %d of %d events meet minimum distance of %s km",
                    n_after, n_before, min_distance_km))
  } else {
    mcc_tabdf$show_event <- TRUE
  }
  
  ## ---- Calculate dispersal intensity and line type ----
  mcc_tabdf$dispersal_intensity <- mcc_tabdf$posterior_prob
  
  mcc_tabdf$line_type <- ifelse(mcc_tabdf$dispersal_intensity >= intensity_threshold, 
                                "solid", "dashed")
  
  mcc_filtered <- mcc_tabdf[mcc_tabdf$show_event, ]
  
  ## ---- Generate HPD polygons ----
  if (!requireNamespace("seraphim", quietly = TRUE)) {
    warning("Package 'seraphim' not available. Skipping HPD polygon generation.")
    polygons <- list()
    polygons_colours <- character(0)
  } else {
    message("Generating HPD polygons for each time slice...")
    
    startDatum <- min(mcc_filtered$startYear, na.rm = TRUE)
    
    if (exists("spreadGraphic2", mode = "function")) {
      message(sprintf("Generating HPD polygons with prob=%.2f, precision=%.2f, startDatum=%.2f", 
                      prob, precision, startDatum))
      
      tryCatch({
        polygons <- suppressWarnings(
          spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision)
        )
        
        message(sprintf("Successfully generated %d HPD polygons", length(polygons)))
        
      }, error = function(e) {
        warning(sprintf("Error generating HPD polygons: %s", e$message))
        polygons <- list()
      })
    } else {
      warning("spreadGraphic2 function not found. Skipping HPD polygons.")
      polygons <- list()
    }
  }
  
  ## ---- Generate color scale ----
  minYear <- min(mcc_filtered$startYear, na.rm = TRUE)
  maxYear <- max(mcc_filtered$endYear, na.rm = TRUE)
  
  endYears_indices <- (((mcc_filtered$endYear - minYear) / (maxYear - minYear)) * 100) + 1
  n_number_colours_needed <- max(round(endYears_indices))
  
  if (color_scheme == "plasma") {
    plasma_full <- viridis::plasma(256)
    colour_scale <- plasma_full[1:230]
      colour_scale <- colorRampPalette(colour_scale)(n_number_colours_needed)
  } else if (color_scheme == "YlOrBr") {
    c1 <- rev(RColorBrewer::brewer.pal(9, "YlOrBr"))
    c2 <- RColorBrewer::brewer.pal(4, "Blues")
    n_repeats_discrete <- 10
    colours <- rev(rep(c(c1, c2), each = n_repeats_discrete))
    colour_scale <- colorRampPalette(colours)(n_number_colours_needed)
  } else if (color_scheme == "PuRd") {
    c1 <- rev(RColorBrewer::brewer.pal(9, "PuRd"))
    c2 <- RColorBrewer::brewer.pal(4, "Blues")
    n_repeats_discrete <- 10
    colours <- rev(rep(c(c1, c2), each = n_repeats_discrete))
    colour_scale <- colorRampPalette(colours)(n_number_colours_needed)
  } else {
    colour_scale <- viridis::viridis(n_number_colours_needed)
  }
  
  # Assign colors to events
  mcc_filtered$year_color <- colour_scale[round(endYears_indices)]
  
  ## ---- Assign colors to HPD polygons (FIXED - no 8-digit hex) ----
  if (length(polygons) > 0) {
    polygons_colours <- rep(NA, length(polygons))
    for (i in seq_along(polygons)) {
      date <- as.numeric(names(polygons[[i]]))
      polygon_index <- round((((date - minYear) / (maxYear - minYear)) * 100) + 1)
      
      # Get base color from colour_scale
      base_color <- colour_scale[polygon_index]
      
      # Convert to proper color with transparency using adjustcolor
      # polygon_alpha = 0.33 corresponds to hex "55" (85/255 ≈ 0.33)
      polygons_colours[i] <- adjustcolor(base_color, alpha.f = polygon_alpha)
    }
  } else {
    polygons_colours <- character(0)
  }
  
  ## ---- Load spatial data ----
  shp_file <- file.path(shp_dsn, paste0(shp_layer, ".shp"))
  if (!file.exists(shp_file)) stop("Shapefile not found: ", shp_file)
  if (!file.exists(raster_file)) stop("Raster file not found: ", raster_file)
  
  my_spdf <- rgdal::readOGR(dsn = shp_dsn, layer = shp_layer, verbose = FALSE)
  
  if (!is.null(borders_fill_col) && borders_fill_col == 'capital') {
    borders <- my_spdf[my_spdf@data$ADM2_EN == 'Kuching', ]
  } else {
    borders <- my_spdf
  }
  
  template_raster <- raster::raster(raster_file)
  borders <- raster::crop(borders, raster::extent(template_raster))
  
  return(list(
    mcc_tabdf = mcc_tabdf,
    mcc_filtered = mcc_filtered,
    colour_scale = colour_scale,
    minYear = minYear,
    maxYear = maxYear,
    borders = borders,
    template_raster = template_raster,
    polygons = polygons,
    polygons_colours = polygons_colours,
    min_distance_km = min_distance_km,
    intensity_threshold = intensity_threshold,
    polygon_alpha = polygon_alpha  # Return alpha value for reference
  ))
}
spatial_data <- process_spatial_data(mcc_tabdf, shp_dsn, shp_layer, raster_file,
                                     intensity_threshold = intensity_threshold)
### Section 3: Plotting Functions
plot_phylogeography <- function(spatial_data, 
                                curvedarrow_lwd = 1.1,
                                arrowhead_length = 0.08,
                                arrowhead_width = 0.06,
                                point_cex = 0.35,
                                intensity_threshold = NULL,
                                solid_lwd_multiplier = 2.0,
                                dashed_lwd_multiplier = 1.0,
                                dashed_lty = 2,
                                show_dashed_lines = FALSE,
                                use_curvedarrow = TRUE,
                                min_distance_km = NULL,
                                pitjit = 0.01,
                                plot_hpd_polygons = TRUE,
                                save_plot = TRUE,
                                filename = "output/rabv_phylogeo_cRRWplot01.png",
                                width = 16,
                                height = 14,
                                dpi = 300,
                                bg = "white") {  # New parameter to control HPD polygon plotting
  # Extract data from spatial_data object
  mcc_filtered <- spatial_data$mcc_filtered
  borders <- spatial_data$borders
  template_raster <- spatial_data$template_raster
  colour_scale <- spatial_data$colour_scale
  polygons <- spatial_data$polygons
  polygons_colours <- spatial_data$polygons_colours
  
  # Use provided parameters or defaults from spatial_data
  if (is.null(min_distance_km)) {
    min_distance_km <- spatial_data$min_distance_km
  }
  
  if (is.null(intensity_threshold)) {
    intensity_threshold <- spatial_data$intensity_threshold
  }
  
  # Set line widths based on intensity
  if (!"line_type" %in% names(mcc_filtered)) {
    mcc_filtered$line_type <- ifelse(mcc_filtered$dispersal_intensity >= intensity_threshold, 
                                     "solid", "dashed")
  }
  
  mcc_filtered$line_width <- ifelse(mcc_filtered$line_type == "solid",
                                    curvedarrow_lwd * solid_lwd_multiplier,
                                    curvedarrow_lwd * dashed_lwd_multiplier)
  
  # Plotting data
  plot_data <- mcc_filtered
  
  message(sprintf("Plotting %d dispersal events (filtered by distance ≥ %s km)", 
                  nrow(plot_data), min_distance_km))
  
  ## ---- Main plot setup ----
  xrange <- c(raster::xmin(template_raster), raster::xmax(template_raster))
  yrange <- c(raster::ymin(template_raster), raster::ymax(template_raster))
  
  par(mar = c(0, 0, 0, 0), oma = c(1.2, 2.4, 1, 0), 
      mgp = c(0, 0.2, 0), lwd = 0.1, bty = "o")
  
  raster::plot(template_raster, col = "white", box = FALSE, 
               axes = FALSE, colNA = "grey90", legend = FALSE,
               xlim = xrange, ylim = yrange)
  
  plot(borders, add = TRUE, lwd = 0.1, border = "gray10")
  ## ---- Plot HPD polygons ----
  convert_hex_to_color <- function(hex_color) {
    if (is.na(hex_color) || is.null(hex_color)) {
      return(adjustcolor("#808080", alpha.f = 0.33))
    }
    if (nchar(hex_color) == 9 && substr(hex_color, 1, 1) == "#") {
      rgb_hex <- substr(hex_color, 1, 7)  # #RRGGBB
      alpha_hex <- substr(hex_color, 8, 9)  # AA
      alpha_decimal <- strtoi(alpha_hex, base = 16)
      return(adjustcolor(rgb_hex, alpha.f = alpha_decimal / 255))
    } else if (nchar(hex_color) == 7 && substr(hex_color, 1, 1) == "#") {
      return(adjustcolor(hex_color, alpha.f = 0.33))
    } else {
      return(adjustcolor(hex_color, alpha.f = 0.33))
    }
  }
  # Use the helper function in the plotting code:
  ## ---- Plot HPD polygons with multiple visualization methods ----
  if (plot_hpd_polygons && !is.null(polygons) && length(polygons) > 0) {
    plot_extent <- c(xrange[1], xrange[2], yrange[1], yrange[2])
    message(sprintf("Plot extent: x=[%.4f, %.4f], y=[%.4f, %.4f]", 
                    plot_extent[1], plot_extent[2], plot_extent[3], plot_extent[4]))
    # Plot polygons in reverse order
    for (i in length(polygons):1) {
      tryCatch({
        poly <- polygons[[i]]
        
        if (!is.null(poly) && (inherits(poly, "SpatialPolygonsDataFrame") || 
                               inherits(poly, "SpatialPolygons"))) {
          
          # Get polygon extent for debugging
          poly_extent <- raster::extent(poly)
          message(sprintf("Polygon %d extent: x=[%.4f, %.4f], y=[%.4f, %.4f]", 
                          i, poly_extent@xmin, poly_extent@xmax, 
                          poly_extent@ymin, poly_extent@ymax))
          
          # Check if polygon is within plot extent
          if (poly_extent@xmin > plot_extent[2] || poly_extent@xmax < plot_extent[1] ||
              poly_extent@ymin > plot_extent[4] || poly_extent@ymax < plot_extent[3]) {
            message(sprintf("Polygon %d is outside plot extent, skipping", i))
            next
          }
          
          # Get polygon color with proper opacity
          if (length(polygons_colours) >= i && !is.na(polygons_colours[i])) {
            poly_color <- polygons_colours[i]
          } else {
            # Use a more visible color for debugging
            poly_color <- adjustcolor("steelblue", alpha.f = 0.6)
          }
          
          # METHOD 1: Direct polygon extraction and plotting (most reliable)
          message(sprintf("Plotting polygon %d using direct polygon() method", i))
          
          # Extract coordinates from SpatialPolygonsDataFrame
          if (inherits(poly, "SpatialPolygonsDataFrame")) {
            # Get the first polygon from the SpatialPolygons
            sp_polygons <- poly@polygons
            
            for (poly_idx in seq_along(sp_polygons)) {
              # Get Polygons object
              polygon_obj <- sp_polygons[[poly_idx]]
              
              # Get all parts of the polygon (including holes)
              for (part_idx in seq_along(polygon_obj@Polygons)) {
                poly_part <- polygon_obj@Polygons[[part_idx]]
                
                # Get coordinates
                coords <- poly_part@coords
                
                # Check if this is a hole
                if (poly_part@hole) {
                  # Plot hole with different color for debugging
                  polygon(coords[,1], coords[,2], 
                          col = "white",  # Fill holes with white
                          border = NA, 
                          lwd = 0)
                } else {
                  # Plot main polygon with color
                  polygon(coords[,1], coords[,2], 
                          col = poly_color, 
                          border = adjustcolor("gray30", alpha.f = 0.5), 
                          lwd = 0.2)
                }
              }
            }
          } else if (inherits(poly, "SpatialPolygons")) {
            # For SpatialPolygons objects
            for (poly_idx in seq_along(poly@polygons)) {
              polygon_obj <- poly@polygons[[poly_idx]]
              
              for (part_idx in seq_along(polygon_obj@Polygons)) {
                poly_part <- polygon_obj@Polygons[[part_idx]]
                coords <- poly_part@coords
                
                if (!poly_part@hole) {
                  polygon(coords[,1], coords[,2], 
                          col = poly_color, 
                          border = adjustcolor("gray30", alpha.f = 0.5), 
                          lwd = 0.2)
                }
              }
            }
          }
        }
      }, error = function(e) {
        message(sprintf("Error processing polygon %d: %s", i, e$message))
      })
    }
  }
  ## Keep the draw_dashed_curve helper function
  draw_dashed_curve <- function(xs, ys, xe, ye, lwd, col, lty = 1, 
                                arrow_length = 0.08, arrow_width = 0.06,
                                curve = 0.3, n_points = 100) {
    # Calculate control point for Bezier curve
    mid_x <- (xs + xe) / 2
    mid_y <- (ys + ye) / 2
    
    # Perpendicular offset for curvature
    dx <- xe - xs
    dy <- ye - ys
    perp_x <- -dy * curve
    perp_y <- dx * curve
    
    ctrl_x <- mid_x + perp_x
    ctrl_y <- mid_y + perp_y
    
    # Generate Bezier curve points
    t <- seq(0, 1, length.out = n_points)
    curve_x <- (1-t)^2 * xs + 2*(1-t)*t * ctrl_x + t^2 * xe
    curve_y <- (1-t)^2 * ys + 2*(1-t)*t * ctrl_y + t^2 * ye
    
    # Draw the line with proper dashing
    lines(curve_x, curve_y, 
          col = col, 
          lwd = lwd, 
          lty = lty,
          type = "l")
    
    # Add arrowhead (solid, not dashed)
    arrow_idx <- floor(n_points * 0.85)  # Arrow at 85% along the curve
    
    if (arrow_idx < n_points - 1) {
      dx_arrow <- curve_x[arrow_idx + 1] - curve_x[arrow_idx]
      dy_arrow <- curve_y[arrow_idx + 1] - curve_y[arrow_idx]
    } else {
      dx_arrow <- curve_x[arrow_idx] - curve_x[arrow_idx - 1]
      dy_arrow <- curve_y[arrow_idx] - curve_y[arrow_idx - 1]
    }
    
    # Normalize and scale arrow
    arrow_norm <- sqrt(dx_arrow^2 + dy_arrow^2)
    if (arrow_norm > 0) {
      dx_arrow <- dx_arrow / arrow_norm * arrow_length * 2
      dy_arrow <- dy_arrow / arrow_norm * arrow_length * 2
    }
    
    # Draw solid arrowhead
    arrows(curve_x[arrow_idx], curve_y[arrow_idx],
           curve_x[arrow_idx] + dx_arrow, curve_y[arrow_idx] + dy_arrow,
           length = arrow_length,
           angle = 20,
           col = col,
           lwd = lwd,
           lty = 1)  # Arrowhead always solid
    
    return(list(x = curve_x, y = curve_y))
  }
  
  ## ---- Plot dispersal arrows ----
  if (use_curvedarrow && requireNamespace("diagram", quietly = TRUE)) {
    message("Using curvedarrow for dispersal lines with arrowheads")
    
    # Separate solid and dashed lines
    solid_data <- plot_data[plot_data$line_type == "solid", ]
    dashed_data <- plot_data[plot_data$line_type == "dashed", ]
    
    # Helper function to calculate midpoint coordinates for curved arrows
    calculate_curved_midpoint <- function(xs, ys, xe, ye, curve = 0.3) {
      # Calculate control point for Bezier curve (same as diagram::curvedarrow uses)
      mid_x <- (xs + xe) / 2
      mid_y <- (ys + ye) / 2
      dx <- xe - xs
      dy <- ye - ys
      
      # Perpendicular offset for curvature
      perp_x <- -dy * curve
      perp_y <- dx * curve
      
      ctrl_x <- mid_x + perp_x
      ctrl_y <- mid_y + perp_y
      
      # Calculate midpoint at t=0.5 on Bezier curve
      t <- 0.5
      mid_curve_x <- (1-t)^2 * xs + 2*(1-t)*t * ctrl_x + t^2 * xe
      mid_curve_y <- (1-t)^2 * ys + 2*(1-t)*t * ctrl_y + t^2 * ye
      
      # Calculate tangent vector at midpoint for arrow direction
      t_dx <- 2*(1-t)*(ctrl_x - xs) + 2*t*(xe - ctrl_x)
      t_dy <- 2*(1-t)*(ctrl_y - ys) + 2*t*(ye - ctrl_y)
      
      # Normalize tangent vector
      length <- sqrt(t_dx^2 + t_dy^2)
      if (length > 0) {
        t_dx <- t_dx / length
        t_dy <- t_dy / length
      }
      
      return(list(x = mid_curve_x, y = mid_curve_y, 
                  dx = t_dx, dy = t_dy, 
                  angle = atan2(t_dy, t_dx) * 180 / pi))
    }
    
    # Draw solid lines with arrowheads
    if (nrow(solid_data) > 0) {
      # First draw all the curved lines
      for (i in 1:nrow(solid_data)) {
        diagram::curvedarrow(
          cbind(solid_data$startLon[i], solid_data$startLat[i]),
          cbind(solid_data$endLon[i], solid_data$endLat[i]),
          arr.length = 0,  # No arrow at end
          arr.width = 0.00,
          lwd = solid_data$line_width[i],
          lty = 1,
          lcol = solid_data$year_color[i],
          arr.col = NA,
          arr.pos = FALSE,
          curve = 0.3,
          dr = NA,
          endhead = FALSE
        )
      }
    }
    
    # Draw dashed lines with arrowheads (if enabled)
    if (show_dashed_lines && nrow(dashed_data) > 0) {
      # First draw all the dashed curved lines
      for (i in 1:nrow(dashed_data)) {
        diagram::curvedarrow(
          cbind(dashed_data$startLon[i], dashed_data$startLat[i]),
          cbind(dashed_data$endLon[i], dashed_data$endLat[i]),
          arr.length = 0,
          arr.width = 0.00,
          lwd = dashed_data$line_width[i],
          lty = dashed_lty,
          lcol = dashed_data$year_color[i],
          arr.col = NA,
          arr.pos = FALSE,
          curve = 0.3,
          dr = NA,
          endhead = FALSE
        )
      }
      
      # Then add arrowheads at midpoints
#      for (i in 1:nrow(dashed_data)) {
#        mid_point <- calculate_curved_midpoint(
#          dashed_data$startLon[i], dashed_data$startLat[i],
#          dashed_data$endLon[i], dashed_data$endLat[i],
#          curve = 0.3
#        )
#        
#        # Draw arrowhead at midpoint (solid arrowhead on dashed line)
#        arrows(
#          x0 = mid_point$x - mid_point$dx * arrowhead_length * 0.5,
#          y0 = mid_point$y - mid_point$dy * arrowhead_length * 0.5,
#          x1 = mid_point$x + mid_point$dx * arrowhead_length * 0.5,
#          y1 = mid_point$y + mid_point$dy * arrowhead_length * 0.5,
#          length = arrowhead_length,
#          angle = 20,
#          col = dashed_data$year_color[i],
#          lwd = dashed_data$line_width[i],
#          lty = 1  # Arrowhead always solid
#        )
#      }
    }
  } else {
    message("Using custom curved lines (diagram package not available)")
    
    # Separate solid and dashed lines
    solid_data <- plot_data[plot_data$line_type == "solid", ]
    dashed_data <- plot_data[plot_data$line_type == "dashed", ]
    
    # Draw solid lines
    if (nrow(solid_data) > 0) {
      for (i in 1:nrow(solid_data)) {
        curve_result <- draw_dashed_curve(solid_data$startLon[i], solid_data$startLat[i],
                                          solid_data$endLon[i], solid_data$endLat[i],
                                          lwd = solid_data$line_width[i],
                                          col = solid_data$year_color[i],
                                          lty = 1,
                                          arrow_length = arrowhead_length)
      }
    }
    
    # Draw dashed lines (if enabled)
    if (show_dashed_lines && nrow(dashed_data) > 0) {
      for (i in 1:nrow(dashed_data)) {
        curve_result <- draw_dashed_curve(dashed_data$startLon[i], dashed_data$startLat[i],
                                          dashed_data$endLon[i], dashed_data$endLat[i],
                                          lwd = dashed_data$line_width[i],
                                          col = dashed_data$year_color[i],
                                          lty = dashed_lty,
                                          arrow_length = arrowhead_length)
      }
    }
  }
  ## ---- Plot points ----
  ptsize <- point_cex
  
  for (i in nrow(plot_data):1) {
    # Add jitter to end points
    xs <- plot_data$startLon[i]
    ys <- plot_data$startLat[i]
    xe <- jitter(plot_data$endLon[i], pitjit)
    ye <- jitter(plot_data$endLat[i], pitjit)
    
    # Plot start point for first event only
    if (i == nrow(plot_data)) {
      points(xs, ys, pch = 16, col = colour_scale[1], cex = ptsize)
      points(xs, ys, pch = 1, col = "gray10", cex = ptsize)
    }
    
    # Plot end points
    points(xe, ye, pch = 16, col = plot_data$year_color[i], cex = ptsize)
    points(xe, ye, pch = 1, col = "gray10", cex = ptsize)
  }
  
  ## ---- Map boundaries ----
  rect(xrange[1], yrange[1], xrange[2], yrange[2], xpd = TRUE, lwd = 0.2)
  
  axis(1, c(ceiling(xrange[1]), floor(xrange[2])), 
       pos = yrange[1], mgp = c(0, 0.2, 0), cex.axis = 0.5, 
       lwd = 0, lwd.tick = 0.2, padj = -0.8, tck = -0.01, col.axis = "gray30")
  
  axis(2, c(ceiling(yrange[1]), floor(yrange[2])), 
       pos = xrange[1], mgp = c(0, 0.5, 0), cex.axis = 0.5, 
       lwd = 0, lwd.tick = 0.2, padj = 1, tck = -0.01, col.axis = "gray30")
  
  ## ---- Line type legend ----
  legend("topright", 
         legend = c(paste("Intensity ≥", intensity_threshold), 
                    paste("Intensity <", intensity_threshold)),
         lty = c(1, if(show_dashed_lines) dashed_lty else NULL),
         lwd = c(curvedarrow_lwd * solid_lwd_multiplier * 0.7,
                 if(show_dashed_lines) curvedarrow_lwd * dashed_lwd_multiplier * 0.7 else NULL),
         col = c("black", "black"),
         cex = 0.5,
         bty = "n", 
         title = "Dispersal Intensity",
         title.cex = 0.6,
         xjust = 0.5)
  
  ## ---- Discrete horizontal year legend ----
  rast <- raster::raster(matrix(nrow = 1, ncol = 2))
  rast[1] <- floor(spatial_data$minYear)
  rast[2] <- ceiling(spatial_data$maxYear)
  
  # Position legend at bottom
  raster::plot(rast, legend.only = TRUE, add = TRUE, col = colour_scale, 
               legend.width = 0.5,
               legend.shrink = 0.3,
               smallplot = c(0.10, 0.40, 0.6, 0.62),
               legend.args = list(text = "", cex = 0.1, line = 0.3, col = "gray30"),
               horizontal = TRUE,
               axis.args = list(cex.axis = 0.6,
                                lwd = 0, 
                                lwd.tick = 0.2,
                                tck = -0.5,
                                col.axis = "gray30", 
                                line = 0, 
                                mgp = c(0, -0.08, 0)))
  text(x = mean(c(0.10, 0.40)), y = 0.65, "Year", cex = 0.6, col = "gray30")
  ## ---- Add filter information text ----
  if (min_distance_km > 0) {
    filter_text <- sprintf("Distance filter: ≥%s km", min_distance_km)
    text(xrange[1] + 0.02 * diff(xrange), yrange[2] - 0.02 * diff(yrange),
         filter_text, cex = 0.5, adj = 0, col = "gray30")
  }
  # ========== SAVING FILE ==========
  if (save_plot) {
    # Close any stray devices first
    while (dev.cur() > 1) {
      dev.off()
    }
    
    # Ensure filename has .png extension
    if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".png")
    }
    
    # Create directory if it doesn't exist
    output_dir <- dirname(filename)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    message(sprintf("\nSaving plot to: %s (%.0fx%.0f inches, %d dpi)", 
                    filename, width, height, dpi))
    
    # Save the plot - FIXED APPROACH
    # Option 1: Re-plot to PNG device (Recommended)
    png(filename = filename, 
        width = width, 
        height = height, 
        units = "in", 
        res = dpi, 
        bg = bg)
    
    # Re-create the plot with all parameters
    # We need to call plot_phylogeography recursively but without saving
    plot_phylogeography_internal <- function(...) {
      # Call the function again but with save_plot = FALSE
      # This is a recursive call but with different parameters
      args <- list(...)
      args$save_plot <- FALSE
      do.call(plot_phylogeography, args)
    }
    
    # Call internal plotting function
    plot_phylogeography_internal(
      spatial_data = spatial_data,
      curvedarrow_lwd = curvedarrow_lwd,
      arrowhead_length = arrowhead_length,
      arrowhead_width = arrowhead_width,
      point_cex = point_cex,
      intensity_threshold = intensity_threshold,
      solid_lwd_multiplier = solid_lwd_multiplier,
      dashed_lwd_multiplier = dashed_lwd_multiplier,
      dashed_lty = dashed_lty,
      show_dashed_lines = show_dashed_lines,
      use_curvedarrow = use_curvedarrow,
      min_distance_km = min_distance_km,
      pitjit = pitjit,
      plot_hpd_polygons = plot_hpd_polygons
    )
    
    # Close device
    dev.off()
    
    message(sprintf("Plot saved successfully!"))
  }
  # =====================================
}
plot_phylogeography(spatial_data)
### Section 4a: Temporal Analysis and Plots: Calculate temporal metrics
calculate_temporal_metrics <- function(mcc_tabdf, origin = NULL) {
  # Load required package
  if (!requireNamespace("geosphere", quietly = TRUE)) stop("Package 'geosphere' required")
  
  # If origin not provided, use earliest event
  if (is.null(origin)) {
    origin <- mcc_tabdf |> 
      dplyr::slice_min(startYear, n = 1) |>
      dplyr::slice(1)
  }
  
  message(sprintf("Using origin at (%.4f, %.4f) from year %.1f", 
                  origin$startLon, origin$startLat, origin$startYear))
  
  # Calculate distances from origin (Haversine is fine, but distGeo is more accurate)
  mcc_tabdf$dist_from_origin_km <- geosphere::distGeo(
    cbind(origin$startLon, origin$startLat),
    cbind(mcc_tabdf$endLon, mcc_tabdf$endLat)
  ) / 1000
  
  # Also calculate total dispersal distance if not already present
  if (!"dispersal_distance_km" %in% names(mcc_tabdf)) {
    mcc_tabdf$dispersal_distance_km <- geosphere::distGeo(
      cbind(mcc_tabdf$startLon, mcc_tabdf$startLat),
      cbind(mcc_tabdf$endLon, mcc_tabdf$endLat)
    ) / 1000
  }
  
  return(list(
    mcc_tabdf = mcc_tabdf,
    origin = origin
  ))
}

### Section 4b: Temporal Analysis and Plots: Calculate monthly statistics
calculate_monthly_stats <- function(mcc_tabdf, origin, bandwidth_months = 2) {
  
  # Check for required columns
  required_cols <- c("endYear", "dist_from_origin_km")
  missing_cols <- setdiff(required_cols, names(mcc_tabdf))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s. Run calculate_temporal_metrics() first.", 
                 paste(missing_cols, collapse = ", ")))
  }
  
  # Get origin start year
  origin_start_year <- if (is.list(origin) && "startYear" %in% names(origin)) {
    origin$startYear
  } else if (is.numeric(origin)) {
    origin
  } else if ("startYear" %in% names(mcc_tabdf)) {
    min(mcc_tabdf$startYear, na.rm = TRUE)
  } else {
    warning("Cannot determine origin year. Using 0 as default.")
    0
  }
  
  message(sprintf("Using origin year: %.2f", origin_start_year))
  
  # Use existing monthly calculations if available, otherwise calculate
  if (!"monthly_index" %in% names(mcc_tabdf)) {
    message("Calculating monthly indices...")
    mcc_tabdf <- mcc_tabdf |>
      dplyr::filter(!is.na(endYear), !is.na(dist_from_origin_km)) |>
      dplyr::mutate(
        year = floor(endYear),
        frac_year = endYear - year,
        month_approx = floor(frac_year * 12) + 1,
        month_approx = pmin(pmax(month_approx, 1), 12),
        monthly_index = year + (month_approx - 1)/12,
        date_label = sprintf("%04d-%02d", year, month_approx)
      )
  } else {
    # Use existing monthly indices
    mcc_tabdf <- mcc_tabdf |>
      dplyr::filter(!is.na(endYear), !is.na(dist_from_origin_km), !is.na(monthly_index))
  }
  
  if (nrow(mcc_tabdf) == 0) {
    warning("No valid data after filtering.")
    return(data.frame())
  }
  
  # Generate monthly sequence
  min_idx <- min(mcc_tabdf$monthly_index, na.rm = TRUE)
  max_idx <- max(mcc_tabdf$monthly_index, na.rm = TRUE)
  monthly_seq <- seq(
    from = floor(min_idx * 12) / 12,
    to = ceiling(max_idx * 12) / 12,
    by = 1/12
  )
  
  bandwidth_years <- bandwidth_months / 12
  
  # Process each month
  results_list <- list()
  
  for (i in seq_along(monthly_seq)) {
    month_time <- monthly_seq[i]
    
    # Filter data within bandwidth
    window_data <- mcc_tabdf |>
      dplyr::filter(
        monthly_index >= month_time - bandwidth_years,
        monthly_index <= month_time + bandwidth_years
      )
    
    if (nrow(window_data) < 2) next
    
    dist_vec <- window_data$dist_from_origin_km
    time_from_origin <- window_data$endYear - origin_start_year
    
    # Calculate basic statistics
    stats <- data.frame(
      month_time = month_time,
      year = floor(month_time),
      month = round((month_time - floor(month_time)) * 12) + 1,
      n_events = length(dist_vec),
      max_dist = max(dist_vec, na.rm = TRUE),
      min_dist = min(dist_vec, na.rm = TRUE),
      mean_dist = mean(dist_vec, na.rm = TRUE),
      median_dist = median(dist_vec, na.rm = TRUE),
      sd_dist = sd(dist_vec, na.rm = TRUE),
      mean_time_from_origin = mean(time_from_origin, na.rm = TRUE),
      median_time_from_origin = median(time_from_origin, na.rm = TRUE),
      mean_sq_disp = mean(dist_vec^2, na.rm = TRUE)
    )
    
    # Diffusion coefficient (only if time from origin > 0)
    if (stats$mean_time_from_origin > 0) {
      stats$diffusion_coeff <- stats$mean_sq_disp / (4 * stats$mean_time_from_origin)
      
      # Bootstrap confidence intervals if enough data
      if (stats$n_events > 2) {
        n_boot <- min(1000, stats$n_events * 10)
        boot_coeffs <- numeric(n_boot)
        
        for (b in 1:n_boot) {
          bs_sample <- sample(dist_vec, replace = TRUE)
          mean_sq <- mean(bs_sample^2, na.rm = TRUE)
          boot_coeffs[b] <- mean_sq / (4 * stats$mean_time_from_origin)
        }
        
        stats$diffusion_lower <- quantile(boot_coeffs, 0.025, na.rm = TRUE)
        stats$diffusion_upper <- quantile(boot_coeffs, 0.975, na.rm = TRUE)
        stats$diffusion_se <- sd(boot_coeffs, na.rm = TRUE)
      } else {
        stats$diffusion_lower <- NA
        stats$diffusion_upper <- NA
        stats$diffusion_se <- NA
      }
    } else {
      stats$diffusion_coeff <- NA
      stats$diffusion_lower <- NA
      stats$diffusion_upper <- NA
      stats$diffusion_se <- NA
    }
    
    results_list[[i]] <- stats
  }
  
  # Combine results
  if (length(results_list) == 0) {
    message("No monthly windows with sufficient data found.")
    return(data.frame())
  }
  
  results <- dplyr::bind_rows(results_list) |>
    dplyr::filter(
      n_events >= 2,
      !is.na(mean_dist)
    ) |>
    dplyr::mutate(
      date_label = sprintf("%04d-%02d", year, month),
      date_plot = as.Date(paste(year, month, "15", sep = "-")),
      time_label = format(date_plot, "%b %Y")
    ) |>
    dplyr::arrange(month_time)
  
  # Summary message
  message(sprintf("Calculated statistics for %d monthly time points", nrow(results)))
  message(sprintf("Time range: %s to %s", 
                  results$time_label[1], 
                  results$time_label[nrow(results)]))
  
  return(results)
}

### Section 4c: Temporal Analysis and Plots: Plot monthly statistics
plot_temporal_analysis <- function(monthly_stats, distance_ribbon_alpha = 0.25) {
  library(ggplot2)
  library(scales)
  
  if (nrow(monthly_stats) == 0) {
    warning("Insufficient data for temporal analysis")
    return(list(
      wavefront_plot = NULL,
      diffusion_plot = NULL
    ))
  }
  # Ensure date_plot exists
  if (!"date_plot" %in% names(monthly_stats)) {
    monthly_stats <- monthly_stats |>
      dplyr::mutate(
        date_plot = as.Date(paste(year, month, "15", sep = "-"))
      )
  }
  # ---- Plot 1: Monthly smoothed maximal wavefront distance ----
  wavefront_plot <- ggplot(monthly_stats, aes(x = date_plot)) + theme_bw()+
    # Uncertainty ribbons
    geom_ribbon(aes(ymin = min_dist, ymax = max_dist),
                fill = "steelblue3", alpha = distance_ribbon_alpha * 0.6) +
    geom_ribbon(aes(ymin = mean_dist - sd_dist,
                    ymax = mean_dist + sd_dist),
                fill = "steelblue2", alpha = distance_ribbon_alpha) +
    # Central tendency lines
    geom_line(aes(y = mean_dist),
              linewidth = 1, color = "steelblue4") +
    geom_line(aes(y = median_dist),
              linewidth = 0.8, color = "navy", linetype = "dashed") +
    # Labels
    labs(
      x = "Time",
      y = "Distance from Origin (km)",
      title = "Wavefront Distance Over Time"
    ) +
    # Scales
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%Y",
                 expand = expansion(mult = 0.05)) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = 0.1)) +
    
    # Theme (NO GRID LINES)
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 10),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # ---- Plot 2: Monthly smoothed diffusion coefficient ----
  # Filter for valid diffusion coefficients
  diffusion_monthly <- monthly_stats |>
    dplyr::filter(!is.na(diffusion_coeff), diffusion_coeff > 0)
  
  if (nrow(diffusion_monthly) > 3) {
    diffusion_plot <- ggplot(diffusion_monthly, aes(x = date_plot)) + theme_bw()+
      # Confidence ribbon (bootstrap CI)
      geom_ribbon(
        aes(ymin = diffusion_lower, ymax = diffusion_upper),
        alpha = distance_ribbon_alpha,
        fill = "darkorange2"
      ) +
      # Diffusion coefficient line
      geom_line(
        aes(y = diffusion_coeff),
        linewidth = 1.2,
        color = "darkorange4"
      ) +
      # Add point markers
#      geom_point(
#        aes(y = diffusion_coeff, size = n_events),
#        color = "darkorange4",
#        alpha = 0.5,
#        show.legend = TRUE
#      ) +
      labs(
        x = "Time",
        y = expression(paste("Diffusion Coefficient (km"^"2", "/year)")),
        title = "Spatial Spread Rate Over Time",
        size = "Events\nper month"
      ) +
      scale_x_date(
        date_breaks = "1 year",
        date_labels = "%Y",
        expand = expansion(mult = 0.05)
      ) +
      scale_y_continuous(
        expand = expansion(mult = 0.1)
      ) +
      scale_size_continuous(
        range = c(1, 6),
        breaks = pretty(diffusion_monthly$n_events, n = 5)
      ) +
      theme(      panel.grid = element_blank(),
                  plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                  axis.title = element_text(face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  axis.text = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 15, 15, 15)
      )
    
    # Add log scale if values span multiple orders of magnitude
    if (max(diffusion_monthly$diffusion_coeff, na.rm = TRUE) / 
        min(diffusion_monthly$diffusion_coeff[diffusion_monthly$diffusion_coeff > 0], na.rm = TRUE) > 100) {
      diffusion_plot <- diffusion_plot + 
        scale_y_log10() +
        labs(y = expression(paste("Diffusion Coefficient (km"^"2", "/year, log scale)")))
    }
    
  } else {
    diffusion_plot <- NULL
    message("Insufficient data for diffusion coefficient calculation (need > 3 time points).")
  }
  
  return(list(
    wavefront_plot = wavefront_plot,
    diffusion_plot = diffusion_plot
  ))
}

# Usage pipeline
temp_metrics <- calculate_temporal_metrics(spatial_data$mcc_tabdf)
monthly_stats <- calculate_monthly_stats(temp_metrics$mcc_tabdf, temp_metrics$origin)
temporal_plots <- plot_temporal_analysis(monthly_stats)

# Display plots
if (!is.null(temporal_plots$wavefront_plot)) {
  print(temporal_plots$wavefront_plot)
}

if (!is.null(temporal_plots$diffusion_plot)) {
  print(temporal_plots$diffusion_plot)
}

# Save plots
if (!is.null(temporal_plots$wavefront_plot)) {
  ggsave("output/rabv_wavefront_distance_plot.png", temporal_plots$wavefront_plot, 
         width = 10, height = 6, dpi = 300)
}

if (!is.null(temporal_plots$diffusion_plot)) {
  ggsave("output/rabv_diffusion_coefficient_plot.png", temporal_plots$diffusion_plot, 
         width = 10, height = 6, dpi = 300)
}


# Access monthly data
if (!is.null(result$monthly_stats)) {
  print(head(result$monthly_stats))
  
  # Create custom monthly plot
  library(ggplot2)
  custom_monthly <- ggplot(result$monthly_stats, aes(x = date_plot)) +
    geom_line(aes(y = max_dist), color = "red") +
    geom_line(aes(y = median_dist), color = "blue") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    theme_minimal()
  
  ggsave("custom_monthly.pdf", custom_monthly, width = 12, height = 6)
}
############################### 6. RABV READ DEPTH QC ###############################
### 6A. MAIN SOURCE CODES: DEFINE PATHS TO FILES, SOURCES  ###
suppressPackageStartupMessages({
library(tidyr)
library(ggplot2)
require(data.table)
library(tidyverse)
library(patchwork)
library(wesanderson)
library(RColorBrewer)
library(ggridges)
rabv_depth <- read.table("01_src/readcount/RABV/ngs_all_depthstat.txt",
                         header=TRUE, sep=" ", na.strings="NA", strip.white=TRUE)
#rabv_read <- read_excel("data/rabv_epidstat.xlsx")
#rabv_stat <- as.data.frame(rabv_read)
#rabv_stat <- melt(rabv_stat, id="epi_week")
features1 <- tribble(~"feature", ~"start", ~"end",
                     "", 1, 1,
                     "G", 3296, 5363, "",11932,11932)
features2 <- tribble(~"feature", ~"start", ~"end",
                     "3'UTR", 1, 58,
                     "N", 59, 1483,
                     "", 1484, 1485,
                     "P", 1486, 2475,
                     "", 2476, 2480,
                     "M", 2481, 3290,
                     "", 3291, 5387,
                     "L", 5388, 11862,
                     "", 11863, 11865,
                     "5'UTR", 11866, 11932)
})
### 6B. MAIN ANALYSIS CODES: PLOTTING AND SAVING (READDEPTH)  ###
{
p1 <-ggplot() +
  geom_line(rabv_depth, mapping=aes(x=pos, y=median), size=0.9, col="#304b56", size=2,alpha=0.9)+
  geom_ribbon(rabv_depth, mapping =aes(ymin=lower_quartile, ymax=upper_quartile, x=pos), fill="#83ba94", alpha=0.75)+
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0.02),breaks = seq(0, 11950, by = 2000))+
  scale_y_continuous(trans = log2_trans(), breaks = c(5,100,1000,2000,5000))+
  theme(legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.direction="vertical", 
        legend.position = c(0.85,0.85), legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.1,'cm'),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(color="black", size=12,vjust=0.5),
        axis.title.y = element_text(color="black", size=18, face="bold"),
        axis.text.y = element_text(color="black", size=12))+
  labs(#title = "Daily new case incidence in Sarawak",
    x = "genome position",
    y = "per base read count")

p2 <- features1 %>%
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end, 
                ymin = 0, ymax = 1,
                fill = feature)) +
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"))+
  scale_x_continuous(expand = c(0.01,0.02),breaks = seq(0, 11950, by = 1000))+
  scale_y_continuous(expand = c(-0.01,-0.01))+
  geom_text(aes(x = (start + end) / 2, y = 0.1,vjust=-1,label = feature),size =3.8,fontface="bold") +
  theme_void() +
  theme(legend.position = "none")

p3 <- features2 %>%
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end, 
                ymin = 0, ymax = 1,
                fill = feature)) +
  scale_fill_brewer(palette = "Set3")+
  scale_x_continuous(expand = c(0.01,0.02),breaks = seq(0, 11950, by = 1000))+
  scale_y_continuous(expand = c(-0.01,-0.01))+
  geom_text(aes(x = (start + end) / 2, y = 0.1,vjust=-1.2,label = feature),hjust=0.7,size =3.8,fontface="bold") +
  theme_void() +
  theme(legend.position = "none")

seqdepth <- p1 / p2 / p3 + plot_layout(nrow = 3, heights = c(1, 0.08, 0.08))
#ggsave(seqdepth, filename = "rabv_seqdepth.png", width = 32, height = 17, units = "cm",device = "png")
all_paths <-
  list.files(path = "./01_src/readcount/RABV/",
             pattern = "*.readdepth.txt",
             full.names = TRUE)
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")
all_filenames <- all_paths %>%
  basename() %>%
  as.list()
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
all_result <- as.data.frame(rbindlist(all_lists, fill = F))
rabv_readcount<-all_result[c(1:3)]
rabv_readcount <- as.data.frame(apply(rabv_readcount, 2, function(x) gsub('.readdepth.txt', '',x)))
colnames(rabv_readcount) <- c('count','length','id')
rabv_readcount <- rabv_readcount[ , c("id", "length", "count")]
rabv_readcount <- rabv_readcount[grepl('187', rabv_readcount$id), ]
rabv_readcount <- rabv_readcount[!grepl('187B63', rabv_readcount$id), ]

#rabv_readcount$source <- ifelse(rabv_readcount$id %in% c("183B45", "183B46", "183B47", "183B48", "183B49", "183B50", "183B51", "183B52"), "primary", "primary")
rabv_readcount$length <- as.numeric(rabv_readcount$length)
rabv_readcount <- rabv_readcount %>% 
  filter((length > 900 & length < 1900)) 
#rabv_rcc_exp <- rabv_readcount[rep(row.names(rabv_readcount), rabv_readcount$count), 1:3]
#rabv_rcc_exp <-rabv_rcc_exp[-c(1)]
#rabv_rcc_exp <- as.data.frame(rabv_rcc_exp[1:19642,])
#rabv_rcc_exp$length <- as.numeric(rabv_rcc_exp$length)
rabv_rccplot <- ggplot(rabv_readcount, aes(x = length, y = id, fill=id)) +
  geom_density_ridges(scale = 3, size = 0.25, rel_min_height = 0.03, alpha=0.6) +
  theme_ridges() +
  scale_x_continuous(limits = c(800, 1800), breaks = seq(800, 1800,by = 200), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(margin=margin(r=50), hjust=0.5),
        axis.title.y = element_text(vjust=-1, hjust=0.5),
        #legend.text = element_text(size = 11, colour = "black"), 
        axis.text.y = element_blank())+xlab("Read length (bp)")+ylab("Read count density")+coord_cartesian(clip = "off")+ labs(fill="sequence id", size = 11)

pclab3<- ggarrange(seqdepth,labels = "A",font.label=list(color="black",size=16),nrow = 1)
pclab4<- ggarrange(rabv_rccplot,labels = "B",font.label=list(color="black",size=16),nrow = 1)
rabv_read_plots <- pclab3 / pclab4 + plot_layout(nrow = 1, ncol=2, heights = c(1, 0, 0))
ggsave(rabv_read_plots, filename = 'rabv_seqdepth_rccplot_2024.png', width = 47, height = 21, units = "cm",limitsize = FALSE,device = "png")
}
### 6C. MAIN ANALYSIS CODES: PLOTTING AND SAVING (READSTATS- PENDING MODIFY) ###
{
rabv_read <- read.table("data/rabv_all_readstats.txt",header=TRUE, sep=",", na.strings="NA", strip.white=TRUE)
rabv_read <- as.data.frame(rabv_read[-c(3:4)])
rabv_read <-rabv_read[order(rabv_read$total_reads, decreasing = FALSE), ]
rabv_read <-rabv_read[order(rabv_read$extr_type, decreasing = FALSE), ]
rabv_readstat<-rabv_read %>% 
  pivot_longer(cols = mapped_reads:total_reads,
               names_to = "read_type", 
               values_to = "count")
rabv_readstat <-rabv_readstat %>%
  mutate(source_no = case_when(
    endsWith(source, "isolate") ~ "2",
    endsWith(source, "primary") ~ "1"
  ))
rabv_readstat$source_no <- as.integer(rabv_readstat$source_no)
rabv_readplot<-ggplot(rabv_readstat, aes(x=reorder(id,extr_type,min), y=count, fill = read_type)) +
  geom_bar(stat='identity',position = position_stack(reverse = TRUE),color = "black")+theme_classic()+
  facet_wrap(~source,nrow=1,scales = "free_x")+
  scale_fill_manual(values = c("#6A51A3", "#7a74b3", "#a8a8d0","#7a74b3", "#a8a8d0"))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_sqrt(expand = c(0, 0),limits = c(0, 350000), breaks = seq(0, 350000,by = 15000))+
  theme(axis.text.y = element_text(size = 10, colour = "black", hjust = 0.2), 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(margin=margin(r=50), hjust=0.5),
        axis.title.y = element_text(vjust=0, hjust=0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, colour = "black"), 
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))+xlab("Metagenomic samples")+ylab("Read count")+labs(fill = " ")+
  coord_cartesian(clip = "off")+
  guides(fill = guide_legend(reverse = TRUE))
ymax <-120000
scaleRight <- ymax/100
rabv_mapplot<- ggplot(rabv_readstat) +
  geom_line(rabv_readstat, mapping = aes(x =reorder(id,extr_type,min), y = scaleRight*mapped_pct), linetype="dashed",stat='identity',color = "#be6171", size = 0.8) +
  geom_boxplot(aes(x=reorder(id,extr_type,min), y=mean_rdepth, fill=source), stat='identity',position = position_stack(reverse = TRUE),color = "black")+theme_classic()+
  facet_wrap(~source,nrow=1,scales = "free_x")+
  scale_fill_manual(values = c("#6A51A3", "#7a74b3", "#a8a8d0","#7a74b3", "#a8a8d0"))+
  scale_x_discrete(expand = c(0, 0))+
  #scale_y_sqrt(expand = c(0, 0),limits = c(0, ymax), breaks = seq(0, ymax,by = 10000), sec.axis = sec_axis(~ . /scaleRight, labels = scales::label_percent()))+
  scale_y_sqrt(expand = c(0, 0),limits = c(0, ymax), breaks = seq(0, ymax,by = 10000))+
  theme(axis.text.y = element_text(size = 10, colour = "black", hjust = 0.2), 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(margin=margin(r=50), hjust=0.5),
        axis.title.y = element_text(vjust=0, hjust=0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, colour = "black"), 
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))+xlab("Metagenomic samples")+ylab("Mean read depth")+labs(fill = " ")+
  coord_cartesian(clip = "off")+
  guides(fill = guide_legend(reverse = TRUE))
pclab3<- ggarrange(seqdepth,labels = "A",font.label=list(color="black",size=16),nrow = 1)
pclab4<- ggarrange(rabv_readplot,labels = "B",font.label=list(color="black",size=16),nrow = 1)
pclab5<- ggarrange(rabv_rccplot,labels = "c",font.label=list(color="black",size=16),nrow = 1)
pclab6<- ggarrange(rabv_mapplot,labels = "d",font.label=list(color="black",size=16),nrow = 1)

rabv_read_plots <- pclab3 / pclab4 / pclab5 / pclab6 + plot_layout(nrow = 2, ncol=2, heights = c(1, 1, 0))
ggsave(rabv_read_plots, filename = 'rabv_seqdepth_rdplot_rccplot_2022.png', width = 47, height = 31, units = "cm",limitsize = FALSE,device = "png")
}

######################## ALLPLEX PCR - FOR VIRUS TYPING ONLY- NOT RELATED TO RABV #############
rvp_meta <- read.table(file = "01_src/RVP_cohort_qpcr_2022.txt", sep = ',', header = TRUE)
#define the colours to use in the figure
rvp_vtype <- melt(rvp_meta, id = c("Epiweek"))
rvp_vtype$Epiweek <- factor(rvp_vtype$Epiweek,levels=unique(rvp_vtype$Epiweek))
#rvp_vtype <- as.data.frame(apply(rvp_vtype, 2, function(x) gsub('X', '',x)))

mx <- ggplot(rvp_vtype, aes(x = Epiweek, fill = variable, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Epidemic week in 2022", y = "Number of positive detections", fill = "Virus type") + 
  scale_fill_manual(values = type_cols)
ggsave('output/rvp_virtype_allplxpcr_2022.png',width = 23, height = 18, units = "cm",limitsize = FALSE, device = 'png')

################################### HADV phylogeny #########################################
hadv_df<- read.table(file = "01_src/swk_hadv_metadata.csv", sep = ',', header = TRUE)
hadv_df <- as.data.frame(hadv_df)
treeadv <- read.newick("02_bootstraps/hadv_swk_globred_aligned_p.fasta.treefile")

advp <- ggtree(treeadv, size=0.75, show.legend = TRUE) %<+% hadv_df + theme_tree2()+
  geom_tippoint(aes(color=genotype), size=3.5)+
  scale_color_manual(values=type_cols)+
  scale_x_continuous(breaks=c(0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70), labels=c("0.00", "0.10", "0.20","0.30","0.40","0.50","0.60","0.70"))+
  expand_limits(x = c(0, 0.9), y=c(0,540))+
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, hjust=-0.02, size = 2.8,fontface = "bold")+ 
  theme(axis.text = element_text(vjust=0.5,size=10, face = "bold", angle = 90), 
        legend.title = element_text(hjust = 0.5, vjust = -0.5, face = "bold"),
        legend.text = element_text(size = 10),
        legend.direction="vertical", 
        legend.position = "left", legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.1,'cm'),
        legend.margin = margin(2, -3, 10, 1, "cm"),
        plot.margin = unit(c(-5,0,0,3), "lines"))+  
  guides(fill=guide_legend(title="HAdv genotypes"))+
  coord_cartesian(xlim = c(-0.02, 0.95), # This focuses the x-axis on the range of interest
                  clip = 'off')
ggsave('output/global_fN_swkrabv_2022.png',width = 21, height = 30, units = "cm",limitsize = FALSE, device = 'png')

library("readxl")
library(readr)
library(dplyr)
library(tidyr)
library(writexl)
library(openxlsx)
nbsb_ac <- read_excel("01_src/book.xlsx", sheet = 1)


result <- nbsb_ac %>%
  # trim whitespace in invoices column
  mutate(invoices = trimws(invoices)) %>%
  
  # split by comma and expand rows
  separate_rows(invoices, sep = ",") %>%
  
  # trim again after split (important)
  mutate(invoices = trimws(invoices))

write_xlsx(result, "01_src/book_output.xlsx")

result2 <- nbsb_ac %>%
  # trim whitespace in invoices column
  mutate(invoices = trimws(invoices)) %>%
  # split by comma and expand rows
  separate_rows(invoices, sep = ",") %>%
  # trim again after split
  mutate(invoices = trimws(invoices)) %>%
  # keep amt only once per rv_id
  group_by(rv_id) %>%
  mutate(
    amt = if_else(row_number() == 1, amt, NA_real_)
  ) %>%
  ungroup()

write_xlsx(result2, "01_src/book_output2.xlsx")
.libPaths()

################################### AUDIT TRAIL FOR NBSB ################################### 
library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)

library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)

# ===============================
# Read paid_inv.xlsx
# ===============================
paid_iv <- read_excel("01_src/paid_inv.xlsx") %>%
  mutate(
    amount = as.numeric(amount),
    paid_in_YE2025 = as.numeric(paid_in_YE2025),
    remaining_amount = as.numeric(amount)
  )

# ===============================
# Read rv_paid.xlsx
# ===============================
rv_paid_raw <- read_excel("01_src/rv_paid.xlsx") %>%
  mutate(
    date = as.Date(date),
    RV_amt = as.numeric(RV_amt)
  )

# ===============================
# Parse INV_applied
# ===============================
rv_paid_parsed <- rv_paid_raw %>%
  mutate(
    INV = str_extract(INV_applied, "^[^\\(]+"),
    bracket_amt = str_extract(INV_applied, "\\(([^\\)]+)\\)") %>%
      str_remove_all("[()]") %>%
      as.numeric()
  )

# ===============================
# Treat each unique (RV, RV_amt) as separate payment
# ===============================
rv_master <- rv_paid_parsed %>%
  distinct(RV, RV_amt, date) %>%
  arrange(RV_amt, date) %>%
  mutate(RV_payment_id = row_number())

# ===============================
# Allocation Engine (STATEFUL)
# ===============================
audit_trail <- list()

# Initialize remaining RV balance tracker
remaining_rv_balance <- data.frame(
  RV = character(),
  RV_total_amt = numeric(),
  RV_applied_amt = numeric(),
  RV_remaining_balance = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(rv_master))) {
  
  rv <- rv_master$RV[i]
  rv_amt_total <- rv_master$RV_amt[i]
  rv_balance <- rv_amt_total
  rv_date <- rv_master$date[i]
  rv_pid <- rv_master$RV_payment_id[i]
  
  rv_lines <- rv_paid_parsed %>%
    filter(RV == rv, RV_amt == rv_amt_total)
  
  rv_applied_amt <- 0  # Initialize total applied amount for this RV
  
  for (k in seq_len(nrow(rv_lines))) {
    
    inv <- rv_lines$INV[k]
    bracket_limit <- rv_lines$bracket_amt[k]
    
    iv_idx <- which(
      paid_iv$INV == inv & 
        paid_iv$remaining_amount > 0
    )
    
    if (length(iv_idx) == 0) next
    
    for (idx in iv_idx) {
      
      if (rv_balance <= 0) break
      
      inv_before <- paid_iv$remaining_amount[idx]
      
      # Determine max deductible for this RV → INV
      if (!is.na(bracket_limit)) {
        max_deduct <- bracket_limit
      } else {
        max_deduct <- paid_iv$paid_in_YE2025[idx]
      }
      
      applied <- min(
        max_deduct,
        inv_before,
        rv_balance
      )
      
      if (applied <= 0) next
      
      audit_trail[[length(audit_trail) + 1]] <- data.frame(
        RV = rv,
        RV_payment_id = rv_pid,
        RV_date = rv_date,
        INV = inv,
        CLIENT_NAME = paid_iv$CLIENT_NAME[idx],
        CLIENT_ID = paid_iv$CLIENT_ID[idx],
        RV_original_amt = rv_amt_total,
        RV_applied_amt = applied,
        RV_balance_after = rv_balance - applied,
        INV_balance_before = inv_before,
        INV_balance_after = inv_before - applied,
        stringsAsFactors = FALSE
      )
      
      # Update balances (STATEFUL)
      rv_balance <- rv_balance - applied
      paid_iv$remaining_amount[idx] <- inv_before - applied
      
      # Track the amount applied for this RV
      rv_applied_amt <- rv_applied_amt + applied
    }
  }
  
  # Track the remaining RV balance after processing each payment
  remaining_rv_balance <- rbind(
    remaining_rv_balance,
    data.frame(
      RV = rv,
      RV_total_amt = rv_amt_total,
      RV_applied_amt = rv_applied_amt,
      RV_remaining_balance = rv_amt_total - rv_applied_amt,
      stringsAsFactors = FALSE
    )
  )
}

# ===============================
# Combine audit trail
# ===============================
audit_df <- bind_rows(audit_trail)

# ===============================
# Create Excel workbook
# ===============================
wb <- createWorkbook()
addWorksheet(wb, "RV_INV_Allocation_Audit")
addWorksheet(wb, "Remaining_Invoice_Balance")
addWorksheet(wb, "RV_Payment_Master")
addWorksheet(wb, "Remaining_RV_Balance")

writeData(wb, "RV_INV_Allocation_Audit", audit_df)
writeData(wb, "Remaining_Invoice_Balance", paid_iv)
writeData(wb, "RV_Payment_Master", rv_master)
writeData(wb, "Remaining_RV_Balance", remaining_rv_balance)

# ===============================
# Save audit file
# ===============================
saveWorkbook(wb, "01_src/RV_INV_Audit_Trail.xlsx", overwrite = TRUE)
