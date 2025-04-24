library(dplyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(forcats)
library(purrr)
library(furrr)
library(progressr)
plan(multisession)

# === PARAMETERS ===
input_summary_file <- "summary_stats.xlsx"
output_dir <- "olfactory_filtered_results2"
dir.create(output_dir, showWarnings = FALSE)

olfactory_roi_left <- c(41, 42, 43, 44, 45, 47, 53, 57, 58, 65, 66, 67, 70, 71, 118, 124, 125, 166)
olfactory_roi_right <- olfactory_roi_left + 1000
olfactory_roi_ids <- c(olfactory_roi_left, olfactory_roi_right)

all_metrics <- excel_sheets(input_summary_file)
metrics_of_interest <- c(
  "mrtrixfa_mean", "mrtrixad_mean", "mrtrixrd_mean", "mrtrixmd_mean",
  "msd_RAS_060622_mean", "ng_RAS_060622_mean", "qiv_RAS_060622_mean",
  "rtap_RAS_060622_mean", "rtop_RAS_060622_mean", "rtpp_RAS_060622_mean",
  "volume_percentage_mean"
)
metrics <- intersect(all_metrics, metrics_of_interest)

process_metric <- function(metric) {
  df <- read_xlsx(input_summary_file, sheet = metric)
  names(df) <- trimws(names(df))
  df <- df %>% rename(p.value = `Pr(>F)`)
  
  cat(glue::glue("\U0001F4C4 Loaded {metric}: {nrow(df)} total rows\n"))
  
  matching_rois <- df %>% filter(ROI %in% olfactory_roi_ids)
  cat(glue::glue("\U0001F9E0 {metric}: {nrow(matching_rois)} rows match olfactory ROIs\n"))
  
  df_filtered <- matching_rois %>%
    filter(!is.na(cohen_f)) %>%
    mutate(p_fdr = p.adjust(p.value, method = "fdr")) %>%
    filter(Factor != "Residuals")
  
  cat(glue::glue("\U0001F50D {metric}: {nrow(df_filtered)} rows after filtering for cohen_f\n"))
  
  if (nrow(df_filtered) == 0) {
    cat(glue::glue("⚠️ No valid olfactory ROI results for {metric}\n"))
    return(NULL)
  }
  
  walk(unique(df_filtered$Factor), function(factor) {
    df_factor <- df_filtered %>% filter(Factor == factor) %>% 
      arrange(desc(cohen_f)) %>% 
      mutate(label = paste0("ROI ", ROI, ": ", structure)) %>% 
      slice_max(cohen_f, n = 10) %>%
      mutate(across(c(f_CI_low, f_CI_high), ~ suppressWarnings(as.numeric(.)))) %>%
      filter(!is.na(f_CI_low), !is.na(f_CI_high)) %>%
      mutate(significance = ifelse(p_fdr < 0.05, "FDR < 0.05", "ns"))
    
    p <- ggplot(df_factor, aes(x = cohen_f, y = fct_reorder(label, cohen_f))) +
      geom_point(aes(color = significance), size = 3) +
      geom_errorbarh(aes(xmin = f_CI_low, xmax = f_CI_high), height = 0.3) +
      scale_color_manual(values = c("ns" = "gray40", "FDR < 0.05" = "firebrick")) +
      labs(
        title = paste("Top 10 for", factor, "in", metric),
        x = "Cohen's f (Effect Size)",
        y = NULL,
        color = "Significance"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.background = element_rect(fill = "white"))
    
    ggsave(filename = file.path(output_dir, paste0("Forest_", factor, "_", metric, ".png")),
           plot = p, width = 8, height = 6, dpi = 300)
  })
  
  return(df_filtered)
}

# === Run in parallel and collect results ===
handlers(global = TRUE)
sig_results_list <- with_progress({
  p <- progressor(along = metrics)
  future_map(metrics, function(metric) {
    p(sprintf("Processing %s", metric))
    process_metric(metric)
  })
})

# Assign names
names(sig_results_list) <- metrics

# === Save all valid non-null data frames into Excel ===
sig_results_wb <- createWorkbook()
for (i in seq_along(sig_results_list)) {
  df <- sig_results_list[[i]]
  if (is.null(df)) next
  sheet_name <- gsub("[^A-Za-z0-9_]", "_", names(sig_results_list)[i])
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(sig_results_wb, sheetName = sheet_name)
  writeData(sig_results_wb, sheet = sheet_name, df)
}
saveWorkbook(sig_results_wb, file = file.path(output_dir, "Olfactory_ROI_Stats.xlsx"), overwrite = TRUE)

cat("✅ All forest plots and Excel workbook saved to:", output_dir, "\n")



###Aggregates
top_by_factor <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  group_by(Factor) %>%
  slice_max(cohen_f, n = 5) %>%
  ungroup()

top_by_roi <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  group_by(ROI, structure) %>%
  slice_max(cohen_f, n = 5) %>%
  ungroup()

# Save as Excel workbook
agg_wb <- createWorkbook()
addWorksheet(agg_wb, "TopByFactor")
writeData(agg_wb, "TopByFactor", top_by_factor)

addWorksheet(agg_wb, "TopByROI")
writeData(agg_wb, "TopByROI", top_by_roi)

saveWorkbook(agg_wb, file = file.path(output_dir, "Olfactory_Aggregates.xlsx"), overwrite = TRUE)

# Save as TSVs
write_tsv(top_by_factor, file.path(output_dir, "TopByFactor.tsv"))
write_tsv(top_by_roi, file.path(output_dir, "TopByROI.tsv"))

####visualize top results****#######
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)

# Define output directory and resolution
output_dir <- "olfactory_filtered_results2"
dpi_val <- 300

# Load TSVs
top_by_factor <- read_tsv(file.path(output_dir, "TopByFactor.tsv"))
top_by_roi <- read_tsv(file.path(output_dir, "TopByROI.tsv"))



# Simplify factor names by extracting only the main effects
top_by_factor <- top_by_factor %>%
  mutate(Factor = gsub(":.*", "", Factor)) %>%
  filter(Factor %in% c("APOE", "Age_Months", "Sex", "Diet", "HN"))

top_by_roi <- read_tsv(file.path(output_dir, "TopByROI.tsv")) %>%
  mutate(Factor = gsub(":.*", "", Factor)) %>%
  filter(Factor %in% c("APOE", "Age_Months", "Sex", "Diet", "HN"))

# === Plot 1: Top ROIs by Factor ===
p1 <- ggplot(top_by_factor, aes(x = cohen_f, y = fct_reorder(paste0("ROI ", ROI, ": ", structure), cohen_f))) +
  geom_point(aes(color = Factor), size = 3) +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(
    title = "Top 5 ROIs by Factor (Cohen's f)",
    x = "Cohen's f (Effect Size)",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "Top_ROIs_by_Factor.png"), plot = p1, width = 12, height = 8, dpi = dpi_val)


# === Plot 2: Top Factors by ROI ===
p2 <- ggplot(top_by_roi, aes(x = cohen_f, y = fct_reorder(paste0("ROI ", ROI, ": ", structure), cohen_f))) +
  geom_point(aes(color = Factor), size = 3) +
  facet_wrap(~ paste("ROI", ROI), scales = "free_y") +
  labs(
    title = "Top Factors per ROI (Cohen's f)",
    x = "Cohen's f (Effect Size)",
    y = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "Top_Factors_per_ROI.png"), plot = p2, width = 12, height = 10, dpi = dpi_val)


# === Plot 3: Frequency of Most Affected ROIs ===
roi_freq_plot <- top_by_roi %>%
  count(structure, sort = TRUE) %>%
  top_n(10) %>%
  ggplot(aes(x = reorder(structure, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Most Frequently Affected ROIs", x = NULL, y = "Count") +
  theme_minimal()

ggsave(file.path(output_dir, "Most_Frequent_ROIs.png"), plot = roi_freq_plot, width = 8, height = 6, dpi = dpi_val)


# === Plot 4: Frequency of Most Influential Factors ===
factor_freq_plot <- top_by_factor %>%
  count(Factor, sort = TRUE) %>%
  ggplot(aes(x = reorder(Factor, n), y = n)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Most Influential Factors", x = NULL, y = "Count") +
  theme_minimal()

ggsave(file.path(output_dir, "Most_Influential_Factors.png"), plot = factor_freq_plot, width = 8, height = 6, dpi = dpi_val)



###end visualize top results
###
library(viridis)
library(tidyr)

# Prepare wide format for heatmap
heatmap_data <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  select(Metric, ROI, structure, cohen_f) %>%
  group_by(Metric, ROI, structure) %>%
  summarize(cohen_f = mean(cohen_f, na.rm = TRUE), .groups = "drop") %>%
  unite(ROI_label, ROI, structure, sep = ": ") %>%
  pivot_wider(names_from = Metric, values_from = cohen_f)

# Replace NA with 0 or leave as NA
heatmap_mat <- as.matrix(heatmap_data[,-1])
rownames(heatmap_mat) <- heatmap_data$ROI_label

# Plot heatmap
pheatmap::pheatmap(
  heatmap_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Cohen's f across Metrics and ROIs"
)

# Save the heatmap
ggsave(file.path(output_dir, "Heatmap_CohenF_ROI_Metric.png"), width = 10, height = 8, dpi = 300)

library(tidyr)
library(dplyr)
library(pheatmap)
library(viridis)

# Create a matrix of ROI vs Metric with average Cohen's f for significant results
roi_metric_matrix <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  select(Metric, ROI, structure, cohen_f) %>%
  group_by(ROI, structure, Metric) %>%
  summarize(mean_f = mean(cohen_f, na.rm = TRUE), .groups = "drop") %>%
  unite(ROI_label, ROI, structure, sep = ": ") %>%
  pivot_wider(names_from = Metric, values_from = mean_f)

# Convert to matrix
roi_mat <- as.matrix(roi_metric_matrix[,-1])
rownames(roi_mat) <- roi_metric_matrix$ROI_label
#Clustered heatmap
pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of ROIs by Cohen's f (Effect Size)",
  fontsize_row = 6,
  fontsize_col = 8
)

# Save figure
ggsave(file.path(output_dir, "Clustered_ROI_EffectSizes.png"), width = 12, height = 10, dpi = 300)

olfactory_labels <- grep("^\\d+: ", roi_metric_matrix$ROI_label, value = TRUE) %>%
  .[as.numeric(sub(":.*", "", .)) %in% olfactory_roi_ids]

olfactory_mat <- roi_mat[rownames(roi_mat) %in% olfactory_labels, ]

pheatmap::pheatmap(
  olfactory_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of Olfactory ROIs by Cohen's f",
  fontsize_row = 7,
  fontsize_col = 8
)

ggsave(file.path(output_dir, "Clustered_Olfactory_ROI_EffectSizes.png"), width = 10, height = 8, dpi = 300)



###using pretty metric names
pretty_metric_names <- c(
  "mrtrixrd_mean"        = "RD",
  "mrtrixad_mean"        = "AD",
  "mrtrixmd_mean"        = "MD",
  "mrtrixfa_mean"        = "FA",
  "volume_percentage_mean" = "Vol",
  "msd_RAS_060622_mean"  = "MSD",
  "qiv_RAS_060622_mean"  = "QIV",
  "rtpp_RAS_060622_mean" = "RTPP",
  "rtap_RAS_060622_mean" = "RTAP",
  "rtop_RAS_060622_mean" = "RTOP",
  "ng_RAS_060622_mean"   = "NG"
)
colnames(roi_mat) <- pretty_metric_names[colnames(roi_mat)]
colnames(olfactory_mat) <- pretty_metric_names[colnames(olfactory_mat)]

walk2(sig_results_list, names(sig_results_list), function(df_filtered, metric) {
  if (is.null(df_filtered)) return(NULL)
  
  walk(unique(df_filtered$Factor), function(factor) {
    df_factor <- df_filtered %>% filter(Factor == factor) %>% 
      arrange(desc(cohen_f)) %>% 
      mutate(label = paste0("ROI ", ROI, ": ", structure)) %>% 
      slice_max(cohen_f, n = 10) %>%
      mutate(across(c(f_CI_low, f_CI_high), ~ suppressWarnings(as.numeric(.)))) %>%
      filter(!is.na(f_CI_low), !is.na(f_CI_high)) %>%
      mutate(significance = ifelse(p_fdr < 0.05, "FDR < 0.05", "ns"))
    
    pretty_name <- pretty_metric_names[[metric]]
    
    p <- ggplot(df_factor, aes(x = cohen_f, y = fct_reorder(label, cohen_f))) +
      geom_point(aes(color = significance), size = 3) +
      geom_errorbarh(aes(xmin = f_CI_low, xmax = f_CI_high), height = 0.3) +
      scale_color_manual(values = c("ns" = "gray40", "FDR < 0.05" = "firebrick")) +
      labs(
        title = paste("Top 10 ROIs for", factor, "in", pretty_name),
        x = "Cohen's f (Effect Size)",
        y = NULL,
        color = "Significance"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.background = element_rect(fill = "white"))
    
    ggsave(filename = file.path(output_dir, paste0("Forest_", factor, "_", pretty_name, ".png")),
           plot = p, width = 8, height = 6, dpi = 300)
  })
})

# Matrix of ROI vs Metric
roi_metric_matrix <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  select(Metric, ROI, structure, cohen_f) %>%
  group_by(Metric, ROI, structure) %>%
  summarize(cohen_f = mean(cohen_f, na.rm = TRUE), .groups = "drop") %>%
  unite(ROI_label, ROI, structure, sep = ": ") %>%
  pivot_wider(names_from = Metric, values_from = cohen_f)

roi_mat <- as.matrix(roi_metric_matrix[,-1])
rownames(roi_mat) <- roi_metric_matrix$ROI_label
colnames(roi_mat) <- pretty_metric_names[colnames(roi_mat)]

# Heatmap
pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of ROIs by Cohen's f (Effect Size)",
  fontsize_row = 6,
  fontsize_col = 8
)

# Save it
ggsave(file.path(output_dir, "Clustered_ROI_EffectSizes_PrettyNames.png"), width = 12, height = 10, dpi = 300)


pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of ROIs by Cohen's f (Effect Size)",
  fontsize_row = 6,
  fontsize_col = 8,
  filename = file.path(output_dir, "Heatmap_Clustered_ROI_EffectSizes.png"),
  width = 12,
  height = 10
)

png(file.path(output_dir, "Heatmap_ROI_EffectSizes_PrettyNames.png"),
    width = 1200, height = 1000, res = 200)

pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of ROIs by Cohen's f (Effect Size)",
  fontsize_row = 6,
  fontsize_col = 8
)

dev.off()
# Hierarchical clustering
row_dendro <- hclust(dist(roi_mat))

# Cut tree into, say, 4 clusters (you can adjust `k`)
cluster_assignments <- cutree(row_dendro, k = 4)

# Create a data frame with ROI labels and cluster IDs
roi_clusters_df <- data.frame(
  ROI_label = names(cluster_assignments),
  Cluster = as.factor(cluster_assignments)
)

# Save to CSV
write.csv(roi_clusters_df, file.path(output_dir, "ROI_Cluster_Assignments.csv"), row.names = FALSE)

# Also save to Excel if needed
wb <- createWorkbook()
addWorksheet(wb, "ROI Clusters")
writeData(wb, "ROI Clusters", roi_clusters_df)
saveWorkbook(wb, file = file.path(output_dir, "ROI_Cluster_Assignments.xlsx"), overwrite = TRUE)

# Matrix of ROI vs Metric
roi_metric_matrix <- bind_rows(sig_results_list, .id = "Metric") %>%
  filter(p_fdr < 0.05) %>%
  select(Metric, ROI, structure, cohen_f) %>%
  group_by(Metric, ROI, structure) %>%
  summarize(cohen_f = mean(cohen_f, na.rm = TRUE), .groups = "drop") %>%
  unite(ROI_label, ROI, structure, sep = ": ") %>%
  pivot_wider(names_from = Metric, values_from = cohen_f)

roi_mat <- as.matrix(roi_metric_matrix[,-1])
rownames(roi_mat) <- roi_metric_matrix$ROI_label
colnames(roi_mat) <- pretty_metric_names[colnames(roi_mat)]

# Heatmap
pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "Clustering of ROIs by Cohen's f (Effect Size)",
  fontsize_row = 6,
  fontsize_col = 8
)

# Save it
ggsave(file.path(output_dir, "Clustered_ROI_EffectSizes_PrettyNames.png"), width = 12, height = 10, dpi = 300)



####save individual heatmaps per factor
# === Loop over selected risk factors and generate separate heatmaps ===

selected_factors <- c("APOE", "Age_Months", "Sex", "Diet", "HN")

walk(selected_factors, function(factor_name) {
  
  # Filter significant results by factor
  roi_metric_matrix <- bind_rows(sig_results_list, .id = "Metric") %>%
    filter(p_fdr < 0.05, Factor == factor_name) %>%
    group_by(Metric, ROI, structure) %>%
    summarize(mean_f = mean(cohen_f, na.rm = TRUE), .groups = "drop") %>%
    unite(ROI_label, ROI, structure, sep = ": ") %>%
    pivot_wider(names_from = Metric, values_from = mean_f)
  
  if (nrow(roi_metric_matrix) == 0) {
    cat(glue::glue("⚠️ No significant ROIs for factor: {factor_name}\n"))
    return(NULL)
  }
  
  # Save CSV
  write.csv(roi_metric_matrix,
            file = file.path(output_dir, paste0("ROI_by_Metric_", factor_name, ".csv")),
            row.names = FALSE)
  
  # Prepare matrix
  roi_mat <- as.matrix(roi_metric_matrix[,-1])
  rownames(roi_mat) <- roi_metric_matrix$ROI_label
  colnames(roi_mat) <- pretty_metric_names[colnames(roi_mat)]
  
  # Save heatmap
  pheatmap::pheatmap(
    roi_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = viridis::viridis(100),
    main = paste("Effect Sizes (Cohen's f):", factor_name),
    fontsize_row = 5,
    fontsize_col = 8,
    filename = file.path(output_dir, paste0("Heatmap_ROI_EffectSizes_", factor_name, ".png")),
    width = 12,
    height = 10
  )
})
