library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
library(tidyr)
library(pheatmap)
library(viridis)

# === PARAMETERS ===
output_dir <- "olfactory_filtered_results2"
dpi_val <- 300
sig_thresh <- 0.05

# Color and shape definitions
factor_colors <- c("APOE" = "#E41A1C", "Age" = "#377EB8",
                   "Sex" = "#4DAF4A", "Diet" = "#984EA3", "HN" = "#FF7F00")
olfactory_shape <- c("Other" = 16, "Olfactory" = 17)

# Olfactory ROI IDs
olfactory_rois <- c(41, 42, 43, 44, 45, 47, 53, 57, 58, 65, 66, 67, 70, 71, 118, 124, 125, 166)
olfactory_rois <- c(olfactory_rois, olfactory_rois + 1000)

# === LOAD DATA ===
top_by_factor <- read_tsv(file.path(output_dir, "TopByFactor.tsv")) %>%
  filter(p_fdr < sig_thresh)

top_by_roi <- read_tsv(file.path(output_dir, "TopByROI.tsv")) %>%
  filter(p_fdr < sig_thresh)

main_factors <- c("APOE", "Age", "Sex", "Diet", "HN")

# === CLEANED DATAFRAMES ===
top_by_factor_clean <- top_by_factor %>%
  mutate(Factor = gsub(":.*", "", Factor)) %>%
  filter(Factor %in% main_factors) %>%
  mutate(Olfactory = if_else(ROI %in% olfactory_rois, "Olfactory", "Other"),
         label = fct_reorder(paste0("ROI ", ROI, ": ", structure), cohen_f))

top_by_roi_clean <- top_by_roi %>%
  mutate(Factor = gsub(":.*", "", Factor)) %>%
  filter(Factor %in% main_factors) %>%
  mutate(Olfactory = if_else(ROI %in% olfactory_rois, "Olfactory", "Other"),
         label = fct_reorder(paste0("ROI ", ROI, ": ", structure), cohen_f))

# === EXPORT CSVs ===
write_csv(top_by_factor_clean, file.path(output_dir, "All_Significant_ROIs_per_Factor.csv"))
write_csv(top_by_roi_clean, file.path(output_dir, "Significant_ROIs_by_Factor.csv"))

# === PLOT 1: All Significant ROIs per Factor ===
p1 <- ggplot(top_by_factor_clean, aes(x = cohen_f, y = label, color = Factor, shape = Olfactory)) +
  geom_point(size = 2.5) +
  facet_wrap(~ Factor, scales = "free") +
  scale_color_manual(values = factor_colors) +
  scale_shape_manual(values = olfactory_shape) +
  labs(title = "All Significant ROIs per Risk Factor (FDR < 0.05)",
       x = "Cohen's f (Effect Size)", y = NULL) +
  xlim(0, max(top_by_factor_clean$cohen_f, na.rm = TRUE) * 1.1) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 7))

ggsave(file.path(output_dir, "All_Significant_ROIs_by_Factor.png"), plot = p1, width = 14, height = 12, dpi = dpi_val)

# === PLOT 2: Most Frequently Affected ROIs (Top 15) ===
roi_counts <- top_by_roi_clean %>%
  count(structure, sort = TRUE) %>%
  slice_max(n, n = 15)

roi_freq_plot <- ggplot(roi_counts, aes(x = reorder(structure, n), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), hjust = -0.1, size = 3) +
  expand_limits(y = max(roi_counts$n) * 1.2) +
  coord_flip() +
  labs(title = "Most Frequently Affected ROIs (FDR < 0.05)", x = NULL, y = "Count") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "Most_Frequent_ROIs_CLEANED.png"), plot = roi_freq_plot, width = 10, height = 6, dpi = dpi_val)

# === PLOT 3: Influential Factors by Count of Significant ROIs ===
factor_counts <- top_by_factor_clean %>%
  count(Factor, sort = TRUE)

factor_freq_plot <- ggplot(factor_counts, aes(x = reorder(Factor, n), y = n, fill = Factor)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.1, size = 3) +
  scale_fill_manual(values = factor_colors) +
  expand_limits(y = max(factor_counts$n) * 1.2) +
  coord_flip() +
  labs(title = "Number of Significant ROIs per Factor (FDR < 0.05)", x = NULL, y = "Count") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "Most_Influential_Factors_CLEANED.png"), plot = factor_freq_plot, width = 8, height = 5, dpi = dpi_val)

# === PLOT 4: Heatmap of ROI × Factor Effect Sizes ===
roi_factor_avg <- top_by_roi_clean %>%
  group_by(Factor, structure) %>%
  summarize(mean_f = mean(cohen_f, na.rm = TRUE), .groups = "drop")

heatmap_data <- roi_factor_avg %>%
  pivot_wider(names_from = Factor, values_from = mean_f, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$structure

pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Mean Cohen's f by ROI × Factor (FDR < 0.05)",
         fontsize_row = 8,
         fontsize_col = 10,
         color = viridis(100),
         filename = file.path(output_dir, "Heatmap_CohensF_ROI_by_Factor_ALLROIs.png"),
         width = 10,
         height = 8)

# === PLOT 5: Metric-Level Breakdown by Factor ===
metric_breakdown <- top_by_factor_clean %>%
  group_by(Factor, Metric) %>%
  summarise(count = n(), mean_f = mean(cohen_f, na.rm = TRUE), .groups = "drop")

metric_plot <- ggplot(metric_breakdown, aes(x = mean_f, y = reorder(Metric, mean_f), fill = Factor)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Factor, scales = "free_x") +
  scale_fill_manual(values = factor_colors) +
  labs(title = "Average Cohen's f by Metric and Factor (FDR < 0.05)",
       x = "Mean Cohen's f", y = "Metric") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 12, face = "bold"))

ggsave(file.path(output_dir, "Metric_Breakdown_by_Factor.png"), plot = metric_plot, width = 14, height = 8, dpi = dpi_val)
