library(pheatmap)
library(viridis)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(tibble)

# --- Define output path ---
output_dir <- "olfactory_filtered_results2/clusters/"
dir.create(output_dir, showWarnings = FALSE)

# --- Assuming roi_mat is already defined with rownames like "41: SomeRegion" ---
# And olfactory_roi_ids is available from earlier:
olfactory_roi_ids <- c(41, 42, 43, 44, 45, 47, 53, 57, 58, 65, 66, 67, 70, 71, 118, 124, 125, 166,
                       1041, 1042, 1043, 1044, 1045, 1047, 1053, 1057, 1058, 1065, 1066, 1067, 1070, 1071, 1118, 1124, 1125, 1166)

# Extract ROI numeric codes from rownames
roi_labels <- rownames(roi_mat)
roi_ids <- as.numeric(sub(":.*", "", roi_labels))
olfactory_flags <- ifelse(roi_ids %in% olfactory_roi_ids, "Olfactory", "Other")

# --- Clustering ---
row_dendro <- hclust(dist(roi_mat))
roi_clusters <- cutree(row_dendro, k = 4)

# --- Annotations ---
annotation_row <- data.frame(
  Cluster = factor(roi_clusters),
  Olfactory = factor(olfactory_flags, levels = c("Olfactory", "Other"))
)
rownames(annotation_row) <- roi_labels

# Define colors
cluster_colors <- RColorBrewer::brewer.pal(length(unique(roi_clusters)), "Set2")
names(cluster_colors) <- levels(annotation_row$Cluster)

olfactory_colors <- c("Olfactory" = "#fc8d62", "Other" = "gray90")

ann_colors <- list(
  Cluster = cluster_colors,
  Olfactory = olfactory_colors
)

# --- Heatmap with annotation ---
pheatmap::pheatmap(
  roi_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis::viridis(100),
  main = "ROI Clusters Highlighting Olfactory Regions",
  fontsize_row = 6,
  fontsize_col = 8,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  filename = file.path(output_dir, "Clustered_ROI_Heatmap_Annotated_Olfactory.png"),
  width = 12,
  height = 10
)

# --- Save annotation table ---
annotated_df <- tibble(ROI_label = roi_labels, ROI_ID = roi_ids, Cluster = roi_clusters, Olfactory = olfactory_flags)
write.csv(annotated_df, file.path(output_dir, "ROI_Clusters_With_Olfactory_Flags.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "ROI Clusters Annotated")
writeData(wb, "ROI Clusters Annotated", annotated_df)
saveWorkbook(wb, file = file.path(output_dir, "ROI_Clusters_With_Olfactory_Flags.xlsx"), overwrite = TRUE)
