# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(broom)
library(purrr)

# === 1. Read the large merged dataframe ===
path_data='/Users/alex/AlexBadea_MyCode/olfaction/data/'
path_out='/Users/alex/AlexBadea_MyCode/olfaction/results/'
merged_df_long <- read_csv(paste0(path_data, "merged_df_long.csv")) # Update with real path

# === 2. Save it as RDS for speed next time ===
saveRDS(merged_df_long, paste0(path_data,"merged_df_long.rds"))

# === 3. Subset to olfactory regions ===
# Define a list of olfactory-related structures (example, adjust as needed!)
olfactory_rois <- read_csv("/Users/alex/AlexBadea_MyCode/olfaction/data/olfactory_ROI_abbreviations.csv")


olfactory_structure_names <- olfactory_rois$structure

# === 3. Read merged_df_long ===
merged_df_long <- readRDS(paste0(path_data,"merged_df_long.rds"))  # update path if needed

# === 4. Filter merged_df_long to retain olfactory structures only (preserve duplicates) ===
olfactory_data <- merged_df_long %>%
  filter(structure %in% olfactory_structure_names)

# === 4. Further subset to metrics of interest ===
olfactory_data_selected <- olfactory_data %>%
  filter(grepl("mean", Metric, ignore.case = TRUE) | Metric == "volume_mm3")

# === Optional: see what we kept ===
olfactory_data_selected %>%
  distinct(Metric)




# === 5. Read eigengenes ===
eigengenes <- read_csv(paste0(path_data, "/eigengenes.csv")) # Your upload

# === 6. Merge with eigengenes ===
olfactory_data_merged <- olfactory_data_selected %>%
  inner_join(eigengenes, by = c("MouseID" = "Badea_ID"))  # carefully matching different names


# Save the merged olfactory metrics + eigengenes
saveRDS(olfactory_data_merged, paste0(path_data,"olfactory_data_with_eigengenes.rds"))

# Optional: save also as CSV if you want it viewable in Excel
write_csv(olfactory_data_merged, paste0(path_data,"olfactory_data_with_eigengenes.csv"))


# === Run the association models ===
association_results <- olfactory_data_merged %>%
  pivot_longer(cols = starts_with("eigengene"), names_to = "Eigengene", values_to = "Eigengene_value") %>%
  group_by(structure, Metric, Eigengene) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Eigengene_value ~ Value, data = .x)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Value") %>%
  select(structure, Metric, Eigengene, estimate, std.error, statistic, p.value)

association_results <- association_results %>%
  mutate(p_fdr = p.adjust(p.value, method = "fdr"))

# === Save results ===
write_csv(association_results, paste0(path_out, "olfactory_metric_eigengene_associations.csv"))

# Filter to significant results
significant_results <- association_results %>%
  filter(p_fdr < 0.05)

# Save significant results separately if you want
write_csv(significant_results, paste0(path_out, "olfactory_metric_eigengene_significant_FDR05.csv"))


library(ggplot2)

# Pick top N significant hits per eigengene
top_hits <- significant_results %>%
  group_by(Eigengene) %>%
  slice_min(order_by = p_fdr, n = 5) %>%
  ungroup()

# Plot
plot1 <- ggplot(top_hits, aes(x = reorder(structure, estimate), y = estimate, fill = Metric)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Eigengene, scales = "free_y") +
  labs(
    title = "Top Olfactory Region-Metric Associations with Eigengenes",
    x = "Brain Region",
    y = "Effect Size (Estimate)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

# === 3. Save plot at 200 dpi ===
ggsave(
  filename = paste0(path_out, "top_olfactory_regions.png"),
  plot = plot1,
  width = 12, height = 8, dpi = 200
)


library(ggplot2)

# === 1. Compute -log10(FDR p-value) ===
top_hits <- significant_results %>%
  mutate(log10_p_fdr = -log10(p_fdr)) %>%
  group_by(Eigengene) %>%
  slice_max(order_by = log10_p_fdr, n = 5) %>%
  ungroup()

# === 2. Plot -log10(FDR p-value) ===
p_logp <- ggplot(top_hits, aes(x = reorder(structure, log10_p_fdr), y = log10_p_fdr, fill = Metric)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Eigengene, scales = "free_y") +
  labs(
    title = "-log10(FDR-corrected p-values) for Top Olfactory Region-Metric Associations",
    x = "Brain Region",
    y = "-log10(FDR p-value)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

# === 3. Save plot at 200 dpi ===
ggsave(
  filename = paste0(path_out, "top_olfactory_regions_logp.png"),
  plot = p_logp,
  width = 12, height = 8, dpi = 200
)


