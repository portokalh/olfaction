# === LOAD LIBRARIES ===
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(readr)

# === PARAMETERS ===
path_out <- "/Users/alex/AlexBadea_MyCode/olfaction/results/"
dir.create(path_out, showWarnings = FALSE)

# === 1. LOAD MERGED DATA ===
olfactory_data_merged <- readRDS(
  "/Users/alex/AlexBadea_MyCode/olfaction/data/olfactory_data_with_eigengenes.rds"
)

# === 2. GET LIST OF METRICS ===
metrics_list <- unique(olfactory_data_merged$Metric)

# === 3. INITIALIZE EMPTY LIST FOR RESULTS ===
all_results <- list()

# === 4. LOOP THROUGH METRICS ===
for (metric_name in metrics_list) {
  cat("Processing Metric:", metric_name, "\n")
  
  data_subset <- olfactory_data_merged %>%
    filter(Metric == metric_name)
  
  result <- data_subset %>%
    pivot_longer(
      cols = starts_with("eigengene"),
      names_to  = "Eigengene",
      values_to = "Eigengene_value"
    ) %>%
    group_by(structure, Eigengene) %>%
    nest() %>%
    mutate(
      model   = map(data, ~ lm(Eigengene_value ~ Value, data = .x)),
      tidied  = map(model, tidy)
    ) %>%
    unnest(tidied) %>%
    filter(term == "Value") %>%
    mutate(Metric = metric_name) %>%
    select(structure, Metric, Eigengene, estimate, std.error, statistic, p.value) %>%
    mutate(p_fdr = p.adjust(p.value, method = "fdr"))
  
  all_results[[metric_name]] <- result
}

# === 5. COMBINE ALL RESULTS ===
association_results_by_metric <- bind_rows(all_results)

# === 6. SAVE FULL RESULTS ===
write_csv(
  association_results_by_metric,
  file.path(path_out, "olfactory_metric_by_metric_eigengene_associations.csv")
)

# === 7. FILTER SIGNIFICANT RESULTS (FDR < 0.05) ===
significant_results_by_metric <- association_results_by_metric %>%
  filter(p_fdr < 0.05)

write_csv(
  significant_results_by_metric,
  file.path(path_out, "olfactory_metric_by_metric_eigengene_significant_FDR05.csv")
)

# === 8. SELECT TOP HITS, CLEAN METRIC NAMES & COMPUTE -log10(p_fdr) ===
top_hits <- significant_results_by_metric %>%
  group_by(Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 5) %>%
  ungroup() %>%
  mutate(
    # clean metric labels
    Metric_clean = Metric,
    Metric_clean = gsub("mrtrix", "", Metric_clean),
    Metric_clean = gsub("_",     "", Metric_clean),
    Metric_clean = gsub("mean",  "", Metric_clean),
    Metric_clean = gsub("RAS_060622",  "", Metric_clean),
    Metric_clean = ifelse(Metric_clean == "volumemm3", "vol", Metric_clean),
    Metric_clean = ifelse(Metric_clean == "qiv",       "qiv", Metric_clean),
    # compute log‐p for plotting
    log10_p_fdr   = -log10(p_fdr)
  )

# === 9. PLOT EFFECT SIZE ===
p_effect_size <- ggplot(top_hits,
                        aes(x = reorder(structure, estimate),
                            y = estimate,
                            fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "free_y") +
  labs(
    title = "Top Olfactory Region Associations by Metric (Effect Size)",
    x = "Brain Region",
    y = "Effect Size (Estimate)"
  ) +
  theme_minimal() +
  theme(
    strip.text   = element_text(size = 7),
    axis.text.y  = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_olfactory_regions_effectsize_by_metric_clean.png"),
  plot     = p_effect_size,
  width    = 16, height = 10, dpi = 200
)

# === 10. PLOT –log10(FDR p-value) ===
p_logp <- ggplot(top_hits,
                 aes(x = reorder(structure, log10_p_fdr),
                     y = log10_p_fdr,
                     fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "free_y") +
  labs(
    title = "-log10(FDR p-value) for Top Olfactory Regions by Metric",
    x = "Brain Region",
    y = "-log10(FDR p-value)"
  ) +
  theme_minimal() +
  theme(
    strip.text   = element_text(size = 7),
    axis.text.y  = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_olfactory_regions_logp_by_metric_clean.png"),
  plot     = p_logp,
  width    = 16, height = 10, dpi = 200
)

# === 11. BONUS: MAKE SUMMARY TABLE ===
top_region_summary <- significant_results_by_metric %>%
  group_by(Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 1) %>%
  select(Metric, Eigengene, structure, estimate, p_fdr) %>%
  rename(
    Top_Region   = structure,
    Effect_Size  = estimate,
    FDR_p_value  = p_fdr
  ) %>%
  arrange(Metric, Eigengene)

write_csv(
  top_region_summary,
  file.path(path_out, "top_region_per_metric_eigengene.csv")
)

cat("🎯 All analyses, plots, and summary tables completed successfully!\n")
