# === LOAD LIBRARIES ===
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(readr)

#the MAP-based mean squared displacement (MSD), q-space inverse variance (QIV), non-Gaussianity (NG) and return-to-origin probability (RTOP)

# === PARAMETERS ===
path_out <- "/Users/alex/AlexBadea_MyCode/olfaction/results/"
dir.create(path_out, showWarnings = FALSE)

# === 1. LOAD MERGED DATA ===
olfactory_data_merged <- readRDS(
  "/Users/alex/AlexBadea_MyCode/olfaction/data/olfactory_data_with_eigengenes.rds"
)

# === 2. GET LIST OF METRICS ===
metrics_list <- unique(olfactory_data_merged$Metric)
metrics_list <- metrics_list[!grepl("b0|dwi|mm3|rtap|rtpp", metrics_list, ignore.case = TRUE)]

# === 3. INITIALIZE EMPTY LIST FOR RESULTS ===
all_results <- list()

# === 4. LOOP THROUGH METRICS & FIT BOTH RAW + STD MODELS ===
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
    #  scale both x and y
    mutate(
      zValue     = as.numeric(scale(Value)),
      zEigengene = as.numeric(scale(Eigengene_value))
    ) %>%
    group_by(structure, Eigengene) %>%
    nest() %>%
    mutate(
      # raw & standardized models
      model_raw  = map(data, ~ lm(Eigengene_value ~ Value,      data = .x)),
      model_std  = map(data, ~ lm(zEigengene      ~ zValue,     data = .x)),
      # extract coefficients
      raw_beta   = map_dbl(model_raw, ~ coef(.x)[["Value"]]),
      std_beta   = map_dbl(model_std, ~ coef(.x)[["zValue"]]),
      # extract R² (same for raw vs. std)
      r_squared  = map_dbl(model_std, ~ summary(.x)$r.squared),
      # p-value, std.error & t-stat for the standardized β
      std_tidy   = map(model_std, ~ tidy(.x) %>% filter(term == "zValue")),
      std_se     = map_dbl(std_tidy, ~ .x$std.error),
      std_stat   = map_dbl(std_tidy, ~ .x$statistic),
      p_value    = map_dbl(std_tidy, ~ .x$p.value)
    ) %>%
    ungroup() %>%
    mutate(
      Metric = metric_name,
      p_fdr  = p.adjust(p_value, method = "fdr")
    ) %>%
    select(
      structure, Metric, Eigengene,
      raw_beta, std_beta,
      std_se, std_stat,
      p_value, p_fdr, r_squared
    )
  
  all_results[[metric_name]] <- result
}

# === 5. COMBINE ALL RESULTS ===
association_results_by_metric <- bind_rows(all_results)

# === 6. SAVE FULL RESULTS ===
write_csv(
  association_results_by_metric,
  file.path(path_out, "olfactory_metric_eigengene_associations_with_betas_R2.csv")
)

# === 7. FILTER SIGNIFICANT RESULTS (FDR < 0.05) ===
significant_results_by_metric <- association_results_by_metric %>%
  filter(p_fdr < 0.05)

write_csv(
  significant_results_by_metric,
  file.path(path_out, "olfactory_metric_eigengene_significant_FDR05.csv")
)

# === 8. SELECT TOP HITS & CLEAN METRIC NAMES ===
top_hits <- significant_results_by_metric %>%
  group_by(Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 10) %>%
  ungroup() %>%
  mutate(
    Metric_clean = Metric %>%
      gsub("mrtrix", "", .) %>%
      gsub("_",     "", .) %>%
      gsub("mean",  "", .) %>%
      gsub("RAS", "", .) %>%
      gsub("060622", "", .) %>%
      gsub("RAS_060622", "", .) %>%
      recode(volumemm3 = "vol"),
    log10_p_fdr = -log10(p_fdr)
  )

# === 9. PLOT STANDARDIZED β ===
p_effect_size <- ggplot(top_hits,
                        aes(x = reorder(structure, std_beta),
                            y = std_beta,
                            fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "free_y") +
  labs(
    title = "Top Olfactory Region Associations by Metric\n(Standardized β)",
    x     = "Brain Region",
    y     = "Standardized β"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_regions_std_beta.png"),
  plot     = p_effect_size,
  width    = 16, height = 16, dpi = 200
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
    x     = "Brain Region",
    y     = "-log10(FDR p-value)"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_regions_logp.png"),
  plot     = p_logp,
  width    = 16, height = 16, dpi = 200
)

# === 11. SUMMARY TABLE WITH RAW & STD β, R² ===
top_region_summary <- significant_results_by_metric %>%
  group_by(Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 10) %>%
  select(
    Metric, Eigengene, Top_Region = structure,
    raw_beta, std_beta,
    FDR_p_value = p_fdr,
    R2 = r_squared
  ) %>%
  arrange(Metric, Eigengene)

write_csv(
  top_region_summary,
  file.path(path_out, "top_region_per_metric_eigengene_summary.csv")
)

# === 11. PLOT │Standardized β│ (Magnitude Only) ===
# first compute abs‐beta
top_hits_abs <- top_hits %>%
  mutate(abs_std_beta = abs(std_beta))

p_abs_beta <- ggplot(top_hits_abs,
                     aes(x = reorder(structure, abs_std_beta),
                         y = abs_std_beta,
                         fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "free_y") +
  labs(
    title = "Top Olfactory Region Associations by Metric\n(|Standardized β|)",
    x     = "Brain Region",
    y     = "|Standardized β|"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_regions_abs_std_beta.png"),
  plot     = p_abs_beta,
  width    = 16, height = 16, dpi = 200
)



cat("✅ Done! Raw β, Standardized β, p-values, FDR, and R² are all in your outputs.\n")
