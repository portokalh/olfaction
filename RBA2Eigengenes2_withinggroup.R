# === LOAD LIBRARIES ===
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(readr)

# the MAP-based mean squared displacement (MSD), q-space inverse variance (QIV),
# non-Gaussianity (NG) and return-to-origin probability (RTOP)

# === PARAMETERS ===
path_out <- "/Users/alex/AlexBadea_MyCode/olfaction/results/"
dir.create(path_out, showWarnings = FALSE)

# === 1. LOAD MERGED DATA ===
olfactory_data_merged <- readRDS(
  "/Users/alex/AlexBadea_MyCode/olfaction/data/olfactory_data_with_eigengenes.rds"
)

# === 2. GET LIST OF METRICS ===
metrics_list <- unique(olfactory_data_merged$Metric)
metrics_list <- metrics_list[!grepl("b0|dwi|mm3|rtap|rtpp",
                                    metrics_list, ignore.case = TRUE)]

# === 3. INITIALIZE EMPTY LIST FOR RESULTS ===
all_results <- list()

# === 4. LOOP THROUGH METRICS & FIT STANDARDIZED (WITHIN‐GROUP) MODELS ===
for (metric_name in metrics_list) {
  cat("Processing Metric:", metric_name, "\n")
  
  data_subset <- olfactory_data_merged %>%
    filter(Metric == metric_name)
  
  result <- data_subset %>%
    pivot_longer(
      cols      = starts_with("eigengene"),
      names_to  = "Eigengene",
      values_to = "Eigengene_value"
    ) %>%
    #  1) group, 2) z‐score within each (structure, Eigengene), 3) nest
    group_by(structure, Eigengene) %>%
    mutate(
      zValue     = as.numeric(scale(Value)),
      zEigengene = as.numeric(scale(Eigengene_value))
    ) %>%
    nest() %>%
    mutate(
      # fit the standardized model only (slope == Pearson r)
      model_std  = map(data, ~ lm(zEigengene ~ zValue, data = .x)),
      std_beta   = map_dbl(model_std, ~ coef(.x)[["zValue"]]),
      r_squared  = map_dbl(model_std, ~ summary(.x)$r.squared),
      tidy_std   = map(model_std, ~ tidy(.x) %>% filter(term == "zValue")),
      std_se     = map_dbl(tidy_std, ~ .x$std.error),
      std_stat   = map_dbl(tidy_std, ~ .x$statistic),
      p_value    = map_dbl(tidy_std, ~ .x$p.value)
    ) %>%
    ungroup() %>%
    mutate(
      Metric = metric_name,
      p_fdr  = p.adjust(p_value, method = "fdr")
    ) %>%
    select(
      structure, Metric, Eigengene,
      std_beta, std_se, std_stat,
      p_value, p_fdr, r_squared
    )
  
  all_results[[metric_name]] <- result
}

# === 5. COMBINE ALL RESULTS ===
association_results_by_metric <- bind_rows(all_results)

# === 6. SAVE FULL RESULTS ===
write_csv(
  association_results_by_metric,
  file.path(path_out, "olfactory_metric_eigengene_std_associations_R2.csv")
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
      gsub("RAS",   "", .) %>%
      gsub("060622","", .) %>%
      recode(volumemm3 = "vol"),
    log10_p_fdr = -log10(p_fdr)
  )

# === 9. COMPUTE DYNAMIC R-LIMIT FOR PLOT ===
beta_max <- max(abs(top_hits$std_beta), na.rm = TRUE)

# === 10. PLOT SIGNED r (Pearson correlation) ===
p_effect_size <- ggplot(top_hits,
                        aes(x = reorder(structure, std_beta),
                            y = std_beta,
                            fill = Metric_clean)) +
  geom_col() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
  scale_y_continuous(limits = c(-beta_max, beta_max), expand = c(0,0)) +
  coord_flip() +
  labs(
    title = "Top Olfactory Region Associations by Metric\n(Standardized β = Pearson r)",
    x     = "Brain Region",
    y     = "r"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_regions_r_fixed.png"),
  plot     = p_effect_size,
  width    = 16, height = 16, dpi = 200
)

# === 11. PLOT –log10(FDR p-value) ===
p_logp <- ggplot(top_hits,
                 aes(x = reorder(structure, log10_p_fdr),
                     y = log10_p_fdr,
                     fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
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
  filename = file.path(path_out, "top_regions_logp_fixed.png"),
  plot     = p_logp,
  width    = 16, height = 16, dpi = 200
)

# === 12. PLOT │r│ (Magnitude Only) ===
top_hits_abs <- top_hits %>%
  mutate(abs_r = abs(std_beta))

p_abs_beta <- ggplot(top_hits_abs,
                     aes(x = reorder(structure, abs_r),
                         y = abs_r,
                         fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
  labs(
    title = "Top Olfactory Region Associations by Metric\n(|r|)",
    x     = "Brain Region",
    y     = "|r|"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(path_out, "top_regions_abs_r_fixed.png"),
  plot     = p_abs_beta,
  width    = 16, height = 16, dpi = 200
)

cat("✅ Done! Within-group standardized correlations (r), R², and FDR results saved, plus fixed-axis plots.\n")
