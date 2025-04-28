# === LOAD LIBRARIES ===
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(readr)

# === PARAMETERS ===
path_out        <- "/Users/alex/AlexBadea_MyCode/olfaction/results/"
path_out_filter <- file.path(path_out, "results_filtered")
dir.create(path_out,        showWarnings = FALSE)
dir.create(path_out_filter, showWarnings = FALSE)

# === 1. LOAD MERGED DATA ===
olfactory_data_merged <- readRDS(
  "/Users/alex/AlexBadea_MyCode/olfaction/data/olfactory_data_with_eigengenes.rds"
)

# === 2. GET LIST OF METRICS ===
metrics_list <- unique(olfactory_data_merged$Metric)
metrics_list <- metrics_list[!grepl("b0|dwi|mm3|rtap|rtpp",
                                    metrics_list, ignore.case = TRUE)]

# === 3. DEFINE YOUR OLFACTORY‐MEMORY REGIONS ===
olf_mem_regions <- c(
  "Piriform_Cortex",
  "Amygdalopiriform_Transition_Area",
  "Posterolateral_Cortical_Amygdaloid_Area",
  "Perirhinal_Cortex",
  "Hippocampus",
  "Postsubiculum",
  "Parasubiculum",
  "Ventral_Orbital_Cortex"
)

# === 4. LOOP THROUGH METRICS & FIT WITHIN‐GROUP STANDARDIZED MODELS ===
all_results <- list()
for (metric_name in metrics_list) {
  message("Processing Metric: ", metric_name)
  
  tmp <- olfactory_data_merged %>%
    filter(Metric == metric_name) %>%
    pivot_longer(
      cols      = starts_with("eigengene"),
      names_to  = "Eigengene",
      values_to = "Eigengene_value"
    ) %>%
    group_by(structure, Eigengene) %>%
    mutate(
      zValue     = as.numeric(scale(Value)),
      zEigengene = as.numeric(scale(Eigengene_value))
    ) %>%
    nest() %>%
    mutate(
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
  
  all_results[[metric_name]] <- tmp
}

# === 5. COMBINE & FILTER SIGNIFICANT ===
assoc <- bind_rows(all_results) %>%
  filter(p_fdr < 0.05)

# === 6. FILTER TO OLFACTORY‐MEMORY REGIONS ===
assoc_olfmem <- assoc %>%
  filter(structure %in% olf_mem_regions)

# === 7. SAVE SIGNIFICANT & FILTERED CSV ===
write_csv(
  assoc_olfmem,
  file.path(path_out_filter, "olfactory_memory_std_associations_FDR05.csv")
)

# === 8. PREPARE TOP HITS FOR PLOTTING ===
top_hits_mem <- assoc_olfmem %>%
  group_by(Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 5) %>%
  ungroup() %>%
  mutate(
    Metric_clean = Metric %>%
      gsub("mrtrix|_|mean|RAS|060622", "", .),
    log10_p_fdr  = -log10(p_fdr)
  )

beta_max_mem <- max(abs(top_hits_mem$std_beta), na.rm = TRUE)

# === 9. PLOT SIGNED r FOR OLFACTORY‐MEMORY ROIs ===
p_mem_r <- ggplot(top_hits_mem,
                  aes(x = reorder(structure, std_beta),
                      y = std_beta,
                      fill = Metric_clean)) +
  geom_col() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
  scale_y_continuous(limits = c(-beta_max_mem, beta_max_mem), expand = c(0,0)) +
  coord_flip() +
  labs(
    title = "Olfactory‐Memory Region Correlations (r)",
    x     = "Region",
    y     = "Pearson r"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(path_out_filter, "olfactory_memory_r_fixed.png"),
  plot     = p_mem_r,
  width    = 12, height = 8, dpi = 200
)

# === 10. PLOT –log10(FDR) FOR OLFACTORY‐MEMORY ROIs ===
p_mem_logp <- ggplot(top_hits_mem,
                     aes(x = reorder(structure, log10_p_fdr),
                         y = log10_p_fdr,
                         fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
  labs(
    title = "Olfactory‐Memory Region Significance (–log10 FDR)",
    x     = "Region",
    y     = "-log10(FDR p-value)"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(path_out_filter, "olfactory_memory_logp_fixed.png"),
  plot     = p_mem_logp,
  width    = 12, height = 8, dpi = 200
)

# === 11. PLOT |r| FOR OLFACTORY‐MEMORY ROIs ===
top_hits_mem_abs <- top_hits_mem %>%
  mutate(abs_r = abs(std_beta))

p_mem_abs <- ggplot(top_hits_mem_abs,
                    aes(x = reorder(structure, abs_r),
                        y = abs_r,
                        fill = Metric_clean)) +
  geom_col() +
  coord_flip() +
  facet_grid(Metric_clean ~ Eigengene, scales = "fixed") +
  labs(
    title = "Olfactory‐Memory Region Association Magnitude (|r|)",
    x     = "Region",
    y     = "|r|"
  ) +
  theme_minimal() +
  theme(
    strip.text  = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(path_out_filter, "olfactory_memory_abs_r_fixed.png"),
  plot     = p_mem_abs,
  width    = 12, height = 8, dpi = 200
)


assoc_grouped <- assoc_olfmem %>%
  mutate(
    RegionGroup = recode(
      structure,
      Hippocampus                     = "Hippocampal_Formation",
      Postsubiculum                   = "Hippocampal_Formation",
      Parasubiculum                   = "Hippocampal_Formation",
      Amygdalopiriform_Transition_Area  = "Olfactory_Amygdala",
      Posterolateral_Cortical_Amygdaloid_Area = "Olfactory_Amygdala",
      Piriform_Cortex                 = "Primary_Olfactory",
      Perirhinal_Cortex               = "Perirhinal_Orbitofrontal",
      Ventral_Orbital_Cortex          = "Perirhinal_Orbitofrontal"
    )
  )

# === B) For each group, pick the single best (lowest p_fdr) hit ===
condensed_results <- assoc_grouped %>%
  group_by(RegionGroup, Metric, Eigengene) %>%
  slice_min(order_by = p_fdr, n = 1) %>%
  ungroup() %>%
  arrange(RegionGroup, Metric, Eigengene)

# === C) Save a condensed CSV ===
write_csv(
  condensed_results,
  file.path(path_out_filter, "olfactory_memory_condensed_by_group.csv")
)

# === D) Quick printed table preview ===
print(condensed_results %>% 
        select(RegionGroup, Metric, Eigengene, std_beta, p_fdr, r_squared))

# === E) (Optional) Plot condensed data ===
beta_max_c <- max(abs(condensed_results$std_beta))
p_condensed <- ggplot(condensed_results,
                      aes(x = reorder(RegionGroup, std_beta),
                          y = std_beta,
                          fill = Metric)) +
  geom_col() +
  facet_grid(Metric ~ Eigengene, scales = "fixed") +
  scale_y_continuous(limits = c(-beta_max_c, beta_max_c), expand = c(0,0)) +
  coord_flip() +
  labs(x = "Region Group", y = "r") +
  theme_minimal()
ggsave(file.path(path_out_filter, "condensed_olf_mem_r.png"),
       plot = p_condensed, width = 12, height = 8, dpi = 200)


cat("✅ Filtered results for olfactory-memory regions are saved in:\n", path_out_filter, "\n")
