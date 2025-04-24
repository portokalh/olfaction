library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(writexl)
library(stringr)
library(dplyr)
library(broom)
library(openxlsx)
library(furrr)
#install.packages("progressr")  # if not already installed
library(progressr)
library(effectsize) 

# Load required packages
library(dplyr)
library(tibble)
library(broom)
library(furrr)
library(progressr)
library(openxlsx)

# --- Define paths ---
stats_dir <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/data/mouse_region_study/corected_stats/"
#metadata_path <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/data/MouseMetaData092024AB.xlsx"
metadata_path <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/data/MouseMetaData042325.xlsx"

output_file <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/results/RBA_Corrected_ANOVA_by_metric.xlsx"
output_path <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/results/"

# --- Load metadata ---
metadata <- read_excel(metadata_path)
colnames(metadata)[1] <- "MouseID"
metadata$MouseID <- as.character(metadata$MouseID)

# --- List all *_stats.xlsx files ---
stat_files <- list.files(stats_dir, pattern = "_stats\\.xlsx$", full.names = TRUE)
cat("📁 Number of stat files found:", length(stat_files), "\n")
# --- Function to load one stats file and tag MouseID ---
load_stats <- function(file) {
  df <- read_xlsx(file)  # ✅ FIXED: removed show_col_types
  df$MouseID <- str_split(basename(file), "_", simplify = TRUE)[1]
  return(df)
}

# --- Combine all files into one dataframe ---
all_stats <- map_dfr(stat_files, load_stats)

long_df <- all_stats %>%
  rename(DWI = MouseID) %>%
  pivot_longer(cols = c(ad_mean, fa_mean, md_mean), names_to = "Metric", values_to = "Value")

cat("📊 Number of unique DWI entries in long_df:", n_distinct(long_df$DWI), "\n")


long_df <- long_df %>%
  mutate(DWI = substr(DWI, 1, 6))

head(unique(long_df$DWI), 209)


metadata <- metadata %>%
  mutate(
    APOE = case_when(
      str_detect(Genotype, "APOE2") ~ "APOE2",
      str_detect(Genotype, "APOE3") ~ "APOE3",
      str_detect(Genotype, "APOE4") ~ "APOE4",
      TRUE ~ NA_character_
    ),
    HN = if_else(str_detect(Genotype, "HN"), 1, -1)
  )

table(is.na(metadata$APOE))
table(is.na(metadata$HN))
metadata %>%
  filter(is.na(APOE) | is.na(HN)) %>%
  select(MouseID, Genotype, APOE, HN) %>%
  print(n = 20)

dropped_rows <- metadata %>%
  filter(is.na(APOE) | is.na(HN)) %>%
  select(MouseID, Genotype, APOE, HN)

print(dropped_rows, n = Inf)  # Shows all rows

metadata <- metadata %>%
  filter(str_starts(DWI, "N")) %>%
  mutate(
    APOE = factor(APOE),
    HN = factor(HN),
    Diet = factor(Diet),
    Sex = factor(Sex)
  )

metadata <- metadata %>%
  filter(!is.na(APOE), !is.na(HN))

print(unique(metadata$HN))
table(metadata$HN, useNA = "always")
print(unique(metadata$Diet))
table(metadata$Diet, useNA = "always")

###sanity check 
# 1. Create a small dataframe with unique MouseID and relevant columns
mouse_summary <- metadata %>%
  filter(!is.na(DWI)) %>%
  select(MouseID, Diet, HN) %>%
  distinct()  # ensures one row per MouseID

# 2. Count how many mice per Diet and HN
mouse_counts <- mouse_summary %>%
  count(Diet, HN) %>%
  arrange(desc(n))

# 3. View the result
print(mouse_counts)
# Total count
total_mice <- nrow(mouse_summary)

cat("🐭 Total unique mice with DWI data:", total_mice, "\n")
###SANITY CHECK


print(unique(metadata$HN))
print(unique(metadata$APOE))
print(unique(metadata$Diet))

summary_table <- metadata %>%
  filter(!is.na(APOE)) %>%
  count(APOE, Sex, Diet, HN, name = "Count")

# View the result
summary_table

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)

addWorksheet(wb, "summary")
writeData(wb, "summary", summary_table)

saveWorkbook(wb, paste0(output_path,"qialmetadata_APOE.xlsx"), overwrite = TRUE)



# --- Merge with metadata ---
metadata_unique <- metadata %>%
  distinct(DWI, .keep_all = TRUE)

merged_df <- inner_join(long_df, metadata_unique, by = "DWI")
cat("🔢 Number of unique DWI entries in merged_df:", n_distinct(merged_df$DWI), "\n")

# DWI values in long_df but not in metadata
anti_join(long_df, metadata_unique, by = "DWI") %>%
  count(DWI)  # Or View() it

# merged_df <- merged_df %>%
#   filter(
#     !is.na(APOE),
#     !is.na(HN),
#     !str_detect(Genotype, "KO")
#   ) %>%
#   mutate(
#     APOE = factor(APOE),
#     HN = factor(HN),
#     Diet = factor(Diet),
#     Sex = factor(Sex)
#   )

#write_xlsx(merged_df, paste0(output_path, "merged_df.xlsx"))
write.csv(merged_df, file = paste0(output_path, "merged_df.csv"), row.names = FALSE)

# --- List of metrics of interest (simplified) ---
metrics_of_interest <- c(
  "volume_mm3", "dwi_mean", "fa_mean", 
  "ad_mean", "rd_mean", "md_mean", "mrtrixfa_mean", "mrtrixad_mean",
  "mrtrixrd_mean", "mrtrixmd_mean", "msd_RAS_060622_mean",
  "ng_RAS_060622_mean", "qiv_RAS_060622_mean", "rtap_RAS_060622_mean",
  "rtop_RAS_060622_mean", "rtpp_RAS_060622_mean", "volume_percentage_mean"
)

cat("🔢 Number of unique DWI values in merged_df:", n_distinct(merged_df$DWI), "\n")
unique_dwids <- unique(merged_df$DWI)
print(unique_dwids)





###ANOVA NOW#####
# Setup parallel backend
future::plan(multisession)  # Use multicore on Linux/macOS
handlers(global = TRUE)     # Enable progress bar in console


# Filter for metrics of interest
filtered_df <- merged_df %>%
  filter(Metric %in% metrics_of_interest)

# Get combinations of Metric × ROI from filtered data
metric_roi_combos <- filtered_df %>%
  distinct(Metric, ROI) %>%
  filter(!is.na(Metric), !is.na(ROI))

# Define the ANOVA function
run_anova <- function(metric, roi) {
  df_sub <- filtered_df %>%
    filter(ROI == roi, Metric == metric) %>%
    filter(!is.na(Value), !is.na(APOE)) %>%
    mutate(APOE = as.factor(APOE))
  
  if (nrow(df_sub) < 10 || length(unique(df_sub$APOE)) < 2) return(NULL)
  
  model <- tryCatch(lm(Value ~ APOE, data = df_sub), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  anova_table <- tryCatch(anova(model), error = function(e) NULL)
  if (is.null(anova_table)) return(NULL)
  
  if ("APOE" %in% rownames(anova_table)) {
    return(tibble(
      ROI = roi,
      Metric = metric,
      df = anova_table["APOE", "Df"],
      F_value = anova_table["APOE", "F value"],
      p_value = anova_table["APOE", "Pr(>F)"]
    ))
  } else {
    return(NULL)
  }
}

# Run in parallel with progress bar
anova_results <- with_progress({
  p <- progressor(along = 1:nrow(metric_roi_combos))
  
  future_pmap_dfr(metric_roi_combos, function(Metric, ROI) {
    p()
    run_anova(Metric, ROI)
  })
})

# Add FDR-corrected p-values
anova_results <- anova_results %>%
  mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"))

# Merge region names (structure) from metadata
roi_names <- merged_df %>%
  select(ROI, structure) %>%
  distinct() %>%
  filter(!is.na(structure))

anova_results_named <- anova_results %>%
  left_join(roi_names, by = "ROI") %>%
  relocate(structure, .after = ROI)

# Write hybrid-format Excel file
wb <- createWorkbook()
addWorksheet(wb, "All_Results")
writeDataTable(wb, "All_Results", anova_results_named)

for (metric in unique(anova_results_named$Metric)) {
  df_metric <- filter(anova_results_named, Metric == metric)
  addWorksheet(wb, sheetName = metric)
  writeDataTable(wb, metric, df_metric)
}

saveWorkbook(wb, paste0(output_path, "anova_lm_results_named_by_metric_and_roi.xlsx"), overwrite = TRUE)

# Reset to sequential plan
future::plan(sequential)

cat("✅ ANOVA complete. Results saved to: anova_lm_results_named_by_metric_and_roi.xlsx\n")


##################################
####make a more complex model#####
##################################
# === Setup ===
future::plan(multisession)
handlers(global = TRUE)

anova_results <- list() 

metrics_of_interest <- c(
  "volume_mm3", "dwi_mean", "fa_mean", 
  "ad_mean", "rd_mean", "md_mean", "mrtrixfa_mean", "mrtrixad_mean",
  "mrtrixrd_mean", "mrtrixmd_mean", "msd_RAS_060622_mean",
  "ng_RAS_060622_mean", "qiv_RAS_060622_mean", "rtap_RAS_060622_mean",
  "rtop_RAS_060622_mean", "rtpp_RAS_060622_mean"
)

output_path <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/results/"
output_file <- paste0(output_path, "anova_fullmodel_effectsizes_tstats.xlsx")

# Filter for metrics of interest
filtered_df <- merged_df %>%
  filter(Metric %in% metrics_of_interest)

# Get combinations of Metric × ROI from filtered data
metric_roi_combos <- filtered_df %>%
  distinct(Metric, ROI) %>%
  filter(!is.na(Metric), !is.na(ROI))



# === ANOVA function with effect sizes and t-ratios ===
run_anova <- function(metric, roi) {
  df_sub <- filtered_df %>%
    filter(ROI == roi, Metric == metric) %>%
    filter(!is.na(Value), !is.na(APOE), !is.na(Age_Months),
           !is.na(Sex), !is.na(Diet), !is.na(HN)) %>%
    mutate(
      APOE = factor(APOE),
      HN = factor(HN),
      Diet = factor(Diet),
      Sex = factor(Sex),
      Age_Months = as.numeric(Age_Months)
    )
  
  if (nrow(df_sub) < 10 || length(unique(df_sub$APOE)) < 2) return(NULL)
  
  model <- tryCatch(lm(Value ~ APOE + Age_Months + Sex + Diet + HN, data = df_sub), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  
  # ANOVA table
  anova_tbl <- tryCatch(anova(model), error = function(e) NULL)
  if (is.null(anova_tbl)) return(NULL)
  
  # Effect sizes
  eta <- tryCatch(eta_squared(model, partial = TRUE, ci = NULL), error = function(e) NULL)
  if (is.null(eta)) return(NULL)
  
  # t-stats from model coefficients
  tidy_tbl <- tryCatch(tidy(model), error = function(e) NULL)
  if (is.null(tidy_tbl)) return(NULL)
  
  # Merge all by term name
  merged_terms <- full_join(
    tibble(term = rownames(anova_tbl), df = anova_tbl$Df, F_value = anova_tbl$`F value`, p_value = anova_tbl$`Pr(>F)`),
    eta %>% rename(term = Parameter, eta2_partial = Eta2_partial),
    by = "term"
  ) %>%
    left_join(
      tidy_tbl %>% select(term, estimate, std.error, statistic, p.value) %>% 
        rename(t_statistic = statistic, t_p_value = p.value),
      by = "term"
    ) %>%
    mutate(ROI = roi, Metric = metric) %>%
    select(ROI, Metric, term, df, F_value, p_value, eta2_partial, estimate, std.error, t_statistic, t_p_value)
  
  return(merged_terms)
}

# === Run in parallel with progress ===
anova_results <- with_progress({
  p <- progressor(along = 1:nrow(metric_roi_combos))
  
  future_pmap_dfr(metric_roi_combos, function(Metric, ROI) {
    p()
    run_anova(Metric, ROI)
  })
})

# === FDR correction ===
anova_results <- anova_results %>%
  mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"))

# === Merge structure names ===
roi_names <- merged_df %>%
  select(ROI, structure) %>%
  distinct() %>%
  filter(!is.na(structure))

anova_results_named <- anova_results %>%
  left_join(roi_names, by = "ROI") %>%
  relocate(structure, .after = ROI)

# === Write to Excel ===
wb <- createWorkbook()
addWorksheet(wb, "All_Results")
writeDataTable(wb, "All_Results", anova_results_named)

for (metric in unique(anova_results_named$Metric)) {
  df_metric <- filter(anova_results_named, Metric == metric)
  addWorksheet(wb, sheetName = metric)
  writeDataTable(wb, metric, df_metric)
}

saveWorkbook(wb, paste0(output_path, "anova_complex.xlsx"), overwrite = TRUE)
future::plan(sequential)

cat("✅ ANOVA + effect sizes + t-stats complete.\n")
cat("Results saved to:", paste0(output_path, "anova_complex.xlsx"), "\n")




##################################
###make plots library(ggplot2)
##################################

# Create plot output folder
plot_dir <- file.path(output_path, "anova_top_term_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# How many top effects per metric to show
top_n <- 10

# Loop over each metric and save top plots
unique_metrics <- unique(anova_results_named$Metric)

for (m in unique_metrics) {
  df_metric <- anova_results_named %>%
    filter(Metric == m, !is.na(p_adj_fdr)) %>%
    arrange(p_adj_fdr) %>%
    filter(p_adj_fdr < 0.05) %>%
    group_by(term) %>%
    slice_min(p_adj_fdr, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(df_metric) == 0) next  # Skip if no significant results
  
  p <- ggplot(df_metric, aes(x = reorder(structure, -F_value), y = F_value, fill = term)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ term, scales = "free_y") +
    coord_flip() +
    labs(
      title = paste("Top", top_n, "Effects by Region for", m),
      x = "Brain Region",
      y = "F Statistic"
    ) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(
    filename = file.path(plot_dir, paste0("Top_Regions_", m, ".pdf")),
    plot = p, width = 10, height = 7
  )
}

