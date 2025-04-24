library(tidyr)
library(dplyr)
library(readr)

# Define metrics of interest
metrics_of_interest <- c(
  "volume_mm3", "dwi_mean", "fa_mean", "ad_mean", "rd_mean", "md_mean",
  "mrtrixfa_mean", "mrtrixad_mean", "mrtrixrd_mean", "mrtrixmd_mean",
  "msd_RAS_060622_mean", "ng_RAS_060622_mean", "qiv_RAS_060622_mean",
  "rtap_RAS_060622_mean", "rtop_RAS_060622_mean", "rtpp_RAS_060622_mean", "volume_percentage_mean", "b0_mean"
)

# Subset only existing columns
existing_metrics <- intersect(metrics_of_interest, colnames(merged_df))

# Drop any existing conflicting columns and pivot
merged_df_long <- merged_df %>%
  select(-any_of(c("Metric", "Value"))) %>%
  pivot_longer(
    cols = all_of(existing_metrics),
    names_to = "Metric",
    values_to = "Value"
  )

# Save result
write_csv(merged_df_long, "merged_df_long.csv")

cat("✅ Long-format data saved with", length(existing_metrics), "metrics → merged_df_long.csv\n")
