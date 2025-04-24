library(openxlsx)
library(dplyr)
library(effectsize)
library(purrr)
library(furrr)
plan(multisession)

# === Get unique metrics ===
all_metrics <- unique(merged_df_long$Metric)

# === Master Excel workbook ===
wb <- createWorkbook()

# === Function to analyze one ROI for a given metric ===
analyze_one_roi <- function(df_metric, roi_id, target_metric) {
  df <- df_metric %>%
    filter(ROI == roi_id) %>%
    filter(!is.na(APOE), !is.na(HN), !is.na(Age_Months), !is.na(Sex), !is.na(Diet))
  
  if (nrow(df) < 10) return(NULL)
  
  tryCatch({
    lm_model <- lm(Value ~ APOE * HN * Age_Months * Sex * Diet, data = df)
    
    anova_df <- as.data.frame(anova(lm_model))
    anova_df$Factor <- rownames(anova_df)
    anova_df <- anova_df %>% select(Factor, everything())
    rownames(anova_df) <- NULL
    
    cf_tbl <- tryCatch({
      cohens_f(lm_model, partial = TRUE, ci = 0.95) %>%
        rename(Factor = Parameter) %>%
        select(Factor, cohen_f = Cohens_f_partial, f_CI_low = CI_low, f_CI_high = CI_high)
    }, error = function(e) {
      tibble(Factor = setdiff(anova_df$Factor, "Residuals"),
             cohen_f = NA_real_, f_CI_low = NA_real_, f_CI_high = NA_real_)
    })
    cf_tbl <- bind_rows(cf_tbl, tibble(Factor = "Residuals", cohen_f = NA, f_CI_low = NA, f_CI_high = NA))
    
    eta_tbl <- tryCatch({
      eta_squared(lm_model, partial = TRUE, ci = 0.95) %>%
        rename(Factor = Parameter) %>%
        select(Factor, eta2 = Eta2_partial, eta2_CI_low = CI_low, eta2_CI_high = CI_high)
    }, error = function(e) {
      tibble(Factor = setdiff(anova_df$Factor, "Residuals"),
             eta2 = NA_real_, eta2_CI_low = NA_real_, eta2_CI_high = NA_real_)
    })
    eta_tbl <- bind_rows(eta_tbl, tibble(Factor = "Residuals", eta2 = NA, eta2_CI_low = NA, eta2_CI_high = NA))
    
    results_df <- anova_df %>%
      left_join(cf_tbl, by = "Factor") %>%
      left_join(eta_tbl, by = "Factor") %>%
      mutate(
        ROI = roi_id,
        structure = df$structure[1],
        Metric = target_metric
      ) %>%
      select(Metric, ROI, structure, everything()) %>%
      mutate(across(contains("CI"), ~ ifelse(is.infinite(.), "Inf", .)))
    
    return(results_df)
  }, error = function(e) {
    cat("❌ Error in ROI", roi_id, "for metric", target_metric, ":", e$message, "\n")
    return(NULL)
  })
}

# === LOOP over each metric ===
walk(all_metrics, function(target_metric) {
  cat("\n📊 Processing metric:", target_metric, "\n")
  
  df_metric <- merged_df_long %>% filter(Metric == target_metric)
  roi_list <- unique(df_metric$ROI)
  
  all_results <- future_map(roi_list, analyze_one_roi, df_metric = df_metric, target_metric = target_metric, .progress = TRUE) %>%
    compact() %>%
    bind_rows()
  
  if (nrow(all_results) > 0) {
    # === Add FDR correction per metric over all ROI × term
    all_results <- all_results %>%
      group_by(Factor) %>%
      mutate(p_fdr = p.adjust(`Pr(>F)`, method = "fdr")) %>%
      ungroup()
    
    # === Write to workbook
    safe_name <- make.names(target_metric)
    addWorksheet(wb, sheetName = safe_name)
    writeData(wb, sheet = safe_name, all_results)
    
    cat("✅ Metric", target_metric, "written to sheet:", safe_name, "\n")
  } else {
    cat("⚠️ No valid results for metric:", target_metric, "\n")
  }
})

# === Save final summary workbook ===
saveWorkbook(wb, "summary_stats.xlsx", overwrite = TRUE)
cat("📁 All metrics saved to: summary_stats.xlsx\n")
