# === LOAD LIBRARIES ===
library(readxl)
library(dplyr)
library(stringr)
library(gprofiler2)
library(ggplot2)
library(tidyr)

# === PARAMETERS ===
input_file <- "eigengene_top_genes.xlsx"   # Your Excel file name
output_folder <- "eigengene_enrichment_results"  # Folder to save results

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# === READ ALL SHEETS ===
sheet_names <- excel_sheets(input_file)

# === FUNCTION TO PLOT ENRICHMENT ===
plot_enrichment <- function(enrichment_df, module_name) {
  top_terms <- enrichment_df %>%
    arrange(p_value) %>%
    slice(1:10)  # top 10 terms
  
  ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top GO Biological Processes - ", module_name),
         x = "GO Term",
         y = "-log10(p-value)") +
    theme_minimal()
}

# === PROCESS EACH MODULE ===
# === PROCESS EACH MODULE ===
for (sheet in sheet_names) {
  # Read each sheet
  df <- read_excel(input_file, sheet = sheet)
  
  # Extract top 100 genes (corrected!)
  top_genes <- df[[1]] %>% head(200)
  
  # Perform enrichment with gprofiler2
  # gost_res <- gost(query = top_genes,
  #                  organism = "mmusculus",  # Mouse
  #                  sources = c("GO:BP"),
  #                  significant = TRUE)
  
  gost_res <- gost(query = top_genes,
                   organism = "hsapiens",  # HUMAN, not mouse
                   sources = c("GO:BP"),
                   significant = FALSE)    # don't restrict
  
  # If there are results, save and plot
  if (!is.null(gost_res) && nrow(gost_res$result) > 0) {
    # Save results table
    # Save results table (fix: flatten list columns first)
    enrichment_fixed <- gost_res$result %>%
      mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))  # collapse lists into comma-separated strings
    
    output_file <- file.path(output_folder, paste0(sheet, "_enrichment_results.tsv"))
    write.table(enrichment_fixed, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Plot top GO terms
    plot <- plot_enrichment(gost_res$result, sheet)
    plot_file <- file.path(output_folder, paste0(sheet, "_top_GO_terms.png"))
    ggsave(plot_file, plot, width = 8, height = 6)
  }
}
