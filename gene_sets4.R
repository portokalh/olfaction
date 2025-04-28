# === LOAD LIBRARIES ===
library(readxl)
library(dplyr)
library(stringr)
library(gprofiler2)
library(ggplot2)
library(tidyr)
library(igraph)
library(ggraph)

# === PARAMETERS ===
input_file <- "eigengene_top_genes.xlsx"
output_folder <- "eigengene_enrichment_results"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# === READ ALL SHEETS ===
sheet_names <- excel_sheets(input_file)

# === FUNCTION TO PLOT ENRICHMENT ===
plot_enrichment <- function(enrichment_df, module_name) {
  top_terms <- enrichment_df %>%
    arrange(p_value) %>%
    slice(1:10)
  
  ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top GO Biological Processes - ", module_name),
         x = "GO Term",
         y = "-log10(p-value)") +
    theme_minimal()
}

# === PROCESS EACH MODULE ===
for (sheet in sheet_names) {
  df <- read_excel(input_file, sheet = sheet)
  top_genes <- df[[1]] %>% head(500)
  
  gost_res <- gost(query = top_genes,
                   organism = "hsapiens",
                   sources = c("GO:BP"),
                   significant = FALSE)
  
  if (!is.null(gost_res) && nrow(gost_res$result) > 0) {
    enrichment_fixed <- gost_res$result %>%
      mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))
    
    output_file <- file.path(output_folder, paste0(sheet, "_enrichment_results.tsv"))
    write.table(enrichment_fixed, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    plot <- plot_enrichment(gost_res$result, sheet)
    plot_file <- file.path(output_folder, paste0(sheet, "_top_GO_terms.png"))
    ggsave(plot_file, plot, width = 8, height = 6)
  }
}

cat("\nAll modules processed! Enrichment results and plots are saved in:", output_folder, "\n")

# === OPTIONAL: Search for olfactory or memory-related pathways ===
keywords <- c("olfactory", "odor", "smell", "memory", "learning", "hippocampus", "cognition")
olf_memory_hits <- list()

for (sheet in sheet_names) {
  result_file <- file.path(output_folder, paste0(sheet, "_enrichment_results.tsv"))
  
  if (file.exists(result_file)) {
    enrichment_data <- read.delim(result_file, stringsAsFactors = FALSE)
    
    hits <- enrichment_data %>%
      filter(grepl(paste(keywords, collapse = "|"), term_name, ignore.case = TRUE)) %>%
      select(term_id, term_name, p_value)
    
    if (nrow(hits) > 0) {
      olf_memory_hits[[sheet]] <- hits
    }
  }
}

if (length(olf_memory_hits) > 0) {
  combined_hits <- bind_rows(olf_memory_hits, .id = "Eigengene_Module")
  
  write.table(combined_hits, file = file.path(output_folder, "olfactory_memory_related_terms.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nOlfactory or memory-related pathways found and saved to olfactory_memory_related_terms.tsv\n")
  
  # === CREATE A CLEAN RANKED BARPLOT ===
  graph_data <- combined_hits %>%
    filter(p_value < 0.05) %>%
    mutate(log_p = -log10(p_value)) %>%
    arrange(p_value) %>%
    slice(1:30) %>%
    distinct(term_name, .keep_all = TRUE)
  
  barplot <- ggplot(graph_data, aes(x = reorder(term_name, log_p), y = log_p)) +
    geom_col(fill = "orange") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top Olfactory/Memory Related Pathways",
         x = "Pathway",
         y = "-log10(p-value)")
  
  ggsave(file.path(output_folder, "olfactory_memory_barplot.png"), barplot, width = 10, height = 8)
  
  cat("\nBarplot of top olfactory/memory pathways saved as olfactory_memory_barplot.png\n")
  
} else {
  cat("\nNo olfactory or memory-related pathways found.\n")
}