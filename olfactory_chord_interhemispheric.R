# ===============================
# Symmetric Olfactory Chord Diagram with Inter-Hemispheric Highlighting and Padding
# ===============================

# Load required libraries
if (!require("circlize")) install.packages("circlize", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("readr")) install.packages("readr", dependencies = TRUE)
library(circlize)
library(dplyr)
library(readr)

pathin <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/data/"
pathout <- "/Users/alex/AlexBadea_MyPapers/olfaction/HSM_v2_041525/results/"
if (!dir.exists(pathout)) dir.create(pathout)

# Load input files
connectivity <- read_csv(paste0(pathin, "olfaction_connectivity_weight.csv"))
roi_annot <- read_csv(paste0(pathin, "olfactory_connectivity_roi_annotations.csv"))
roi_abbrev <- read_csv(paste0(pathin, "olfactory_ROI_abbreviations.csv"))

# Annotate with hemispheres and abbreviations
connectivity <- connectivity %>%
  mutate(
    Hemi1 = ifelse(Hemi1 == "Left", "L", "R"),
    Hemi2 = ifelse(Hemi2 == "Left", "L", "R")
  ) %>%
  left_join(roi_abbrev, by = c("ROI1" = "FullName")) %>%
  rename(ROI1_Abbr = Abbreviation) %>%
  left_join(roi_abbrev, by = c("ROI2" = "FullName")) %>%
  rename(ROI2_Abbr = Abbreviation) %>%
  mutate(
    ROI1_Label = paste(ROI1_Abbr, Hemi1, sep = "_"),
    ROI2_Label = paste(ROI2_Abbr, Hemi2, sep = "_"),
    InterHemispheric = Hemi1 != Hemi2,
    col = case_when(
      InterHemispheric ~ "black",
      Weight < 0 ~ "red",
      TRUE ~ "blue"
    )
  )

# Sort all unique sectors by hemisphere for symmetry
all_nodes <- unique(c(connectivity$ROI1_Label, connectivity$ROI2_Label))
left_nodes <- sort(all_nodes[grepl("_L$", all_nodes)])
right_nodes <- sort(all_nodes[grepl("_R$", all_nodes)])
sector_order <- c(left_nodes, right_nodes)

# Assign colors to nodes by hemisphere
sector_colors <- setNames(
  ifelse(grepl("_L$", sector_order), "#9b59b6", "#1abc9c"),
  sector_order
)

# Filter top 50 strongest connections
top_df <- connectivity %>%
  arrange(desc(abs(Weight))) %>%
  slice(1:50)

# Save plot to PNG with larger margins and text
png(paste0(pathout, "olfactory_chord_interhemispheric.png"), width = 1200, height = 1200, res = 150)

circos.clear()
circos.clear()
circos.par(
  track.margin = c(0.02, 0.02),   # less padding between tracks
  canvas.xlim = c(-1.2, 1.2),     # reduce from 1.5
  canvas.ylim = c(-1.2, 1.2)
)

chordDiagram(
  x = top_df[, c("ROI1_Label", "ROI2_Label", "Weight")],
  order = sector_order,
  grid.col = sector_colors,
  col = top_df$col,
  annotationTrack = "grid",
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  transparency = 0.3
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = get.cell.meta.data("ylim")[1] + 0.5,
      labels = get.cell.meta.data("sector.index"),
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.9
    )
  },
  bg.border = NA
)

#title("Symmetric Olfactory Chord Diagram with Inter-Hemispheric Highlighting")
dev.off()

cat("Diagram saved with more margin and larger text.\n")
