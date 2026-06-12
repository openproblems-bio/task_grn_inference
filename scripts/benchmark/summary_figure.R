library(tidyverse)
library(funkyheatmap)

library(rprojroot)

# Get the directory of the current script
script_dir <- dirname(rprojroot::thisfile())


# Source the helpers.R from the same directory
source(file.path(script_dir, "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)

# Check if a file path is provided
if (length(args) == 0) {
  stop("Please provide a file path to the summary file.")
}

# Read the provided file path from the command line
file_path <- args[1] #summary_all.tsv 
to_save <- args[2] #"figure.pdf"

# Read the TSV file using the provided file path
summary_all <- read_tsv(file_path)

##################################################
# FIGURE 3a
##################################################

# Add the new column for method types
# Using exact surrogate names from Python config:
# METRICS: r_precision, r_recall, vc, sem, ws_precision, ws_recall, t_rec_precision, t_rec_recall, replicate_consistency, tfb_f1, gs_f1
# Their surrogate names: R (precision), R (recall), Virtual cell, SEM, WS (precision), WS (recall), TF recovery (precision), TF recovery (recall), Replica consistency, TF binding, Gene sets
column_info <- bind_rows(
  tribble(
    ~id, ~id_color, ~name, ~group, ~geom, ~palette, ~options,
    "method_name", NA_character_, "Name", "method", "text", NA_character_, list(width = 6, hjust = 0),
    "method_type", NA_character_, "Modality", "method", "text", NA_character_, list(width = 2, hjust = 0),
    "overall_score", "overall_score", "Score", "overall", "bar", "overall", list(width = 4),
    # METRICS with surrogate names
    "Regression (precision)", "Regression (precision)", "Regression (precision)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "Regression (recall)", "Regression (recall)", "Regression (recall)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "Virtual cell", "Virtual cell", "Virtual cell", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "SEM", "SEM", "SEM", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "WS (precision)", "WS (precision)", "WS distance (precision)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "WS (recall)", "WS (recall)", "WS distance (recall)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "TF recovery (precision)", "TF recovery (precision)", "TF recovery (precision)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "TF recovery (recall)", "TF recovery (recall)", "TF recovery (recall)", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "Replicate consistency", "Replicate consistency", "Replicate consistency", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "TF binding", "TF binding", "TF binding", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    "Gene sets recovery", "Gene sets recovery", "Gene sets recovery", "metric_1", "funkyrect",  "metric_1", list(width = 1.5),
    # Datasets
    "OPSCA", "OPSCA", "OPSCA", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "MSCIC", "MSCIC", "MSCIC", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "ParseBioscience", "ParseBioscience", "ParseBioscience", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "300BCG", "300BCG", "300BCG", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "IBD:UC", "IBD:UC", "IBD:UC", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "IBD:CD", "IBD:CD", "IBD:CD", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "Replogle", "Replogle", "Replogle", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "Xaira:HEK293T", "Xaira:HEK293T", "Xaira:HEK293T", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "Xaira:HCT116", "Xaira:HCT116", "Xaira:HCT116", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "Nakatake", "Nakatake", "Nakatake", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "Norman", "Norman", "Norman", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "SoundLife", "SoundLife", "SoundLife", "dataset", "funkyrect", "dataset", list(width = 1.5),
    "SoundLife:Vaccine", "SoundLife:Vaccine", "SoundLife:Vaccine", "dataset", "funkyrect", "dataset", list(width = 1.5),
  ),
  tribble(
    ~id, ~name, ~geom,
    "memory_log", "Peak memory (GB)", "rect",
    "memory_str", "", "text",
    "duration_log", "Duration (hour)", "rect",
    "duration_str", "", "text",
    "complexity_log", "Complexity", "rect",
    "Complexity", "", "text"
  ) %>% mutate(
    group = "resources",
    palette = ifelse(geom == "text", NA_character_, "resources"),
    options = map(geom, function(geom) {
      if (geom == "text") {
        list(overlay = TRUE, size = 2.5)
      } else {
        list(width = 2)
      }
    })
  )
)

# Define the method type mapping
method_type_mapping <- tribble(
  ~method_name, ~method_type,
  "GRNBoost2", "S",
  "GENIE3", "S",
  "PPCOR", "S",
  "Scenic", "S",
  "Portia", "S",
  "scGLUE", "M",
  "CellOracle", "M",
  "Scenic+", "M",
  "FigR", "M",
  "GRaNIE", "M",
  "Positive Ctrl", "C",
  "Negative Ctrl", "C",
  "Pearson Corr.", "C",
  "Spearman Corr.", "S",
  "scPRINT", "F",
  "Geneformer", "F",
  "scGPT", "F"
)


# Include the method types in the summary_all DataFrame
summary_all <- summary_all %>%
  left_join(method_type_mapping, by = "method_name")
print(summary_all)
# Update column groups to include the new "Type" column
column_groups <- tribble(
  ~group, ~palette, ~level1,
  "method", NA_character_, "",
  "overall", "overall", "Overall",
  "metric_1", "metric_1", "Metrics",
  "resources", "resources", "Resources", 
  "dataset", "dataset", "Datasets"
)

# Add palettes for the new column
palettes <- list(
  overall = "Greys",
  metric_1 = "Reds",
  # metric_2 = "Reds",
  resources = "YlOrBr"
)

# Update legends if necessary
legends <- list(
  list(
    title = "Rank",
    palette = "overall",
    geom = "rect",
    labels = c("", "", "worst", "", "", "", "best", ""),
    label_hjust = rep(.5, 8),
    size = c(0, 0, 1, 1, 1, 1, 1, 0)
  ),
  list(
    title = "Scaled score",
    palette = "overall",
    geom = "funkyrect",
    labels = c("", "0", "", "", "", "0.4", "", "0.6", "", "0.8", "", "1"),
    size = c(0, seq(0, 1, by = .1)),
    label_hjust = rep(.5, 12)
  ),
  list(palette = "metric_1", enabled = FALSE),
  list(palette = "dataset", enabled = FALSE),
  # list(palette = "metric_2", enabled = FALSE),
  list(
    title = "Resources",
    palette = "resources",
    geom = "rect",
    labels = c("min", "", "", "", "max"),
    label_hjust = c(0, .5, .5, .5, 1),
    color = colorRampPalette(rev(funkyheatmap:::default_palettes$numerical$YlOrBr))(5),
    size = c(1, 1, 1, 1, 1)
  )
)

# Create the funkyheatmap
g3 <- funky_heatmap(
  data = summary_all,
  column_info = column_info %>% filter(id %in% colnames(summary_all)),
  column_groups = column_groups,
  palettes = palettes,
  position_args = position_arguments(
    expand_xmax = 2,
    col_annot_offset = max(str_length(column_info$name)) / 5
  ),
  add_abc = TRUE,
  scale_column = FALSE,
  legends = legends
) +
  theme(
    text = element_text(family = "Liberation Sans", size = 10), # Change font family and size
    plot.title = element_text(family = "Liberation Sans", face = "bold", size = 12), # Title font customization
    axis.text = element_text(size = 10),  # Axis text customization
    legend.text = element_text(size = 10) # Legend text customization
  )
ggsave(
  paste0(to_save,'.pdf'),
  g3,
  width = g3$width +2,
  height = g3$height + 2
)

ggsave(
  paste0(to_save,'.png'),
  g3,
  width = g3$width +2,
  height = g3$height + 2,
  dpi=300
)
