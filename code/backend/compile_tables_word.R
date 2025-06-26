#!/usr/env/bin Rscript

library(gt)
library(officer)
library(gto)
library(jsonlite)
library(here)

# Path to metadata
meta_path <- here("tables/metadata.json")
# Directory containing .rds gt tables
tables_dir <- here("tables")

# Load metadata
metadata <- fromJSON(meta_path)

clean <- latex <- function(x) {
  # Remove math delimiters: $...$, \(...\), \[...\]
  x <- gsub("\\$([^$]+)\\$", "\\1", x)
  x <- gsub("\\\\\\(([^)]+)\\\\\\)", "\\1", x)
  x <- gsub("\\\\\\[([^\\]]+)\\\\\\]", "\\1", x)

  # Remove * delimiters: *...*
  x <- gsub("\\*([^$]+)\\*", "\\1", x)

  # Remove \textit{}, \textbf{}, \text{}, \emph{}, etc.
  x <- gsub("\\\\text(it|bf|normal|sc|tt)?\\{([^}]+)\\}", "\\2", x)
  x <- gsub("\\\\emph\\{([^}]+)\\}", "\\1", x)

  # Remove inline math operators like ~\text{CI}
  x <- gsub("~\\\\text\\{([^}]+)\\}", " \\1", x)
  x <- gsub("\\\\text\\{([^}]+)\\}", "\\1", x)

  # Remove bold/italic markup: \mathbf{}, \mathit{}, etc.
  x <- gsub("\\\\math(it|bf)\\{([^}]+)\\}", "\\2", x)

 # Remove LaTeX special characters: \%, \ <- , \&, etc.
  x <- gsub("\\\\%", "%", x)
  x <- gsub("\\\\_", "_", x)
  x <- gsub("\\\\&", "&", x)
  x <- gsub("\\\\#", "#", x)
  x <- gsub("\\\\,", "", x)
  x <- gsub("\\\\ ", " ", x)

  # Remove unnecessary backslashes
  x <- gsub("\\\\", "", x)

  # Clean up multiple spaces
  x <- gsub(" +", " ", x)

  trimws(x)
}

# Initialize a new Word document
doc1 <- read_docx()
doc2 <- read_docx()

for (i in seq_along(metadata)) {
    fname <- paste0(metadata[[i]]$filename, ".rds")
    gt_path <- file.path(tables_dir, fname)
  if (!file.exists(gt_path)) stop("Missing table: ", gt <- path)

    gt_tbl <- readRDS(gt_path) # This should be a gt object

  if (i %in% 1:2) {
    doc1 <- doc1 |>
      body_add_par(metadata[[i]]$title, style = "heading 2") |>
      body_add_gt(gt_tbl) |>
      body_add_par(metadata[[i]]$caption, style = "Normal") |>
      body_add_par("", style = "Normal") # blank line for spacing
  } else {
    doc2 <- doc2 |>
      body_add_par(metadata[[i]]$title, style = "heading 2") |>
      body_add_gt(gt_tbl) |>
      body_add_par(metadata[[i]]$caption, style = "Normal") |>
      body_add_par("", style = "Normal") # blank line for spacing
  }
}

# Save the final documents
print(doc1, target = here("tables/tables.docx"))
print(doc2, target = here("tables/tables_supp.docx"))
