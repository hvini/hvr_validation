#!/usr/env/bin Rscript

library(officer)
library(jsonlite)
library(here)

# Path to metadata
meta_path <- here("plots/metadata.json")
# Directory containing figures
plots_dir <- here("plots")

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
    plot_path <- file.path(plots_dir, metadata[[i]]$fname)
  if (!file.exists(plot_path)) stop("Missing figure: ", plot_path)

  if (i %in% 1:6) {
    doc1 <- doc1 |>
      body_add_par(
        "Plot titles and captions",
        style = "heading 1"
      ) |>
      body_add_par(
        paste0("Fig. ", i, ": ", metadata[[i]]$title),
        style = "heading 2"
      ) |>
      #body_add_img(
        #plot_path,
        #width = metadata[[i]]$width,
        #height = metadata[[i]]$height
      #) |>
      body_add_par(metadata[[i]]$caption, style = "Normal") |>
      body_add_par("", style = "Normal") # blank line for spacing
  } else {
    doc2 <- doc2 |>
      body_add_par(
        "Supplementary plot titles and captions",
        style = "heading 1"
      ) |>
      body_add_par(
        paste0("Fig. S", i - 6, ": ", metadata[[i]]$title),
        style = "heading 2"
      ) |>
      #body_add_img(
        #plot_path,
        #width = metadata[[i]]$width,
        #height = metadata[[i]]$height
      #) |>
      body_add_par(metadata[[i]]$caption, style = "Normal") |>
      body_add_par("", style = "Normal") # blank line for spacing
  }
}

# Save the final documents
print(doc1, target = here("plots/plots_captions.docx"))
print(doc2, target = here("plots/plots_supp_captions.docx"))
