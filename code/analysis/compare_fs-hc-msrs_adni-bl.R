#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
#library(stringr)
#library(bootES)
library(ggplot2)
library(GGally)
#library(ggtext)
#library(gridExtra)
#library(rlang)
#library(gtsummary)
#library(dunn.test)

## Remake plots
ReDoPlots   <- TRUE

### Read RDS objects
## ADNIMERGE
fpath       <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_adnimerge-bl.R') |> source()
}

# ICC volume and ScaleFactors
fpath           <- here("data/derivatives/adni_icc_scale.csv")
if (!file.exists(fpath)) {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}
icc_scales      <- fread(fpath)
icc_scales[, SCANDATE := lubridate::ymd(SCANDATE)]

## HCvols measures head-size adjusted
fpath       <- here("data/rds/adni-bl_fs-hc-msrs_icv-adjusted.rds")
if (file.exists(fpath)) {
  volumes   <- read_rds(fpath)
} else {
  here('code/analysis/adjust_fs-hc-msrs_adni-bl.R') |> source()
}

## Effect sizes
fpath       <- here("data/rds/adni-bl_effect-sizes_fs-hc-msrs_dx_labs.rds")
if (file.exists(fpath)) {
  efflabs   <- read_rds(fpath)
} else {
  here('code/analysis/effect-sizes_fs-hc-msrs_adni-bl.R') |> source()
}

## Get clinical label and average two sides
DT          <- icc_scales[, .(PTID, SCANDATE, ICC = ICC / 1000)
                          ][adnimerge[, .(PTID, SCANDATE, DX)
                                      ], on = .(PTID, SCANDATE), .(PTID, DX, ICC)
                          ][volumes[, .(AVG = mean(VAL)),
                                    .(PTID, MEASURE, ADJ)], on = "PTID"]

DT[, ADJ := factor(ADJ, labels = c("Unadj.", "STX", "Proport.",
                                   "Power-corr.", "Residuals"))]

rm(icc_scales, adnimerge, volumes)

### Plots
# Palette
cbPalette   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

diag_fun  <- function(data, mapping, var, labels.dt, measure,...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(alpha = 1) +
    stat_summary(aes(xintercept = ..x.., y = 0), fun = mean,
                 geom = "vline", orientation = "y",
                 linetype = "dashed", alpha = 1) +
    geom_richtext(data = labels.dt[labels.dt$MSR %in% c("ICC", measure) &
                                   labels.dt$ADJ == as_label(mapping$x)],
                  aes(label = LABEL), inherit.aes = FALSE,
                  colour = "Black", fill = "White",
                  size = 3.5, alpha = .8,
                  x = -Inf, y = -Inf,
                  hjust = -0.1, vjust = -0.25)
}

for (measure in DT[, levels(MEASURE)]) {
  fnames      <- sprintf("plots/adni-bl_similarity_fs-%s-adj_dx.%s",
                         stringr::str_to_lower(measure),
                         c("png", "tiff")) |> here()

  if (!file.exists(fnames[1]) || !file.exists(fnames[2]) || ReDoPlots) {
    dt        <- DT[MEASURE == measure, -"MEASURE"] |>
                dcast(... ~ ADJ, value.var = "AVG")
    g <- ggpairs(dt, columns = 3:length(dt),
                 aes(colour = DX, alpha = 0.7),
                 upper = list(continuous = wrap("cor", method = "spearman")),
                 diag = list(continuous = wrap(diag_fun, labels.dt = efflabs,
                                               measure = measure))) +
      theme_classic(base_size = 12) +
      theme(text = element_text(size = 14)) +
      scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
      scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
      labs(title = sprintf("Head-size adjustment: %s (FS)", measure),
           caption = "* p < 0.05; ** p < 0.01; *** p < 0.001")
  }

  if (measure == "HC") {
    width   <- 14
    height  <- 11
  } else {
    width   <- 11
    height  <- 8
  }

  if (!file.exists(fnames[1]) || ReDoPlots) {
    png(fnames[1], width = width, height = height, units = "in", res = 600)
    print(g)
    dev.off()
  }

  #if (!file.exists(fnames[2]) || ReDoPlots) {
    #tiff(fnames[2], width = width, height = height, units = "in", res = 600)
    #print(g)
    #dev.off()
  #}
}
