#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
#library(stringr)
library(bootES)
library(ggplot2)
#library(GGally)
library(ggtext)
library(gridExtra)
#library(rlang)
#library(gtsummary)
#library(dunn.test)

## Remake plots
ReDoPlots   <- FALSE
ReRunSims   <- FALSE

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
fpath       <- here("data/rds/adni-bl_hc-msrs_icv-adjusted.rds")
if (file.exists(fpath)) {
  volumes   <- read_rds(fpath)
} else {
  here('code/analysis/adjust_hc-msrs_adni-bl.R') |> source()
}

## Get clinical label and average two sides
DT          <- icc_scales[, .(PTID, SCANDATE, ICC)
                          ][adnimerge[, .(PTID, SCANDATE, DX)
                                      ], on = .(PTID, SCANDATE), .(PTID, DX, ICC)
                          ][volumes[, .(AVG = mean(VAL)),
                                    .(PTID, MEASURE, ADJ)], on = "PTID"]

#volumes     <- adnimerge[, .(PTID, DX)
                         #][volumes[, .(AVG = mean(VAL)),
                                   #.(PTID, MEASURE, ADJ)], on = "PTID"]
rm(icc_scales, adnimerge, volumes)

### HC volume
#hcv.dt        <- dcast(volumes[!is.na(METHOD), -"HVR"],
                       #... ~ METHOD, value.var = "HC")

#setnames(hcv.dt,
         #c("cnn", "malf", "nlpb", "fs6"),
         #c("CNN", "MALF", "NLPB", "FS_V6"))

#setcolorder(hcv.dt,
            #c("PTID", "DX", "CNN", "NLPB", "MALF", "FS_V4_V5", "FS_V6"))

#hcv.dt.long   <- melt(hcv.dt, id.vars = c(1:2),
                      #variable.name = "METHOD", value.name = "HCV")
#hcv.dt.long   <- hcv.dt.long[!is.na(HCV)]

### HVR
#hvr.dt        <- dcast(volumes[!is.na(METHOD), -"HC"],
                       #... ~ METHOD, value.var = "HVR")

#hvr.dt[, FS_V4_V5 := NULL]
#setnames(hvr.dt,
         #c("cnn", "malf", "nlpb", "fs6"),
         #c("CNN", "MALF", "NLPB", "FS_V6"))

#setcolorder(hvr.dt,
            #c("PTID", "DX", "CNN", "NLPB", "MALF", "FS_V6"))

#hvr.dt.long   <- melt(hvr.dt, id.vars = c(1:2),
                      #variable.name = "METHOD", value.name = "HVR")
#hvr.dt.long   <- hvr.dt.long[!is.na(HVR)]
#hvr.dt.long[, METHOD := factor(METHOD)]

### Effect sizes
msrs        <- DT[, levels(MEASURE)]
adjs        <- DT[, levels(ADJ)]
dxs         <- DT[, levels(DX)][-2] # Focus on CH-AD difference

## Glass' delta (CH sd only)
## CH vs AD
fnames      <- sprintf("data/rds/adni-bl_effect-sizes_hc-msrs_dx%s.rds",
                       c("", "_labs", "_sims")) |> here()
if (all(file.exists(fnames), !ReRunSims)) {
  effvals <- read_rds(fnames[1])
  efflabs <- read_rds(fnames[2])
  effsims <- read_rds(fnames[3])
} else {
  effvals <- expand.grid(msrs, adjs) |> as.data.table()
  setnames(effvals, c("MSR", "ADJ"))
  #effs  <- bounds_l <- bounds_h <- vector()
  sims  <- vector("list", effvals[, .N + 1])
  names(sims) <- c(effvals[, paste(MSR, ADJ, sep = "_")], "ICC")
  # Adjustment-methods by HC-measure
  for (i in 1:effvals[, .N]) {
    dt  <- DT[MEASURE == effvals[i, MSR] &
              ADJ == effvals[i, ADJ] &
              DX %in% dxs]

    if (dt[, .N == 0]) next

    effect <- bootES(dt,
                     data.col       = "AVG",
                     group.col      = "DX",
                     contrast       = c("CH", "AD"),
                     effect.type    = "cohens.d",
                     glass.control  = "CH")

    sims[[i]] <- effect$t
    effvals[i, `:=`(EFFECT  = round(effect$t0, 2),
                    BOUNDS_l = round(effect$bounds[1], 2),
                    BOUNDS_h = round(effect$bounds[2], 2))]
  }
  # ICC for comparison
  effect <- DT[, 1:3] |>
              unique() |>
              bootES(data.col       = "ICC",
                     group.col      = "DX",
                     contrast       = c("CH", "AD"),
                     effect.type    = "cohens.d",
                     glass.control  = "CH")

  sims[[effvals[, .N + 1]]] <- effect$t
  effvals <- rbind(effvals,
                   data.table(MSR       = "ICC",
                              ADJ       = "icc",
                              EFFECT    = round(effect$t0, 2),
                              BOUNDS_l  = round(effect$bounds[1], 2),
                              BOUNDS_h  = round(effect$bounds[2], 2)))
  rm(effect)

  effsims <- as.data.table(sims)
  na_cols <- names(effsims)[sapply(effsims, \(x) all(is.na(x)))]
  effsims <- effsims[, -..na_cols]
  rm(na_cols)

  effvals <- effvals[!is.na(EFFECT)]
  efflabs <- effvals[, .(DX     = NA,
                         MSR,
                         ADJ    = factor(ADJ, labels = c("Unadj.", "STX",
                                                         "Proport.",
                                                         "Power-corr.",
                                                         "Residuals", "ICC")),
                         LABEL  = paste0("&Delta; = ", EFFECT,
                                         " [", BOUNDS_l, ", ", BOUNDS_h, "]"))]

  write_rds(effvals, fnames[1])
  write_rds(efflabs, fnames[2])
  write_rds(effsims, fnames[3])
}

### Plots
## Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Effect-size only plots
DT[, ADJ := factor(ADJ, labels = c("Unadj.", "STX", "Proport.",
                                   "Power-corr.", "Residuals"))]

p1 <- ggplot(aes(x = AVG, colour = DX), data = DT[MEASURE == "HC"]) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = DT[MEASURE == "HC", mean(AVG), .(DX, ADJ)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .9) +
  facet_wrap(facets = vars(ADJ), scales = "free", nrow = 1) +
  geom_richtext(data = efflabs[MSR == "HC"],
                aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 3, alpha = .8,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL, colour = NULL, subtitle = "Hippocampal volume")

p2 <- ggplot(aes(x = AVG, colour = DX), data = DT[MEASURE == "HVR"]) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = DT[MEASURE == "HVR", mean(AVG), .(DX, ADJ)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .9) +
  facet_wrap(facets = vars(ADJ), scales = "free", nrow = 1) +
  geom_richtext(data = efflabs[MSR == "HVR"],
                aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 3, alpha = .8,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL, colour = NULL,
       subtitle = "Hippocampal-to-Ventricle ratio")

p3 <- ggplot(aes(x = AVG, colour = DX), data = DT[MEASURE == "HAVR"]) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = DT[MEASURE == "HAVR", mean(AVG), .(DX, ADJ)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .9) +
  facet_wrap(facets = vars(ADJ), scales = "free", nrow = 1) +
  geom_richtext(data = efflabs[MSR == "HAVR"],
                aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 3, alpha = .8,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL, colour = NULL,
       subtitle = "Hippocampal&Amygdala-to-Ventricle ratio")

p4 <- ggplot(aes(x = AVG, colour = DX), data = DT[MEASURE == "HAVAS"]) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = DT[MEASURE == "HAVAS", mean(AVG), .(DX, ADJ)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .9) +
  facet_wrap(facets = vars(ADJ), scales = "free", nrow = 1) +
  geom_richtext(data = efflabs[MSR == "HAVAS"],
                aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 3, alpha = .8,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL, colour = NULL,
       subtitle = "Hippocampal & Amygdala & Ventricle (~HAVAS)")

p <- grid.arrange(p1, p2, p3, p4, nrow = 4)
fnames <- here(paste("plots/adni-bl_hc-msrs_effsizes-dx",
                     c("png", "tiff"), sep = "."))

if (!file.exists(fnames[1]) || ReDoPlots) {
  fnames[1] |>
  ggsave(p, width = 13, height = 10, units = "in", dpi = 600)
}

#if (!file.exists(fnames[2]) || ReDoPlots) {
  #fnames[2] |>
  #ggsave(p, width = 13, height = 10, units = "in", device = "tiff", dpi = 600)
#}
