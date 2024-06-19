#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(bootES)
library(ggplot2)
library(GGally)
library(ggtext)
library(rlang)
library(gtsummary)
library(dunn.test)

## Remake plots
ReDoPlots     <- FALSE
ReRunSims     <- FALSE

## Read RDS objects
# ADNIMERGE
fpath    <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge     <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_adnimerge-bl.R') |> source()
}

fpath    <- here("data/rds/adni-bl_volumes_freesurfer.rds")
if (file.exists(fpath)) {
  fs_vols       <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_freesurfer-vols.R') |> source()
}

fpath    <- here("data/rds/adni-bl_volumes_icv-adjusted.rds")
if (file.exists(fpath)) {
  volumes       <- read_rds(fpath)
} else {
  here('code/analysis/adjust_hc-hvr_adni-bl.R') |> source()
}


## Merge
adni          <- fs_vols[adnimerge, on = .(PTID, SCANDATE)]
volumes       <- adni[volumes, on = "PTID",
                         .(PTID, METHOD, DX,
                           FS_V4_V5 = FS_ucsf * SCALEFACTOR / 2000,
                           #HC = HC_stx_l + HC_stx_r,
                           HC = HC_stx_mean,
                           HVR = HVR_mean)]
rm(adnimerge, fs_vols, adni)

## HC volume
hcv.dt        <- dcast(volumes[!is.na(METHOD), -"HVR"],
                       ... ~ METHOD, value.var = "HC")

setnames(hcv.dt,
         c("cnn", "malf", "nlpb", "fs6"),
         c("CNN", "MALF", "NLPB", "FS_V6"))

setcolorder(hcv.dt,
            c("PTID", "DX", "CNN", "NLPB", "MALF", "FS_V4_V5", "FS_V6"))

hcv.dt.long   <- melt(hcv.dt, id.vars = c(1:2),
                      variable.name = "METHOD", value.name = "HCV")
hcv.dt.long   <- hcv.dt.long[!is.na(HCV)]

# Table
hcv.dt[, -"PTID"] |>
  tbl_summary(by = DX,
              label = list(FS_V4_V5 ~ "FreeSurfer (v4.3 & v5.1)",
                           FS_V6 ~ "FreeSurfer (v6.0)"),
              statistic = all_continuous() ~ "{mean} ({sd})",
              missing_text = "Failures") |>
  modify_header(label ~ "**Automatic Method**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
  add_n() |>
  as_flex_table() |>
  flextable::save_as_docx(path = "data/derivatives/adni-bl_hcv-table.docx")

# Post-hoc differences between Methods by DX
#hcv_dunn.dt <- hcv.dt.long[, dunn.test(HCV, METHOD,
                                       #method = "bonferroni", list = TRUE), DX]

#hcv_dunn.dt[order(DX, comparisons)] |>
  #flextable() |>
  #separate_header() |>
  #colformat_double(digits = 3) |>
  #labelizor(part = "header", labels = c("chi2" = "Chi-squared",
                                        #"P" = "p-value",
                                        #"P.adjusted" = "Adjusted p-value",
                                        #"comparisons" = "Comparisons")) |>
  #bold(part = "header") |>
  ##fix_border_issues() |>
  ##bold(~ P < 0.05, j = 4) |>
  ##bold(~ P.adjusted < 0.05, j = 5) |>
  #autofit() |>
  #save_as_docx(path = "data/derivatives/adni-bl_hcv-posthoc.docx")

## HVR
hvr.dt        <- dcast(volumes[!is.na(METHOD), -"HC"],
                       ... ~ METHOD, value.var = "HVR")

hvr.dt[, FS_V4_V5 := NULL]
setnames(hvr.dt,
         c("cnn", "malf", "nlpb", "fs6"),
         c("CNN", "MALF", "NLPB", "FS_V6"))

setcolorder(hvr.dt,
            c("PTID", "DX", "CNN", "NLPB", "MALF", "FS_V6"))

hvr.dt.long   <- melt(hvr.dt, id.vars = c(1:2),
                      variable.name = "METHOD", value.name = "HVR")
hvr.dt.long   <- hvr.dt.long[!is.na(HVR)]
hvr.dt.long[, METHOD := factor(METHOD)]

# Table
hvr.dt[, -"PTID"] |>
  tbl_summary(by = DX,
              label = FS_V6 ~ "FreeSurfer (v6.0)",
              statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no") |>
  modify_header(label ~ "**Automatic Method**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
  add_n() |>
  as_flex_table() |>
  flextable::save_as_docx(path = "data/derivatives/adni-bl_hvr-table.docx")


## Effect sizes
# HC volume (sum of sides)
mtds  <- hcv.dt.long[, levels(METHOD)]
dxs   <- hcv.dt.long[, levels(DX)][-2] # Focus on CH-AD difference

# CH vs AD
# Glass' delta (CH sd only)
fnames <- here(paste0("data/rds/adni-bl_effect-sizes_hcv_dx",
                     c(".rds", "_sims.rds")))
if (all(file.exists(fnames), !ReRunSims)) {
  effvals_hcv_dx <- read_rds(fnames[1])
  effsims_hcv_dx <- read_rds(fnames[2])
} else {
  effs <- bounds_l <- bounds_h <- vector()
  sims <- vector("list", length(mtds))
  names(sims) <- mtds
  for (mtd in mtds) {
    effect <- bootES(hcv.dt.long[METHOD == mtd & DX %in% dxs],
                     data.col       = "HCV",
                     group.col      = "DX",
                     contrast       = c("CH", "AD"),
                     effect.type    = "cohens.d",
                     glass.control  = "CH")
    sims[[mtd]] <- effect$t
    effs        <- c(effs, effect$t0)
    bounds_l    <- c(bounds_l, effect$bounds[1])
    bounds_h    <- c(bounds_h, effect$bounds[2])
  }

  effsims_hcv_dx <- as.data.table(sims)

  effvals_hcv_dx <- data.table(METHOD    = mtds,
                               EFFECT    = round(effs, 2),
                               BOUNDS_l  = round(bounds_l, 2),
                               BOUNDS_h  = round(bounds_h, 2))

  effvals_hcv_dx[, `:=`(PTGENDER = NA, DX = NA,
                        LABEL = paste0("&Delta; = ", EFFECT,
                                       " [", BOUNDS_l,
                                       ", ", BOUNDS_h, "]"))]
  write_rds(effvals_hcv_dx, fnames[1])
  write_rds(effsims_hcv_dx, fnames[2])
}

# HVR (average of sides)
mtds  <- hvr.dt.long[, levels(METHOD)]
dxs   <- hvr.dt.long[, levels(DX)][-2] # Focus on CH-AD difference

# CH vs AD
# Glass' delta (CH sd only)
fnames <- here(paste0("data/rds/adni-bl_effect-sizes_hvr_dx",
                     c(".rds", "_sims.rds")))
if (all(file.exists(fnames), !ReRunSims)) {
  effvals_hvr_dx <- read_rds(fnames[1])
  effsims_hvr_dx <- read_rds(fnames[2])
} else {
  effs <- bounds_l <- bounds_h <- vector()
  sims <- vector("list", length(mtds))
  names(sims) <- mtds
  for (mtd in mtds) {
    effect <- bootES(hvr.dt.long[METHOD == mtd & DX %in% dxs],
                     data.col       = "HVR",
                     group.col      = "DX",
                     contrast       = c("CH", "AD"),
                     effect.type    = "cohens.d",
                     glass.control  = "CH")
    sims[[mtd]] <- effect$t
    effs        <- c(effs, effect$t0)
    bounds_l    <- c(bounds_l, effect$bounds[1])
    bounds_h    <- c(bounds_h, effect$bounds[2])
  }

  effsims_hvr_dx <- as.data.table(sims)
  effvals_hvr_dx <- data.table(METHOD    = mtds,
                               EFFECT    = round(effs, 2),
                               BOUNDS_l  = round(bounds_l, 2),
                               BOUNDS_h  = round(bounds_h, 2))

  effvals_hvr_dx[, `:=`(PTGENDER = NA, DX = NA,
                        LABEL = paste0("&Delta; = ", EFFECT,
                                       " [", BOUNDS_l,
                                       ", ", BOUNDS_h, "]"))]
  write_rds(effvals_hvr_dx, fnames[1])
  write_rds(effsims_hvr_dx, fnames[2])
}

## Plots
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Effect-size only plots
g <- ggplot(aes(x = HCV, colour = DX), data = hcv.dt.long) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = hcv.dt.long[, mean(HCV), .(DX, METHOD)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .8) +
  facet_wrap(facets = vars(METHOD), scales = "free", nrow = 1) +
  geom_richtext(data = effvals_hcv_dx, aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 4.5,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL)

fnames <- here(paste("plots/adni-bl_hcv_effsizes-dx",
                     c("png", "tiff"), sep = "."))
if (!file.exists(fnames[1]) || ReDoPlots) {
  png(fnames[1], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

if (!file.exists(fnames[2]) || ReDoPlots) {
  tiff(fnames[2], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

g <- ggplot(aes(x = HVR, colour = DX), data = hvr.dt.long) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = hvr.dt.long[, mean(HVR), .(DX, METHOD)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .8) +
  facet_wrap(facets = vars(METHOD), scales = "free", nrow = 1) +
  geom_richtext(data = effvals_hvr_dx, aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 4.5,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL)

fnames <- here(paste("plots/adni-bl_hvr_effsizes-dx",
                     c("png", "tiff"), sep = "."))
if (!file.exists(fnames[1]) || ReDoPlots) {
  png(fnames[1], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

if (!file.exists(fnames[2]) || ReDoPlots) {
  tiff(fnames[2], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

# Both HCV and HVR
hcv_hvr.dt.long <- rbindlist(list(hcv.dt.long[, .(MSR = "HCV", DX, METHOD,
                                                  VAL = HCV)],
                                  hvr.dt.long[, .(MSR = "HVR", DX, METHOD,
                                                  VAL = HVR)]))

effvals_dx <- rbindlist(list(effvals_hcv_dx[, .(MSR = "HCV", METHOD, LABEL)],
                             effvals_hvr_dx[, .(MSR = "HVR", METHOD, LABEL)]))

g <- ggplot(aes(x = VAL, colour = DX), data = hcv_hvr.dt.long) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  geom_histogram(fill = "transparent", bins = 45) +
  geom_vline(data = hcv_hvr.dt.long[, mean(VAL), .(MSR, DX, METHOD)],
             aes(xintercept = V1, colour = DX), linetype = "dashed",
             alpha = .8) +
  facet_wrap(facets = vars(MSR, METHOD), scales = "free") +
  geom_richtext(data = effvals_dx, aes(label = LABEL), inherit.aes = FALSE,
                colour = "Black", fill = "White", size = 3,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25) +
  labs(x = NULL, y = NULL)

fnames <- here(paste("plots/adni-bl_effsizes-dx",
                     c("png", "tiff"), sep = "."))
if (!file.exists(fnames[1]) || ReDoPlots) {
  png(fnames[1], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

if (!file.exists(fnames[2]) || ReDoPlots) {
  tiff(fnames[2], width = 13, height = 5, units = "in", res = 600)
  print(g)
  dev.off()
}

# Diagonal plots: add mean vertical lines and effect sizes
diag_fun  <- function(data, mapping, var, labels.dt,...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(alpha = 1) +
    stat_summary(aes(xintercept = ..x.., y = 0), fun = mean,
                 geom = "vline", orientation = "y",
                 linetype = "dashed", alpha = 1) +
    geom_richtext(data = labels.dt[labels.dt$METHOD == as_label(mapping$x)],
                  aes(label = LABEL), inherit.aes = FALSE,
                  colour = "Black", fill = "White",
                  size = 3.5, x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25)
}

# HC By DX
fnames <- here(paste("plots/adni-bl_similarity_hcv_dx",
                     c("png", "tiff"), sep = "."))
if (!file.exists(fnames[1]) || !file.exists(fnames[2]) || ReDoPlots) {
  g <- ggpairs(hcv.dt, columns = 3:7,
               aes(colour = DX, alpha = 0.7),
               upper = list(continuous = wrap("cor", method = "spearman")),
               diag = list(continuous = wrap(diag_fun,
                                             labels.dt = effvals_hcv_dx))) +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 14)) +
    scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
    scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
    labs(caption = "* p < 0.05; ** p < 0.01; *** p < 0.001")
}

if (!file.exists(fnames[1]) || ReDoPlots) {
  png(fnames[1], width = 13, height = 7, units = "in", res = 600)
  print(g)
  dev.off()
}

if (!file.exists(fnames[2]) || ReDoPlots) {
  tiff(fnames[2], width = 13, height = 7, units = "in", res = 600)
  print(g)
  dev.off()
}

# HVR By DX
fnames <- here(paste("plots/adni-bl_similarity_hvr_dx",
                     c("png", "tiff"), sep = "."))
if (!file.exists(fnames[1]) || !file.exists(fnames[2]) || ReDoPlots) {
  g <- ggpairs(hvr.dt, columns = 3:6,
               aes(colour = DX, alpha = 0.7),
               upper = list(continuous = wrap("cor", method = "spearman")),
               diag = list(continuous = wrap(diag_fun,
                                             labels.dt = effvals_hvr_dx))) +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 14)) +
    scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
    scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
    labs(caption = "* p < 0.05; ** p < 0.01; *** p < 0.001")
}

if (!file.exists(fnames[1]) || ReDoPlots) {
  png(fnames[1], width = 13, height = 7, units = "in", res = 600)
  print(g)
  dev.off()
}

if (!file.exists(fnames[2]) || ReDoPlots) {
  tiff(fnames[2], width = 13, height = 7, units = "in", res = 600)
  print(g)
  dev.off()
}
