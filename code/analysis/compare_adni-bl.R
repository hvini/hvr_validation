#!/usr/env/bin Rscript

library(here)
library(data.table)
library(bootES)
library(ggplot2)
library(GGally)
library(ggtext)
library(gt)
#library(rlang)
#library(dunn.test)

### CONSTANTS
REDOTABLE <- FALSE
REDOPLOTS <- FALSE
RERUNSIMS <- FALSE

### INPUT
fpaths <- list(
  RDS = c(
    "adnimerge_baseline",
    "adni-bl_volumes_freesurfer",
    "adni-bl_volumes_icv-adjusted"
  ) |> sprintf(fmt = "data/rds/%s.rds") |> here(),
  SRC = c(
    "data_parsing/parse_adnimerge-bl",
    "data_parsing/parse_freesurfer-vols",
    "analysis/adjust_hc-hvr_adni-bl"
  ) |> sprintf(fmt = "code/%s.R") |> here()
)

## Parse DATA
## Adnimerge
if (file.exists(fpaths$RDS[1])) {
  adnimerge <- readRDS(fpaths$RDS[1])
} else {
  source(fpaths$SRC[1])
  adnimerge <- data.lst$ADNIMERGE
  rm(data.lst)
}

## FreeSurfer volumes
if (file.exists(fpaths$RDS[2])) {
  fs_vols <- readRDS(fpaths$RDS[2])
} else {
  source(fpaths$SRC[2])
  fs_vols <- data.lst$FS
  rm(data.lst)
}

# Adjusted HC & HVR volumes
if (file.exists(fpaths$RDS[3])) {
  vols.lst <- readRDS(fpaths$RDS[3])
} else {
  source(fpaths$SRC[1])
  vols.lst <- data.lst$ADJ
  rm(data.lst)
}

data.lst <- list(
  ADNIMERGE = adnimerge[, DX, .(PTID, SCANDATE)],
  FS = fs_vols[!is.na(FS_ucsf), FS_ucsf, .(PTID, SCANDATE)],
  HC = vols.lst$STX[, .(SCALEFACTOR, HC_mean), .(PTID, METHOD)],
  HVR = vols.lst$HVR[, HVR_mean, .(PTID, METHOD)]
)

rm(fpaths, adnimerge, fs_vols, vols.lst)

### Data CLEANING
data.dt <- data.lst$ADNIMERGE |>
  merge(data.lst$FS, all = TRUE) |>
  merge(data.lst$HC) |>
  merge(data.lst$HVR, by = c("PTID", "METHOD")) |>
  (
    \(DT)
    DT[, let(
      FS_ucsf = FS_ucsf * SCALEFACTOR / 2000,
      METHOD = factor(
        METHOD,
        levels = c("cnn", "nlpb", "malf", "fs6"),
        labels = c("CNN", "NLPB", "MALF", "FS_V6")
      ),
      SCALEFACTOR = NULL)
    ]
  )() |>
  setnames(c("HC_mean", "HVR_mean"), c("HC", "HVR"))

data.dt <- data.dt[, .(METHOD = "FS_V4_V5", HC = FS_ucsf), .(PTID, DX)] |>
  unique() |>
  rbind(data.dt[, -"FS_ucsf"], fill = TRUE) |>
  setorder(DX, METHOD)

### TABLE summary
fname <- "adni-bl_table-2.tex"
fpath <- here("tables")
if (!file.exists(here(fpath, fname)) | REDOTABLE) {
  data.dt |>
    melt(measure = c("HC", "HVR")) |>
    na.omit() |>
    {function(DT)
      DT[
        ,
        .(value = sprintf( "%.2f (%.2f)", mean(value), sd(value))),
        by = DX:variable
      ]
    }() |>
    rbind(
      data.dt[
        is.na(HC),
        .(variable = "Fail", value = .N),
        by = .(METHOD, DX)
      ]
    ) |>
    {function(DT) {
      DT[
        ,
        variable := fcase(
          variable == "HC", "HC",
          variable == "HVR", "HVR",
          METHOD == "FS_V4_V5", "M",
          default = "F"
        ) |>
        factor(
          levels = c("HC", "HVR", "M", "F"),
          labels = c("Hippocampus", "HVR", "Missing", "Failures")
        )
      ]
      DT[
        data.dt[, .(N = .N - .SD[is.na(HC), .N]), METHOD],
        on = "METHOD",
        METHOD := sprintf(
          "**%s**, N: %s",
          fcase(
            METHOD == "FS_V4_V5", "FreeSurfer v4.3 & v5.1",
            METHOD == "FS_V6", "FreeSurfer v6.0",
            default = as.character(METHOD)
          ),
          format(N, big.mark = ","))
      ]
    }}() |>
    dcast(METHOD + variable ~ DX, value.var = "value") |>
    gt(
      rowname_col = "variable",
      groupname_col = "METHOD",
      process_md = TRUE
    ) |>
    tab_spanner(
      label = "Clinical Label",
      columns = c("CH", "MCI", "AD")
    ) |>
    tab_options(
      latex.tbl.pos = "h",
      footnotes.multiline = FALSE
    ) |>
    cols_align("center", columns = c("CH", "MCI", "AD")) |>
    cols_label(
      CH = "**CH**, N: %i" |>
        sprintf(data.dt[!duplicated(PTID)]["CH", on = "DX", .N]) |>
        md(),
      MCI = "**MCI**, N: %i" |>
        sprintf(data.dt[!duplicated(PTID)]["MCI", on = "DX", .N]) |>
        md(),
      AD = "**AD**, N: %i" |>
        sprintf(data.dt[!duplicated(PTID)]["AD", on = "DX", .N]) |>
        md()
    ) |>
    tab_stub_indent(rows = everything(), indent = 2) |>
    tab_footnote(
      footnote = "Mean (standard deviation).",
      locations = cells_stub(rows = contains("H"))
    ) |>
    tab_footnote(
      footnote = "Number of failed segmentations.",
      locations = cells_stub(rows = matches("Failures"))
    ) |>
    tab_footnote(
      footnote = "Number of unreported cases from ADNI.",
      locations = cells_stub(rows = matches("Missing"))
    ) |>
    gtsave(filename = fname, path = fpath)
}
rm(fname, fpath)


### EFFECT sizes
# HC volume (sum of sides??)
# HVR (average of sides)

# CH vs AD
# Glass' delta (CH sd only)
mtds  <- data.dt[, levels(METHOD)] # FS_V4_V5 does not have HVR
dxs   <- data.dt[, levels(DX)][-2] # Focus on CH-AD difference
effsizes.lst <- list()
for (roi in c("HC", "HVR")) {
  fnames <- "%s/adni-bl_effect-sizes_%s_dx%s.rds" |>
    sprintf("data/rds", tolower(roi), c("", "_sims")) |>
    here()
  if (all(file.exists(fnames), !RERUNSIMS)) {
    effsizes.lst[[roi]] <- list(
      VALS = readRDS(fnames[1]),
      SIMS = readRDS(fnames[2])
    )
  } else {
    effsizes.lst[[roi]] <- list(
      VALS = list(),
      SIMS = list()
    )

    effs <- ci_l <- ci_h <- vector()
    if (roi == "HVR") mtds <- setdiff(mtds, "FS_V4_V5")
    sims <- vector("list", length(mtds))
    names(sims) <- mtds
    for (mtd in mtds) {
      effect <- bootES(
        data.dt[mtd, on = "METHOD"][dxs, on = "DX"][!is.na(get(roi))],
        data.col       = roi,
        group.col      = "DX",
        contrast       = c("CH", "AD"),
        effect.type    = "cohens.d",
        glass.control  = "CH"
      )
      sims[[mtd]] <- effect$t
      effs        <- c(effs, effect$t0)
      ci_l        <- c(ci_l, effect$bounds[1])
      ci_h        <- c(ci_h, effect$bounds[2])
    }

    effsizes.lst[[roi]][["VALS"]] <- data.table(
      PTGENDER  = NA,
      DX        = NA,
      METHOD    = mtds,
      EFFECT    = effs,
      ci_l      = ci_l,
      ci_h      = ci_h,
      LABEL     = sprintf("&Delta; = %.2f [%.2f, %.2f]", effs, ci_l, ci_h)
    )
    saveRDS(effsizes.lst[[roi]][["VALS"]], fnames[1])

    effsizes.lst[[roi]][["SIMS"]] <- as.data.table(sims)
    saveRDS(effsizes.lst[[roi]][["SIMS"]], fnames[2])

    rm(effs, ci_l, ci_h, sims, effect)
  }
}
rm(mtds, dxs, roi, fnames)

### PLOTS
## Palette
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# Diagonal plots: add mean vertical lines and effect sizes
diag_fun  <- function(data, mapping, var, labels.dt,...) {
  ## Calculate density first to get the y-axis limits
  p <- ggplot(mapping = mapping) + geom_density(data = data, alpha = .1)
  y_limits <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
  ## Add mean lines (complete data) and labels (effs_labs)
  p +
    theme(text = element_text(size = 10)) +
    stat_summary(
    aes(xintercept = ..x.., y = 0),
    data = data,
    fun = mean,
    geom = "vline",
    orientation = "y",
    linetype = "dashed",
    linewidth = .6,
    alpha = 1
  ) + geom_richtext(
    aes(label = LABEL),
    data = labels.dt[labels.dt$METHOD == as_label(mapping$x)],
    size = 2.5,
    colour = "Black",
    fill = "White",
    alpha = 1,
    x = -Inf,
    y = y_limits[2] + y_limits[2] * .17,
    hjust = -0.25
  ) + scale_y_continuous(limits = c(y_limits[1], y_limits[2] * 1.3))
}

# HC By DX
outdir <- here("plots")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
for (roi in c("HC", "HVR")) {
  fnames <- sprintf(
    "%s/adni-bl_similarity_%s.%s",
    outdir,
    tolower(roi),
    c("png", "tiff")
  )

  if (any(!file.exists(fnames), REDOPLOTS)) {
    p <- data.dt |>
      melt(id = c("PTID", "DX", "METHOD")) |>
      (\(DT) DT[roi, on = "variable", -"variable"])() |>
      na.omit() |>
      dcast(... ~ METHOD) |>
      setcolorder(c("MALF", "NLPB", "CNN"), after = "DX") |>
      ggpairs(
        columns = 3:ifelse(roi == "HC", 7, 6),
        mapping = aes(colour = DX, alpha = 0.7),
        upper = list(
          continuous = wrap("cor", method = "spearman", stars = FALSE)
        ),
        diag = list(
          continuous = wrap(diag_fun, labels.dt = effsizes.lst[[roi]]$VALS)
        ),
        lower = list(
          continuous = wrap("points", size = .65, shape = 21)
        )
      ) +
      theme_classic(base_size = 10) +
      theme(text = element_text(size = 10)) +
      scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
      scale_colour_manual(values = cbPalette[c(2:3, 8)])
  }

  if (any(!file.exists(fnames[1]), REDOPLOTS)) {
    png(
      fnames[1],
      width = ifelse(roi == "HC", 8, 7),
      height = ifelse(roi == "HC", 6, 5),
      units = "in",
      res = 600
    )
    print(p)
    dev.off()
  }

  if (any(!file.exists(fnames[2]), REDOPLOTS)) {
    tiff(
      fnames[2],
      width = ifelse(roi == "HC", 8, 7),
      height = ifelse(roi == "HC", 6, 5),
      units = "in",
      res = 600
    )
    print(p)
    dev.off()
  }
}
rm(roi, fnames)
