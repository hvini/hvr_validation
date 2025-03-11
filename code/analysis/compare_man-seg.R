#!/usr/env/bin Rscript

library(here)
library(data.table)
library(WRS2)
library(progress)
library(ggplot2)
library(ggsignif)
library(ggtext)

### CONSTANT
REDOPLOTS <- FALSE

### INPUT
segm.lst <- valid.lst <- list()
## Trained segmentation methods: MALF & NLPB & CNN
for (segm in c("malf", "nlpb", "cnn")) {
  fpath     <- "data/derivatives/man-seg_kappa_hcvc_%s.csv" |>
    sprintf(segm) |>
    here()

  if (!file.exists(fpath)) {
    fpath |>
      sprintf(fmt = "File: %s is required but could not be found.") |>
      stop()
  }
  segm.lst[[toupper(segm)]] <- fread(fpath)
}

## Validation datasets: ICBM & ADNI
for (dataset in c("icbm", "adni")) {
  fpath   <- "data/derivatives/man-seg_kappa_hc_cnn_%s.csv" |>
    sprintf(dataset) |>
    here()

  if (!file.exists(fpath)) {
    fpath |>
      sprintf(fmt = "File: %s is required but could not be found.") |>
      stop()
  }
  valid.lst[[toupper(dataset)]] <- fread(fpath)
}

## ADNI groups — 20/20/20 NC/MCI/AD
fpath <- here("data/adni_validation_groups.csv")
if (!file.exists(fpath)) fpath |>
  sprintf(fmt = "File: %s is required but could not be found.") |>
  stop()

adni_groups.dt <- fread(fpath) |> setnames(c("ID", "GROUP"))

rm(segm, dataset, fpath)

### Data CLEANING
data.lst <- list()
## Trained segmentations
Map(
  \(Name, DT) {
    DT[, let(
      segm  = factor(Name),
      id    = factor(id),
      roi   = roi |> toupper() |> factor(),
      side  = side |> substr(1, 1) |> toupper() |> factor()
    )]
    DT |> setcolorder("segm") |> setnames(toupper)
  },
  names(segm.lst),
  segm.lst
) |> invisible()

data.lst$Training <- segm.lst |>
  rbindlist() |>
  melt(id = 1:4, variable = "MSR", value = "VAL")

## Validation datasets
Map(
  \(Name, DT) {
    DT[, let(
      dataset = factor(Name),
      id      = factor(id),
      side    = hc |> substr(4, 4) |> toupper() |> factor()
    )][, hc := NULL]
    DT |> setcolorder(c("dataset", "id", "side")) |> setnames(toupper)
  },
  names(valid.lst),
  valid.lst
) |> invisible()

## Add results of training variables for comparison
data.lst$Validation <- valid.lst |>
  rbindlist() |>
  rbind(segm.lst$CNN["HC", on = "ROI", -c("SEGM", "ROI")], fill = TRUE) |>
  (\(DT) DT[is.na(DATASET), DATASET := "TRAINING"])() |>
  melt(id = 1:3, variable = "MSR", value = "VAL")
rm(segm.lst, valid.lst)


### ANALYSIS
## 2-way ANOVA w/trimmed means
comparisons.lst <- list()
pb <- progress_bar$new(
  format = "Comparisons | :what [:bar] :current/:total",
  total = (
    # Training
    data.lst$Training$ROI |> levels() |> length() *
    data.lst$Training$MSR |> levels() |> length() +
    # Validation
    data.lst$Validation$MSR |> levels() |> length()
  ),
  clear = FALSE,
  width = 75
)

## Trained segmentations
comparisons.lst$Training          <- list()
comparisons.lst$Training$ANOVA    <- list()
comparisons.lst$Training$PostHoc  <- list()
for (roi in levels(data.lst$Training$ROI)) {
  for (msr in levels(data.lst$Training$MSR)) {
    pb$tick(
      tokens = list(what = sprintf("Training :: %s: %s", roi, msr))
    )

    aov_res <- tryCatch({
      t2way(
        VAL ~ SEGM * SIDE,
        data.lst$Training[roi, on = "ROI"][msr, on = "MSR"]
      )
    },
      error = \(e) {
        cat(sprintf(
          "Error in t2way for %s - %s: %s\n", roi, msr, e$message))
        NULL
    })

    if (!is.null(aov_res)) {
      identifier <- paste(roi, msr, sep = "_")
      comparisons.lst$Training$ANOVA[[identifier]] <- data.table(
        ROI = roi,
        MSR = msr,
        COMP = c("SEGM", "SIDE", "SEGM:SIDE"),
        STAT = as.numeric(aov_res[c("Qa", "Qb", "Qab")]),
        Pval = as.numeric(aov_res[c("A.p.value", "B.p.value", "AB.p.value")])
      )

      # If ANOVA ran, perform the relevant Post-hoc analysis
      # SEGM1 : MALFvNLPB; SEGM2 : MALFvCNN; SEGM3 : NLPBvCNN
      posthoc <- mcp2atm(
        VAL ~ SEGM * SIDE,
        data.lst$Training[roi, on = "ROI"][msr, on = "MSR"]
      )
      comparisons.lst$Training$PostHoc[[identifier]] <- data.table(
        ROI = roi,
        MSR = msr,
        CONTRAST = names(posthoc$contrasts),
        PSIHAT = posthoc$effects |> lapply(\(x) x$psihat) |> unlist(),
        CI_l = posthoc$effects |>
          lapply(\(x) if (is.null(names(x$conf.int))) {
            x$conf.int[,1]
          } else {
            x$conf.int[1]
          }) |> unlist(),
        CI_h = posthoc$effects |>
          lapply(\(x) if (is.null(names(x$conf.int))) {
            x$conf.int[,2]
          } else {
            x$conf.int[2]
          }) |> unlist(),
        Pval = posthoc$effects |> lapply(\(x) x$p.value) |> unlist()
      )
      rm(identifier, posthoc)
    }
    rm(aov_res)
  }
  rm(msr)
}
rm(roi)

## Validation datasets
comparisons.lst$Validation          <- list()
comparisons.lst$Validation$ANOVA    <- list()
comparisons.lst$Validation$PostHoc  <- list()
for (msr in levels(data.lst$Validation$MSR)) {
  pb$tick(tokens = list(what = sprintf("Validation :: %s", msr)))

  aov_res <- tryCatch({
    t2way(VAL ~ DATASET * SIDE, data.lst$Validation[msr, on = "MSR"])
  },
    error = \(e) {
      cat(sprintf("Error in t2way for %s (validation): %s\n", msr, e$message))
      NULL
  })

  if (!is.null(aov_res)) {
    comparisons.lst$Validation$ANOVA[[msr]] <- data.table(
      ROI = "HC",
      MSR = msr,
      COMP = c("DATASET", "SIDE", "DATASET:SIDE"),
      STAT = as.numeric(aov_res[c("Qa", "Qb", "Qab")]),
      Pval = as.numeric(aov_res[c("A.p.value", "B.p.value", "AB.p.value")])
    )

    # If ANOVA ran, perform the relevant Post-hoc analysis
    # DT1 : ICBMvADNI ; DT2 : ICBMvTraining; DT3 : ADNIvTraining
    posthoc <- mcp2atm(
      VAL ~ DATASET * SIDE,
      data.lst$Validation[msr, on = "MSR"]
    )
    comparisons.lst$Validation$PostHoc[[msr]] <- data.table(
      ROI = "HC",
      MSR = msr,
      CONTRAST = names(posthoc$contrasts),
      PSIHAT = posthoc$effects |> lapply(\(x) x$psihat) |> unlist(),
      CI_l = posthoc$effects |>
        lapply(\(x) if (is.null(names(x$conf.int))) {
          x$conf.int[,1]
        } else {
          x$conf.int[1]
        }) |> unlist(),
      CI_h = posthoc$effects |>
        lapply(\(x) if (is.null(names(x$conf.int))) {
          x$conf.int[,2]
        } else {
          x$conf.int[2]
        }) |> unlist(),
      Pval = posthoc$effects |> lapply(\(x) x$p.value) |> unlist()
    )
    rm(posthoc)
  }
  rm(aov_res)
}
rm(pb, msr)

comparisons.lst <- lapply(comparisons.lst, \(level1) {
  Map(
    \(name, level2) {
      if (name == "ANOVA") {
        rbindlist(level2) |>
        (\(DT) DT[, P_adj := p.adjust(Pval, method = "BH", n = 2 * .N)])()
      } else {
        rbindlist(level2) |>
        (\(DT) DT[!CONTRAST %like% "SIDE"])()
      }
    },
    names(level1),
    level1
  )
})

## Create manual significance data.frame for plotting with ggsignif
## Use only PostHoc comparisons,
## i.e. segm/dset vs segm/dset (controlling for side)
annotation.lst <- list(
  Training = comparisons.lst$Training$PostHoc[
    !MSR %like% "ACC|SPEC",
  ][
    Pval < 0.05,
    .(
      ROI = ROI |> factor(labels = c("Hippocampus", "Lateral Ventricles")),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      # SEGM1 : NLPBvMALF; SEGM2 : CNNvMALF; SEGM3 : CNNvNLPB
      x1 = fifelse(CONTRAST %like% 1, "NLPB", "CNN"),
      x2 = fifelse(CONTRAST %like% 3, "NLPB", "MALF"),
      y = fcase(
        CONTRAST %like% 3, .99,
        ROI %like% "H" & CONTRAST %like% 1, 1,
        ROI %like% "H" & CONTRAST %like% 2, 1.015,
        #ROI %like% "V" & CONTRAST %like% 1, 1, ## This is not present
        ROI %like% "V" & CONTRAST %like% 2, 1.03
      ),
      lab = fcase(Pval < 0.001, "***", Pval < 0.01, "**", Pval < 0.05, "*")
    )
  ],
  Validation = comparisons.lst$Validation$PostHoc[
    !MSR %like% "ACC|SPEC",
  ][
    Pval < 0.05,
    .(
      ROI = factor("Hippocampus"),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      # DT1 : ADNIvICBM ; DT2 : TrainingvICBM; DT3 : TrainingvADNI
      x1 = fifelse(CONTRAST %like% 1, "ADNI", "Training"),
      x2 = fifelse(CONTRAST %like% 3, "ADNI", "ICBM"),
      y = fcase(
        CONTRAST %like% 1, .985,
        CONTRAST %like% 2, 1,
        CONTRAST %like% 3, .98
      ),
      lab = fcase(Pval < 0.001, "***", Pval < 0.01, "**", Pval < 0.05, "*")
    )
  ]
)

## Boxplot Kappas HC(VC) Training/Validation
outdir <- here("plots")
fplots <- lapply(
  c("Training", "Validation"),
  \(dset) dset |>
    tolower() |>
    substr(1, 5) |>
    sprintf(fmt = "plots/man-seg_comp_%s.%s", c("png", "tiff")) |>
    here()
)
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Different trained segmentation methods
if (any(REDOPLOTS, !file.exists(fplots[[1]]))) {
  ## Accuracy & Specificity ~ 1.0
  p <- data.lst$Training[
    !MSR %like% "ACC|SPEC",
    .(
      SEGM = SEGM |> factor(levels = c("CNN", "NLPB", "MALF")),
      SIDE = SIDE |> factor(labels = c("Left", "Right")),
      ROI = ROI |> factor(labels = c("Hippocampus", "Lateral Ventricles")),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      VAL
    )
  ] |>
    ggplot(mapping = aes(x = SEGM, y = VAL)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_violin(
      mapping = aes(fill = SIDE),
      position = position_dodge(width = 0.95),
      scale = "width",
      linewidth = 0.3,
      alpha = .02
    ) +
    stat_summary(
      mapping = aes(colour = SIDE),
      position = position_dodge(width = 0.95),
      fun.data = "median_hilow",
      geom = "pointrange",
      size = 0.01,
      linewidth = 0.2
    ) +
    geom_signif(
      data = annotation.lst$Training,
      mapping = aes(xmin = x1, xmax = x2, annotations = lab, y_position = y),
      textsize = 3,
      vjust = .5,
      tip_length = 0,
      manual = TRUE
    ) +
    scale_colour_manual(values = c("darkred", "midnightblue")) +
    scale_fill_manual(values = c("darkred", "midnightblue")) +
    facet_grid(rows = vars(ROI), cols = vars(MSR), scales = "free_y") +
    labs(
      y = "Overlap/Accuracy",
      x = paste(
        "Segmentation method (<span style='color:darkred;'>Left</span> &amp;",
        "<span style='color:midnightblue;'>Right</span>)"
      )
    )

  if (any(REDOPLOTS, !file.exists(fplots[[1]][1]))) {
    ggsave(
      filename = fplots[[1]][1],
      plot = p,
      width = 6,
      height = 5,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[1]][2]))) {
    ggsave(
      filename = fplots[[1]][2],
      plot = p,
      width = 6,
      height = 5,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}

# Validation datasets
if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  ## Accuracy & Specificity > .97
  p <- data.lst$Validation[
    !MSR %like% "ACC|SPEC",
    .(
      DSET = DATASET |> factor(
        levels = c("TRAINING", "ADNI", "ICBM"),
        labels = c("Training", "ADNI", "ICBM")
      ),
      SIDE = SIDE |> factor(labels = c("Left", "Right")),
      ROI = factor("Hippocampus"),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      VAL
    )
  ] |>
    ggplot(mapping = aes(x = DSET, y = VAL)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_violin(
      mapping = aes(fill = SIDE),
      position = position_dodge(width = 0.95),
      scale = "width",
      linewidth = 0.3,
      alpha = .02
    ) +
    stat_summary(
      mapping = aes(colour = SIDE),
      position = position_dodge(width = 0.95),
      fun.data = "median_hilow",
      geom = "pointrange",
      size = 0.01,
      linewidth = 0.2
    ) +
    geom_signif(
      data = annotation.lst$Validation,
      mapping = aes(xmin = x1, xmax = x2, annotations = lab, y_position = y),
      textsize = 3,
      vjust = .5,
      tip_length = 0,
      manual = TRUE
    ) +
    scale_colour_manual(values = c("darkred", "midnightblue")) +
    scale_fill_manual(values = c("darkred", "midnightblue")) +
    facet_grid(rows = vars(ROI), cols = vars(MSR), scales = "free_y") +
    labs(
      y = "Overlap/Accuracy",
      x = paste(
        "Segmentation method (<span style='color:darkred;'>Left</span> &amp;",
        "<span style='color:midnightblue;'>Right</span>)"
      )
    )

  if (any(REDOPLOTS, !file.exists(fplots[[2]][1]))) {
    ggsave(
      filename = fplots[[2]][1],
      plot = p,
      width = 6,
      height = 3,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[2]][2]))) {
    ggsave(
      filename = fplots[[2]][2],
      plot = p,
      width = 6,
      height = 3,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}


### Volumes Correlation
## Read Volume CSVs
vols.lst <- list()

vols.lst["Training"] <- list()
for (segm in c("manual", "malf", "nlpb", "cnn")) {
  fpath     <- "data/derivatives/man-seg_volumes_hcvc_%s.csv" |>
    sprintf(segm) |>
    here()

  if (!file.exists(fpath)) {
    fpath |>
      sprintf(fmt = "File: %s is required but could not be found.") |>
      stop()
  }

  DT <- fpath |>
    fread(drop = c("HC", "CSF")) |>
    melt(measure = patterns(HC = "HC", CSF = "CSF"), variable = "SIDE") |>
    melt(measure = patterns("(HC|CSF)"), variable = "ROI", value = "CC") |>
    (\(DT) DT[
      , let(
        ID = ID |> regexpr(pattern = "\\d{3}") |> regmatches(x = ID),
        SEGM = toupper(segm),
        SIDE = factor(SIDE, labels = c("Left", "Right")),
        ROI = factor(ROI, labels = c("Hippocampus", "Lateral Ventricles")),
        CC = CC / 1000
      )
    ])() |>
    setcolorder(c("ID", "SEGM"))
  vols.lst$Training[[toupper(segm)]] <- copy(DT)
}

vols.lst$Validation <- list()
for (dset in c("adni", "icbm")) {
  fpath     <- "data/derivatives/man-seg_volumes_hc_%s.csv" |>
    sprintf(dset) |>
    here()

  if (!file.exists(fpath)) {
    fpath |>
      sprintf(fmt = "File: %s is required but could not be found.") |>
      stop()
  }

  DT <- fpath |>
    fread(drop = "HC") |>
    melt(measure = patterns("HC"), variable = "SIDE") |>
    (\(DT) {
      DT[
        , let(
          SEGM = fifelse(ID %like% "cnn", "CNN", "MAN"),
          DSET = toupper(dset),
          SIDE = factor(SIDE, labels = c("Left", "Right")),
          ROI = "Hippocampus",
          value = value / 1000
        )
      ]
      DT[, ID := sprintf("hc_%s_(man|cnn)_", dset) |> sub("", ID)]
    })() |>
    dcast(... ~ SEGM, value = "value") |>
    setcolorder(c("ID", "DSET"))
  vols.lst$Validation[[toupper(dset)]] <- copy(DT)
}
rm(fpath, segm, dset, DT)

vols.lst <- lapply(vols.lst, rbindlist)

dt1 <- vols.lst[[1]][SEGM != "MANUAL", .(SEGM, CC), keyby = .(ID, SIDE, ROI)]
dt2 <- vols.lst[[1]][SEGM == "MANUAL", .(MAN = CC), keyby = .(ID, SIDE, ROI)]
vols.lst[[1]] <- dt1[dt2]
rm(dt1, dt2)

## Correlation plots
fplots  <- lapply(
  c("Training", "Validation"),
  \(dset) dset |>
    tolower() |>
    substr(1, 5) |>
    sprintf(fmt = "plots/man-seg_corr_%s.%s", c("png", "tiff")) |>
    here()
)

cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

if (any(REDOPLOTS, !file.exists(fplots[[1]]))) {
  p <- vols.lst$Training |>
    ggplot(aes(x = CC, y = MAN, colour = SEGM)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = cbPalette[1],
      linetype = "dashed"
    ) +
    geom_smooth(method = "lm", alpha = .2, linewidth = .5) +
    geom_point(size = 2, shape = 21) +
    ggpubr::stat_cor(
      aes(label = ..r.label..),
      r.accuracy = 0.001,
      size = 2.7,
      label.x.npc = "right",
      label.y.npc = "bottom",
      hjust = "inward",
      method = "spearman"
    ) +
    facet_grid(rows = vars(ROI), cols = vars(SIDE), scales = "free") +
    scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
    labs(
      y = "Manually segmented volumes",
      x = paste(
        "Volumes segmented by",
        "<span style='color:%s;'>MALF</span>,",
        "<span style='color:%s;'>NLPB</span>, or",
        "<span style='color:%s;'>CNN</span>"
      ) |> sprintf(cbPalette[3], cbPalette[8], cbPalette[2])
    )

  if (any(REDOPLOTS, !file.exists(fplots[[1]][1]))) {
    ggsave(
      filename = fplots[[1]][1],
      plot = p,
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[1]][2]))) {
    ggsave(
      filename = fplots[[1]][2],
      plot = p,
      width = 5,
      height = 5,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}

if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  p <- vols.lst$Validation |>
    ggplot(aes(x = CNN, y = MAN, colour = DSET)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = cbPalette[1],
      linetype = "dashed"
    ) +
    geom_smooth(method = "lm", alpha = .2, linewidth = .5) +
    geom_point(size = 2, shape = 21) +
    ggpubr::stat_cor(
      aes(label = ..r.label..),
      r.accuracy = 0.001,
      size = 2.7,
      label.x.npc = "right",
      label.y.npc = "bottom",
      hjust = "inward",
      method = "spearman"
    ) +
    facet_grid(rows = vars(ROI), cols = vars(SIDE), scales = "free") +
    scale_colour_manual(values = c("darkred", "midnightblue")) +
    labs(
      y = "Manually segmented volumes",
      x = paste(
        "Volumes segmented by CNN",
        "(<span style='color:darkred;'>ADNI</span> &amp;",
        "<span style='color:midnightblue;'>ICBM</span>)"
      )
    )

  if (any(REDOPLOTS, !file.exists(fplots[[2]][1]))) {
    ggsave(
      filename = fplots[[2]][1],
      plot = p,
      width = 5,
      height = 5,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[2]][2]))) {
    ggsave(
      filename = fplots[[2]][2],
      plot = p,
      width = 5,
      height = 5,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}

## Bland-Altman plots
fplots  <- lapply(
  c("Training", "Validation"),
  \(dset) dset |>
    tolower() |>
    substr(1, 5) |>
    sprintf(fmt = "plots/man-seg_blandaltman_%s.%s", c("png", "tiff")) |>
    here()
)

if (any(REDOPLOTS, !file.exists(fplots[[1]]))) {
  DT <- vols.lst[[1]][
    , .(
      AVG = (CC + MAN) / 2,
      DIFF = CC - MAN
    ),
    .(SIDE, ROI, SEGM)
  ][
    , .(
      AVG, DIFF,
      MEAN_DIFF = mean(DIFF),
      CI_low = mean(DIFF) - 1.96 * sd(DIFF),
      CI_high = mean(DIFF) + 1.96 * sd(DIFF)
    ),
    .(ROI, SEGM)
  ]

  p <- ggplot(DT, aes(x = AVG, y = DIFF)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_point(shape = 21, colour = cbPalette[1]) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, SEGM)],
      aes(yintercept = MEAN_DIFF),
      colour = cbPalette[2],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[ROI %like% "Hipp", .SD[1], SEGM],
      aes(x = 1.2, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
      size = 2.5
    ) +
    geom_text(
      data = DT[ROI %like% "Vent", .SD[1], SEGM],
      aes(x = 5.5, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, SEGM)],
      aes(yintercept = CI_low),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[ROI %like% "Hipp", .SD[1], SEGM],
      aes(x = 1.2, y = CI_low, label = round(CI_low, 2)),
      nudge_y = -.100,
      size = 2.5
    ) +
    geom_text(
      data = DT[ROI %like% "Vent", .SD[1], SEGM],
      aes(x = 5.5, y = CI_low, label = round(CI_low, 2)),
      nudge_y = -.15,
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, SEGM)],
      aes(yintercept = CI_high),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[ROI %like% "Hipp", .SD[1], SEGM],
      aes(x = 1.2, y = CI_high, label = round(CI_high, 2)),
      nudge_y = .1,
      size = 2.5
    ) +
    geom_text(
      data = DT[ROI %like% "Vent", .SD[1], SEGM],
      aes(x = 5.5, y = CI_high, label = round(CI_high, 2)),
      nudge_y = .15,
      size = 2.5
    ) +
    xlim(.95, 6.2) +
    facet_grid(rows = vars(ROI), cols = vars(SEGM), scales = "free") +
    labs(
      x = "Mean computed & manual volumes",
      y = "Computed - manual volumes"
    )

  if (any(REDOPLOTS, !file.exists(fplots[[1]][1]))) {
    ggsave(
      filename = fplots[[1]][1],
      plot = p,
      width = 10,
      height = 5,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[1]][2]))) {
    ggsave(
      filename = fplots[[1]][2],
      plot = p,
      width = 10,
      height = 5,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}

if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  DT <- vols.lst[[2]][
    , .(
      AVG = (CNN + MAN) / 2,
      DIFF = CNN - MAN
    ),
    .(SIDE, ROI, DSET)
  ][
    , .(
      AVG, DIFF,
      MEAN_DIFF = mean(DIFF),
      CI_low = mean(DIFF) - 1.96 * sd(DIFF),
      CI_high = mean(DIFF) + 1.96 * sd(DIFF)
    ),
    .(ROI, DSET)
  ]

  p <- ggplot(DT, aes(x = AVG, y = DIFF)) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_markdown(),
      legend.position = "none"
    ) +
    geom_point(shape = 21, colour = cbPalette[1]) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, DSET)],
      aes(yintercept = MEAN_DIFF),
      colour = cbPalette[2],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], DSET],
      aes(x = 1.2, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, DSET)],
      aes(yintercept = CI_low),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], DSET],
      aes(x = 1.2, y = CI_low, label = round(CI_low, 2)),
      nudge_y = -.100,
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, DSET)],
      aes(yintercept = CI_high),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], DSET],
      aes(x = 1.2, y = CI_high, label = round(CI_high, 2)),
      nudge_y = .1,
      size = 2.5
    ) +
    xlim(.95, 6.2) +
    facet_grid(rows = vars(ROI), cols = vars(DSET), scales = "free") +
    labs(
      x = "Mean computed & manual volumes",
      y = "Computed - manual volumes"
    )

  if (any(REDOPLOTS, !file.exists(fplots[[2]][1]))) {
    ggsave(
      filename = fplots[[2]][1],
      plot = p,
      width = 6,
      height = 3,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[2]][2]))) {
    ggsave(
      filename = fplots[[2]][2],
      plot = p,
      width = 6,
      height = 3,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}
