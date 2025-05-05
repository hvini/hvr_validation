#!/usr/env/bin Rscript

library(here)
library(data.table)
library(WRS2)
library(progress)
library(gt)
library(ggplot2)
library(ggsignif)
library(ggtext)

### CONSTANT
REDOPLOTS <- TRUE

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

## Validation dataset: ADNI
fpath   <- here("data/derivatives/man-seg_kappa_hc_cnn_adni.csv")

if (!file.exists(fpath)) {
  fpath |>
    sprintf(fmt = "File: %s is required but could not be found.") |>
    stop()
}
valid.lst$ADNI <- fread(fpath)

## ADNI groups — 20/20/20 NC/MCI/AD
fpath <- here("data/adni_validation_groups.csv")
if (!file.exists(fpath)) fpath |>
  sprintf(fmt = "File: %s is required but could not be found.") |>
  stop()

adni_groups.dt <- fread(fpath) |>
  setnames(c("ID", "GROUP")) |>
  (\(DT) DT[
    ,
    GROUP := factor(
      GROUP,
      levels = c("NC", "MCI", "AD"),
      labels = c("CH", "MCI", "AD")
    )
  ])()

rm(segm, fpath)

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

# Add group data of the validation dataset
data.lst$Validation <- adni_groups.dt[valid.lst[[1]], on = "ID"] |>
  melt(id = 1:4, variable = "MSR", value = "VAL")
rm(adni_groups.dt, segm.lst, valid.lst)


### ANALYSIS
comparisons.lst                   <- list()
## Cross-validation (Training)
# 2-way ANOVA w/trimmed means (rmanova; within-sample)
comparisons.lst$Training          <- list()
comparisons.lst$Training$ANOVA    <- list()
comparisons.lst$Training$PostHoc  <- list()

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

for (roi in levels(data.lst$Training$ROI)) {
  for (msr in levels(data.lst$Training$MSR)) {
    pb$tick(
      tokens = list(what = sprintf("Training :: %s: %s", roi, msr))
    )
    aov_res <- tryCatch({
      data.lst$Training[
        roi, on = "ROI"
      ][
        msr, on = "MSR"
      ][
        order(SEGM, SIDE),
        rmanova(VAL, interaction(SEGM, SIDE), ID)
      ] |> suppressWarnings()
    },
      error = \(e) {
        cat(sprintf(
          "Error in rmanova for %s - %s: %s\n", roi, msr, e$message))
        NULL
    })

    if (!is.null(aov_res)) {
      identifier <- paste(roi, msr, sep = "_")
      comparisons.lst$Training$ANOVA[[identifier]] <- data.table(
        ROI = roi,
        MSR = msr,
        Fstat = aov_res$test,
        DF = sprintf("(%.2f, %.2f)", aov_res$df1, aov_res$df2),
        Pval = aov_res$p.value,
        Padj = p.adjust(aov_res$p.value, method = "fdr", n = 10)
      )

      # If ANOVA ran, perform the relevant Post-hoc analysis
      # No need to use if; all comparisons are significant.
      posthoc <- data.lst$Training[
        roi, on = "ROI"
      ][
        msr, on = "MSR"
      ][
        order(SEGM, SIDE),
        rmmcp(VAL, interaction(SEGM, SIDE), ID)
      ] |> suppressWarnings()
      comparisons.lst$Training$PostHoc[[identifier]] <- posthoc$comp |>
        data.table() |>
        setnames(1:2, paste0("G", 1:2)) |>
        (\(DT) DT[, let(
          G1 = factor(G1, labels = posthoc$fnames[-6]),
          G2 = factor(G2, labels = posthoc$fnames[-1]),
          ROI = roi,
          MSR = msr
        )])()
      rm(identifier, posthoc)
    }
    rm(aov_res)
  }
  rm(msr)
}
rm(roi)

## Out-of-sample validation
# Robust mixed ANOVA
comparisons.lst$Validation          <- list()
comparisons.lst$Validation$ANOVA    <- list()
comparisons.lst$Validation$PostHoc  <- list()
for (msr in levels(data.lst$Validation$MSR)) {
  pb$tick(tokens = list(what = sprintf("Validation :: %s", msr)))

  aov_res <- tryCatch({
    bwtrim(VAL ~ GROUP * SIDE, ID, data.lst$Validation[msr[1], on = "MSR"])
  },
    error = \(e) {
      cat(sprintf("Error in bwtrim for %s (validation): %s\n", msr, e$message))
      NULL
  })

  if (!is.null(aov_res)) {
    aov.dt <- data.table(
      ROI = "HC",
      MSR = msr,
      COMP = c("GROUP", "SIDE", "GROUP:SIDE"),
      Fstat = sprintf("%.2f", as.numeric(aov_res[c("Qa", "Qb", "Qab")])),
      DF = sprintf(
        "(%i, %.2f)",
        aov_res[c("A.df", "B.df", "AB.df")] |> lapply(\(x) x[[1]]) |> unlist(),
        aov_res[c("A.df", "B.df", "AB.df")] |> lapply(\(x) x[[2]]) |> unlist()
      ),
      Pval = as.numeric(aov_res[c("A.p.value", "B.p.value", "AB.p.value")])
    )

    aov.dt[, Padj := p.adjust(Pval, method = "fdr", n = 5)]

    comparisons.lst$Validation$ANOVA[[msr]] <- aov.dt

    # If ANOVA ran, perform the relevant Post-hoc analysis
    if (aov.dt["GROUP", on = "COMP", Pval <= 0.05]) {
      posthoc <- lincon(VAL ~ GROUP, data.lst$Validation[msr, on = "MSR"])
      comparisons.lst$Validation$PostHoc[[msr]] <- posthoc$comp |>
        data.table() |>
        setnames(1:2, paste0("G", 1:2)) |>
        (\(DT) DT[, let(
          G1 = factor(G1, labels = posthoc$fnames[-3]),
          G2 = factor(G2, labels = posthoc$fnames[-1]),
          MSR = msr
        )])()
      rm(posthoc)
    }
  }
  rm(aov_res)
}
rm(pb, msr)

## Unlist
comparisons.lst <- lapply(
  comparisons.lst,
  \(level1) lapply(level1, rbindlist)
)

## Annotation data.table for plotting below
segms <- levels(data.lst$Training$SEGM)
annotation.dt <- rbind(
  comparisons.lst$Training$PostHoc[
    G1 %like% segms[1] & G2 %like% segms[2],
    .(sum(p.value < p.crit), x1 = segms[1], x2 = segms[2]),
    .(ROI, MSR)
  ][
    V1 == 4, -"V1"
  ],
  comparisons.lst$Training$PostHoc[
    G1 %like% segms[1] & G2 %like% segms[3],
    .(sum(p.value < p.crit), x1 = segms[1], x2 = segms[3]),
    .(ROI, MSR)
  ][
    V1 == 4, -"V1"
  ],
  comparisons.lst$Training$PostHoc[
    G1 %like% segms[2] & G2 %like% segms[3],
    .(sum(p.value < p.crit), x1 = segms[2], x2 = segms[3]),
    .(ROI, MSR)
  ][
    V1 == 4, -"V1"
  ]
) |>
 (\(DT) DT[
    !MSR %like% "ACC|SPEC",
    .(
      ROI = ROI |> factor(labels = c("Hippocampus", "Lateral Ventricles")),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      x1,
      x2,
      # SEGM1 : NLPBvMALF; SEGM2 : CNNvMALF; SEGM3 : CNNvNLPB
      y = fcase(
        x1 == "NLPB", .99,
        ROI %like% "H" & x2 == "NLPB", 1,
        ROI %like% "H" & x2 == "CNN", 1.015,
        #ROI %like% "V" & CONTRAST %like% 1, 1, ## This is not present
        ROI %like% "V" & x2 == "CNN", 1.03
      ),
      lab = "***"
    )
  ])()


### TABLES
fpath <- here("tables")
if (!file.exists(fpath)) dir.create(fpath)
## Cross-validation
# Median & SD - Hippocampus & Ventricles
fname <- "table-s1_man-seg_hcv-hvr_cv.tex"
data.lst$Training[
  ,
  .(VAL = sprintf("%.2f (%.2f)", median(VAL), sd(VAL))),
  .(
    SEGM = factor(SEGM, labels = sprintf("**%s**", c("MALF", "NLPB", "CNN"))),
    MSR = stringr::str_to_title(MSR),
    ROI_SIDE = paste(ROI, SIDE, sep = "_")
  )
] |>
  dcast(... ~ ROI_SIDE, value.var = "VAL") |>
  gt(rowname_col = "MSR", groupname_col = "SEGM", process_md = TRUE) |>
  tab_options(latex.tbl.pos = "h") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = c("HC_L", "HC_R", "VC_L", "VC_R")) |>
  tab_spanner(label = "Hippocampus", columns = starts_with("H"), level = 2) |>
  tab_spanner(label = "Ventricles", columns = starts_with("V"), level = 2) |>
  cols_label(ends_with("L") ~ "Left", ends_with("R") ~ "Right") |>
  tab_footnote(footnote = "Median (SD).") |>
  gtsave(filename = fname, path = fpath)

# Repeated measures ANOVA
fname <- "table-s2_man-seg_anova_cv.tex"
comparisons.lst$Training$ANOVA[
  order(MSR),
  .(
    Fstat,
    DF,
    Pval = sprintf("%.3f", Padj)
  ),
  .(
    ROI = fifelse(ROI == "HC", "**Hippocampus**", "**Ventricles**"),
    MSR = stringr::str_to_title(MSR)
  )
] |>
  gt(rowname_col = "MSR", groupname_col = "ROI", process_md = TRUE) |>
  tab_options(latex.tbl.pos = "h") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = c("Fstat", "DF", "Pval")) |>
  cols_label(
    "Fstat" ~ md("*F*"), "DF" ~ "df", "Pval" ~ md("*p*")
  ) |>
  tab_footnote(footnote = paste(
    "Repeated-measures robust two-way ANOVA with trimmed means.",
    "*p* values controlled for false discovery rate based on",
    "Benjamini & Hochberg."
  ) |> md()) |>
  gtsave(filename = fname, path = fpath)

# Posthoc (Sides)
fname <- "table-s3_man-seg_posthoc_cv_side.tex"
comparisons.lst$Training$PostHoc[
  substr(G1, 1, 1) == substr(G2, 1, 1),
  .(
    PSI = fcase(
      psihat == 0, "0",
      abs(psihat) < 0.01, sprintf("%.1e", psihat),
      abs(psihat) >= 0.01, sprintf("%.2f", psihat)
    ),
    CI = sprintf(
      "(%s, %s)",
      fcase(
        ci.lower == 0, "0",
        abs(ci.lower) < 0.01, sprintf("%.1e", ci.lower),
        abs(ci.lower) >= 0.01, sprintf("%.2f", ci.lower)
      ),
      fcase(
        ci.upper == 0, "0",
        abs(ci.upper) < 0.01, sprintf("%.1e", ci.upper),
        abs(ci.upper) >= 0.01, sprintf("%.2f", ci.upper)
      )
    ),
    Pval = fifelse(is.na(p.value), NA, sprintf("%.3f", p.value)),
    SIGN = p.value < p.crit
  ),
  .(
    ROI,
    SEGM = G1 |>
      gsub(pattern = "\\.(L|R)", replacement = "") |>
      factor(
        levels = c("MALF", "NLPB", "CNN"),
        labels = sprintf("**%s**", c("MALF", "NLPB", "CNN"))
      ),
    MSR = stringr::str_to_title(MSR)
  )
] |>
  dcast(... ~ ROI, value.var = c("PSI", "CI", "Pval", "SIGN")) |>
  setorder(MSR) |>
  setcolorder(c(1:3, 5, 7, 9)) |>
  gt(
    rowname_col = "MSR",
    groupname_col = "SEGM",
    #row_group_as_column = T,
    process_md = TRUE
  ) |>
  tab_options(latex.tbl.pos = "h") |>
  #tab_stubhead(label = "Contrasts") |>
  #tab_style(style = cell_text(size = "small"), locations = cells_body()) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_HC", "CI_HC", "Pval_HC"),
      rows = SIGN_HC == TRUE
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_VC", "CI_VC", "Pval_VC"),
      rows = SIGN_VC == TRUE
    )
  ) |>
  cols_hide(columns = contains("SIGN")) |>
  sub_missing(missing_text = "-") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = everything()) |>
  opt_horizontal_padding(scale = .5) |>
  tab_spanner(label = "Hippocampus", columns = ends_with("HC")) |>
  tab_spanner(label = "Ventricles", columns = ends_with("VC")) |>
  cols_label(
    starts_with("PSI") ~ md("$\\psi$"),
    starts_with("CI") ~ md("$95\\%~\\text{CI}$" ),
    starts_with("Pval") ~ md("$p$")
  ) |>
  tab_footnote(
    footnote = paste(
      "Robust pair-wise posthoc comparisons using trimmed means.",
      "Bold cells show significance after multiple comparison correction."
    )
  ) |>
  gtsave(filename = fname, path = fpath)

# Posthoc (Segmentations)
subDT <- comparisons.lst$Training$PostHoc[substr(G1,1,1) != substr(G2,1,1)] |>
  (\(DT){
    DT[, c("G1", "G1s") := DT[, tstrsplit(G1, "\\.")]]
    DT[, c("G2", "G2s") := DT[, tstrsplit(G2, "\\.")]]
    DT[, COMP := sprintf("%s-%s — %s-%s", G1, G1s, G2, G2s)]
    DT[MSR %like% "ACC" & psihat == 0, COMP := "MALF — NLPB"]
    DT[, COMP := factor(COMP, levels = DT[order(G2), unique(COMP)])]
    DT[
      ,
      .(
        CONTRAST = COMP,
        PSI = fcase(
          psihat == 0, "0",
          abs(psihat) < 0.01, sprintf("%.1e", psihat),
          abs(psihat) >= 0.01, sprintf("%.2f", psihat)
        ),
        CI = sprintf(
          "(%s, %s)",
          fcase(
            ci.lower == 0, "0",
            abs(ci.lower) < 0.01, sprintf("%.1e", ci.lower),
            abs(ci.lower) >= 0.01, sprintf("%.2f", ci.lower)
          ),
          fcase(
            ci.upper == 0, "0",
            abs(ci.upper) < 0.01, sprintf("%.1e", ci.upper),
            abs(ci.upper) >= 0.01, sprintf("%.2f", ci.upper)
          )
        ),
        Pval = fifelse(is.na(p.value), NA, sprintf("%.3f", p.value)),
        SIGN = p.value < p.crit
      ),
      .(
        ROI,
        MSR = MSR |> stringr::str_to_title() |> sprintf(fmt = "**%s**")
      )
    ]
  })() |>
  unique() |>
  dcast(... ~ ROI, value.var = c("PSI", "CI", "Pval", "SIGN")) |>
  setcolorder(c(1:3, 5, 7, 9))

# Overall agreement (Accuracy, Dice & Kappa)
fname <- "table-s4_man-seg_posthoc_cv_segm1.tex"
subDT[!MSR %like% "S"] |>
  gt(
    rowname_col = "CONTRAST",
    groupname_col = "MSR",
    #row_group_as_column = T,
    process_md = TRUE
  ) |>
  tab_options(
    latex.tbl.pos = "h"
  ) |>
  #tab_stubhead(label = "Contrasts") |>
  tab_style(style = cell_text(size = "small"), locations = cells_body()) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_HC", "CI_HC", "Pval_HC"),
      rows = SIGN_HC == TRUE
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_VC", "CI_VC", "Pval_VC"),
      rows = SIGN_VC == TRUE
    )
  ) |>
  cols_hide(columns = contains("SIGN")) |>
  sub_missing(missing_text = "-") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = everything()) |>
  opt_horizontal_padding(scale = .5) |>
  tab_spanner(label = "Hippocampus", columns = ends_with("HC")) |>
  tab_spanner(label = "Ventricles", columns = ends_with("VC")) |>
  cols_label(
    starts_with("PSI") ~ md("$\\psi$"),
    starts_with("CI") ~ md("$95\\%~\\text{CI}$" ),
    starts_with("Pval") ~ md("$p$")
  ) |>
  tab_footnote(
    footnote = paste(
      "Robust pair-wise posthoc comparisons using trimmed means.",
      "Bold cells show significance after multiple comparison correction."
    )
  ) |>
  gtsave(filename = fname, path = fpath)

# Class-wise performance (Sensitivity & Specificity)
fname <- "table-s5_man-seg_posthoc_cv_segm2.tex"
subDT[MSR %like% "S"] |>
  gt(
    rowname_col = "CONTRAST",
    groupname_col = "MSR",
    #row_group_as_column = T,
    process_md = TRUE
  ) |>
  tab_options(latex.tbl.pos = "h") |>
  #tab_stubhead(label = "Contrasts") |>
  tab_style(style = cell_text(size = "small"), locations = cells_body()) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_HC", "CI_HC", "Pval_HC"),
      rows = SIGN_HC == TRUE
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c("PSI_VC", "CI_VC", "Pval_VC"),
      rows = SIGN_VC == TRUE
    )
  ) |>
  cols_hide(columns = contains("SIGN")) |>
  sub_missing(missing_text = "-") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = everything()) |>
  opt_horizontal_padding(scale = .5) |>
  tab_spanner(label = "Hippocampus", columns = ends_with("HC")) |>
  tab_spanner(label = "Ventricles", columns = ends_with("VC")) |>
  cols_label(
    starts_with("PSI") ~ md("$\\psi$"),
    starts_with("CI") ~ md("$95\\%~\\text{CI}$" ),
    starts_with("Pval") ~ md("$p$")
  ) |>
  tab_footnote(
    footnote = paste(
      "Robust pair-wise posthoc comparisons using trimmed means.",
      "Bold cells show significance after multiple comparison correction."
    )
  ) |>
  gtsave(filename = fname, path = fpath)


## Out-of-sample Validation
# Mean & SD of Hippocampus by Group
fname <- "table-s6_man-seg_hcv-hvr_val.tex"
data.lst$Validation[
  order(MSR),
  .(VAL = sprintf("%.2f (%.2f)", median(VAL), sd(VAL))),
  .(
    GROUP = factor(
      GROUP,
      labels = sprintf("**%s**, N: 20", c("CH", "MCI", "AD"))
    ),
    MSR = stringr::str_to_title(MSR), SIDE
  )
] |>
  dcast(... ~ SIDE, value.var = "VAL") |>
  gt(rowname_col = "MSR", groupname_col = "GROUP", process_md = TRUE) |>
  tab_options(latex.tbl.pos = "h") |>
  cols_align("center", columns = c("L", "R")) |>
  tab_spanner(label = "Hippocampus", columns = c("L", "R")) |>
  cols_label(ends_with("L") ~ "Left", ends_with("R") ~ "Right") |>
  tab_footnote(footnote = "Median (SD).") |>
  gtsave(filename = fname, path = fpath)

# Mixed ANOVA of Group & Side
fname <- "table-s7_man-seg_anova_val.tex"
comparisons.lst$Validation$ANOVA[
  order(MSR),
  .(
    Fstat,
    DF,
    Pval = sprintf("%.3f", Padj)
  ),
  .(
    MSR = stringr::str_to_title(MSR),
    COMP = factor(
      COMP,
      levels = c("GROUP", "SIDE", "GROUP:SIDE"),
      labels = c("**Group**", "**Side**", "**Group:Side**")
    )
  )
] |>
  gt(rowname_col = "MSR", groupname_col = "COMP", process_md = TRUE) |>
  tab_options(latex.tbl.pos = "h") |>
  tab_stub_indent(rows = everything(), indent = 1) |>
  cols_align("center", columns = c("Fstat", "DF", "Pval")) |>
  cols_label(
    "Fstat" ~ md("*F*"), "DF" ~ "df", "Pval" ~ md("*p*")
  ) |>
  tab_footnote(footnote = paste(
    "Robust two-way mixed ANOVA with trimmed means.",
    "*p* values controlled for false discovery rate based on",
    "Benjamini & Hochberg."
  ) |> md()) |>
  gtsave(filename = fname, path = fpath)
rm(fname)

### PLOTS
## Boxplot Kappas HC(VC)
outdir <- here("plots")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
fplots <- Map(
  \(dset, fig) {
    dset <- dset |> tolower() |> substr(1, 5)
    sprintf(
      "%s/%s_man-seg_comp_%s.%s",
      outdir,
      fig,
      dset,
      c("png", "tiff")
    )
  },
  c("Training", "Validation"),
  c("fig-p1", "fig-s1")
)

# Cross-validation
if (any(REDOPLOTS, !file.exists(fplots[[1]]))) {
  ## Accuracy & Specificity ~ 1.0
  p <- data.lst$Training[
    !MSR %like% "ACC|SPEC",
    .(
      SEGM = SEGM |> factor(levels = c("MALF", "NLPB", "CNN")),
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
      data = annotation.dt,
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

# Out-of-sample validation
if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  ## Accuracy & Specificity > .97
  p <- data.lst$Validation[
    !MSR %like% "ACC|SPEC",
    .(
      SIDE = SIDE |> factor(labels = c("Left", "Right")),
      ROI = factor("Hippocampus"),
      MSR = MSR |> stringr::str_to_title() |> factor(),
      GROUP,
      VAL
    )
  ] |>
    ggplot(mapping = aes(x = GROUP, y = VAL)) +
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
    scale_colour_manual(values = c("darkred", "midnightblue")) +
    scale_fill_manual(values = c("darkred", "midnightblue")) +
    facet_grid(rows = vars(ROI), cols = vars(MSR), scales = "free_y") +
    labs(
      y = "Overlap/Accuracy",
      x = paste(
        "Clinical group (<span style='color:darkred;'>Left</span> &amp;",
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

## Volumes Correlation
# Read Volume CSVs
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
fpath     <- here("data/derivatives/man-seg_volumes_hc_adni.csv")

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
        SIDE = factor(SIDE, labels = c("Left", "Right")),
        ROI = "Hippocampus",
        value = value / 1000
      )
    ]
    DT[, ID := sub("hc_adni_(man|cnn)_", "", ID)]
    data.lst$Validation[!duplicated(ID), GROUP, keyby = ID][DT]
  })() |>
  dcast(... ~ SEGM, value = "value") |>
  setcolorder(c("ID"))
vols.lst$Validation[["ADNI"]] <- copy(DT)
rm(fpath, segm, DT)

vols.lst <- lapply(vols.lst, rbindlist)

dt1 <- vols.lst[[1]][SEGM != "MANUAL", .(SEGM, CC), keyby = .(ID, SIDE, ROI)]
dt2 <- vols.lst[[1]][SEGM == "MANUAL", .(MAN = CC), keyby = .(ID, SIDE, ROI)]
vols.lst[[1]] <- dt1[dt2]
rm(dt1, dt2)

## Correlation plots
fplots <- Map(
  \(dset, fig) {
    dset <- dset |> tolower() |> substr(1, 5)
    sprintf(
      "%s/%s_man-seg_corr_%s.%s",
      outdir,
      fig,
      dset,
      c("png", "tiff")
    )
  },
  c("Training", "Validation"),
  c("fig-2", "fig-s2")
)

cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# Cross-validation
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
      cor.coef.name = "rho",
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

# Out-of-sample validation
if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  p <- vols.lst$Validation |>
    ggplot(aes(x = CNN, y = MAN, colour = GROUP)) +
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
      cor.coef.name = "rho",
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
        "Volumes segmented by CNN for: ",
        "<span style='color:%s;'>CH</span>,",
        "<span style='color:%s;'>MCI</span>, and",
        "<span style='color:%s;'>AD</span>"
      ) |> sprintf(cbPalette[2], cbPalette[3], cbPalette[8])
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
fplots <- Map(
  \(dset, fig) {
    dset <- dset |> tolower() |> substr(1, 5)
    sprintf(
      "%s/%s_man-seg_blandaltman_%s.%s",
      outdir,
      fig,
      dset,
      c("png", "tiff")
    )
  },
  c("Training", "Validation"),
  c("fig-3", "fig-s3")
)

# Cross-validation
if (any(REDOPLOTS, !file.exists(fplots[[1]]))) {
  DT <- vols.lst[[1]][
    , SEGM := factor(SEGM, levels = c("MALF", "NLPB", "CNN"))
  ][
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
      width = 7,
      height = 4,
      units = "in",
      dpi = 600
    )
  }

  if (any(REDOPLOTS, !file.exists(fplots[[1]][2]))) {
    ggsave(
      filename = fplots[[1]][2],
      plot = p,
      width = 7,
      height = 4,
      units = "in",
      device = "tiff",
      dpi = 600
    )
  }
}

# Out-of-sample validation
if (any(REDOPLOTS, !file.exists(fplots[[2]]))) {
  DT <- vols.lst[[2]][
    , .(
      AVG = (CNN + MAN) / 2,
      DIFF = CNN - MAN
    ),
    .(SIDE, ROI, GROUP)
  ][
    , .(
      AVG, DIFF,
      MEAN_DIFF = mean(DIFF),
      CI_low = mean(DIFF) - 1.96 * sd(DIFF),
      CI_high = mean(DIFF) + 1.96 * sd(DIFF)
    ),
    .(ROI, GROUP)
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
      data = DT[, .SD[1], .(ROI, GROUP)],
      aes(yintercept = MEAN_DIFF),
      colour = cbPalette[2],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], GROUP],
      aes(x = 1.2, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, GROUP)],
      aes(yintercept = CI_low),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], GROUP],
      aes(x = 1.2, y = CI_low, label = round(CI_low, 2)),
      nudge_y = -.100,
      size = 2.5
    ) +
    geom_hline(
      data = DT[, .SD[1], .(ROI, GROUP)],
      aes(yintercept = CI_high),
      colour = cbPalette[3],
      linetype = "dashed",
      alpha = .7
    ) +
    geom_text(
      data = DT[, .SD[1], GROUP],
      aes(x = 1.2, y = CI_high, label = round(CI_high, 2)),
      nudge_y = .1,
      size = 2.5
    ) +
    xlim(.95, 6.2) +
    facet_grid(rows = vars(ROI), cols = vars(GROUP), scales = "free") +
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
