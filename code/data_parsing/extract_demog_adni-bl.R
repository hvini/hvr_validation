#!/usr/bin/env Rscript

library(here)
library(data.table)
library(gt)
library(dunn.test)

### INPUT
fpath         <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge   <- readRDS(fpath)
} else {
  here("code/data_parsing/parse_adnimerge-bl.R") |> source()
}

fpath         <- here("data/rds/adni-bl_volumes_hcvc.rds")
if (file.exists(fpath)) {
  volumes     <- readRDS(fpath)
} else {
  here("code/data_parsing/qc_segmentations_adni-bl.R") |> source()
}


### Data CLEANING
demog.dt      <- adnimerge[
  volumes[METHOD == "cnn"],
  on = "PTID",
  .(
    PTID, DX, PTGENDER, AGE, PTEDUCAT, ADAS13,
    RAVLT_learning = as.numeric(RAVLT_learning)
  )
]

### Data EXPLORATION
## N
demog.dt[, .N, DX]
#        DX     N
#    <fctr> <int>
# 1:     CH   501
# 2:     AD   321
# 3:    MCI   819

## Females (N; percentage)
demog.dt[
  ,
  .(
    N = sum(PTGENDER %like% "F"),
    Perc = sum(PTGENDER %like% "F") / .N * 100
  ),
  DX
]
#        DX     N     Perc
#    <fctr> <int>    <num>
# 1:     CH   260 51.89621
# 2:     AD   146 45.48287
# 3:    MCI   336 41.02564

## Education
demog.dt[, .(M = mean(PTEDUCAT), SD = sd(PTEDUCAT)), DX]
#        DX        M       SD
#    <fctr>    <num>    <num>
# 1:     CH 16.34930 2.704024
# 2:     AD 15.18692 2.980554
# 3:    MCI 15.91819 2.827891

## Age
demog.dt[, .(M = mean(AGE), SD = sd(AGE)), DX]
#        DX        M       SD
#    <fctr>    <num>    <num>
# 1:     CH 74.25349 5.749318
# 2:     AD 74.84953 7.810430
# 3:    MCI 72.89426 7.620770

## Cognition
demog.dt[!is.na(ADAS13), .(M = mean(ADAS13), SD = sd(ADAS13)), DX]
#        DX         M       SD
#    <fctr>     <num>    <num>
# 1:     CH  9.255772 4.405977
# 2:     AD 29.933419 7.977899
# 3:    MCI 16.633741 6.750174

# Memory
demog.dt[
  !is.na(RAVLT_learning),
  .(M = mean(RAVLT_learning), SD = sd(RAVLT_learning)),
  DX
]
#        DX        M       SD
#    <fctr>    <num>    <num>
# 1:     CH 5.845070 2.302291
# 2:     AD 1.810726 1.790131
# 3:    MCI 4.067485 2.593530

## Demographics table
fname <- "adni-bl_table-1.tex"
fpath <- here("tables")
demog.dt[, .SD, DX, .SDcols = AGE:RAVLT_learning] |>
  melt(id = "DX", variable = "VAR") |>
  suppressWarnings() |>
  na.omit() |>
  (
    \(DT) DT[
      ,
      sprintf("%.1f (%.1f)", mean(value), sd(value)),
      .(DX, VAR)
    ]
  )() |>
  rbind(
    demog.dt[
      , lapply(.SD, \(col) sum(is.na(col))), DX, .SDcols = AGE:RAVLT_learning
    ] |>
      melt(id = "DX", variable = "VAR", value = "V1") |>
      (
        \(DT)
        DT[V1 > 0, .(VAR = sprintf("%s_NA", substr(VAR, 1, 4)), V1), DX]
      )(),
    demog.dt[
      , .N, keyby = DX
    ][
      demog.dt[, .N, keyby = .(DX, PTGENDER)]
    ][
      "Female", on = "PTGENDER",
      .(
        sprintf("%i (%.0f%%)", i.N, 100 * i.N /N),
        VAR = "SEXF"
      ),
      DX
    ],
    use.names = TRUE
  ) |>
  dcast(VAR ~ DX, value.var = "V1") |>
  (\(DT) DT[
    ,
    VAR := factor(
      VAR,
      levels = c(
        "SEXF",
        "AGE",
        "PTEDUCAT",
        "ADAS13",
        "ADAS_NA",
        "RAVLT_learning",
        "RAVL_NA"
      ),
      labels = c(
        "Sex (F)",
        "Age (years)",
        "Education (years)",
        "ADAS-13",
        "NA_ADAS13",
        "RAVLT-learning",
        "NA_RAVLT"
      )
    )][order(VAR)]
  )() |>
  gt(rowname_col = "VAR", process_md = TRUE) |>
  tab_spanner(label = "Clinical Label", columns = c("CH", "MCI", "AD")) |>
  tab_options(
    latex.tbl.pos = "h",
    footnotes.multiline = FALSE
  ) |>
  cols_align("center", columns = c("CH", "MCI", "AD")) |>
  cols_label(
    CH = sprintf("**CH**, N: %i", demog.dt["CH", on = "DX", .N]) |> md(),
    MCI = sprintf("**MCI**, N: %i", demog.dt["MCI", on = "DX", .N]) |> md(),
    AD = sprintf("**AD**, N: %i", demog.dt["AD", on = "DX", .N]) |> md()
  ) |>
  tab_stub_indent(starts_with("NA"), indent = 3) |>
  sub_values(values = c("NA_ADAS13", "NA_RAVLT"), replacement = "Missing") |>
  tab_footnote(
    footnote = "N (%).",
    locations = cells_stub(rows = contains("Sex"))
  ) |>
  tab_footnote(
    footnote = "N of subjects without cognitive data.",
    locations = cells_stub(rows = contains("NA"))
  ) |>
  tab_footnote(
    footnote = "Mean (SD).",
    locations = cells_stub(rows = c(2:4, 6))
  ) |>
  gtsave(filename = fname, path = fpath)


### Post-hoc analyses
## Chi2: DX v Sex
demog.dt[, chisq.posthoc.test::chisq.posthoc.test(table(DX, PTGENDER))]
#   Dimension     Value     Female       Male
# 1        CH Residuals  3.6042584 -3.6042584
# 2        CH  p values  0.0018780  0.0018780
# 3       MCI Residuals -3.4046518  3.4046518
# 4       MCI  p values  0.0039750  0.0039750
# 5        AD Residuals  0.1069794 -0.1069794
# 6        AD  p values  1.0000000  1.0000000

## Dunn.tests
demog.dt[, dunn.test(AGE, DX, method = "bonferroni")]
#   Kruskal-Wallis rank sum test
#
# data: AGE and DX
# Kruskal-Wallis chi-squared = 18.0925, df = 2, p-value = 0
#
#
#                             Comparison of AGE by DX
#                                  (Bonferroni)
# Col Mean-|
# Row Mean |         AD         CH
# ---------+----------------------
#       CH |   1.645034
#          |     0.1499
#          |
#      MCI |   4.029487   2.604685
#          |    0.0001*    0.0138*
#
# alpha = 0.05
# Reject Ho if p <= alpha/2
#        chi2        Z            P   P.adjusted comparisons
#       <num>    <num>        <num>        <num>      <char>
# 1: 18.09248 1.645035 4.998130e-02 1.499439e-01     AD - CH
# 2: 18.09248 4.029488 2.794928e-05 8.384784e-05    AD - MCI
# 3: 18.09248 2.604686 4.597931e-03 1.379379e-02    CH - MCI

demog.dt[, dunn.test(PTEDUCAT, DX, method = "bonferroni")]
#   Kruskal-Wallis rank sum test
#
# data: PTEDUCAT and DX
# Kruskal-Wallis chi-squared = 31.3186, df = 2, p-value = 0
#
#
#                          Comparison of PTEDUCAT by DX
#                                  (Bonferroni)
# Col Mean-|
# Row Mean |         AD         CH
# ---------+----------------------
#       CH |  -5.594794
#          |    0.0000*
#          |
#      MCI |  -3.799494   2.640956
#          |    0.0002*    0.0124*
#
# alpha = 0.05
# Reject Ho if p <= alpha/2
#        chi2         Z            P   P.adjusted comparisons
#       <num>     <num>        <num>        <num>      <char>
# 1: 31.31857 -5.594794 1.104418e-08 3.313254e-08     AD - CH
# 2: 31.31857 -3.799494 7.249584e-05 2.174875e-04    AD - MCI
# 3: 31.31857  2.640957 4.133612e-03 1.240084e-02    CH - MCI
demog.dt[, dunn.test(ADAS13, DX, method = "bonferroni")]
#   Kruskal-Wallis rank sum test
#
# data: ADAS13 and DX
# Kruskal-Wallis chi-squared = 850.4675, df = 2, p-value = 0
#
#
#                           Comparison of ADAS13 by DX
#                                  (Bonferroni)
# Col Mean-|
# Row Mean |         AD         CH
# ---------+----------------------
#       CH |   29.02225
#          |    0.0000*
#          |
#      MCI |   17.25578  -16.62969
#          |    0.0000*    0.0000*
#
# alpha = 0.05
# Reject Ho if p <= alpha/2
#        chi2         Z             P    P.adjusted comparisons
#       <num>     <num>         <num>         <num>      <char>
# 1: 850.4675  29.02226 1.723423e-185 5.170269e-185     AD - CH
# 2: 850.4675  17.25578  5.061886e-67  1.518566e-66    AD - MCI
# 3: 850.4675 -16.62969  2.123926e-62  6.371777e-62    CH - MCI

demog.dt[, dunn.test(RAVLT_learning, DX, method = "bonferroni")]
#   Kruskal-Wallis rank sum test
#
# data: RAVLT_learning and DX
# Kruskal-Wallis chi-squared = 446.2512, df = 2, p-value = 0
#
#
#                       Comparison of RAVLT_learning by DX
#                                  (Bonferroni)
# Col Mean-|
# Row Mean |         AD         CH
# ---------+----------------------
#       CH |  -21.07070
#          |    0.0000*
#          |
#      MCI |  -12.84059   11.67727
#          |    0.0000*    0.0000*
#
# alpha = 0.05
# Reject Ho if p <= alpha/2
#        chi2         Z            P   P.adjusted comparisons
#       <num>     <num>        <num>        <num>      <char>
# 1: 446.2512 -21.07070 7.385824e-99 2.215747e-98     AD - CH
# 2: 446.2512 -12.84059 4.856327e-38 1.456898e-37    AD - MCI
# 3: 446.2512  11.67728 8.326994e-32 2.498098e-31    CH - MCI
