#!/usr/bin/env Rscript

library(here)
library(data.table)
library(progress)
library(DescTools)
library(ggplot2)
library(gridExtra)
library(ggtext)
#library(ggridges)
#library(ggsignif)
#library(ggpubr)

## Calculate and compare correlations of HC & Age | Memory | Cognition
## ADNI data CN|MCI|AD

RERUNPERMS <- FALSE

# Read RDS objects
fpaths <- list(
  rds = c(
    "adnimerge_baseline.rds",
    #"adni-bl_volumes_freesurfer.rds",
    "adni-bl_volumes_icv-adjusted.rds"
  ),
  scripts = c(
    "data_parsing/parse_adnimerge-bl.R",
    #"data_parsing/parse_freesurfer-vols.R",
    "analysis/adjust_hc-hvr_adni-bl.R"
  )
)

data.lst <- vector("list", 2)
for (i in seq_along(data.lst)){
  rds <- fpaths[["rds"]][i] |> sprintf(fmt = "data/rds/%s") |> here()
  if (file.exists(rds)) {
    data.lst[[i]] <- readRDS(rds)
  } else {
    fpaths[["scripts"]][i] |> sprintf(fmt = "code/%s") |> here() |> source()
  }
}

rm(i, rds, fpaths)

# Merge
data.dt <- merge(
  data.lst[[1]][
    !is.na(ADAS13)
  ][
    !is.na(RAVLT_learning),
    .(
      DX,
      AGE,
      ADAS13,
      RAVLT_immediate = as.numeric(RAVLT_immediate),
      RAVLT_perc_forgetting = as.numeric(RAVLT_perc_forgetting),
      RAVLT_learning = as.numeric(RAVLT_learning)
    ),
    keyby = PTID
  ],
  data.lst[[2]][
    !is.na(HVR_mean),
    .(
      METHOD,
      HCv_l = HC_stx_l,
      HCv_r = HC_stx_r,
      HCv = HC_stx_mean, # Head-size normalized
      HVR_l,
      HVR_r,
      HVR = HVR_mean
    ),
    keyby = "PTID"
  ]
)
rm(data.lst)


## CORRELATION | Permutation tests
fpath <- here("data/rds/adni-bl_hcv-hvr_corrs_non-parametric.rds")
if (all(file.exists(fpath), !RERUNPERMS)) {
  corr.lst <- readRDS(fpath)
} else {
  n_perms <- 10000

  params.lst <- list(
    DX   = list(),
    COV  = list(),
    MTD  = list(),
    HC   = list(),
    Rval = list(),
    Tval = list(),
    Pval = list(),
    CIl  = list(),
    CIh  = list()
  )

  # HCv vs HVR
  perms.lst <- list(
    HCvHVR = list(
      DX      = list(),
      COVAR   = list(),
      METHOD  = list(),
      P_diff  = list()
    ),
    CNNvFS = list(
      METHOD  = "CNN-FS6",
      DX      = list(),
      COVAR   = list(),
      P_diff  = list()
    )
  )

  set.seed(1618)
  # 3 DXs ; 3 Covars ; 4 segmentation methods
  ticks1 <- 3 * 3 * 4 * n_perms
  ticks2 <- 3 * 3 * n_perms
  pb <- progress_bar$new(
    format = "Permutations | :what [:bar] :current/:total",
    total = ticks1 + ticks2,
    clear = FALSE,
    width = 75
  )
  rm(ticks1, ticks2)

  for (dx in data.dt[, unique(DX)]) {
    for (covar in c("AGE", "RAVLT_learning", "ADAS13")) {
      ## CNN vs FS6 (HVR)
      subDT <- data.dt[dx, on = "DX"]
      for (p in 1:n_perms) {
        pb$tick(tokens = list(what = sprintf("%s: CNN vs FS6\n", dx)))
        perms.lst$CNNvFS$DX     <- c(perms.lst$CNNvFS$DX, dx)
        perms.lst$CNNvFS$COVAR  <- c(perms.lst$CNNvFS$COVAR, covar)
        perms.lst$CNNvFS$P_diff <- c(
          perms.lst$CNNvFS$P_diff,
          subDT[
            c("cnn", "fs6"),
            on = "METHOD",
            .(HVR, COVAR = get(covar), METHOD = sample(METHOD))
          ][
            ,
            cor(HVR, COVAR, method = "spearman", use = "complete.obs"),
            by = "METHOD"
          ][
            # CNN before FS, Positive values are CNN > FS6
            order(METHOD), diff(V1)
          ]
        )
      }
      for (mtd in c("cnn", "fs6", "malf", "nlpb")) {
        subDT <- data.dt[dx, on = "DX"][mtd, on = "METHOD"]
        ## HCv vs HVR
        for (p in 1:n_perms) {
          pb$tick(tokens = list(what = sprintf("%s: %s\n", dx, mtd)))
          perms.lst$HCvHVR$DX     <- c(perms.lst$HCvHVR$DX, dx)
          perms.lst$HCvHVR$COVAR  <- c(perms.lst$HCvHVR$COVAR, covar)
          perms.lst$HCvHVR$METHOD <- c(perms.lst$HCvHVR$METHOD, mtd)
          perms.lst$HCvHVR$P_diff <- c(
            perms.lst$HCvHVR$P_diff,
            subDT[, .(PTID, COVAR = get(covar), HCv, HVR)] |>
            melt(measure = c("HCv", "HVR"), value = "MSR") |>
            (
              \(DT) DT[
                , .(MSR, COVAR, HC = sample(variable))
              ][
                ,
                cor(MSR, COVAR, method = "spearman", use = "complete.obs"),
                by = "HC"
              ][
                # HCv first; Positive values are HCv > HVR
                order(HC), diff(V1)
              ]
            )()
          )
        }
        for (msr in c("HCv", "HVR")) {   # HCv,HVR
          corr <- list(
            RHO = subDT[
              ,
              cor.test(get(covar), get(msr), method = "spearman") |>
                suppressWarnings()
            ],
            CI = subDT[, DescTools::SpearmanRho(
              get(covar), get(msr), conf.level = .95
              )
            ]
          )
          params.lst$DX   <- c(params.lst$DX, dx)
          params.lst$COV  <- c(params.lst$COV, covar)
          params.lst$MTD  <- c(params.lst$MTD, mtd)
          params.lst$HC   <- c(params.lst$HC, msr)
          params.lst$Rval <- c(params.lst$Rval, corr$RHO$estimate)
          params.lst$Tval <- c(params.lst$Tval, corr$RHO$statistic)
          params.lst$Pval <- c(params.lst$Pval, corr$RHO$p.value)
          params.lst$CIl  <- c(params.lst$CIl,  corr$CI["lwr.ci"])
          params.lst$CIh  <- c(params.lst$CIh,  corr$CI["upr.ci"])
          rm(corr)
        }
      }
    }
  }
  rm(n_perms, dx, covar, p, mtd, msr, subDT)

  params.dt <- as.data.table(params.lst) |>
    (\(DT) DT[
      ,
      c("DX", "COV", "MTD", "HC") := lapply(.SD, as.character),
      .SDcols = DX:HC
    ][
      , DX := factor(DX, levels = c("CH", "MCI", "AD"))
    ][
      , COV := factor(
        COV,
        levels = c("AGE", "RAVLT_learning", "ADAS13"),
        labels = c("Age", "Memory", "Cognition")
      )
    ][
      ,
      c("Rval", "Tval", "Pval", "CIl", "CIh") := lapply(.SD, as.numeric),
      .SDcols = Rval:CIh
    ][
      , Pa := p.adjust(Pval, method = "bonferroni")
    ][
      , SIGN := fcase(
        Pa < 0.001, "***",
        Pa < 0.01, "**",
        Pa < 0.05, "*",
        default = ""
      )
    ])() |>
    setnames(c("COV", "MTD", "Pa"), c("COVAR", "METHOD", "Pval_adj"))

  perms.dt <- perms.lst |>
    lapply(as.data.table) |>
    rbindlist(use.names = TRUE) |>
    (\(DT) DT[
      ,
      c("DX", "COVAR", "METHOD") := lapply(.SD, as.character),
      .SDcols = DX:METHOD
    ][
      , COVAR := factor(
        COVAR,
        levels = c("AGE", "RAVLT_learning", "ADAS13"),
        labels = c("Age", "Memory", "Cognition")
      )
    ][
      , P_diff := as.numeric(P_diff)
    ][
      , COMP := fifelse(
        METHOD == "CNN-FS6",
        "CNN-FS6_HVR",
        sprintf("%s_HCv-HVR", METHOD)
      )
    ])()

  corr.lst <- list(
    COEFS = params.dt,
    PERMS = perms.dt,
    CONTR = rbind(
      params.dt[order(HC), .(DIFF = diff(Rval)), by = .(DX, COVAR, METHOD)],
      params.dt[
        "HVR", on = "HC"
      ][
        c("cnn", "fs6"), on = "METHOD"
      ][
        order(METHOD),
        .(METHOD = "CNN-FS6", DIFF = diff(Rval)),
        by = .(DX, COVAR)
      ],
      use.names = TRUE
    ) |>
    merge(perms.dt, by = c("DX", "COVAR", "METHOD")) |>
    (\(DT) DT[
      ,
      .(
        Pval = fifelse(
          COVAR == "Memory",
          sum(P_diff >= DIFF) / .N,
          sum(P_diff <= DIFF) / .N
        )
      ),
      .(DX, COVAR, METHOD)
      ][
        Pval < 0.05,
        LABEL := fcase(Pval < 0.001, "***", Pval < 0.01, "**", default = "*")
      ]
    )() |>
    merge(params.dt[, .(Y = max(CIh)), keyby = DX:METHOD], all.x = T)
  )
  rm(params.lst, perms.lst, params.dt, perms.dt)
  saveRDS(corr.lst, fpath)
}
rm(fpath)

### PLOTS
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

plots.lst <- plot_params.lst <- list()
plot_params.lst[["TITLE"]] <- c("CNN", "FS v6.0", "MALF", "NLPB")
plot_params.lst[["N"]] <- data.dt[order(METHOD), .N, METHOD]$N
plot_params.lst[["COLOUR"]] <- cbPalette[c(2:3, 8, 6)]
plot_params.lst |> lapply(
  setattr,
  "names",
  corr.lst$COEFS$METHOD |> unique() |> sort()
) |> invisible()
plot_params.lst[["X"]] <- c(
  "HCv" = "<span style='color: darkred;'>HCv</span>",
  "HVR" = "<span style='color: midnightblue;'>HVR</span>"
)

# Helper function to create layers for each data subset
plot_layers.fn <- function(dt, hc, side) {
  plot_params.lst <- list(
    L = list(
      colour = "darkred",
      pos_nudge = .2,
      hjust = "right"
    ),
    R = list(
      colour = "midnightblue",
      pos_nudge = -.2,
      hjust = "left"
    )
  )
  DT <- dt[hc, on = "HC"]
  list(
    geom_errorbar(
      data = DT,
      aes(ymin = CIl, ymax = CIh),
      colour = plot_params.lst[[side]]$colour,
      width = 0.2,
      position = position_nudge(x = plot_params.lst[[side]]$pos_nudge)
    ),
    geom_point(
      data = DT,
      shape = 21,
      colour = plot_params.lst[[side]]$colour,
      fill = "white",
      size = .9,
      stroke = 0.3,
      position = position_nudge(x = plot_params.lst[[side]]$pos_nudge)
    ),
    geom_text(
      data = DT,
      aes(label = SIGN, y = CIh),
      colour = plot_params.lst[[side]]$colour,
      size = 2.5,
      vjust = 0.1,
      position = position_nudge(x = plot_params.lst[[side]]$pos_nudge)
    ),
    geom_text(
      data = DT,
      aes(label = round(Rval, 2)),
      colour = plot_params.lst[[side]]$colour,
      size = 2.5,
      hjust = plot_params.lst[[side]]$hjust
    )
  )
}


for (mtd in unique(corr.lst$COEFS$METHOD)) {
  subDT <- corr.lst$COEFS[mtd, on = "METHOD"]
  plots.lst[[mtd]] <- ggplot(subDT, aes(HC, Rval)) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(size = 11),
      axis.text.x = element_markdown(),
      axis.title.x = element_blank(),
      plot.caption = element_markdown(size = 10),
      legend.position = "none"
    ) +
    facet_grid(rows = vars(DX), cols = vars(COVAR)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      alpha = .5,
      colour = cbPalette[1]
    ) +
    plot_layers.fn(subDT, "HCv", "L") +
    plot_layers.fn(subDT, "HVR", "R") +
    ylim(-.65, .5) +
    scale_x_discrete(labels = plot_params.lst$X) +
    labs(
      y = "Spearman's rho",
      caption = sprintf(
        "**%s:** N = %i;<br>* p < 0.05; ** p < 0.01; *** p < 0.001",
        plot_params.lst$TITLE[[mtd]],
        plot_params.lst$N[[mtd]]
      )
    )
  }


#GGSIGNIF
#CNN
geom_signif(
  data = corr.perm.sign.dt[METHOD == "CNN"],
  aes(xmin = HCv, xmax = HVR, annotations = LABEL, y_position = Y + .13),
  manual = TRUE,
  colour = cbPalette[1],
  textsize = 3,
  inherit.aes = FALSE
)

# MALF
geom_signif(data = corr.perm.sign.dt[METHOD == "MALF"],
              aes(xmin = HCv, xmax = HVR, annotations = LABEL,
                                    y <- position = Y + .1), manual = TRUE, colour = cbPalette[1],
              textsize = 3, inherit.aes = FALSE)

# NLPB
geom_signif(data = corr.perm.sign.dt[METHOD == "NLPB"],
              aes(xmin = HCv, xmax = HVR, annotations = LABEL,
                                    y <- position = Y + .13), manual = TRUE, colour = cbPalette[1],
              textsize = 3, inherit.aes = FALSE)

## Also, figure out correct order
p <- grid.arrange(
  plots.lst[["cnn"]],
  plots.lst[["malf"]],
  plots.lst[["nlpb"]],
  plots.lst[["fs6"]],
  nrow = 2
)

fpaths <- c("png", "tiff") |>
  sprintf(fmt = "plots/adni-bl_hcv-hvr_corrs.%s") |>
  here()

ggsave(fpaths[[1]], p, width = 7, height = 7, units = "in", dpi = 600)
ggsave(fpaths[[2]], p, width = 7, height = 7, units = "in", device = "tiff", dpi = 600)
