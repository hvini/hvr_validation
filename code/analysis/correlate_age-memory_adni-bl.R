#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(progress)
library(DescTools)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggtext)
library(gridExtra)
library(ggpubr)
library(glue)

## Calculate and compare correlations of HC & Age | Memory | Cognition
## ADNI data CN|MCI|AD

ReRunPerms      <- FALSE

# Read RDS objects
# ADNIMERGE
fpath    <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge     <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_adnimerge-bl.R') |> source()
}

fpath    <- here("data/rds/adni-bl_volumes_icv-adjusted.rds")
if (file.exists(fpath)) {
  volumes       <- read_rds(fpath)
} else {
  here('code/analysis/adjust_hc-hvr_adni.R') |> source()
}

# Merge
DT            <- adnimerge[volumes, on = "PTID",
                         .(PTID, METHOD, DX, AGE, ADAS13,
                           #RAVLT_immediate, RAVLT_learning, RAVLT_forgetting,
                           RAVLT_immediate = as.numeric(RAVLT_immediate),
                           RAVLT_perc_forgetting = as.numeric(RAVLT_perc_forgetting),
                           RAVLT_learning = as.numeric(RAVLT_learning),
                           #HC = HC_mean,     # Native space
                           HCv_l = HC_stx_l,
                           HCv_r = HC_stx_r,
                           HCv = HC_stx_mean, # Head-size normalized
                           HVR_l, HVR_r,
                           HVR = HVR_mean)
                         ][!is.na(HVR)    &
                           !is.na(ADAS13) &
                           !is.na(RAVLT_learning)]



## Correlation | Permutation tests
fpaths <- here("data/rds", paste0("adni-bl_hcv-hvr_corrs_",
                                  c("", "permutations_"),
                                  "non-parametric.rds"))
if (all(file.exists(fpaths), !ReRunPerms)) {
  corr.dt       <- read_rds(fpaths[1])
  perms.dif.dt  <- read_rds(fpaths[2])
} else {
  dxs           <- DT[, unique(DX)]
  mtds1         <- c("cnn", "fs6", "malf", "nlpb")
  mtds2         <- c("CNN", "FS_V6", "MALF", "NLPB")
  msrs          <- c("HCv", "HVR")
  covs          <- c("AGE", "RAVLT_learning", "ADAS13")

  n_perms       <- 10000

  r <- t <- dfs <- pval <- cil <- cih <-
    rep(NA, length(dxs) * length(covs) * length(mtds1) * length(msrs))

  # HCv vs HVR
  cor.difs1 <- rep(NA, length(dxs) * length(covs) * length(mtds1) * n_perms)
  # CNN vs FS
  cor.difs2 <- rep(NA, length(dxs) * length(covs))

  set.seed(1618)
  i <- j <- k <- 0
  pb <- progress_bar$new(format = "Permutations | :what [:bar] :current/:total",
                         total = length(dxs) * length(covs) * length(mtds1) *
                           n_perms + length(dxs) * length(covs) * n_perms,
                         clear = FALSE, width = 75)
  for (dx in dxs) {
    for (cov in covs) {
      for (mtd in mtds1) {
        wideDT      <- DT[DX == dx & METHOD == mtd]

        longDT      <- wideDT |>
                      melt(measure.vars   = msrs,
                           variable.name  = "MEASURE",
                           value.name     = "VAL1") |>
                      melt(measure.vars   = covs,
                           variable.name  = "COVAR",
                           value.name     = "VAL2")

        for (msr in msrs) {
          i <- i + 1
          corr    <- longDT[COVAR == cov & MEASURE == msr,
                            cor.test(VAL1, VAL2, method = "spearman")]
          conf    <- longDT[COVAR == cov & MEASURE == msr,
                            SpearmanRho(VAL1, VAL2, conf.level = 0.95)]
          r[i]    <- corr$estimate
          t[i]    <- corr$statistic
          pval[i] <- corr$p.value
          cil[i]  <- conf[2]
          cih[i]  <- conf[3]
        }

        for (p in 1:n_perms) {
          pb$tick(tokens = list(what = paste(dx, mtd, sep = ":")))
          longDT[, SHUFFLE := sample(MEASURE)]
          cor.1   <- longDT[COVAR == cov & SHUFFLE == "HVR",
                            cor(VAL1, VAL2, method = "spearman")]
          cor.2   <- longDT[COVAR == cov & SHUFFLE == "HCv",
                            cor(VAL1, VAL2, method = "spearman")]
          cor.difs1[p + j * n_perms] <- cor.1 - cor.2
        }
        j <- j + 1 # Increase counter
      }
      wideDT      <- DT[DX == dx & METHOD %in% c("cnn", "fs6")]

      longDT      <- wideDT |>
                    melt(measure.vars   = covs,
                         variable.name  = "COVAR",
                         value.name     = "VAL2")

      for (p in 1:n_perms) {
        pb$tick(tokens = list(what = paste(dx, "CNN vs FS", sep = ":")))
        longDT[, SHUFFLE := sample(METHOD)]
        cor.1   <- longDT[COVAR == cov & SHUFFLE == "cnn",
                          cor(HVR, VAL2, method = "spearman")]
        cor.2   <- longDT[COVAR == cov & SHUFFLE == "fs6",
                          cor(HVR, VAL2, method = "spearman")]
        cor.difs2[p + k * n_perms] <- cor.1 - cor.2
      }
      k <- k + 1
    }
  }

  corr.dt       <- data.table(DX      = rep(dxs, each = 3 * 4 * 2),
                              COVAR   = rep(covs, times = 3, each = 4 * 2),
                              METHOD  = rep(mtds2, times = 3 * 3, each = 2),
                              HC      = rep(msrs, times = 3 * 4 * 3),
                              R       = r,
                              Tstat   = t,
                              Pval    = pval,
                              DF      = dfs,
                              CIhigh  = cih,
                              CIlow   = cil)

  write_rds(corr.dt, fpaths[1])

  perms.dif1.dt <- data.table(DX      = rep(dxs, each = 3 * 4 * n_perms),
                              COVAR   = rep(covs, times = 3,
                                            each = 4 * n_perms),
                              METHOD  = rep(mtds2, times = 3 * 3,
                                            each = n_perms),
                              DIFF_p  = cor.difs1)

  perms.dif2.dt <- data.table(DX      = rep(dxs, each = 3 * n_perms),
                              COVAR   = rep(covs, times = 3, each = n_perms),
                              DIFF_p  = cor.difs2)

  perms.dif2.dt[, METHOD := "CNN-FS_V6"]

  perms.dif.dt  <- rbindlist(list(perms.dif1.dt, perms.dif2.dt),
                             use.names = TRUE)

  write_rds(perms.dif.dt, fpaths[2])
}

# R differences DTs
corr.hc.dt    <- corr.dt[, .(DX, COVAR, METHOD, HC, R)] |>
                dcast(... ~ HC, value.var = "R")

corr.seg.dt   <- corr.hc.dt[METHOD %in% c("CNN", "FS_V6"),
                            .(DX, COVAR, METHOD, HVR)] |>
                dcast(... ~ METHOD, value.var = "HVR")

corr.hc.dt[, DIFF := HVR - HCv][, c("HCv", "HVR") := NULL]
corr.seg.dt[, DIFF := CNN - FS_V6][, c("CNN", "FS_V6") := NULL]

corr.seg.dt[, METHOD := "CNN-FS_V6"]
corr.dif.dt   <- rbindlist(list(corr.hc.dt, corr.seg.dt), use.names = TRUE)
rm(corr.hc.dt, corr.seg.dt)


# Add Comparison labels
perms.dif.dt[, COMP := paste0(METHOD, ": HVR - HCv")]
perms.dif.dt[METHOD == "CNN-FS_V6", COMP := "HVR: CNN - FS_V6"]


## Plots
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Correlations
corr.dt[, COVAR := factor(COVAR, levels = c("AGE", "RAVLT_learning", "ADAS13"),
                          labels = c("Age", "Memory", "Cognition"))]
corr.dt[, Pval_adj := p.adjust(Pval, method = "bonferroni")]

corr.dt[Pval_adj < 0.05, SIGN := "*"]
corr.dt[Pval_adj < 0.01, SIGN := "**"]
corr.dt[Pval_adj < 0.001, SIGN := "***"]

# Factor DX
corr.perm.dt1 <- corr.dif.dt[perms.dif.dt,
                            on = .(DX, COVAR, METHOD)
                            ][COVAR != "RAVLT_learning",
                            .(Pval = sum(DIFF_p <= DIFF) / .N),
                            .(DX, COVAR, METHOD)]

corr.perm.dt2 <- corr.dif.dt[perms.dif.dt,
                            on = .(DX, COVAR, METHOD)
                            ][COVAR == "RAVLT_learning",
                            .(Pval = sum(DIFF_p <= DIFF) / .N),
                            .(DX, COVAR, METHOD)]

corr.perm.dt <- rbindlist(list(corr.perm.dt1, corr.perm.dt2))
rm(corr.perm.dt1, corr.perm.dt2)

# Create dt for permutation significance (only MCI&AD | CNN | AGE)
corr.perm.sign.dt <- corr.perm.dt[Pval < 0.05,
                                  .(DX, COVAR, METHOD, Pval,
                                    HCv = "HCv", HVR = "HVR", LABEL = "*")]

corr.perm.sign.dt[Pval < 0.01, LABEL := "**"]
corr.perm.sign.dt[Pval < 0.001, LABEL := "***"]

corr.perm.sign.dt[, COVAR := factor(COVAR,
                                    levels = c("AGE", "RAVLT_learning", "ADAS13"),
                                    labels = c("Age", "Memory", "Cognition"))]

corr.perm.sign.dt <- corr.dt[, .(Y = max(CIhigh)), .(DX, COVAR, METHOD)
                             ][corr.perm.sign.dt, on = .(DX, COVAR, METHOD)]


# FS_V6
p1  <- ggplot(corr.dt[METHOD == "FS_V6"], aes(HC, R, colour = HC)) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "none") +
  facet_grid(rows = vars(DX), cols = vars(COVAR)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             alpha = .5, colour = cbPalette[1]) +
  geom_errorbar(data = corr.dt[METHOD == "FS_V6" & HC == "HCv"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = .2)) +
  geom_point(data = corr.dt[METHOD == "FS_V6" & HC == "HCv"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "FS_V6" & HC == "HCv"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "FS_V6" & HC == "HCv"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = -.03, hjust = "right") +
  geom_errorbar(data = corr.dt[METHOD == "FS_V6" & HC == "HVR"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = -.2)) +
  geom_point(data = corr.dt[METHOD == "FS_V6" & HC == "HVR"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "FS_V6" & HC == "HVR"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "FS_V6" & HC == "HVR"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = .03, hjust = "left") +
  scale_colour_manual(values = cbPalette[2:3]) +
  ylim(-.6, .5) +
  labs(title = "FS_V6", x = "HC measure", y = "Spearman's rho",
       caption = glue("N = {DT[METHOD == 'fs6', .N]}",
                      "\n* p < 0.05; ** p < 0.01; *** p < 0.001"))

#here("plots/adni-bl_hcv-hvr_corrs_fs6.png") |>
  #ggsave(width = 4, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_fs6.tiff") |>
  #ggsave(width = 4, height = 7, units = "in",
           #device = "tiff", dpi = 600)

# CNN
p2  <- ggplot(corr.dt[METHOD == "CNN"], aes(HC, R, colour = HC)) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "none") +
  facet_grid(rows = vars(DX), cols = vars(COVAR)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             alpha = .5, colour = cbPalette[1]) +
  geom_errorbar(data = corr.dt[METHOD == "CNN" & HC == "HCv"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = .2)) +
  geom_point(data = corr.dt[METHOD == "CNN" & HC == "HCv"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "CNN" & HC == "HCv"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "CNN" & HC == "HCv"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = -.03, hjust = "right") +
  geom_errorbar(data = corr.dt[METHOD == "CNN" & HC == "HVR"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = -.2)) +
  geom_point(data = corr.dt[METHOD == "CNN" & HC == "HVR"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "CNN" & HC == "HVR"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "CNN" & HC == "HVR"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = .03, hjust = "left") +
  geom_signif(data = corr.perm.sign.dt[METHOD == "CNN"],
              aes(xmin = HCv, xmax = HVR, annotations = LABEL,
                  y_position = Y + .13), manual = TRUE, colour = cbPalette[1],
              textsize = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = cbPalette[2:3]) +
  ylim(-.6, .5) +
  labs(title = "CNN", x = "HC measure", y = "Spearman's rho",
       caption = glue("N = {DT[METHOD == 'cnn', .N]}",
                      "\n* p < 0.05; ** p < 0.01; *** p < 0.001"))

#here("plots/adni-bl_hcv-hvr_corrs_cnn.png") |>
  #ggsave(width = 4, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_cnn.tiff") |>
  #ggsave(width = 4, height = 7, units = "in",
           #device = "tiff", dpi = 600)

# MALF
p3  <- ggplot(corr.dt[METHOD == "MALF"], aes(HC, R, colour = HC)) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "none") +
  facet_grid(rows = vars(DX), cols = vars(COVAR)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             alpha = .5, colour = cbPalette[1]) +
  geom_errorbar(data = corr.dt[METHOD == "MALF" & HC == "HCv"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = .2)) +
  geom_point(data = corr.dt[METHOD == "MALF" & HC == "HCv"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "MALF" & HC == "HCv"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "MALF" & HC == "HCv"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = -.03, hjust = "right") +
  geom_errorbar(data = corr.dt[METHOD == "MALF" & HC == "HVR"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = -.2)) +
  geom_point(data = corr.dt[METHOD == "MALF" & HC == "HVR"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "MALF" & HC == "HVR"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "MALF" & HC == "HVR"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = .03, hjust = "left") +
  geom_signif(data = corr.perm.sign.dt[METHOD == "MALF"],
              aes(xmin = HCv, xmax = HVR, annotations = LABEL,
                  y_position = Y + .1), manual = TRUE, colour = cbPalette[1],
              textsize = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = cbPalette[2:3]) +
  ylim(-.6, .4) +
  labs(title = "MALF", x = "HC measure", y = "Spearman's rho",
       caption = glue("N = {DT[METHOD == 'malf', .N]}",
                      "\n* p < 0.05; ** p < 0.01; *** p < 0.001"))

#here("plots/adni-bl_hcv-hvr_corrs_malf.png") |>
  #ggsave(width = 4, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_malf.tiff") |>
  #ggsave(width = 4, height = 7, units = "in",
           #device = "tiff", dpi = 600)

# NLPB
p4  <- ggplot(corr.dt[METHOD == "NLPB"], aes(HC, R, colour = HC)) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "none") +
  facet_grid(rows = vars(DX), cols = vars(COVAR)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             alpha = .5, colour = cbPalette[1]) +
  geom_errorbar(data = corr.dt[METHOD == "NLPB" & HC == "HCv"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = .2)) +
  geom_point(data = corr.dt[METHOD == "NLPB" & HC == "HCv"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "NLPB" & HC == "HCv"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = .2)) +
  geom_text(data = corr.dt[METHOD == "NLPB" & HC == "HCv"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = -.03, hjust = "right") +
  geom_errorbar(data = corr.dt[METHOD == "NLPB" & HC == "HVR"],
                aes(ymin = CIlow, ymax = CIhigh), width = 0.2,
                position = position_nudge(x = -.2)) +
  geom_point(data = corr.dt[METHOD == "NLPB" & HC == "HVR"],
             shape = 21, fill = "white", size = 1.5, stroke = .5,
             position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "NLPB" & HC == "HVR"],
            aes(label = SIGN, y = CIhigh), size = 3, vjust = .1,
            position = position_nudge(x = -.2)) +
  geom_text(data = corr.dt[METHOD == "NLPB" & HC == "HVR"],
            aes(label = round(R, 2)),
            size = 3, nudge_x = .03, hjust = "left") +
  geom_signif(data = corr.perm.sign.dt[METHOD == "NLPB"],
              aes(xmin = HCv, xmax = HVR, annotations = LABEL,
                  y_position = Y + .13), manual = TRUE, colour = cbPalette[1],
              textsize = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = cbPalette[2:3]) +
  ylim(-.6, .5) +
  labs(title = "NLPB",
       x = "HC measure", y = "Spearman's rho",
       caption = glue("N = {DT[METHOD == 'nlpb', .N]}",
                      "\n* p < 0.05; ** p < 0.01; *** p < 0.001"))

#here("plots/adni-bl_hcv-hvr_corrs_nlpb.png") |>
  #ggsave(width = 4, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_nlpb.tiff") |>
  #ggsave(width = 4, height = 7, units = "in",
           #device = "tiff", dpi = 600)

p <- grid.arrange(p2, p3, p4, p1, nrow = 2)
here("plots/adni-bl_hcv-hvr_corrs.png") |>
  ggsave(p, width = 9, height = 9, units = "in", dpi = 600)

here("plots/adni-bl_hcv-hvr_corrs.tiff") |>
  ggsave(p, width = 9, height = 9, units = "in", device = "tiff", dpi = 600)

### Permutation tests
## Age
#plot.dt <- corr.dif.dt[perms.dif.dt[COVAR == "AGE" & METHOD != "CNN-FS_V6"],
                       #on = .(DX, COVAR, METHOD)]
#plot.dt <- corr.perm.dt[plot.dt, on = .(DX, COVAR, METHOD)]

#ggplot(plot.dt, aes(x = DIFF_p, y = DX, fill = factor(after_stat(quantile)))) +
  #theme_classic(base_size = 12) +
  #theme(text = element_text(size = 12), axis.text.y = element_blank(),
        #axis.ticks.y = element_blank()) +
  #facet_grid(rows = vars(DX), cols = vars(COMP), scales = "free_y") +
  #stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      #quantiles = 0.05, scale = 1, alpha = .3) +
  #geom_vline(aes(xintercept = DIFF), colour = cbPalette[3]) +
  #scale_fill_manual(values = cbPalette[2:1], name = "One-sided\nhypothesis",
                    #labels = c("lower 5%", "upper 95%")) +
  #geom_richtext(data = unique(plot.dt[, .(DX, COMP, Pval)]),
                #aes(label = paste0("<i>p</i> = ", Pval)),
                #inherit.aes = F, colour = "Black", fill = "White",
                #size = 2.5, x = 0, y = -Inf, vjust = -0.25) +
  #labs(title = "Permutation tests: Difference in correlation with Age",
       #x = "Difference in r", y = NULL)

#here("plots/adni-bl_hcv-hvr_corrs_perms_age.png") |>
  #ggsave(width = 13, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_perms_age.tiff") |>
  #ggsave(width = 13, height = 7, units = "in",
           #device = "tiff", dpi = 600)

## Memory (RAVLT_learning)
## Correlation is in opposite direction
#plot.dt <- corr.dif.dt[perms.dif.dt[COVAR == "RAVLT_learning" & METHOD != "CNN-FS_V6"],
                       #on = .(DX, COVAR, METHOD)]
#plot.dt <- corr.perm.dt[plot.dt, on = .(DX, COVAR, METHOD)]

#ggplot(plot.dt, aes(x = DIFF_p, y = DX, fill = factor(after_stat(quantile)))) +
  #theme_classic(base_size = 12) +
  #theme(text = element_text(size = 12), axis.text.y = element_blank(),
        #axis.ticks.y = element_blank()) +
  #facet_grid(rows = vars(DX), cols = vars(COMP), scales = "free_y") +
  #stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      #quantiles = 0.95, scale = 1, alpha = .3) +
  #geom_vline(aes(xintercept = DIFF), colour = cbPalette[3]) +
  #scale_fill_manual(values = cbPalette[1:2], name = "One-sided\nhypothesis",
                    #labels = c("lower 95%", "upper 5%")) +
  #geom_richtext(data = unique(plot.dt[, .(DX, COMP, Pval)]),
                #aes(label = paste0("<i>p</i> = ", Pval)),
                #inherit.aes = F, colour = "Black", fill = "White",
                #size = 2.5, x = 0, y = -Inf, vjust = -0.25) +
  #labs(title = "Permutation tests: Difference in correlation with Memory",
       #x = "Difference in r", y = NULL)

#here("plots/adni-bl_hcv-hvr_corrs_perms_mem.png") |>
  #ggsave(width = 13, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_perms_mem.tiff") |>
  #ggsave(width = 13, height = 7, units = "in",
           #device = "tiff", dpi = 600)

## Cognition
#plot.dt <- corr.dif.dt[perms.dif.dt[COVAR == "ADAS13" & METHOD != "CNN-FS_V6"],
                       #on = .(DX, COVAR, METHOD)]
#plot.dt <- corr.perm.dt[plot.dt, on = .(DX, COVAR, METHOD)]

#ggplot(plot.dt, aes(x = DIFF_p, y = DX, fill = factor(after_stat(quantile)))) +
  #theme_classic(base_size = 12) +
  #theme(text = element_text(size = 12), axis.text.y = element_blank(),
        #axis.ticks.y = element_blank()) +
  #facet_grid(rows = vars(DX), cols = vars(COMP), scales = "free_y") +
  #stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      #quantiles = 0.05, scale = 1, alpha = .3) +
  #geom_vline(aes(xintercept = DIFF), colour = cbPalette[3]) +
  #scale_fill_manual(values = cbPalette[2:1], name = "One-sided\nhypothesis",
                    #labels = c("lower 5%", "upper 95%")) +
  #geom_richtext(data = unique(plot.dt[, .(DX, COMP, Pval)]),
                #aes(label = paste0("<i>p</i> = ", Pval)),
                #inherit.aes = F, colour = "Black", fill = "White",
                #size = 2.5, x = 0, y = -Inf, vjust = -0.25) +
  #labs(title = "Permutation tests: Difference in correlation with Cognition",
       #x = "Difference in r", y = NULL)

#here("plots/adni-bl_hcv-hvr_corrs_perms_cog.png") |>
  #ggsave(width = 13, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_perms_cog.tiff") |>
  #ggsave(width = 13, height = 7, units = "in",
           #device = "tiff", dpi = 600)

## HVR: CNN - FS6
## Correlation is in opposite direction
#plot.dt <- corr.dif.dt[perms.dif.dt[METHOD == "CNN-FS_V6"],
                       #on = .(DX, COVAR, METHOD)]
#plot.dt <- corr.perm.dt[plot.dt, on = .(DX, COVAR, METHOD)]

#ggplot(plot.dt, aes(x = DIFF_p, y = DX, fill = factor(after_stat(quantile)))) +
  #theme_classic(base_size = 12) +
  #theme(text = element_text(size = 12), axis.text.y = element_blank(),
        #axis.ticks.y = element_blank()) +
  #facet_grid(rows = vars(DX), cols = vars(COVAR), scales = "free_y") +
  #stat_density_ridges(data = plot.dt[COVAR == "RAVLT_learning"],
                      #geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      #quantiles = 0.95, scale = 1, alpha = .3) +
  #stat_density_ridges(data = plot.dt[COVAR != "RAVLT_learning"],
                      #geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      #quantiles = 0.95, scale = 1, alpha = .3) +
  #geom_vline(aes(xintercept = DIFF), colour = cbPalette[3]) +
  #scale_fill_manual(values = cbPalette[1:2], name = "One-sided\nhypothesis",
                    #labels = c("lower 95%", "upper 5%")) +
  #geom_richtext(data = unique(plot.dt[, .(DX, COMP, Pval)]),
                #aes(label = paste0("<i>p</i> = ", Pval)),
                #inherit.aes = F, colour = "Black", fill = "White",
                #size = 2.5, x = 0, y = -Inf, vjust = -0.25) +
  #labs(title = "Permutation tests: Difference in correlation with Memory",
       #x = "Difference in r", y = NULL)

#here("plots/adni-bl_hcv-hvr_corrs_perms_mem.png") |>
  #ggsave(width = 13, height = 7, units = "in", dpi = 600)

#here("plots/adni-bl_hcv-hvr_corrs_perms_mem.tiff") |>
  #ggsave(width = 13, height = 7, units = "in",
           #device = "tiff", dpi = 600)


## Correlations with Memory
DT_long <- melt(DT,
                measure.vars = patterns("^H"),
                variable.name = "HC_msr",
                value.name = "VOL")

mem_hc_cnn <- DT_long[METHOD == "cnn",
                      .(DX, HC_msr, VOL, RAVLT_learning)]

ggplot(mem_hc_cnn, aes(x = RAVLT_learning, y = VOL, colour = DX)) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "bottom") +
  geom_point(size = 2, shape = 21) +
  geom_abline(intercept = 0, slope = 1,
              colour = cbPalette[1], linetype = "dashed") +
  geom_smooth(method = "lm", alpha = .2) +
  stat_cor(size = 2.7, label.x.npc = "right", label.y.npc = "bottom",
           hjust = "inward") +
  facet_wrap(vars(HC_msr, DX), ncol = 3, scales = "free") +
  scale_colour_manual(values = cbPalette[-1]) +
  labs(x = "Memory", y = "Volume",
       colour = "Clinical label")

ggsave("plots/adni-bl_hcv_hvr_memory_corrs.png", width = 12, height = 24,
       units = "in", dpi = 600)

