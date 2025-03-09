#!/usr/bin/env Rscript

library(here)
library(readr)
library(glue)
library(data.table)
library(lubridate)
library(stringr)

## Dependencies
# ICC volume and ScaleFactors
fpath           <- here("data/derivatives/adni_icc_scale.csv")
if (!file.exists(fpath)) {
  sprintf("File: %s is required but could not be found.", fpath) |> stop()
}
icc_scales      <- fread(fpath)
icc_scales[, SCANDATE := ymd(SCANDATE)]

# Curated HC CSF volumes
fpath       <- here("data/rds/adni-bl_volumes_hcvc-ag.rds")
if (file.exists(fpath)) {
  volumes_hcvcag <- read_rds(fpath)
} else {
  here("code/data_parsing/qc_segmentations_adni-bl.R") |> source()
  rm(volumes) # Unused
}

## Controls
fpath    <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge     <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_adnimerge-bl.R') |> source()
}

controls    <- adnimerge[DX == "CH", PTID]
rm(adnimerge, fpath)

# Parse volumes from CNN
volumes     <- volumes_hcvcag[QC != "Fail",
                              .(PTID, SCANDATE,
                                LHC   = L_HC_T + L_HC_B + L_HC_H,
                                RHC   = R_HC_T + R_HC_B + R_HC_H,
                                LCSF  = L_VC_T + L_VC_B + L_VC_H,
                                RCSF  = R_VC_T + R_VC_B + R_VC_H,
                                LAMY  = L_AMY,
                                RAMY  = R_AMY)]

rm(volumes_hcvcag)

# Merge ICC and volumes
volumes     <- icc_scales[volumes, on = .(PTID, SCANDATE)]
rm(icc_scales)
volumes[, ICC := ICC / 1000]
icc_cn      <- volumes[PTID %in% controls, mean(ICC)]

## Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
# Scale everything to cm^3
volumes[, `:=`(HC_l_raw   = LHC   / (SCALEFACTOR * 1000),
               HC_r_raw   = RHC   / (SCALEFACTOR * 1000),
               CSF_l_raw  = LCSF  / (SCALEFACTOR * 1000),
               CSF_r_raw  = RCSF  / (SCALEFACTOR * 1000),
               AMY_l_raw  = LAMY  / (SCALEFACTOR * 1000),
               AMY_r_raw  = RAMY  / (SCALEFACTOR * 1000),
               HC_l_stx   = LHC   / 1000,
               HC_r_stx   = RHC   / 1000)]

volumes[, c("LHC", "RHC", "LCSF", "RCSF", "LAMY", "RAMY") := NULL]

## Regression slopes for PCP & Residual normalizations
# Use only Controls for the models
cols        <- str_subset(names(volumes), "ICC|_raw")
volumes_lng <- volumes[PTID %in% controls, ..cols] |>
  melt(id = "ICC", variable.name = "ROI")
volumes_lng[, ROI := str_sub(ROI, end = -5)]
rm(controls, cols)

## Power-corrected proportion:
# VOL_adj = VOL / ICC ** b
# b: slope of log(VOL) ~ log(ICC)
b_pcp       <-
  volumes_lng[value != 0,
              .(b_pcp = summary(lm(log(value) ~ log(ICC)))$coefficients[2]),
              ROI]

## Residuals
# Remove the residuals from VOL ~ ICC regression
# VOL_adj = VOL - b(ICC - ICC_cn)
b_res       <-
  volumes_lng[, .(b_res = summary(lm(value ~ ICC))$coefficients[2]), ROI]

volumes_lng <- volumes |>
  melt(measure = patterns("_l|_r"), variable = "ROI")
volumes_lng[, ADJ := str_sub(ROI, -3)]
volumes_lng[, ROI := str_sub(ROI, end = -5)]
volumes_lng <- volumes_lng |>
  dcast(... ~ ADJ, value.var = "value")

## Apply adjustment methods
volumes_lng <- b_pcp[b_res, on = "ROI"][volumes_lng, on = "ROI"]
volumes_lng[, SIDE := ifelse(grepl("_l$", ROI), "L", "R")]
volumes_lng[, ROI := str_sub(ROI, end = -3)]

volumes_lng[, `:=`(prop = raw / ICC,
                   pcp = raw / ICC ** b_pcp,
                   res = raw - b_res * (ICC - icc_cn))]

rm(b_pcp, b_res, icc_cn)

volumes     <- volumes_lng[, -c("b_pcp", "b_res", "ICC", "SCALEFACTOR")] |>
  melt(measure = c("raw", "stx", "prop", "pcp", "res"), variable = "ADJ") |>
  dcast(PTID + SCANDATE + SIDE + ADJ ~ ROI)

rm(volumes_lng)

## HVR + HAVR + HAVAS
volumes[!ADJ %in% c("stx", "prop"),
        `:=`(HVR    = HC / (HC + CSF),
             HAVR   = (HC + AMY) / (HC + AMY + CSF),
             HAVAS  = HC + AMY + CSF)]

volumes     <- volumes[, -c("AMY", "CSF")] |>
  melt(measure = patterns("^H"), variable = "MEASURE", value = "VAL")

volumes     <- volumes[!is.na(VAL)]

## Export RDS
write_rds(volumes, here("data/rds/adni-bl_hc-msrs_icv-adjusted.rds"))
