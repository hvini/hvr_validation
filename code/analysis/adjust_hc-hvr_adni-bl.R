#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)


### INPUT
fpaths <- list(
  CSV = here("data/derivatives/adni_icc_scale.csv"),
  RDS = c(
    "adni-bl_volumes_hcvc", "adni-bl_volumes_freesurfer", "adnimerge_baseline"
  ) |> sprintf(fmt = "data/rds/%s.rds") |> here(),
  SRC = c(
    "qc_segmentations_adni-bl", "parse_freesurfer-vols", "parse_adnimerge-bl"
  ) |> sprintf(fmt = "code/data_parsing/%s.R") |> here()
)

## Dependency
if (!file.exists(fpaths$CSV)) fpaths$CSV |>
  sprintf(fmt = "Required file(s) could not be found: %s.") |>
  stop()

## Parse DATA
# Curated HC CSF volumes
if (file.exists(fpaths$RDS[1])) {
  volumes <- readRDS(fpaths$RDS[1])
} else {
  source(fpaths$SRC[1])
  volumes <- data.lst$VOL
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

## Adnimerge
if (file.exists(fpaths$RDS[3])) {
  adnimerge <- readRDS(fpaths$RDS[3])
} else {
  source(fpaths$SRC[3])
  adnimerge <- data.lst$ADNIMERGE
  rm(data.lst)
}

### Parse DATA
data.lst <- list(
  ICC_SCL =  fread(fpaths$CSV) |> (\(DT) DT[, SCANDATE := ymd(SCANDATE)])(),
  SEGMS = volumes,
  FS = fs_vols[
    , c("BRAIN", "UCSFFS", "FS_house", "FS_ucsf") := NULL
  ][
    , METHOD := "fs6"
  ],
  CTRLS = adnimerge[DX == "CH", PTID]
)

rm(fpaths, volumes, fs_vols, adnimerge)

## Merge SEGMS & FSICC & SCALE
## Filter QC
data.lst$VOLS <- data.lst$SEGMS["Pass", on = "QC", -"QC"] |>
  rbind(data.lst$FS, use.names = TRUE) |>
  merge(data.lst$ICC_SCL, by = c("PTID", "SCANDATE")) |>
  na.omit() |>
  (\(DT) DT[, ICC := ICC / 1000])() # Bring ICC to CC

### Head-size Adjustments
data.lst$ADJ <- list()

## Unadjusted
# Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
# Scale everything to CC
data.lst$ADJ$NAT <- rbind(
  data.lst$VOLS[
    "fs6", on = "METHOD",
    .(
      HC_l  = LHC  / 1000,
      HC_r  = RHC  / 1000,
      CSF_l = LCSF / 1000,
      CSF_r = RCSF / 1000
    ),
    .(PTID, METHOD, ICC, SCALEFACTOR) #TODO: Do I need all of these?
  ],
  data.lst$VOLS[
    !"fs6", on = "METHOD",
    .(
      HC_l  = LHC  / (SCALEFACTOR * 1000),
      HC_r  = RHC  / (SCALEFACTOR * 1000),
      CSF_l = LCSF / (SCALEFACTOR * 1000),
      CSF_r = RCSF / (SCALEFACTOR * 1000)
    ),
    .(PTID, METHOD, ICC, SCALEFACTOR) #TODO: Do I need all of these?
  ]
)

## STX
data.lst$ADJ$STX <- rbind(
  data.lst$VOLS[
    "fs6", on = "METHOD",
    .(HC_l = LHC * SCALEFACTOR / 1000, HC_r = RHC * SCALEFACTOR / 1000),
    .(PTID, METHOD, ICC, SCALEFACTOR) #TODO: Do I need all of these?
  ],
  data.lst$VOLS[
    !"fs6", on = "METHOD",
    .(HC_l = LHC / 1000, HC_r = RHC / 1000),
    .(PTID, METHOD, ICC, SCALEFACTOR)
  ]
)

## HVR
data.lst$ADJ$HVR <- data.lst$ADJ$NAT[
  ,
  .(HVR_l = HC_l  / (HC_l + CSF_l), HVR_r = HC_r / (HC_r + CSF_r)),
  .(PTID, METHOD, ICC, SCALEFACTOR)
]

## Robust methods
## Residuals:
## VOL_adj = VOL - b(ICC - ICC_cn)
data.lst$ADJ$RES <- data.lst$ADJ$NAT |>
  melt(measure = patterns("_(l|r)$"), variable = "ROI", value = "CC") |>
  # Mean ICC (Cognitively Healthy only)
  merge(
    data.lst$VOLS[
      data.lst$CTRLS,
      on = "PTID",
      nomatch = NULL,
      .(ICC_ch = mean(ICC)),
      METHOD
    ],
    by = "METHOD"
  ) |>
  # Residual slope
  merge(
    data.lst$ADJ$NAT[
      data.lst$CTRLS, on = "PTID",
      .(METHOD, ICC, HC_l, HC_r, CSF_l, CSF_r)
    ] |>
      melt(id = 1:2, variable = "ROI", value = "CC") |>
      na.omit() |>
      (\(DT) DT[
        ,
        .(b = summary(lm(CC ~ ICC))$coefficients[2]),
        .(ROI, METHOD)
      ])(),
    by = c("METHOD", "ROI")
  ) |>
  # Adjustment
  (\(DT) DT[, CC := CC - b * (ICC - ICC_ch)][, c("ICC_ch", "b") := NULL])() |>
  dcast(... ~ ROI, value.var = "CC") |>
  setcolorder(c("PTID", "METHOD", "ICC", "SCALEFACTOR"))


## Power-corrected proportion:
# VOL_adj = VOL / ICC ** b
# b: slope of log(VOL) ~ log(ICC)
data.lst$ADJ$PCP <- data.lst$ADJ$NAT |>
  melt(measure = patterns("_(l|r)$"), variable = "ROI", value = "CC") |>
  # PCP slope
  merge(
    data.lst$ADJ$NAT[
      data.lst$CTRLS, on = "PTID",
      .(METHOD, ICC, HC_l, HC_r, CSF_l, CSF_r)
    ] |>
      melt(id = 1:2, variable = "ROI", value = "CC") |>
      na.omit() |>
      (\(DT) DT[
        ,
        .(b = summary(lm(log(CC) ~ log(ICC)))$coefficients[2]),
        .(ROI, METHOD)
      ])(),
    by = c("METHOD", "ROI")
  ) |>
  # Adjustment
  (\(DT) DT[, CC := CC / ICC ** b][, b := NULL])() |>
  dcast(... ~ ROI, value.var = "CC") |>
  setcolorder(c("PTID", "METHOD", "ICC", "SCALEFACTOR"))

# HVR (Robust)
data.lst$ADJ[c("HVR_RES", "HVR_PCP")] <- lapply(
  data.lst$ADJ[c("RES", "PCP")],
  \(DT) DT[
    ,
    .(HVR_l = HC_l  / (HC_l + CSF_l), HVR_r = HC_r / (HC_r + CSF_r)),
    .(PTID, METHOD, ICC, SCALEFACTOR)
  ]
)

## Mean by sides
lapply(
  data.lst$ADJ,
  \(DT) {
    # HC
    if ("HC_l" %in% names(DT)) DT[, HC_mean := (HC_l + HC_r) / 2]
    if ("CSF_L" %in% names(DT)) DT[, CSF_mean := (CSF_l + CSF_r) / 2]
    if ("HVR_l" %in% names(DT)) DT[, HVR_mean := (HVR_l + HVR_r) / 2]
    setkey(DT, PTID, METHOD)
  }
) |> invisible()

### OUTPUT
outdir <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
outdir |>
  here("adni-bl_volumes_icv-adjusted.rds") |>
  saveRDS(object = data.lst$ADJ)
rm(outdir)
