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
    unique(volumes$PTID),
    on = "PTID"
  ][
    , let(METHOD = fifelse(FSv == 6, "fs6", paste0("fs", FSv)), FSv = NULL)
  ],
  CTRLS = adnimerge[DX == "CH", PTID]
)

rm(fpaths, volumes, fs_vols, adnimerge)

## Merge SEGMS & FSICC & SCALE
## Filter QC
#data.lst$VOLS <- data.lst$SEGMS["Pass", on = "QC", -"QC"] |>
data.lst$VOLS <- data.lst$SEGMS |>
  rbind(data.lst$FS[, QC := "Pass"], use.names = TRUE) |>
  merge(data.lst$ICC_SCL, by = c("PTID", "SCANDATE")) |>
  # Bring ICC to CC
  (\(DT) DT[, ICC := ICC / 1000])()

## Add Failures/missing to FS
data.lst$VOLS <- data.lst$VOLS |>
  # Failures on FS are missing values
  rbind(
    data.lst$VOLS[
      "cnn", on = "METHOD"
    ][
      !data.lst$VOLS["fs6", on = "METHOD"],
      .(PTID, SCANDATE, ICC, SCALEFACTOR, METHOD = "fs6", QC = "Fail")
    ],
    use.names = TRUE, fill = TRUE
  ) |>
  # Missing values on ADNI are just missing, but still keep
  rbind(
    data.lst$VOLS[
      "cnn", on = "METHOD"
    ][
      !data.lst$VOLS[c("fs4.3", "fs5.1"), on = "METHOD"],
      .(PTID, SCANDATE, ICC, SCALEFACTOR, METHOD = "fs", QC = "Missing")
    ],
    use.names = TRUE, fill = TRUE
  )


### Head-size Adjustments
data.lst$ADJ <- list()

## Unadjusted
# Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
# Scale everything to CC
data.lst$ADJ$NAT <- rbind(
  data.lst$VOLS[
    METHOD %like% "fs",
    #c("fs4.3", "fs5.1", "fs6"), on = "METHOD",
    .(
      QC,
      HC_l  = LHC  / 1000,
      HC_r  = RHC  / 1000,
      CSF_l = LCSF / 1000,
      CSF_r = RCSF / 1000
    ),
    .(PTID, METHOD) #TODO: Do I need all of these?
  ],
  data.lst$VOLS[
    !METHOD %like% "fs",
    .(
      QC,
      HC_l  = LHC  / (SCALEFACTOR * 1000),
      HC_r  = RHC  / (SCALEFACTOR * 1000),
      CSF_l = LCSF / (SCALEFACTOR * 1000),
      CSF_r = RCSF / (SCALEFACTOR * 1000)
    ),
    .(PTID, METHOD) #TODO: Do I need all of these?
  ]
)

## STX method
data.lst$ADJ$HC <- rbind(
  data.lst$VOLS[
    METHOD %like% "fs",
    .(QC, HC_l = LHC * SCALEFACTOR / 1000, HC_r = RHC * SCALEFACTOR / 1000),
    .(PTID, METHOD) #TODO: Do I need all of these?
  ],
  data.lst$VOLS[
    !METHOD %like% "fs",
    .(QC, HC_l = LHC / 1000, HC_r = RHC / 1000),
    .(PTID, METHOD)
  ]
)

## HVR
data.lst$ADJ$HVR <- data.lst$ADJ$NAT[
  !is.na(CSF_l) & !is.na(CSF_r),
  .(HVR_l = HC_l  / (HC_l + CSF_l), HVR_r = HC_r / (HC_r + CSF_r)),
  .(PTID, METHOD)
]

## Mean by sides
lapply(
  data.lst$ADJ,
  \(DT) {
    # HC
    if ("HC_l" %in% names(DT)) DT[, HC_mean := (HC_l + HC_r) / 2]
    if ("CSF_l" %in% names(DT)) DT[, CSF_mean := (CSF_l + CSF_r) / 2]
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
