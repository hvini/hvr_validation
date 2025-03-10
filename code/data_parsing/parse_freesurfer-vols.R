#!/usr/bin/env Rscript
library(here)
library(data.table)

### INPUT
fpaths <- list(
  RDS = c("adni-bl_volumes_hcvc", "adnimerge_baseline") |>
    sprintf(fmt = "data/rds/%s.rds") |> here(),
  CSV = c(
    "UCSFFSX_11_02_15_20Nov2023",
    "UCSFFSX51_11_08_19_20Nov2023",
    "ADNI_FS_hc",
    "ADNI_FS_hc_vc"
  ) |> sprintf(fmt = "data/%s.csv") |> here(),
  SRC = c("qc_segmentations_adni-bl", "parse_adnimerge-bl") |>
    sprintf(fmt = "code/data_parsing/%s.R") |> here()
)

## Missing CSV files
missing_csv <- fpaths$CSV[!file.exists(fpaths$CSV)]
if (length(missing_csv) > 0) {
  missing_csv |>
    paste(collapse = ", ") |>
    sprintf(fmt = "Required file(s) could not be found: %s.") |>
    stop()
}
rm(missing_csv)

## Subjects that passed CNN segmentation
if (file.exists(fpaths$RDS[1])) {
  volumes <- readRDS(fpaths$RDS[1])
} else {
  source(fpaths$SRC[1])
}

## ADNI FSvols
if (file.exists(fpaths$RDS[2])) {
  adnimerge <- readRDS(fpaths$RDS[2])
} else {
  source(fpaths$SRC[2])
}

### Parse DATA
data.lst <- list(
  CNN_SUBS = volumes["cnn", on = "METHOD", .(PTID, SCANDATE)],
  ADNIMERGE = adnimerge[, .(PTID, RID, SCANDATE, Hippocampus, FSVERSION)],
  FS4 = fread(
    fpaths$CSV[1],
    # LeftHC: ST29SV; RightHC; ST88SV
    select = c("RID", "OVERALLQC", "ST29SV", "ST88SV")
  ) |> setnames(c("ST29SV", "ST88SV"), c("LHC", "RHC")),
  FS5 = fread(
    fpaths$CSV[2],
    # LeftHC: ST29SV; RightHC; ST88SV
    select = c("RID", "OVERALLQC", "LHIPQC", "RHIPQC", "ST29SV", "ST88SV")
  ) |> setnames(c("ST29SV", "ST88SV"), c("LHC", "RHC")),
  FS6 = merge(
    fread(
      fpaths$CSV[3],
      select = c(1, 3, 5:7),
      col.names = c("PTID", "DATE", "LHC", "RHC", "BRAIN"),
      key = c("PTID", "DATE")
    ),
    fread(
      fpaths$CSV[4],
      select = c(1, 3, 8:9),
      col.names = c("PTID", "DATE", "LCSF", "RCSF"),
      key = c("PTID", "DATE")
    )
  )
)

rm(fpaths, adnimerge, volumes)

## UCSF HCvols
# Calculate whole HC
# Join with ADNIMERGE
# Set column for FS version
data.lst[["UCSF"]] <- Map(
  \(DT, Version_str, Version_number) {
    DT[, let(Hippocampus = LHC + RHC, UCSFFS = Version_number)]
    DT[
      data.lst$ADNIMERGE[Version_str, on = "FSVERSION"],
      on = .(RID, Hippocampus)
    ][
      !is.na(Hippocampus), .(PTID, SCANDATE, UCSFFS, Hippocampus)
    ]
  },
  data.lst[c("FS4", "FS5")],
  c(
    "Cross-Sectional FreeSurfer (FreeSurfer Version 4.3)",
    "Cross-Sectional FreeSurfer (5.1)"
  ),
  c(4.3, 5.1)
) |> rbindlist()


## House FSvols
# Calculate whole HC
# Filter subjects who have CNN segmentation
data.lst$FS6 <- (
  \(DT)
  DT[
    , FS_house := LHC + RHC
  ][
    data.lst$CNN_SUBS, on = .(PTID, DATE = SCANDATE)
  ]
)(data.lst$FS6)

## Merge FS
data.lst$FS <- data.lst$UCSF[
  data.lst$FS6, on = "PTID",
  .(
    PTID,
    SCANDATE = DATE,
    LHC,
    RHC,
    HC = LHC + RHC,
    LCSF,
    RCSF,
    CSF = LCSF + RCSF,
    BRAIN,
    UCSFFS,
    FS_house,
    FS_ucsf = Hippocampus
  )
] |> unique()

## Versions
data.lst$FS[, .N, UCSFFS]
#    UCSFFS     N
#     <num> <int>
# 1:    4.3   635
# 2:     NA   227
# 3:    5.1   779

### OUTPUT
outdir <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
here(outdir, "adni-bl_volumes_freesurfer.rds") |>
  saveRDS(object = data.lst$FS)
rm(outdir)
