#!/usr/bin/env Rscript
library(here)
library(readr)
library(data.table)
library(glue)

## Parse volumes
fpath         <- here("data/rds/adni-bl_volumes_hcvc.rds")
if (file.exists(fpath)) {
  volumes     <- fpath |> read_rds()
} else {
  here("code/data_parsing/qc_segmentations_adni-bl.R") |> source()
}

# ADNI FSvols
fpath         <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge   <- fpath |> read_rds()
} else {
  here("code/data_parsing/parse_adnimerge-bl.R") |> source()
}

adni_vols     <- adnimerge[, .(PTID, RID, SCANDATE, Hippocampus, FSVERSION)]

# UCSF HCvols
fpath         <- here("data/UCSFFSX_11_02_15_20Nov2023.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()
fs4_vols      <- fread(fpath, select = c("RID", "OVERALLQC",
                                         "ST29SV", # LeftHC
                                         "ST88SV"))# RightHC

fsv4          <- "Cross-Sectional FreeSurfer (FreeSurfer Version 4.3)"
fs4_vols      <- fs4_vols[, .(RID, OVERALLQC, LHC = ST29SV, RHC = ST88SV,
                              Hippocampus = ST29SV + ST88SV)
                          ][adni_vols[FSVERSION == fsv4],
                            on = .(RID, Hippocampus)]
fs4_vols[, UCSFFS := 4.3]


fpath         <- here("data/UCSFFSX51_11_08_19_20Nov2023.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()
fs5_vols      <- fread(fpath, select = c("RID", "OVERALLQC",
                                         "LHIPQC", "RHIPQC",
                                         "ST29SV", # LeftHC
                                         "ST88SV"))# RightHC

fsv5          <- "Cross-Sectional FreeSurfer (5.1)"
fs5_vols      <- fs5_vols[, .(RID, OVERALLQC, LHIPQC, RHIPQC,
                              LHC = ST29SV, RHC = ST88SV,
                              Hippocampus = ST29SV + ST88SV)
                          ][adni_vols[FSVERSION == fsv5],
                            on = .(RID, Hippocampus)]
fs5_vols[, UCSFFS := 5.1]

ucsf_vols     <- rbindlist(list(fs4_vols[!is.na(Hippocampus),
                                .(PTID, SCANDATE, UCSFFS, Hippocampus)],
                                fs5_vols[!is.na(Hippocampus),
                                .(PTID, SCANDATE, UCSFFS, Hippocampus)]))

# House FSvols
fpath         <- here("data/ADNI_FS_hc.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()
fs6_vols1     <- fread(fpath, select = c(1, 3, 5:7),
                       col.names = c("PTID", "DATE", "LHC", "RHC", "BRAIN"))

fpath         <- here("data/ADNI_FS_hc_vc.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()
fs6_vols2     <- fread(fpath, select = c(1, 3, 8:9),
                      col.names = c("PTID", "DATE", "LCSF", "RCSF"))

fs6_vols      <- fs6_vols1[fs6_vols2, on = .(PTID, DATE)]
fs6_vols[, FS_house := LHC + RHC]

# Keep only subjects that went through the segmentation algorithms
fs6_vols      <- fs6_vols[volumes[METHOD == "cnn", .(PTID, SCANDATE)],
                          on = .(PTID, DATE = SCANDATE)]

# Merge
fs_vols       <- ucsf_vols[fs6_vols, on = "PTID",
                          .(PTID, SCANDATE = DATE,
                            LHC, RHC, HC = LHC + RHC,
                            LCSF, RCSF, CSF = LCSF + RCSF,
                            BRAIN, UCSFFS, FS_house, FS_ucsf = Hippocampus)]

fs_vols       <- unique(fs_vols)

# Versions
fs_vols[, .N, UCSFFS]
## 5.1: 701
## 4.3: 545

## Save RDS
outdir            <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
write_rds(fs_vols, here(outdir, "adni-bl_volumes_freesurfer.rds"))

## Clean workspace
rm(adnimerge, adni_vols, discarded, fpath, fsv4, fs4_vols, fsv5, fs5_vols,
   fs6_vols1, fs6_vols2, fs6_vols, outdir, ucsf_vols)
