#!/usr/bin/env Rscript
library(here)
library(data.table)
library(stringr)
library(lubridate)

### INPUT
fpaths <- list(
  RDS = "adnimerge_baseline" |>
    sprintf(fmt = "data/rds/%s.rds") |> here(),
  LST = c(
    "adni_acquisition_failures.lst",
    "qrater_malf_2022-12-20.csv",
    "qrater_malf_reg_fails_2022-12-20.csv",
    "qrater_nlpb_2022-12-20.csv",
    "qrater_nlpb_reg_fails_2022-12-20.csv",
    "qrater_cnn1_2022-12-20.csv",
    "qrater_cnn1_reg_fails_2022-12-20.csv"
  ) |> sprintf(fmt = "lists/%s") |> here(),
  VOL = c(
    "adni-bl_volumes_hcvc_malf",
    "adni-bl_volumes_hcvc_nlpb",
    "adni-bl_volumes_hcvc_cnn"
  ) |> sprintf(fmt = "data/derivatives/%s.csv") |> here(),
  SRC = "parse_adnimerge-bl" |>
    sprintf(fmt = "code/data_parsing/%s.R") |> here()
)

## Missing files
dependencies <- c(fpaths$LST, fpaths$VOL)
missing_files <- dependencies[!file.exists(dependencies)]
if (length(missing_files) > 0) {
  missing_files |>
  paste(collapse = "\n") |>
  sprintf(fmt = "Required file(s) could not be found:\n%s.") |>
  stop()
}
rm(missing_files, dependencies)

## ADNIMERGE
if (file.exists(fpaths$RDS)) {
  adnimerge_subs <- fpaths$RDS |> readRDS() |> (\(DT) DT$PTID)()
} else {
  source(fpaths$SRC)
  adnimerge_subs <- data.lst$ADNIMERGE$PTID
  rm(data.lst)
}

### Parse DATA
data.lst <- list(QC = list(), VOL = list(), DROPPED = list())

## Acquisition Fails
data.lst$QC$ACQ <- fpaths$LST[1] |>
  fread(header = FALSE, col.names = "FAILS") |>
  (
    \(DT)
    DT[
      , let(
        PTID = str_extract(FAILS, "\\d{3}_S_\\d{4}"),
        SCANDATE = ymd(str_extract(FAILS, "(?<=\\d{4}_)S*\\d+"))
      )
    ]
  )()

## Qrater QC
data.lst$QC[c("MALF", "NLPB", "CNN")] <- list(
  fpaths$LST[2:3], fpaths$LST[4:5], fpaths$LST[6:7]
) |> lapply(
  \(L) lapply(L, fread, header = F, col.names = c("ID", "QC")) |>
  rbindlist()|>
  (
    \(DT) DT[
      ,
      let(
        PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
        SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+"))
      )
    ][
      , ID := NULL
    ][
      "Fail", on = "QC"
    ]
  )()
) |> invisible()

## Volumes
data.lst$VOL[c("MALF", "NLPB", "CNN")] <- fpaths$VOL |>
  lapply(fread) |>
  lapply(
    \(DT) DT[
      ,
      let(
        PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
        SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+"))
      )
    ]
  )

## Dropped subjects by ADNI -> 70
data.lst$DROPPED <- lapply(
  data.lst$VOL,
  \(DT) DT[!adnimerge_subs, on = "PTID", PTID]
)
rm(fpaths, adnimerge_subs)


### Parse QC
data.lst$VOL <- Map(
  \(DT, QC, dropped_subs, name) {
    DT[
      !dropped_subs, on = "PTID" # This leaves 1641 subjects
    ][
      !data.lst$QC$ACQ, on = .(PTID, SCANDATE) # This varies by segmentation
    ][
      , let(ID = NULL, METHOD = tolower(name))
    ][
      QC, on = .(PTID, SCANDATE), QC := "Fail"
    ][
      is.na(QC), QC := "Pass"
    ] |> setcolorder(c("PTID", "SCANDATE", "METHOD", "QC"))
  },
  data.lst$VOL,
  data.lst$QC[-1],
  data.lst$DROPPED,
  names(data.lst$VOL)
) |> rbindlist()


### OUTPUT
outdir <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
outdir |>
  sprintf(
    fmt = "%s/%s.rds",
    c("ptid_qc_discarded", "adni-bl_volumes_hchvc")
  ) |>
  here() |>
  Map(
    f = \(Outfile, Object) saveRDS(Object, Outfile),
    list(data.lst$QC$ACQ, data.lst$VOL)
  ) |>
  invisible()
rm(outdir)
