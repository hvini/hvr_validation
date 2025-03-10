#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)


### INPUTS
fpaths <- list(
  LST = here("lists/adni_baseline.lst"),
  CSV = c("MRILIST", "ADNIMERGE") |> sprintf(fmt = "data/%s.csv") |> here()
)

## Missing CSV files
dependencies <- unlist(fpaths, use.names = F)
missing_files <- dependencies[!file.exists(dependencies)]
if (length(missing_files) > 0) {
  missing_files |>
    paste(collapse = ", ") |>
    sprintf(fmt = "Required file(s) could not be found: %s.") |>
    stop()
}
rm(missing_files, dependencies)


### Parse DATA
## ADNIMERGE columns of interest
cols <- c(
  "PTID", "RID",
  "DX_bl", "DX",
  "AGE", "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "APOE4",
  "ADAS13", "CDRSB", "MMSE",
  "RAVLT_immediate", "RAVLT_learning", "RAVLT_forgetting", "RAVLT_perc_forgetting",
  "EXAMDATE", "Month",
  "ABETA", "PIB", "AV45", "FBB",
  "TAU", "PTAU",
  "Ventricles", "Hippocampus",
  "ICV", "FSVERSION"
)

data.lst <- list(
  BASELINE = fread(fpaths$LST, col.names = c("PTID", "SCANDATE")),
  MRILIST = fpaths$CSV[1] |>
    fread(select = c("SUBJECT", "VISIT", "SCANDATE")) |>
    setnames("SUBJECT", "PTID") |>
    unique(),
  ADNI_DATES = fread(fpaths$CSV[2], select = c("PTID", "VISCODE", "EXAMDATE")),
  ADNIMERGE = fread(fpaths$CSV[2], select = cols)
)
rm(fpaths, cols)

## Parse DATE
lapply(
  data.lst[c("BASELINE", "MRILIST")],
  \(DT) DT[, SCANDATE := ymd(SCANDATE)]
) |> invisible()

## "Baseline" segmented images are not necessarily "Baseline" in ADNIMERGE
## Need to identify visit with SCANDATE of MRILIST
data.lst$BL_MRI <- data.lst$BASELINE |>
  merge(data.lst$MRILIST, by = c("PTID", "SCANDATE")) |>
  merge(data.lst$ADNI_DATES, by = "PTID") |>
  na.omit() |>
  ## Pick the closest ADNIMERGE EXAMDATE for each PTID's SCANDATE
  (\(DT) DT[, .SD[which.min(abs(ymd(EXAMDATE) - SCANDATE))], PTID])()

## Parse ADNIMERGE
data.lst$ADNIMERGE <- merge(
  data.lst$ADNIMERGE,
  data.lst$BL_MRI[, .(PTID, EXAMDATE, SCANDATE)],
  by = c("PTID", "EXAMDATE")
) |>
  (\(DT) {
    DT["", on = "DX", DX := DX_bl]
    DT[DX_bl == "CN", DX := "CN"]
    DT[DX_bl %like% "MCI", DX := "MCI"]
    DT[DX == "AD", DX := "Dementia"]
    DT[
      ,
      DX := factor(
        DX, levels = c("CN", "MCI", "Dementia"), labels = c("CH", "MCI", "AD")
      )
    ]
  })()

### OUTPUT
outdir <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
outdir |>
  # TODO: Check this works
  here("adnimerge_baseline.rds") |>
  saveRDS(object = data.lst$ADNIMERGE)
rm(outdir)
