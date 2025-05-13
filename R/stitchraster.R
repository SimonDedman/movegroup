#' Stitch together movegroup data individuals core and general use areas
#'
#' If over-large datasets cause RAM crashes for movegroup, one can run batches
#' of individuals in movegroup then join the individual saved area.ct csv files.
#'
#' See www.GitHub.com/SimonDedman/movegroup for issues, feedback, and
#' development suggestions. Install 'move' development version with:
#' remotes::install_git('https://gitlab.com/bartk/move.git')
#'
#' @import utils
#' @importFrom dplyr mutate rename select bind_rows pull
#' @importFrom rlang .data
#' @importFrom tidyselect all_of
#'
#' @export stitchraster
#'
#' @param data Data frame object containing the data. Requires columns Lat Lon DateTime ID and
#' potentially a grouping column (not currently implemented, email to request). Column names
#' specified in later parameters.
#' @param ID Name of animal tag ID column in data. "Character".
#' @param absVolumeAreaSaveName File name plus extension where UD estimates are saved. Default
#' "VolumeArea_AbsoluteScale.csv".
#' @param savedir Save outputs to a temporary directory (default) else change
#' to desired directory e.g. "/home/me/folder". Do not use getwd() for this.
#' Do NOT include terminal slash. Directory must exist. Default tempdir().
#'
#' @return Calculated volume area estimates for 50 and 95pct contours csv.
#' @details
#' Parameters values should match those used in movegroup.
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
stitchraster <- function(
  data = NULL, # Data frame object containing the data. Requires columns Lat Lon DateTime ID and optionally a grouping column.
  ID = NULL, # Name of animal tag ID column in data.
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = tempdir() # save outputs to a temporary directory (default) else change to current
  # directory e.g. "/home/me/folder". Do not use getwd() here.
) {
  data <- dplyr::rename(
    # Rename user entry to "ID", Ditto Datetime Lat & Lon
    .data = data,
    ID = tidyselect::all_of(ID)
  ) |>
    dplyr::mutate(ID = make.names(ID)) |>
    # remove all extraneous columns, massively reducing computational need
    dplyr::select(tidyselect::all_of(c("ID", "Datetime", "Lat", "Lon")))

  ID <- unique(data$ID)
  rm(data)
  bb.list <- list()
  counter <- 0

  for (i in ID) {
    counter <- counter + 1
    # read in area cts
    area.ct <- read.csv(file = file.path(savedir, paste0("areact_", i, ".csv")))
    # Put in list
    bb.list[[counter]] <- area.ct
  } # close for i in ID

  # Put everything in a data.frame
  md <- dplyr::bind_rows(bb.list, .id = "column_label") |>
    dplyr::select(!"column_label") # remove column_label column
  # 2023-08-30 quoted column_label to hopefully address gbm.factorplot: no visible binding for global variable ‘column_label’

  # read resterres from file
  rasterres <- read.csv(file = file.path(savedir, "Resolutions.csv")) |>
    dplyr::pull(rasterres)
  md$core.use <- (rasterres * md$core.use.new) / 1000000 # convert from cells/pixels to metres squared area based on cell size, then to kilometres squared area
  md$general.use <- (rasterres * md$general.use.new) / 1000000

  write.csv(
    md,
    file = file.path(savedir, absVolumeAreaSaveName),
    row.names = FALSE
  )

  # 2023-10-04 Vital memory bug warning
  if (all(md$core.use == md$core.use[1]))
    message(
      "All core UDs identical. Maybe insufficient memory for raster calcs - check rasterResolution"
    )
  if (all(md$general.use == md$general.use[1]))
    message(
      "All general UDs identical. Maybe insufficient memory for raster calcs - check rasterResolution"
    )
  if (length(which(md$core.use == max(md$core.use, na.rm = TRUE))) > 1)
    message(
      "More than 1 individual share exactly the same max value for core use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
  if (length(which(md$general.use == max(md$general.use, na.rm = TRUE))) > 1)
    message(
      "More than 1 individual share exactly the same max value for general use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
  if (length(which((md$core.use - md$general.use) == 0)) > 0)
    message(
      "1 or more individuals have exactly the same value for core and general use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
} # close function
