#' Automates dynamic Brownian bridge movement model construction across individuals
#'
#' Automates dynamic Brownian bridge movement model calculation for utilization distribution (UD)
#' estimation for multiple individuals simultaneously, using functions in the 'move' package. The
#' authors are indebted to the move package authors Bart Kraunstauber, Marco Smolla, and Anne K
#' Scharf, and to Sarah Becker for seed code which inspired the development of the
#' movegroup::movegroup function.
#'
#' Step 1. Filter individuals.
#' Remove those individuals for which there are insufficient data i.e. number of re-locations is
#' smaller than the window size parameter value (default = 31).
#'
#' Step 2. Generate universal raster.
#' Based on all remaining data, a universal raster is generated where the calculated UDs are plotted
#'  into.
#'
#' Step 3. Loop through individuals.
#' Individuals are looped through to construct individual-level movement models (on an absolute
#' scale).
#'
#' See www.GitHub.com/SimonDedman/movegroup for issues, feedback, and development
#' suggestions.
#'
#' install_git('https://gitlab.com/bartk/move.git') #Installs 'move' development version
#' @import utils
#' @importFrom lubridate is.POSIXt
#' @importFrom sp CRS SpatialPoints spTransform proj4string
#' @importFrom beepr beep
#' @importFrom dplyr mutate rename group_by summarise filter semi_join distinct ungroup arrange across bind_cols pull
#' @importFrom tidyr drop_na
#' @importFrom raster raster projectExtent res ncell setValues calc values writeRaster
#' @importFrom move move timeLag burst brownian.bridge.dyn brownian.motion.variance.dyn getVolumeUD
#' @importFrom rlang .data
#' @importFrom grDevices graphics.off
#' @importFrom graphics par
#' @importFrom methods new
#' @importFrom stats setNames
#'
#' @export movegroup
#'
#' @param data Data frame object containing the data. Requires columns Lat Lon DateTime ID and
#' potentially a grouping column (not currently implemented, email to request). Column names
#' specified in later parameters.
#' @param ID Name of animal tag ID column in data. "Character".
#' @param Datetime Column name in data that contains date/time stamps for each recorded detection.
#' Must be in POSIXct format. "Character".
#' @param Lat Name of latitude column in data. "Character".
#' @param Lon Name of longitude column in data. "Character".
# #' @param Group Name of grouping column in data. Unused 2023-08-14.
#' @param dat.TZ Timezone of data for as.POSIXct. Default "US/Eastern".
#' @param proj CRS for move function. Default sp::CRS("+proj=longlat +datum=WGS84").
#' @param sensor Sensor for move function. Single character or vector with length of the number of
#' coordinates. Optional. Default "VR2W".
#' @param moveLocError Location error (m) in the 'brownian.bridge.dyn' function in the 'move'
#' package. Numeric. Either single or a vector of length nrow data. If using passive acoustic data
#' this is the detection range of the receiver(s). Default 1. See MoveLocErrorCalc function for
#' satellite data with state space modelled locations with 95% confidence intervals for latlon i.e.
#' lat and lon025 and 975.
#' @param timeDiffLong Single numeric value. Threshold value in timeDiffUnits designating the length
#'  of long breaks in re-locations. Used for bursting a movement track into segments, thereby
#'  removing long breaks from the movement track. See ?move::bursted for details. Default 2 hours is
#'   arbitrary, looping through 18, 24, and 36 hours for satellite data on great hammerhead sharks
#'  revealed volume areas for core and general use gradually rise with timeDiffLong increases,
#'  multiple small dots of presence get blobbed together, and therefore sometimes this covers land.
#'  Ideally one would not discard any data, in which case one should choose a value higher than the
#'  largest between-detections gap in their dataset (or just pick a very large number). This
#'  parameter is useful when the model would otherwise get stuck trying to calculate a UD for an
#'  individual with a very large home range that is inadequately captured by a receiver array.
#'  Default 2.
#' @param timeDiffUnits Character. Unit for timeDiffLong. Default "hours".
#' @param center US English alternate to centre. Do you want to center the move object within
#' extent? See spTransform. Default TRUE.
#' @param centre British English alternate to center. Do you want to centre the move object within
#' extent? See spTransform. Default NULL.
#' @param buffpct Buffer extent for raster creation, proportion of 1. Default 0.3, can try e.g. 3
#' for a large buffer to avoid clipping, at the cost of file size, but later cropping in
#' plotraster.R will remove extraneous blank space.
#' @param rasterExtent Extent of raster created around data. If NULL (default), calculated from
#' data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon
#' degrees. Don't go to 90 (degrees) north or south for ymax or ymin. Doesn't prevent constraint to
#' data limits (in plot anyway), but prevents raster clipping crash.
#' @param rasterCRS CRS for raster creation. Default sp::CRS("+proj=utm +zone=17 +datum=WGS84").
#' @param rasterResolution Single numeric value to set raster resolution - cell size (width and
#' height) in metres. 111000: 1 degree lat = 111km. Trade-off between small res = big file & processing
#' time. Should be a function of the spatial resolution of your receivers or positioning tags.
#' Higher resolution will lead to more precision in the volume areas calculations. Try using
#' 2*dbblocationerror, if dbblocationerror is a single value. Default 50 = 50m = 50m² = 0.00005 km²
#' (divide by 1000000) = 0.00045 degrees. Try around the median of your moveLocError.
# #' @param dbblocationerror Location.error param in 'brownian.bridge.dyn' function in the 'move'
# #' package. single numeric value or vector of the length of coordinates that describes the error of
# #' the location (sender/receiver) system in map units. Or a character string with the name of the
# #' column containing the location error can be provided. Default is moveLocError. See
# #' MoveLocErrorCalc function for satellite data with state space modelled locations with 95%
# #' confidence intervals for latlon i.e. lat and lon025 and 975.
#' @param movemargin Margin size for variance calc in move::brownian.motion.variance.dyn and
#' behavioural change point analysis in move::brownian.bridge.dyn. Must be an odd number. Default 11.
#' Motion variance based on only the middle section of the trajectory; the ends
#' of the movement trajectory where no changes are allowed because at some stage
#' you want to have a few locations to base your estimation of the variance on
#' and how many locations in either side of the window we use for this, is
#' called the margin. Smaller values for window size and margin is expected to
#' give a higher frequency of behavioural changes; make these large for looking
#' at migrations.
#' @param dbbext Ext param in the 'brownian.bridge.dyn' function in the 'move' package. Extends
#' bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax)
#' . Default 0.3 - extends bounding box by 30 percent. Excessive buffering will get cropped automatically.
#' @param dbbwindowsize The window.size param in the 'brownian.bridge.dyn' function
#' in the 'move' package. The size of the moving window along the track. Larger windows provide more
#' stable/accurate estimates of the brownian motion variance but are less well able to capture more
#' frequent changes in behaviour. Number must be odd. Individuals with fewer
#' detections than the window size (default 23) will be removed. A window size
#' of 23 means the sliding window is moved every 23 locations or every 23 hours
#'  (via timeDiffUnits). Must be >= 2*movemargin (default 11).
#' @param writeRasterFormat Character. Output file type for raster::writeRaster param format.
#' Default "ascii". TO DEPRECIATE.
#' @param writeRasterExtension Character. Output file extension for raster::writeRaster param
#' extension. Default ".asc". TO DEPRECIATE.
#' @param writeRasterDatatype Character. Data type for writing values to disk for
#' raster::writeRaster param Datatype. Default "FLT4S". TO DEPRECIATE.
#' @param absVolumeAreaSaveName File name plus extension where UD estimates are saved. Default
#' "VolumeArea_AbsoluteScale.csv".
#' @param savedir Save outputs to a temporary directory (default) else change
#' to desired directory e.g. "/home/me/folder". Do not use getwd() for this.
#' Do NOT include terminal slash. Directory must exist. Default tempdir().
#' @param saveAreaCT Save tiny individual core and general use areas tables to
#' disk. These are the only things retained in the per-individual loop, so if
#' your large dataset causes memory crashes, you can run it in chunks and stitch
#'  the results together later with stitchraster. Default FALSE.
#' @param alerts Audio warning for failures. Default TRUE.
#'
#' @return Individual-level utilization distributions, saved as rasters, as well
#'  as calculated volume area estimates for 50 and 95pct contours, saved in a
#'  .csv file. Motion variance csvs per individual ("MotionVariance.csv"), see
#'  move::brownian.motion.variance.dyn. No processed object is returned, i.e. bad: "objectname <-
#'  movegroup()", good: "movegroup()".
#' @details
#' When used together, the order of functions would be: movegroup, (stitchraster
#'  if required), scaleraster, (alignraster if required then stitchraster again)
#' , plotraster.
#'
#' ## Errors and their origins:
#'
#' 1. Error in .local(object, raster, location.error = location.error, ext = ext: Higher y or x grid not
#' large enough, consider extending the raster in that direction or enlarging the ext argument.
#' Increase buffpct, e.g. to 3.
#'
#' 2. Error in .data\[\[dttm\]\]: Must subset the data pronoun with a string, not a <POSIXct/POSIXt>
#' object. Use "ColName" not dataframe$ColName syntax for Datetime, ID, Lat, Lon.
#'
#' 3. Error in splice(dot_call(capture_dots, frame_env = frame_env, named = named,: object
#' 'DateTime' not found. Use "ColName" not ColName syntax for Datetime, ID, Lat, Lon.
#'
#' 4. Error in .local(object, raster, location.error = location.error, ext = ext: Higher x grid not
#' large enough, consider extending the raster in that direction or enlarging the ext argument. Try
#' "buffpct = 1," , then larger e.g. 3, if still getting the error.
#'
#' 5. cannot allocate vector of size (BIG) Gb: Increase rasterResolution value.
#'
#' 6. In min/max: No non-missing arguments to min; returning Inf: likely not enough memory, increase
#'  rasterResolution value.
#'
#' 7. Error in tmp\[\[i\]\]: subscript out of bounds. dbbmmwindowsize may be too large relative to nrow
#' of that individual. Try lowering movemargin (default 11, has to be odd) and then lowering
#' dbbmmwindowsize (default 23, has to be >=2*movemargin, has to be odd).
#'
#' 8. Error in validityMethod(as(object, superClass)): The used raster is not a
#' UD (sum unequal to 1), sum is NaN. Potentially from memory overrun from large
#' datasets. Close all other programs, restart your session, remove other
#' objects from R, and try again, watching the RAM usage piechart by the
#' Environment tab. If unsuccessful, run function in chunks of individuals until
#' all (of length >= windowsize) are saved to disc as asc files, then produce
#' final summary stats with stitchraster.
#'
#' @examples
#' \donttest{
#' # load data
#' data("TracksCleaned")
#' # run function
#' movegroup(
#'  data = TracksCleaned,
#'  ID = "Shark",
#'  Datetime = "Datetime",
#'  Lat = "Lat",
#'  Lon = "Lon",
#'  savedir = tempdir())
#'  }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}
#'
#' @references Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M. and Safi, K. (2012), A
#' dynamic Brownian bridge movement model to estimate utilization distributions for heterogeneous
#' animal movement. Journal of Animal Ecology. doi: 10.1111/j.1365-2656.2012.01955.x
#'
#' Kranstauber, B., M. Smolla & A. K. Scharf. 2019. Move: visualizing and analyzing animal track
#' data. R package version 4.2.4 (at 2023-08-15). https://CRAN.R-project.org/package=move.

movegroup <- function(
  data = NULL, # Data frame object containing the data. Requires columns Lat Lon DateTime ID and optionally a grouping column.
  ID = NULL, # Name of animal tag ID column in data.
  Datetime = NULL, # Column name in data that contains date/time stamps for each recorded detection. Must be in POSIXct format.
  Lat = NULL, # name of Lat & Lon columns in data.
  Lon = NULL,
  # Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  dat.TZ = "US/Eastern", # timezone for as.POSIXct.
  proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
  sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
  moveLocError = 1, # location error in metres for move function. Numeric. Either single or a vector of length nrow data.
  timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  timeDiffUnits = "hours", # units for time difference for move function.
  center = TRUE, # center move object within extent? See spTransform.
  centre = NULL, # British English alternate to center.
  buffpct = 0.3, # buffer extent for raster creation, proportion of 1.
  rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
  # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
  rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. This is around Bimini, Bahamas.
  rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres. 111000: 1 degree lat = 111km
  # dbblocationerror = moveLocError, # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
  movemargin = 11, # Margin size for variance calc in move::brownian.motion.variance.dyn and behavioral change point analysis in
  # move::brownian.bridge.dyn. Must be an odd number. Default 11.
  dbbext = 3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 3.
  dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more
  # stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior.
  # This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 23). Must be >=2*margin (movemargin) which is 11 so >=22,
  # but odd so >=23.
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = tempdir(), # save outputs to a temporary directory (default) else change to current
  # directory e.g. "/home/me/folder". Do not use getwd() here.
  saveAreaCT = FALSE, # save tiny individual core and general use areas tables
  # to disk. These are the only things retained in the per-individual loop, so
  # if your large dataset causes memory crashes, you can run it in chunks and
  # stitch the results together later with stitchraster.
  alerts = TRUE # audio warning for failures
) {
  # If savedir has a terminal slash, remove it, it's added later
  if (substr(x = savedir, start = nchar(savedir), stop = nchar(savedir)) == "/")
    savedir = substr(x = savedir, start = 1, stop = nchar(savedir) - 1)

  # if savedir doesn't exist, can't save or setwd into it, therefore create it
  if (!file.exists(savedir)) dir.create(savedir)

  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts)
    options(
      error = function() {
        beepr::beep(9) # give warning noise if it fails
        graphics.off() # kill all graphics devices
        setwd(oldwd) # reinstate original working directory. Probably redundant given on.exit
      } # close options subcurly
    ) # close options

  # edit center if centre used
  if (!is.null(centre)) center <- centre

  # Create writeRasterExtension from writeRasterFormat
  if (is.null(writeRasterExtension))
    writeRasterExtension <- switch(
      EXPR = writeRasterFormat,
      "ascii" = ".asc",
      "raster" = ".grd",
      "SAGA" = ".sdat",
      "IDRISI" = ".rst",
      "CDF" = ".nc",
      "GTiff" = ".tif",
      "ENVI" = ".envi",
      "EHdr" = ".bil",
      "HFA" = ".img"
    )

  data <- dplyr::rename(
    # Rename user entry to "ID", Ditto Datetime Lat & Lon
    .data = data,
    Datetime = all_of(Datetime), # updated for tidyselect1.2.0
    ID = all_of(ID),
    Lat = all_of(Lat),
    Lon = all_of(Lon)
  ) |>
    dplyr::mutate(ID = make.names(ID)) |>
    # remove all extraneous columns, massively reducing computational need
    dplyr::select(tidyselect::all_of(c("ID", "Datetime", "Lat", "Lon")))

  if (!any(class(data$Datetime) == "POSIXct")) {
    stop("Datetime column must be POSIXct")
  }

  # Add movelocationerror and dbblocationerror as columns if they're vectors, else their length
  # desyncs from nrow(data) if any IDs < windowsize
  if (nrow(data) == length(moveLocError)) data$moveLocError <- moveLocError
  # if (nrow(data) == length(dbblocationerror)) data$dbblocationerror <- dbblocationerror

  # Convert coordinate sets to SpatialPoints and project
  cord.dec <- sp::SpatialPoints(
    cbind(data$Lon, data$Lat),
    proj4string = sp::CRS("+proj=longlat")
  )

  # Transform to UTM by setting the EPSG to local UTM zone for the data.
  cord.UTM <- as.data.frame(sp::spTransform(cord.dec, rasterCRS))
  colnames(cord.UTM) <- c("NewEastingUTM", "NewNorthingUTM")
  data <- cbind(data, cord.UTM) # 1308 x 7
  rm(list = c("cord.UTM", "cord.dec"))

  # A dBBMM is not run if total detections of individual < window size (default value, 31).
  # Below code checks and filters out individuals with insufficient data.
  enoughrelocs <- data |>
    group_by(ID) |>
    summarise(n = n()) |>
    filter(n >= dbbwindowsize) |>
    pull(ID)
  data <- data |> filter(ID %in% enoughrelocs)
  rm(enoughrelocs)

  # Create per-individual move object, project it, construct a dBBMM, calculate the volume area
  # within the 50% and 95% contours, save as ASCII.

  bb.list <- list()

  data <- tidyr::drop_na(data = data, ID) |>
    dplyr::group_by(ID) |>
    dplyr::distinct(Datetime, .keep_all = TRUE) |> # distinct (grouped): removed one row (<1%), 1,286 rows remaining
    # prevents duplicate Datetime crash in move() later
    dplyr::ungroup() # 1253 x 7 after removing as.numeric above

  # Make a raster for the UD to plot into. Start with UTM.
  # These coordinates need to be big enough to cover your data.
  # May need to expand x and y ranges (i.e. the "e" variable below) if you encounter errors when constructing the DBBMM in the next code chunk.
  # NOTE: "xUTM" grid should be larger than the "x" i.e. lower minima and higher maxima.
  # The variable 'e' influences the extent [now done with buffpct]
  # The resolution determines the grid size.
  # Finer resolution/grid size (in m) i.e. lower value means longer computation time. May need/want to play with this value too.
  # e <- 80 * 1000 # make e a function of range of extent as %

  # Error in .local(object, raster, location.error = location.error, ext = ext,  :
  # Lower x grid not large enough, consider extending the raster in that direction or enlarging the ext argument
  # make buffpct larger
  if (!is.null(rasterExtent)) {
    # if xUTM object exists
    xUTM <- sp::SpatialPoints(
      cbind(
        c(rasterExtent[1], rasterExtent[2]),
        c(rasterExtent[3], rasterExtent[4])
      ),
      proj4string = sp::CRS("+proj=longlat")
    )
    xUTM <- as.data.frame(sp::spTransform(xUTM, rasterCRS)) # Transform to UTM
    xUTM <- raster::raster(
      # create a raster
      xmn = xUTM[1, 1] - ((xUTM[2, 1] - xUTM[1, 1]) * buffpct), # with xUTM's values as extents
      xmx = xUTM[2, 1] + ((xUTM[2, 1] - xUTM[1, 1]) * buffpct),
      ymn = xUTM[1, 2] - ((xUTM[2, 2] - xUTM[1, 2]) * buffpct),
      ymx = xUTM[2, 2] + ((xUTM[2, 2] - xUTM[1, 2]) * buffpct),
      crs = rasterCRS, # +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs
      resolution = rasterResolution # 50
    ) # close raster
  } else {
    # exists(xUTM), else create from data
    xUTM <- raster::raster(
      xmn = min(data$NewEastingUTM, na.rm = TRUE) -
        ((max(data$NewEastingUTM, na.rm = TRUE) -
          min(data$NewEastingUTM, na.rm = TRUE)) *
          buffpct),
      xmx = max(data$NewEastingUTM, na.rm = TRUE) +
        ((max(data$NewEastingUTM, na.rm = TRUE) -
          min(data$NewEastingUTM, na.rm = TRUE)) *
          buffpct),
      ymn = min(data$NewNorthingUTM, na.rm = TRUE) -
        ((max(data$NewNorthingUTM, na.rm = TRUE) -
          min(data$NewNorthingUTM, na.rm = TRUE)) *
          buffpct),
      ymx = max(data$NewNorthingUTM, na.rm = TRUE) +
        ((max(data$NewNorthingUTM, na.rm = TRUE) -
          min(data$NewNorthingUTM, na.rm = TRUE)) *
          buffpct),
      crs = rasterCRS, # +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs
      resolution = rasterResolution # 50
    ) # close raster
  } # close else

  # run move on all data together to generate global CRS for raster
  alldata <- data |>
    dplyr::arrange(Datetime) |>
    dplyr::group_by(Datetime) |> # remove duplicates
    dplyr::summarise(
      dplyr::across(
        tidyselect::where(~ is.numeric(.)),
        \(x) mean(x, na.rm = TRUE)
      ), # mean, na.rm = TRUE
      dplyr::across(
        tidyselect::where(~ is.character(.) | lubridate::is.POSIXt(.)),
        dplyr::first
      )
    ) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::where(~ is.numeric(.)),
        ~ ifelse(is.nan(.), NA, .)
      ), #convert NaN to NA. POSIX needs lubridate
      dplyr::across(
        tidyselect::where(~ is.character(.)),
        ~ ifelse(is.nan(.), NA, .)
      )
    ) |> # https://community.rstudio.com/t/why-does-tidyrs-fill-work-with-nas-but-not-nans/25506/5
    dplyr::ungroup()

  alldata <- as.data.frame(alldata)
  moveall <- move::move(
    x = alldata$Lon, # numeric vector with x coordinates if non-Movebank data are provided (e.g. data$x).
    y = alldata$Lat, # numeric vector with y coordinates if non-Movebank data are provided (e.g. data$y).
    time = alldata$Datetime, # vector of time stamps with POSIXct conversion if non-Movebank data are provided,
    # i.e. as.POSIXct(data$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    data = alldata, # extra data associated with the relocations
    proj = proj, # projection method for non-Movebank data; requires a valid CRS (see CRS-class) object, e.g. CRS("+proj=longlat +ellps=WGS84")
    # animal = data$ID,
    sensor = sensor # Sensor name(s), either single character or a vector with length of the number of coordinates.
  )

  # Convert projection to Azimuthal Equi-Distance projection (aeqd)
  rall <- sp::spTransform(moveall, center = center)
  # sp::proj4string(rall) # "+proj=aeqd +lat_0=44.99489997 +lon_0=-17.48575004 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  rallCRS <- sp::CRS(sp::proj4string(rall))
  # class(rallCRS) # CRS
  write.csv(
    sp::proj4string(rall),
    file.path(savedir, "CRS.csv"),
    row.names = FALSE
  )
  saveRDS(rallCRS, file = file.path(savedir, "CRS.Rds"))
  rm(moveall)
  rm(alldata)

  # Reproject to AEQD. Make a dummy object to get the correct projection.
  x.i <- raster::raster(
    xmn = min(data$NewEastingUTM),
    xmx = max(data$NewEastingUTM),
    ymn = min(data$NewNorthingUTM),
    ymx = max(data$NewNorthingUTM),
    crs = sp::proj4string(rall) # if rasterCRS then same as above! trying proj4string(r.i) which is AEQD from spTransform
    # could add here: resolution = rasterResolution # 50 ####
  )

  # Reproject UTM to AEQD
  # ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), i.e. projectExtent(from, to, ...))
  newTemplate <- raster::projectExtent(xUTM, sp::proj4string(x.i))

  # change newTemplate to xAEQD. Change xAEQD res so x & y match. Kills values.
  raster::res(newTemplate) <- rep(rasterResolution, 2)

  # Give newTemplate some values. Make Rep equal to the ncell dimension
  ones <- rep(1, raster::ncell(newTemplate))
  xAEQD <- raster::setValues(newTemplate, ones)
  rm(newTemplate)
  rm(ones)

  # get resolution from raster in rasterlist, assign it object, squared
  rasterres <- (raster::res(xAEQD)[1])^2
  write.csv(
    x = data.frame(rasterres = rasterres, rasterResolution = rasterResolution),
    file = file.path(savedir, "Resolutions.csv")
  )

  # Loop through all unique tags
  counter <- 0
  for (i in unique(data$ID)) {
    # i <- unique(data$ID)[1]

    # Print individual
    counter <- counter + 1
    message(
      "processing ",
      which(unique(data$ID) %in% i),
      " of ",
      length(unique(data$ID))
    )

    # Filter individual data
    data.i <- as.data.frame(data[data$ID == i, ])

    # create move object
    move.i <- move::move(
      x = data.i$Lon,
      y = data.i$Lat,
      time = data.i$Datetime,
      proj = proj,
      data = data.i,
      animal = data.i$ID,
      sensor = sensor
    )

    # Check the current projection
    # sp::proj4string(move.i) # "+proj=longlat +datum=WGS84 +no_defs"

    # Incorporate uncertainty in the model by including a location error.
    if ("moveLocError" %in% colnames(data.i)) {
      move.i$LocationError <- data.i$moveLocError
    } else {
      if (exists("moveLocError")) {
        if (length(moveLocError) == 1) {
          move.i$LocationError <- moveLocError # single value, replicated down for each relocation
        } else {
          # else if multiple values
          if (length(moveLocError) == nrow(data.i)) {
            # should be same length as full dataset
            move.i$LocationError <- data.i |> # take the full dataset,
              dplyr::bind_cols(moveLocError = moveLocError) |> # cbind the full movelocerror
              dplyr::filter(ID == i) |> # filter for just this ID
              dplyr::pull(moveLocError) # and pull just the movelocerror for this ID
          } else {
            # if not length 1 and not length of nrow(data)
            stop(message(
              "moveLocError must be either length 1 or length(nrow(data))"
            )) # if not stop and tell user
          } # close not length1 not same length as full dataset else
        } # close length1 or else
      } # close if exists movelocerror
    }

    # Convert projection to Azimuthal Equi-Distance projection (aeqd)
    r.i <- sp::spTransform(move.i, CRSobj = rallCRS, center = center)
    rm(move.i) # not needed after this in loop

    # Calculate time lag between consecutive detections.
    # Obviously there is no time lag for the first detection, so we need to tell R this.
    # Knowing time lags in acoustic telemetry is important, as the motion variance is based on the
    # time difference b/w detections at consecutive locations i.e. larger time lag creates a wider
    # bridge where the animal could've been and therefore potentially inflates the motion variance.
    # Below code deals with this issue by identifying the location and number of detections with
    # large time gaps. Below we will use an arbitrary value of 2 h.
    TimeDiff <- move::timeLag(r.i, units = timeDiffUnits)
    r.i$TimeDiff <- append(0, TimeDiff)
    long <- which(r.i$TimeDiff > timeDiffLong)
    # if no timediffs are longer than timedifflong, long is integer(0), a potentially dangerous object.

    # Below we exclude the data points that are 'long' i.e. create a large time gap. move::burst ?
    bursted <- move::burst(
      r.i,
      c("normal", "long")[
        1 + (move::timeLag(r.i, units = timeDiffUnits) > timeDiffLong)
      ]
    )
    rm(r.i)
    # There are 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command.

    if ("moveLocError" %in% colnames(data.i)) {
      moveLocError.i <- data.i$moveLocError
    } else {
      if (exists("moveLocError")) {
        if (length(moveLocError) == 1) {
          moveLocError.i <- moveLocError # single value, replicated down for each relocation
          # Robs: 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate
        } else {
          # else if multiple values
          if (length(moveLocError) == nrow(data)) {
            # should be same length as full dataset
            moveLocError.i <- data |> # take the full dataset,
              dplyr::bind_cols(moveLocError = moveLocError) |> # cbind the full movelocerror
              dplyr::filter(ID == i) |> # filter for just this ID
              dplyr::pull(moveLocError) # and pull just the movelocerror for this ID
          } else {
            # if not length 1 and not length of nrow(data)
            stop(message(
              "moveLocError must be either length 1 or length(nrow(data))"
            )) # if not stop and tell user
          } # close not length1 not same length as full dataset else
        } # close length1 or else
      } # close if exists moveLocError
    } # close ifelse

    rm(data.i) #not needed after this in loop

    # location error needs to be a positive number. Replace zeroes with 0.00001
    moveLocError.i[which(moveLocError.i == 0)] <- 0.00001

    # Construct the model. The time.step should reflect the ping frequency of the tag (in minutes)
    # This is what takes CPU & memory
    bursted_dbbmm <- move::brownian.bridge.dyn(
      bursted,
      burstType = "normal",
      raster = xAEQD, # has to be AEQD for metre-based calculations, presuambly?
      # Error:The projection of the raster and the Move object are not equal.
      # Need bursted, r.i, to be the same projection as xAEQD
      location.error = moveLocError.i,
      ext = dbbext, # dbbext
      margin = movemargin,
      window.size = dbbwindowsize #  must be >=2*margin which is 11 so >=22, but odd so >=23
    )

    # Build dBMvarianceBurst object with brownian.motion.variance.dyn
    bursted_motionvar <- move::brownian.motion.variance.dyn(
      bursted,
      location.error = moveLocError.i,
      window.size = dbbwindowsize, #  must be >=2*margin which is 11 so >=22, but odd so >=23
      margin = movemargin
    )
    mtnvar <- as.data.frame(bursted_motionvar) # extract the data
    mtnvar$motionVariance <- move::getMotionVariance(bursted_motionvar) # add motionvariance
    write.csv(
      mtnvar,
      file = file.path(savedir, paste0(i, "_MotionVariance.csv")),
      row.names = FALSE
    )
    rm(bursted)
    rm(mtnvar)

    # Re-standardize (Dr. Kranstauber's trouble shooting solution).
    # Occurring errors are "due to limits of accuracy during calculations.
    # Given the error is really small (<.000001 %) I would not worry to much and re-standardize".
    # Then aggregate UD segments.
    tmp <- raster::calc(bursted_dbbmm, sum) # raster::calc :Calculate values for a new Raster* object from another Raster* object
    # using a formula. returns a RasterLayer if fun returns a single value (e.g. sum)
    rm(bursted_dbbmm)
    bb <- new(".UD", tmp / sum(raster::values(tmp))) # new() creates object from Class ("UD.")
    rm(tmp)

    # Calculate volume area (m^2) within 50% (core) and 95% (general use) contours. Note: absolute scale
    # Note: Based on 'An introduction to the 'move' package', it's the opposite of what you might expect:
    # "A cell with a very high value in the UD raster will have a very low value in the contour raster"
    # area.50 <- round(sum(raster::values(move::getVolumeUD(bb) <= .50)), 4)
    # area.95 <- round(sum(raster::values(move::getVolumeUD(bb) <= .95)), 4)
    # area.50 <- sum(raster::values(move::getVolumeUD(bb) <= .50)) # 2024-01-08 remove rounding here. Seemingly had no effect
    # area.95 <- sum(raster::values(move::getVolumeUD(bb) <= .95)) # 2024-01-08 remove rounding here. Seemingly had no effect

    # Below calc introduced in 2022-10-08 commit, message =
    #"changed code volume area: instead of using getVolumeUD() from the move package,
    # we just did the calculations on the UD raster and not using the function.
    # these seem to produce very plausible estimates!"
    # This change creates unnaturally small areas.
    # Changed back to previous calcs 2023-02-24. Added rounding
    # area.50.new <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.5)), 4)
    # area.95.new <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.05)), 4)
    area.50.new <- sum(raster::values(bb) >= (max(bb@data@values) * 0.5)) # 2024-01-08 remove rounding here. Seemingly had no effect
    area.95.new <- sum(raster::values(bb) >= (max(bb@data@values) * 0.05)) # 2024-01-08 remove rounding here. Seemingly had no effect

    # Combine in single df
    area.ct <- data.frame(
      core.use.new = area.50.new,
      general.use.new = area.95.new
    ) # 2024-01-08 add rounding here if necessary. Maybe not.

    # Add ID id
    area.ct$ID <- i

    # Save area.ct's if requested due to memory issues
    if (exists("saveAreaCT") && saveAreaCT) {
      write.csv(
        area.ct,
        file = file.path(savedir, paste0("areact_", i, ".csv")),
        row.names = FALSE
      )
    }

    # Put in list
    bb.list[[counter]] <- area.ct

    # need to reproject bb back to UTM else the resolution is unequal in x & y
    # only an issue for ascii raster format
    names(bb) <- i # puts name in slot, means scaleraster can name its output properly with lapply

    # Export aggregated individual UDs to as ascii files for further exploration/spatial analysis in GIS software.
    # Needed for when the aim is to plot population-level and normalized UDs per species.
    raster::writeRaster(
      bb, # x has unequal horizontal and vertical resolutions. Such data cannot be stored in arc-ascii format
      file.path(savedir, paste0(i, writeRasterExtension)),
      format = writeRasterFormat, # "ascii"
      datatype = writeRasterDatatype, # "FLT4S"
      if (writeRasterFormat != "CDF") bylayer = TRUE, # bylayer kills ncdf4
      overwrite = TRUE
    )
    rm(bb)
    gc() # cleanup
  } # close for i in unique data$ID

  # MEMORY OVERRUN ISSUE ####
  # conceptually I could save area.ct
  # (and bb? bb is created as a blank list but then overwritten as a new UD then named (may do nothing) then reoverwritten in the next loop)
  # then load all area.ct objects in a loop
  # so users could break their runs into chunks

  # Put everything in a data.frame
  md <- dplyr::bind_rows(bb.list, .id = "column_label") |>
    dplyr::select(!"column_label") # remove column_label column
  # 2023-08-30 quoted column_label to hopefully address gbm.factorplot: no visible binding for global variable ‘column_label’
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
