#' Automates dynamic Brownian bridge movement model and dynamic bivariate Gaussian bridge construction across individuals
#'
#' Automates dynamic Brownian bridge movement model (dBBMM) and dynamic bivariate Gaussian bridge (dBGB) calculation for 
#' utilization distribution (UD) estimation for multiple individuals simultaneously, using functions in the 'move' package. 
#' The authors are indebted to the move package authors Bart Kraunstauber, Marco Smolla, and Anne K
#' Scharf, and to Sarah Becker for seed code which inspired the development of the
#' movegroup::movegroup function.
#'
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
#' @param model Character. Movement model to use for utilization distribution estimation. Either
#' "dbbmm" (dynamic Brownian bridge movement model, default) or "dbgb" (dynamic bivariate Gaussian
#' bridge). Please familiarize yourself with the different models using the references provided in the details before 
#' making a decision on which model you choose. Remember, in dBBMMs the animal's movement between two locations
#' is equally likely to happen in all directions, i.e. is an isotropic diffusive process. In dBGBs, the diffusion
#' is factorized in two directions, i.e. the direction towards the next location, and the direction perpendicular to it.
#' The uncertainty in your animals whereabouts between locations when the time gap between two locations is large, can 
#' therefore have different extents and patterns depending on your model choice, and depending on the temporal resolution
#' and dynamics, your data may be more suitable for one model or the other. Please also see our suggestion on how to quickly
#' and qualitatively compare the two models and their suitability for your data in the Details below.
#' @param ID Name of animal tag ID column in data. "Character".
#' @param Datetime Column name in data that contains date/time stamps for each recorded detection.
#' Must be in POSIXct format. "Character".
#' @param Lat Name of latitude column in data. "Character".
#' @param Lon Name of longitude column in data. "Character".
# #' @param Group Name of grouping column in data. Unused 2023-08-14.
#' @param dat.TZ Timezone of data for as.POSIXct. Default "US/Eastern".
#' @param proj CRS for move function. Default sp::CRS("+proj=longlat +datum=WGS84").
#' @param sensor Sensor for move function. Single character or vector with length of the number of
#' coordinates. Optional. Default "acoustic-telemetry". See https://www.movebank.org/cms/movebank-content/movebank-attribute-dictionary#event_attributes for
#' possible sensor types. This is just for personal record, and not used in any calculations. 
#' @param integrationStep corresponds to time.step and timeStep in move::brownian.bridge.dyn and move::dynBGB, respectively. 
#' I.e. corresponds to the size of the timer intervals taken for every integration step in minutes. If left NULL, this function uses the default values 
#' from the 'move' package where 15 (for dBBMM) or 20.1 (for dBGB) steps are taken in the shortest time interval 
#' For example, in dBBMM-calculations this means that if there is a location recorded e.g. every 60 minutes, the function divides each
#' segment into 4 minute chunks upon which it does the calculation. However, please be advised, that if you have even just
#' one very short time interval, each segment of the track will be divided in the nr. of chunks based on this one short time interval.
#' This can markedly increase computing time. Check the time intervals in your data first. As quote from the 'move' vignette:
#' "If the track contains time intervals much shorter than the scheduled on the tag, set the time.step e.g. to the scheduled time 
#' interval at which the tag was set (in minutes). Higher values to this argument reduce calculation time." 
#' The here provided function automatically divides your chosen duration by 15 or 20.1 depending on your model choice.
#' @param moveLocError location.error (m) in the 'brownian.bridge.dyn' function and locErr (m) in the 'dynBGB' function
#' in the 'move' package. Numeric. Either single or a vector of length nrow data. If using passive acoustic data
#' this is the detection range of the receiver(s). Default 1. See MoveLocErrorCalc function for
#' satellite data with state space modelled locations with 95% confidence intervals for latlon i.e.
#' lat and lon025 and 975.
#' @param timeDiffLong Single numeric value. Threshold value in units set via 'timeDiffUnits' designating the length
#' of long breaks in re-locations. Allows setting movement segments that correspond to long time gaps to FALSE so they
#' are ignored in subsequent model calculations avoiding model and UD inflation. Default of 2 hours is arbitrary. Looping 
#' through 18, 24, and 36 hours for satellite data on great hammerhead sharks revealed volume areas for core and 
#' general use gradually rise with timeDiffLong increases, multiple small dots of presence get blobbed together, and 
#' therefore sometimes this covers land. Ideally one would not discard any data, in which case one should choose a value 
#' higher than the largest between-detections gap in their dataset (or just pick a very large number). This
#' parameter is useful when the model would otherwise get stuck trying to calculate a UD for an
#' individual with a very large home range that is inadequately captured by a receiver array.
#' Default 2.
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
#' behavioural change point analysis in move::brownian.bridge.dyn. It is also used in move::dynBGBvariance
#' for the behavioral change point analysis and in move::dynBGB. Must be an odd number. Depending on model choice
#' there are different defaults: dbbmm - 11, dbgb - 15.
#' Motion variance based on only the middle section of the trajectory; the ends
#' of the movement trajectory where no changes are allowed because at some stage
#' you want to have a few locations to base your estimation of the variance on
#' and how many locations in either side of the window we use for this, is
#' called the margin. Smaller values for window size and margin is expected to
#' give a higher frequency of behavioural changes; make these large for looking
#' at migrations.
#' @param dbbext Ext param in the 'brownian.bridge.dyn' and 'dynBGB' functions in the 'move' package. Extends
#' bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax)
#' . Default 0.3 - extends bounding box by 30 percent. Excessive buffering will get cropped automatically.
#' @param dbbwindowsize The window.size param in the 'brownian.bridge.dyn' and windowSize param in the 'dynBGB' functions
#' in the 'move' package. The size of the moving window along the track. Larger windows provide more
#' stable/accurate estimates of the brownian motion variance but are less well able to capture more
#' frequent changes in behaviour. For example: a window size of 31 means the sliding window is moved every 31 locations 
#' or every 31 hours (via timeDiffUnits). This number must be odd. Must be >= 2*movemargin. Given the different 'movemargin' defaults for 
#' dBBMM (11) and dBGB (15), defaults for this parameter are set to 23 and 31 depending on your model choice 'dbbmm' and 'dbgb', respectively.
#' Individuals with fewer detections than the window size will be removed. 
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
#' Do NOT include terminal slash. Default tempdir().
#' @param saveAreaCT Save tiny individual core and general use areas tables to
#' disk. These are the only things retained in the per-individual loop, so if
#' your large dataset causes memory crashes, you can run it in chunks and stitch
#'  the results together later with stitchraster. Default FALSE.
#' @param alerts Audio warning for failures. Default TRUE.
#'
#' @return Individual-level utilization distributions, saved as rasters, as well
#' as calculated volume area estimates for 50 and 95pct contours, saved in a
#' .csv file. Motion variance csvs per individual ("MotionVariance.csv"), see
#' move::brownian.motion.variance.dyn or move::dynBGBvariance. No processed object is returned, i.e. bad: "objectname <-
#' movegroup()", good: "movegroup()".
#' 
#' @details We strongly recommend that prior to choosing any of the model options, you familiarize yourself with the underlying
#' concepts of dBBMMs and dBGBs by reading the corresponding publications from Kranstauber et al. (2012) for dBBMMs and
#' Kranstauber et al. (2014) for dBGBs. Please see the full reference details in the 'References' section below.
#' 
#' Here we provide a very short overview of the most important points from the listed references. Dynamic BBMMs and dBGBs are an improvement over traditional UD calculation methods, as they are able to integrate a 
#' temporal component to the movement data and explicitly model the movement trajectory between two consecutive locations. 
#' The dBBMM does this by defining the variance of the Brownian motion between two locations and estimating how irregular the
#' movement path is. Importantly, in dBBMMs the movement between locations is assumed to be diffusive and isotropic. Effectively, this means
#' that the movement between two locations is equally likely to happen in all directions. 
#' Dynamic BGBs allow the incorporation of a directional bias in the movement path analysis of bridging models. Kranstauber et al.'s
#' Bivariate Gaussian Bridges model does this by factorizing the diffusion, i.e. the Brownian variance estimates, into 
#' the direction towards and perpendicular to a straight line to the next location, resulting in two normally distributed
#' probability densities for each direction.
#' 
#' In either model choice, long time gaps between two consecutive detections can therefore present an issue, as these long time
#' gaps create a very high uncertainty in regard to your animal's whereabouts between two locations. However, given the different
#' underlying assumptions in regard to the diffusion process in each model, longer time gaps can cause different problematic outputs.
#' We would like to thank Dr. Anne K. Scharf for the further explanations regarding the model choice. 
#' In short, dBBMM does not distinguish between the animal changing movement in direction or speed between two locations. But, dBGBs make a difference
#' between the animal changing movement in direction and in speed between two locations. Long time gaps in dBGBs can therefore create
#' stronger uncertainty in change in direction so that the outputs of your model runs might show biologically unrealistic, very narrow and long peaks along or perpendicular to the 
#' trajectory between two locations, whereas such a case in dBBMMs would create a uninformative, inflated circular area with a large diameter within which
#' your animal might have been (i.e. pretty much anywhere, which is not informative).
#' 
#' The solution to this problem is in both cases to remove the variance of the segments corresponding to the large time gaps. This can be done using the 'move'
#' package, and is directly implemented in this function already (please see params timeDiffLong and timeDiffUnits above). However, the ideal settings
#' may be different based on your model choice due to the reasons explained above. And, based on your data structure, your data may be inadequate for your model choice.
#' It is therefore imperative that you understand the assumptions of the different models and how they relate to your own data.
#' 
#' Once you have decided, which model is most appropriate for your data, we here provide the movegroup::movegroup() function
#' that automates the following steps in calculating dBBMM/dBGB-based UDs:
#' 
#' Step 1. Filters individuals.
#' Remove those individuals for which there are insufficient data i.e. number of re-locations is
#' smaller than the window size parameter value (defaults: 23 for model "dbbmm" and 31 for model "dBGB").
#'
#' Step 2. Generates universal raster.
#' Based on all remaining data, a universal raster is generated where the calculated UDs are plotted
#' into.
#'
#' Step 3. Loops through individuals.
#' Individuals are looped through to construct individual-level movement models (on an absolute
#' scale).
#'
#' See www.GitHub.com/SimonDedman/movegroup for issues, feedback, and development
#' suggestions.
#'
#' install_git('https://gitlab.com/bartk/move.git') #Installs 'move' development version
#' 
#' The movegroup package provides additional functions to automate and simplify the visualisation and quantification of space
#' use data of a group or populations of animals. When used together, the order of these functions would be: movegroup, (stitchraster
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
#' 9. Error in validityMethod(as(object, superClass)) : The used raster is not a UD (sum unequal to 1), sum is: NaN One possible cause is 
#' loss of accuracy due to writing raster to disk with dataType FLT4S this can be solved preventing disk usage or changing data type. 
#' Same error message as in nr. 8 above. This error can also occur if you are calculating dBGBs and set your timeDiffLong shorter than the shortest
#' segment. Check your time gaps between locations and compare to timeDiffLong, check the timeDiffUnits parameter and make sure you chose the correct unit.
#' Try shorter timeDiffLong options.
#'
#' 10. Function runs but no findable output: The default directory into which outputs are saved (i.e. savedir) is set to 
#' 'tempdir()'. This means, that if you forget to change this parameter, all outputs will be saved to your system's temporary folder.
#' 
#' @examples
#' \donttest{
#' # load data
#' data("TracksCleaned")
#' # run function
#' movegroup(
#'   data = TracksCleaned,
#'   ID = "Shark",
#'   Datetime = "Datetime",
#'   Lat = "Lat",
#'   Lon = "Lon",
#'   savedir = tempdir()
#' )
#' }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}
#' @author Vital Heim, \email{vital.heim@@gmail.com}
#'
#' @references Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M. and Safi, K. (2012), A
#' dynamic Brownian bridge movement model to estimate utilization distributions for heterogeneous
#' animal movement. Journal of Animal Ecology. doi: 10.1111/j.1365-2656.2012.01955.x
#' @references Kranstauber, B., Safi, K., and Bartumeus, F. (2014). Bivariate Gaussian bridges: directional
#' factorization of diffusion in Brownian bridge models. Movement Ecology. https://doi.org/10.1186/2051-3933-2-5
#' @references Kranstauber, B., M. Smolla & A. K. Scharf. 2019. Move: visualizing and analyzing animal track
#' data. R package version 4.2.4 (at 2023-08-15). https://CRAN.R-project.org/package=move.

movegroup <- function(
    data = NULL, # Data frame object containing the data. Requires columns Lat Lon DateTime ID and optionally a grouping column.
    model = "dbbmm", # Movement model choice: "dbbmm" (dynamic Brownian bridge movement model,
    # default) or "dbgb" (dynamic bivariate Gaussian bridge). Outputs are saved into a
    # subdirectory of savedir named after the chosen model (savedir/dbbmm or savedir/dbgb).
    # Note: dBGB uses the non-bursted move object r.i and does not support bursting via
    # timeDiffLong/timeDiffUnits. dBBMM uses the bursted move object as before.
    ID = NULL, # Name of animal tag ID column in data.
    Datetime = NULL, # Column name in data that contains date/time stamps for each recorded detection. Must be in POSIXct format.
    Lat = NULL, # name of Lat & Lon columns in data.
    Lon = NULL,
    # Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
    dat.TZ = "US/Eastern", # timezone for as.POSIXct.
    proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
    sensor = "acoustic-telemetry", # sensor for move function. Single character or vector with length of the number of coordinates. Optional. This is just for personal record, and not used in any calculations. 
    integrationStep = NULL, # in minutes. Defines the size of the time intervals taken for every integration step and thus specifies the temporal resolution of the numerical integration. If acoustic transmitters are used, we suggest that the time.step should reflect the ping frequency of the tag (in minutes). If NULL 15 steps are taken in the shortest time interval for dBBMM, and 20.1 steps in the shortest time interval for dBGB. Optional.
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
    movemargin = NULL, # Margin size for variance calc in move::brownian.motion.variance.dyn and behavioral change point analysis in
    # move::brownian.bridge.dyn. Must be an odd number. Default 11.
    dbbext = 3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 3.
    dbbwindowsize = NULL, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more
    # stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior.
    # This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 23). Must be >=2*margin (movemargin) which is 11 so >=22,
    # but odd so >=23.
    writeRasterFormat = "ascii",
    writeRasterExtension = ".asc",
    writeRasterDatatype = "FLT4S",
    absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
    savedir = tempdir(), # save outputs to a temporary directory (default) else change to current
    # directory e.g. "/home/me/folder". Subdirectory based on model choice is automatically created. Do not use getwd() here. This is the "root" directory. Do not specify potential grouping variables here.
    saveAreaCT = FALSE, # save tiny individual core and general use areas tables
    # to disk. These are the only things retained in the per-individual loop, so
    # if your large dataset causes memory crashes, you can run it in chunks and
    # stitch the results together later with stitchraster.
    alerts = TRUE # audio warning for failures
) {
  
  # first, validate a valid model choice
  model <- match.arg(model, choices = c("dbbmm", "dbgb"))
  
  # based on model choice we set default values for parameters window size and margin
  model_defaults <- list(
    dbbmm = list(windowsize = 23, margin = 11),
    dbgb  = list(windowsize = 31, margin = 15)
  )
  if (is.null(dbbwindowsize)) dbbwindowsize <- model_defaults[[model]]$windowsize
  if (is.null(movemargin)) movemargin <- model_defaults[[model]]$margin
  
  # If savedir has a terminal slash, remove it, it's added later
  if (substr(x = savedir, start = nchar(savedir), stop = nchar(savedir)) == "/") {
    savedir <- substr(x = savedir, start = 1, stop = nchar(savedir) - 1)
  }
  
  # automatically create a subdirectory within savedir based on model choice and potential grouping variables
  # savedir <- file.path(savedir, model)
  
  # if savedir doesn't exist, can't save or setwd into it, therefore create it
  if (!file.exists(savedir)) dir.create(savedir, recursive = T)
  
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts) {
    options(
      error = function() {
        beepr::beep(9) # give warning noise if it fails
        graphics.off() # kill all graphics devices
        setwd(oldwd) # reinstate original working directory. Probably redundant given on.exit
      } # close options subcurly
    )
  } # close options
  
  # edit center if centre used
  if (!is.null(centre)) center <- centre
  
  # Create writeRasterExtension from writeRasterFormat
  if (is.null(writeRasterExtension)) {
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
  }
  
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
  
  # A dBBMM or dBGB is not run if total detections of individual < window size (default value, 31).
  # Below code checks and filters out individuals with insufficient data.
  enoughrelocs <- data |>
    group_by(ID) |>
    summarise(n = n()) |>
    filter(n >= dbbwindowsize) |>
    pull(ID)
  data <- data |> filter(ID %in% enoughrelocs)
  rm(enoughrelocs)
  
  # Create per-individual move object, project it, construct a dBBMM or dBGB, calculate the volume area
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
      ), # convert NaN to NA. POSIX needs lubridate
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
      paste0("processing ", model, " "),
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
    
    # Calculate time lag between consecutive detections. - DEPRECATED
    # Obviously there is no time lag for the first detection, so we need to tell R this. - DEPRECATED
    # Knowing time lags in acoustic telemetry is important, as the motion variance is based on the-  DEPRECATED
    # time difference b/w detections at consecutive locations i.e. larger time lag creates a wider - DEPRECATED
    # bridge where the animal could've been and therefore potentially inflates the motion variance. - DEPRECATED
    # Below code deals with this issue by identifying the location and number of detections with - DEPRECATED
    # large time gaps. Below we will use an arbitrary value of 2 h. - DEPRECATED
    # TimeDiff <- move::timeLag(r.i, units = timeDiffUnits) - DEPRECATED
    # r.i$TimeDiff <- append(0, TimeDiff) - DEPRECATED
    # long <- which(r.i$TimeDiff > timeDiffLong) - DEPRECATED
    # if no timediffs are longer than timedifflong, long is integer(0), a potentially dangerous object. - DEPRECATED
    
    # Below we exclude the data points that are 'long' i.e. create a large time gap. move::burst ? - DEPRECATED
    # bursted <- move::burst( - DEPRECATED
    #   r.i, - DEPRECATED
    #   c("normal", "long")[ - DEPRECATED
    #     1 + (move::timeLag(r.i, units = timeDiffUnits) > timeDiffLong) - DEPRECATED
    #   ] - DEPRECATED
    # ) - DEPRECATED
    # rm(r.i) - DEPRECATED
    # There are 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command. - DEPRECATED
    
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
    
    rm(data.i) # not needed after this in loop
    
    # location error needs to be a positive number. Replace zeroes with 0.00001
    moveLocError.i[which(moveLocError.i == 0)] <- 0.00001
    
    # START of code block for model calculation based on model choice by user 
    # if statement based on your model choice for model fitting and variance calculation
    if (model == "dbbmm") {
      
      # First, calculate the dynamic Brownian motion variance for the movement track
      dynbbmm_var <- move::brownian.motion.variance.dyn(
        r.i,
        location.error = moveLocError.i,
        window.size = dbbwindowsize, #  must be >=2*margin which is default 11 so >=22, but odd so >=23
        margin = movemargin
      )
      mtnvar <- as.data.frame(dynbbmm_var) # extract the data
      mtnvar$motionVariance <- move::getMotionVariance(dynbbmm_var) # add motionvariance
      write.csv(
        mtnvar,
        file = file.path(savedir, paste0(i, "_dynBrownianMotionVariance.csv")),
        row.names = FALSE
      )
      
      # save your dynbbmm variance object
      saveRDS(
        dynbbmm_var,
        file = file.path(savedir, paste0(i, "_dBBMMvariance.Rds"))
      )
      
      # Data that contains large time gaps between two consecutive locations can be problematic, as the uncertainty
      # where the animal might have been during this time is very large. This can result in very large, uninformative UDs following
      # model calculation.
      # To solve this problem, one can remove the variance of the segments corresponding to long time gaps as suggested by Kranstauber et al.
      # The dynbbmm_var object calculated above contains am '@interest' slot, that when set to FALSE for a given segment, 
      # this segment will be ignored in the calculation of the dBBMM. 
      # Here we use the 'timeDiffLong' argument in 'timeDiffUnits' as threshold so that any segment with time lag larger than the specified value will be set to FALSE.
      dynbbmm_var@interest[move::timeLag(r.i, units = timeDiffUnits) > timeDiffLong] <- FALSE
      
      # Now, the adjusted 'dynbbmm_var' object can be used for the dBBMM calculation.
      
      # Construct the model. The time.step should reflect the ping frequency of the tag (in minutes). If specified as NULL, use the function without the timeStep argument.
      # This is what takes CPU & memory
      # bursted_dbbmm <- move::brownian.bridge.dyn( - DEPRECATED
      #   bursted, - DEPRECATED
      #   burstType = "normal", - DEPRECATED
      #   raster = xAEQD, # has to be AEQD for metre-based calculations, presuambly? - DEPRECATED
      #   # Error:The projection of the raster and the Move object are not equal. - DEPRECATED
      #   # Need bursted, r.i, to be the same projection as xAEQD - DEPRECATED
      #   location.error = moveLocError.i, - DEPRECATED
      #   ext = dbbext, # dbbext - DEPRECATED
      #   margin = movemargin, - DEPRECATED
      #   window.size = dbbwindowsize, #  must be >=2*margin which is 11 so >=22, but odd so >=23 - DEPRECATED
      # ) - DEPRECATED
      
      dbbmm_mod <- if (is.null(integrationStep)) {
        move::brownian.bridge.dyn(
          dynbbmm_var,
          raster = xAEQD,
          location.error = moveLocError.i,
          ext = dbbext,
          margin = movemargin,
          window.size = dbbwindowsize # must be >=2*margin which is 11 so >=22, but odd so >=23
        )
      } else {
        move::brownian.bridge.dyn(
          dynbbmm_var,
          raster = xAEQD,
          location.error = moveLocError.i,
          ext = dbbext,
          margin = movemargin,
          window.size = dbbwindowsize, # must be >=2*margin which is default 11 so >=22, but odd so >=23
          time.step = integrationStep/15
        )
      }
      
      # Save your dBBMM object
      saveRDS(
        dbbmm_mod,
        file = file.path(savedir, paste0(i, "_dbbmm.Rds"))
      )
      
      # rm(bursted) - DEPRECATED
      rm(dynbbmm_var)
      rm(mtnvar)
      
      # Re-standardize (Dr. Kranstauber's trouble shooting solution)
      # Occurring errors are "due to limits of accuracy during calculations.
      # Given the error is really small (<.000001 %) I would not worry to much and re-standardize".
      # Then aggregate UD segments.
      # 20260630 - SD and VH discussed removing this block, to investigate in future
      tmp <- raster::calc(dbbmm_mod, sum) # raster::calc :Calculate values for a new Raster* object from another Raster* object
      # using a formula. returns a RasterLayer if fun returns a single value (e.g. sum)
      rm(dbbmm_mod)
      bb <- new(".UD", tmp / sum(raster::values(tmp))) # new() creates object from Class ("UD.")
      rm(tmp)
      
    } else {
      # dBGB branch based on your model choice
      # Similar to the dBBMM branch, first the dynBGB-variance is calculated, then segments of interests within the object are labelled
      # as false and ignored, if time lag is larger than specified in timeDiffLong, and finally the updated dBGBvariance object
      # is used in the dynBGB() calculation.
      
      # calculate the bivariate Gaussian bridge orthogonal and parallel variance along the movement track
      dynbgb_var <- move::dynBGBvariance(
        r.i,
        locErr = moveLocError.i,
        windowSize = dbbwindowsize,
        margin = movemargin
      )
      mtnvar <- as.data.frame(dynbgb_var) # extract the data
      mtnvar$OrthoParallelVariance <- move::getMotionVariance(dynbgb_var) # add motionvariance
      write.csv(
        mtnvar,
        file = file.path(savedir, paste0(i, "_OrthogonalParallelVariance.csv")),
        row.names = FALSE
      )
      
      # save your dynbgb variance object
      saveRDS(
        dynbgb_var,
        file = file.path(savedir, paste0(i, "_dBGBvariance.Rds"))
      )
      
      # now, we ignore all segments that have a larger time lag than timeDiffLong. The 'dBGBvariance' object resulting from the function above,
      # contains the slot '@segInterest' in which those segments marked as FALSE won't be included in the calculation of the dBGB.
      dynbgb_var@segInterest[move::timeLag(r.i, units = timeDiffUnits) > timeDiffLong] <- FALSE
      
      # now we construct the dynamic bivariate Gaussian bridge model using the dynBGBvariance object
      # The time.step should reflect the ping frequency of the tag (in minutes). If specified as NULL, use the function without the timeStep argument.
      # dynbgb_mod <- move::dynBGB(
      #   dynbgb_var,
      #   raster = xAEQD,
      #   locErr = moveLocError.i,
      #   ext = dbbext,
      #   windowSize = dbbwindowsize,
      #   margin = movemargin
      #   )
      
      dynbgb_mod <- if (is.null(integrationStep)) {
        move::dynBGB(
          dynbgb_var,
          raster = xAEQD,
          locErr = moveLocError.i,
          ext = dbbext,
          windowSize = dbbwindowsize,
          margin = movemargin
        )
      } else {
        move::dynBGB(
          dynbgb_var,
          raster = xAEQD,
          locErr = moveLocError.i,
          ext = dbbext,
          windowSize = dbbwindowsize,
          margin = movemargin,
          timeStep = integrationStep/20.1
        )
      }
      
      # save your dBGB object
      saveRDS(
        dynbgb_mod,
        file = file.path(savedir, paste0(i, "_dbgb.Rds"))
      )
      
      # remove files we dont need anymore
      rm(dynbgb_var)
      rm(mtnvar)
      
      # Re-standardize and aggregate UD segments, following the same pattern as dBBMM above
      # 20260630 - SD and VH discussed removing this block, to investigate in future
      tmp <- raster::calc(dynbgb_mod, sum)
      rm(dynbgb_mod)
      bb <- new(".UD", tmp / sum(raster::values(tmp)))
      rm(tmp)
    } # END of branches for model choice
    # END of code block for model calculations
    
    # Calculate volume area (m^2) within 50% (core) and 95% (general use) contours. Note: absolute scale
    # Note: Based on 'An introduction to the 'move' package', it's the opposite of what you might expect:
    # "A cell with a very high value in the UD raster will have a very low value in the contour raster"
    
    # Below calc introduced in 2022-10-08 commit, message =
    # "changed code volume area: instead of using getVolumeUD() from the move package,
    # we just did the calculations on the UD raster and not using the function.
    # these seem to produce very plausible estimates!"
    # This change creates unnaturally small areas.
    # Changed back to previous calcs 2023-02-24. Added rounding - probelmatic as these are raw probabilities - UD should be volumetric UD. Changed back to getVolumeUD()
    # area.50.new <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.5)), 4)
    # area.95.new <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.05)), 4)
    # area.50.new <- sum(raster::values(bb) >= (max(bb@data@values) * 0.5)) # 2024-01-08 remove rounding here. Seemingly had no effect - DEPRECATED
    # area.95.new <- sum(raster::values(bb) >= (max(bb@data@values) * 0.05)) # 2024-01-08 remove rounding here. Seemingly had no effect - DEPRECATED
    
    # Below calculation introduced 20260624: get back to getVolumeUD()
    # getVolumeUD() can be used to calculate the area in which an animal has some specified probability of being located as the
    # cumulative distribution function. This means that rather than using raw probability values as above (2022-10-08) that are
    # probability densities and not probabilities so that their magnitude depend on arbitrary choices such as e.g. rastersize,
    # getVolumeUD() ranks cells from most to least used and adds them up until they accumulate the specified % of all usage.
    area.50.new <- sum(raster::values(move::getVolumeUD(bb) <= .50))
    area.95.new <- sum(raster::values(move::getVolumeUD(bb) <= .95)) 
    
    # Combine in single df
    area.ct <- data.frame(
      core.use.new = area.50.new,
      general.use.new = area.95.new
    ) # 2024-01-08 add rounding here if necessary. Maybe not. No, not needed.
    
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
      if (writeRasterFormat != "CDF") bylayer <- TRUE, # bylayer kills ncdf4
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
  
  # calculate space use areas in km2
  md$core.use.km2 <- (rasterres * md$core.use.new) / 1000000 # convert from cells/pixels to metres squared area based on cell size, then to kilometres squared area
  md$general.use.km2 <- (rasterres * md$general.use.new) / 1000000
  
  md <- md |> 
    dplyr::rename(
      core.use.cells = core.use.new,
      general.use.cells = general.use.new
    )
  
  write.csv(
    md,
    file = file.path(savedir, absVolumeAreaSaveName),
    row.names = FALSE
  )
  
  # 2023-10-04 Vital memory bug warning
  if (all(md$core.use.km2 == md$core.use.km2[1])) {
    message(
      "All core UDs identical. Maybe insufficient memory for raster calcs - check rasterResolution"
    )
  }
  if (all(md$general.use.km2 == md$general.use.km2[1])) {
    message(
      "All general UDs identical. Maybe insufficient memory for raster calcs - check rasterResolution"
    )
  }
  if (length(which(md$core.use.km2 == max(md$core.use.km2, na.rm = TRUE))) > 1) {
    message(
      "More than 1 individual share exactly the same max value for core use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
  }
  if (length(which(md$general.use.km2 == max(md$general.use.km2, na.rm = TRUE))) > 1) {
    message(
      "More than 1 individual share exactly the same max value for general use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
  }
  if (length(which((md$core.use.km2 - md$general.use.km2) == 0)) > 0) {
    message(
      "1 or more individuals have exactly the same value for core and general use, maybe insufficient memory for raster calcs - check rasterResolution"
    )
  }
} # close function