#' Automates dynamic Brownian bridge construction across individuals
#'
#' Automates dynamic Brownian bridge movement model calculation for utilization distribution (UD) estimation for multiple 
#' individuals simultaneously, via the 'brownian.bridge.dyn()' function in the 'move' package. 
#' 
#' Step 1. Filter individuals. 
#' Remove those individuals for which there is insufficient data i.e. number of re-locations is smaller than the window size parameter value (default = 31). 
#' 
#' Step 2. Generate universal raster. 
#' Based on all remaining data, a universal raster is generated where the calculated UDs are plotted into. 
#' 
#' Step 3. Loop through individuals.
#' Individuals are looped through to construct individual-level movement models (on an absolute scale). See www.GitHub.com/SimonDedman/dBBMMhomeRange for issues, feedback, and development suggestions. 
#'
#' @param data Data frame object containing the data. Requires columns Lat Lon 
#' DateTime ID and optionally a grouping column. Names specifide in later 
#' parameters. Grouping not currently implemented 2022-11-15, see Group 
#' parameter below.
#' @param ID Name of animal tag ID column in data.
#' @param Datetime Column name in data that contains date/time stamps for each 
#' recorded detection. Must be in POSIXct format.
#' @param Lat Name of latitude column in data.
#' @param Lon Name of longitude column in data.
#' @param Group Name of grouping column in data. UNUSED 2022-11-15
#' @param dat.TZ Timezone of data for as.POSIXct.
#' @param proj CRS for move function.
#' @param projectedCRS EPSG code for CRS for initial transform of latlon points;
#'  corresponds to rasterCRS zone.
#' @param sensor Sensor for move function. Single character or vector with 
#' length of the number of coordinates. Optional.
#' @param moveLocError Location error (m) in the 'brownian.bridge.dyn' function 
#' in the 'move' package. Numeric. Either single or a vector of length nrow 
#' data. If using passive acoustic data this is the detection range of the 
#' receiver(s).
#' @param timeDiffLong Single numeric value. Threshold value in timeDiffUnits 
#' designating the length of long breaks in re-locations. Used for bursting a 
#' movement track into segments, thereby removing long breaks from the movement 
#' track. See ?move::bursted for details.
#' @param timeDiffUnits Character. Unit for timeDiffLong.
#' @param center Do you want to center the move object within extent? See 
#' spTransform.
#' @param buffpct Buffer extent for raster creation, proportion of 1.
#' @param rasterExtent Extent of raster created around data. If NULL, raster 
#' extent calculated from data, buffpct, rasterResolution. Else length 4 vector,
#'  c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 (degrees) north
#'  or south for ymax or ymin. Doesn't prevent constraint to data limits (in 
#' plot anyway), but prevents raster clipping crash.
#' @param rasterCRS CRS for raster creation.
#' @param rasterResolution Single numeric value to set raster resolution - cell 
#' size in metres? 111000: 1 degree lat = 111km.
#' @param dbblocationerror Location.error param in 'brownian.bridge.dyn' 
#' function in the 'move' package. single numeric value or vector of the length 
#' of coordinates that describes the error of the location (sender/receiver) 
#' system in map units. Or a character string with the name of the column 
#' containing the location error can be provided. Could use the same as 
#' moveLocError?.
#' @param dbbext Ext param in the 'brownian.bridge.dyn' function in the 'move' 
#' package. Extends bounding box around track. Numeric single (all edges), 
#' double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3.
#' @param dbbwindowsize window.size param in the 'brownian.bridge.dyn' function 
#' in the 'move' package. The size of the moving window along the track. Larger 
#' windows provide more stable/accurate estimates of the brownian motion 
#' variance but are less well able to capture more frequent changes in behavior.
#'  This number has to be odd. A dBBMM is not run if total detections of 
#'  individual < window size (default 31).
#' @param writeRasterFormat Character. Output file type. ascii. TO DEPRECIATE.
#' @param writeRasterExtension Character. Output file extension. Should match 
#' character writeRasterFormat. TO DEPRECIATE.
#' @param writeRasterDatatype Character. Data type for writing values to disk. 
#' TO DEPRECIATE.
#' @param absVolumeAreaSaveName File name plus extension where UD estimates are 
#' saved. TO DEPRECIATE.
#' @param savedir Save outputs to a temporary directory (default) else change 
#' to desired directory e.g. "/home/me/folder/". Do not use getwd() for this. 
#' Include terminal slash.
#' @param alerts Audio warning for failures.
#' 
#' @return Individual-level utilization distributions, saved as rasters, as well
#'  as calculated volume area estimates for 50 and 95pct contours, saved in a 
#'  .csv file
#' 
#' saved to disk.
#' @details Errors and their origins:
#' @examples
#' \donttest{
#' # Not run
#' }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}
#'
#' @export

# install_git('https://gitlab.com/bartk/move.git') #Installs 'move' development version
#' @import utils
#' @import magrittr
#' @importFrom lubridate is.POSIXt
#' @importFrom sp CRS SpatialPoints spTransform proj4string
#' @importFrom beepr beep
#' @importFrom dplyr mutate rename group_by summarise filter semi_join distinct ungroup arrange across bind_cols pull
#' @importFrom tidyr drop_na
#' @importFrom raster raster projectExtent res ncell setValues calc values writeRaster
#' @importFrom move move timeLag burst brownian.bridge.dyn getVolumeUD
#' @importFrom rlang .data
#' @importFrom grDevices graphics.off
#' @importFrom graphics par
#' @importFrom methods new
#' @importFrom stats setNames

dBBMMhomeRange <- function(
    data = NULL, # Data frame object containing the data. Requires columns Lat Lon DateTime ID and optionally a grouping column.
    ID = NULL, # Name of animal tag ID column in data.
    Datetime = NULL, # Column name in data that contains date/time stamps for each recorded detection. Must be in POSIXct format.
    Lat = NULL, # name of Lat & Lon columns in data.
    Lon = NULL,
    Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
    dat.TZ = "US/Eastern", # timezone for as.POSIXct.
    proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
    projectedCRS = "+init=epsg:32617", # EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone. This is around Bimini, Bahamas.
    sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
    moveLocError = 1, # location error in metres for move function. Numeric. Either single or a vector of length nrow data.
    timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
    timeDiffUnits = "hours", # units for time difference for move function.
    center = TRUE, # center move object within extent? See spTransform.
    buffpct = 0.3, # buffer extent for raster creation, proportion of 1.
    rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
    # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
    rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. This is around Bimini, Bahamas.
    rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
    dbblocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
    dbbext = 3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
    dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
    writeRasterFormat = "ascii",
    writeRasterExtension = ".asc",
    writeRasterDatatype = "FLT4S",
    absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
    savedir = tempdir(),  # save outputs to a temporary directory (default) else.
    # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
    alerts = TRUE # audio warning for failures
) {
  # TODO:
  # option to use motion variance as the dependent variable, not the UD
  # clean up all notes, package elements, authorship, dependencies etc.
  # Add examples
  
  # Generalised Boosting Model / Boosted Regression Tree process chain automater.
  # Simon Dedman, 2012-6 simondedman@gmail.com GitHub.com/SimonDedman/gbm.auto
  
  # source("scaleraster.R") # might not be needed if included as function in same package
  
  # @import circular
  # @import devtools
  # @import knitr
  
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts) options(error = function() {
    beepr::beep(9)# give warning noise if it fails
    graphics.off()# kill all graphics devices
    setwd(oldwd) # reinstate original working directory. Probably redundant given on.exit
  } # close options subcurly
  ) # close options
  
  # library(sp)
  # library(move)
  # if (writeRasterFormat == "CDF") library(ncdf4)
  # do this for other formats?####
  
  
  # Create writeRasterExtension from writeRasterFormat
  if (is.null(writeRasterExtension)) writeRasterExtension <- switch(EXPR = writeRasterFormat,
                                                                    "ascii" = ".asc",
                                                                    "raster" = ".grd",
                                                                    "SAGA" = ".sdat",
                                                                    "IDRISI" = ".rst",
                                                                    "CDF" = ".nc",
                                                                    "GTiff" = ".tif",
                                                                    "ENVI" = ".envi",
                                                                    "EHdr" = ".bil",
                                                                    "HFA" = ".img")
  
  # Function to automate the many steps required to use boosted regression trees
  # to predict abundances in a delta process, i.e. binary (0/1) proportion
  # prediction coupled with presence-only abundance prediction to give total
  # prediction. Loops through all permutations of parameters provided (learning
  # rate, tree complexity, bag fraction), chooses the best, then tries to simplify
  # that. Generates line, dot and bar plots, and outputs these and the predictions
  # and a report of all variables used, statistics for tests, variable
  # interactions, predictors used and dropped, etc.. If selected, generates
  # predicted abundance maps, and Unrepresentativeness surfaces.
  #
  # Underlying functions are from packages gbm and dismo, functions from Elith
  # et al. 2008 (bundled as gbm.utils.R), mapplots, and my own functions gbm.map,
  # gbm.rsb, gbm.valuemap, gbm.cons, gbm.basemap
  
  ####1. Check packages, start loop####
  # fam1 <- match.arg(fam1) # populate object from function argument in proper way
  # # tibble's don't collapse into a vector, instead an X x 1 df, which breaks various functionality.
  # if ("tbl" %in% class(grids)) grids <- as.data.frame(grids)
  
  # prefixes "X" to numerical-named sharks to avoid issues later
  data %<>% 
    # tmp <- data %>% mutate(.data[[ID]] = make.names(.data[[ID]]))
    # tmp <- data %>% mutate(.data[, ID] = make.names(.data[, ID]))
    # check this works with a dynamic ID name####
  # Unless I can rename the user entry to "ID" # done
  # Same for Lat & Lon
  # And Datetime. And group.
  dplyr::rename(Datetime = .data[[Datetime]],
                ID = .data[[ID]],
                Lat = .data[[Lat]],
                Lon = .data[[Lon]]) %>%
    dplyr::mutate(ID = make.names(ID))
  
  if (!any(class(data$Datetime) == "POSIXct")) stop("Datetime column must be POSIXct")
  
  # Add movelocationerror and dbblocationerror as columns if they're vectors, else their length
  # desyncs from nrow(data) if any IDs < windowsize
  if (nrow(data) == length(moveLocError)) data$moveLocError <- moveLocError
  if (nrow(data) == length(dbblocationerror)) data$dbblocationerror <- dbblocationerror
  
  
  
  
  # title: "Lemon sharks ecology Bimini"
  # author: "Sprinkles (Maurits van Zinnicq Bergmann) & Simon Dedman" 
  # date: "2021-09-07"
  
  # param: working directory
  
  ## Data loading
  
  # Load data set
  
  # ```{r read_files}
  # print(getwd()) # "/home/simon/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R"
  # data <- read.csv("../data/TRACKS1.csv")
  # the script environment location is where the script is i.e. dBBMM_HomeRange/R.
  # the console environment is the project root i.e. dBBMM_HomeRange
  # Error in file(file, "rt") : cannot open the connection
  # works if you paste it into the console. The command is correct. What's going on? Something with markdown?
  # env1 <- Sys.getenv() checked against console, only difference is console width
  # ```
  
  # Check data directory and create it if not present
  
  # ```{r data directory}
  # # set path to general folder where output files will be written
  # if (!dir.exists(data.dir)) {
  #   stop(paste0(
  #     "Directory not found: \n",
  #     "       ", gsub("\\.", getwd(), data.dir)
  #   ))
  # }
  # 
  # out.dir <- "dBBMM ASCII"
  # if (!dir.exists(out.dir)) dir.create(file.path(data.dir, out.dir))
  
  ### would be useful to include in the function a path to output dir that is user-defined
  ### instead of stopping if R detects no data.dir, would be better to create this dir like for out.dir. add later
  # ```
  
  
  # Explore the data set
  # 
  # ```{r explore}
  # names(data)
  # str(data)
  # dim(data)
  # head(data)
  # 
  # table(data$Tidal.Phase) %>%
  #   kbl(caption = "Detection frequency per tidal phase") %>%
  #   kable_classic(full_width = F, html_font = "Cambria")
  # 
  # table(data$ID) %>%
  #   kbl(caption = "Detection frequency per shark") %>%
  #   kable_classic(full_width = F, html_font = "Cambria")
  # 
  # table(data$ID, data$Tidal.Phase) %>%
  #   kbl(caption = "Detection frequency per shark and tidal phase") %>%
  #   kable_classic(full_width = F, html_font = "Cambria")
  # 
  # table(data$ID, data$Tidal.Phase)
  # ```
  
  # Housekeeping and tidy up
  # 
  # ```{r housekeeping}
  # data %<>%
  #   mutate(Datetime = as.POSIXct(Datetime,
  #                                format = "%m/%d/%y %H:%M",
  #                                tz = dat.TZ
  #   ),
  #   ID = make.names(ID) # prefixes "X" to numerical-named sharks to avoid issues later
  #   ) %>%
  #   rename(
  #     Lat = N,
  #     Lon = W,
  #     T.Ph = Tidal.Phase
  #   ) %>%
  #   select(Datetime, ID, T.Ph, Lat, Lon) %>% # dropped 2 variables (Date, Time)
  #   arrange(ID, Datetime) # df size: 1308 x 5
  
  # write.csv(x = data, file = "TracksCleaned.csv", row.names = FALSE) # Could load this directly here
  # ```
  
  # Convert lonlat to UTM
  # 
  # ```{r coord_UTM}
  
  # First convert coordinate sets to SpatialPoints and project
  cord.dec <- sp::SpatialPoints(cbind(data$Lon, data$Lat),
                                proj4string = sp::CRS("+proj=longlat")
  )
  
  # Transform to UTM by setting the EPSG to 32617 for WGS 84, UTM zone 17, northern hemisphere. This is where Bimini is located.
  cord.UTM <- as.data.frame(sp::spTransform(cord.dec, sp::CRS(projectedCRS)))
  colnames(cord.UTM) <- c("NewEastingUTM", "NewNorthingUTM")
  data <- cbind(data, cord.UTM) # 1308 x 7
  # ```
  
  # Construct movement models per individual.
  
  # Some notes regarding the 'move package and construction of movement models: Several arguments need to be used to run the model. 1. Window size: corresponds with number of locations and moves along a given trajectory to estimate the MA parameter within defined subsections of the path. This increases the ability to detect breakpoints where changes in behaviour occur. The window size should relate to what kind of behaviours the model is desired to identify e.g., a window size of 23 means the sliding window is moved every 23 locations or every 23 hours (has to do with sampling interval) 2. Margin: motion variance based on only the middle section of the trajectory; the ends of the movement trajectory where no changes are allowed because at some stage you want to have a few locations to base your estimation of the variance on and how many locations in either side of the window we use for this, is called the margin. Smaller values for window size and margin is expected to give a higher frequency of behavioural changes; make these large for looking at migrations. 3. The raster dictates the grid cell size for the UD to be calculated per grid cell per individual. Create the raster of certain size that matches with coordinates used to make the move object 4. Extent: is incorporated if there are animal locations that border the edges of the raster.
  
  # A dBBMM is not run if total detections of individual < window size (default value, 31).
  # Below code checks and filters individuals with insufficient data.
  # Below we set window size arbitrarily to 23
  
  # ```{r filter_data}
  check1 <- data %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::summarise(relocations = length(.data$Datetime))
  check2 <- dplyr::filter(check1, .data$relocations >= dbbwindowsize) # filter: removed 2 rows (14%), 12 rows remaining
  
  if (length(check1$ID) != length(check2$ID)) {
    data <- dplyr::semi_join(data, check2) # Joining, by = "ID". semi_join: added no columns
    check1 <- data %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(relocations = length(.data$Datetime))
    check2 <- dplyr::filter(check1, .data$relocations >= dbbwindowsize) # filter: no rows removed
    length(check1$ID) == length(check2$ID)
  } # data: 1253 x 7
  # TODO: improve this####
  # also
  # print which IDs get removed ####
  
  
  # In the next R code we create per individual a move object, project it, 
  # construct a dBBMM, calculate the volume area within the 50% and 95% contours and
  # finally save the outcome as ASCII for plotting in an external GIS software.
  
  # ```{r construct_dBMMMM}
  bb <- list()
  bb.list <- list()
  
  data %<>%
    tidyr::drop_na(ID) %>%
    dplyr::group_by(ID) %>%
    dplyr::distinct(Datetime, .keep_all = TRUE) %>% # distinct (grouped): removed one row (<1%), 1,286 rows remaining
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
  if (!is.null(rasterExtent)) { # if xUTM object exists
    xUTM <- sp::SpatialPoints(cbind(c(rasterExtent[1], rasterExtent[2]), c(rasterExtent[3], rasterExtent[4])), proj4string = sp::CRS("+proj=longlat"))
    xUTM <- as.data.frame(sp::spTransform(xUTM, sp::CRS(projectedCRS))) # Transform to UTM
    xUTM <- raster::raster( # create a raster
      xmn = xUTM[1,1] - ((xUTM[2,1] - xUTM[1,1]) * buffpct), # with xUTM's values as extents
      xmx = xUTM[2,1] + ((xUTM[2,1] - xUTM[1,1]) * buffpct),
      ymn = xUTM[1,2] - ((xUTM[2,2] - xUTM[1,2]) * buffpct),
      ymx = xUTM[2,2] + ((xUTM[2,2] - xUTM[1,2]) * buffpct),
      crs = rasterCRS, # +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs
      resolution = rasterResolution # 50
    ) # close raster
  } else { # exists(xUTM), else create from data
    xUTM <- raster::raster(
      xmn = min(data$NewEastingUTM, na.rm = TRUE) - ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
      xmx = max(data$NewEastingUTM, na.rm = TRUE) + ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
      ymn = min(data$NewNorthingUTM, na.rm = TRUE) - ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
      ymx = max(data$NewNorthingUTM, na.rm = TRUE) + ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
      crs = rasterCRS, # +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs
      resolution = rasterResolution # 50
    ) # close raster
  } # close else
  sp::proj4string(xUTM) # "+proj=utm +zone=27 +datum=WGS84 +units=m +no_defs"
  xUTM
  # class      : RasterLayer 
  # dimensions : 97, 148, 14356  (nrow, ncol, ncell)
  # resolution : 111000, 111000  (x, y)
  # extent     : -11061354, 5366646, 462297.3, 11229297  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs 
  
  
  
  # run move on all data together to generate global CRS for raster
  alldata <- data %>%
    dplyr::arrange(Datetime) %>%
    dplyr::group_by(Datetime) %>% # remove duplicates
    dplyr::summarise(dplyr::across(where(~ is.numeric(.)), mean, na.rm = TRUE),
                     dplyr::across(where(~ is.character(.) | lubridate::is.POSIXt(.)), dplyr::first)) %>% # library(lubridate)
    dplyr::mutate(dplyr::across(where(~ is.numeric(.)), ~ ifelse(is.nan(.), NA, .)), #convert NaN to NA. POSIX needs lubridate
                  dplyr::across(where(~ is.character(.)), ~ ifelse(is.nan(.), NA, .))) %>% # https://community.rstudio.com/t/why-does-tidyrs-fill-work-with-nas-but-not-nans/25506/5
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
  # document() error####
  # Error in (function (classes, fdef, mtable): unable to find an inherited method for function ‘move’ for signature ‘"numeric", "numeric", "character", "data.frame", "CRS"’
  
  
  # Convert projection to Azimuthal Equi-Distance projection (aeqd)
  rall <- sp::spTransform(moveall, center = center)
  sp::proj4string(rall) # "+proj=aeqd +lat_0=44.99489997 +lon_0=-17.48575004 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  rallCRS <- sp::CRS(sp::proj4string(rall))
  class(rallCRS) # CRS
  write.csv(sp::proj4string(rall), paste0(savedir, "CRS.csv"), row.names = FALSE)
  saveRDS(rallCRS, file = paste0(savedir, "CRS.Rds"))
  rm(moveall)
  rm(alldata)
  
  
  # We now need to reproject this into AEQD. Make a dummy object to get the correct projection.
  x.i <- raster::raster(
    xmn = min(data$NewEastingUTM),
    xmx = max(data$NewEastingUTM),
    ymn = min(data$NewNorthingUTM),
    ymx = max(data$NewNorthingUTM),
    crs = sp::proj4string(rall) # if rasterCRS then same as above! trying proj4string(r.i) which is AEQD from spTransform
    # could add here: resolution = rasterResolution # 50 ####
  )
  
  # before loop, r.i doesn't exist
  # in loop, r.i is aeqd'd from move.i which creates different lon & lat centres each time, hence the rasters don't stack
  # need to create aeqd from data####
  
  # +proj=aeqd +lat_0=24.75136 +lon_0=-36.01593 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  sp::proj4string(x.i) # "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  x.i
  # class      : RasterLayer 
  # dimensions : 180, 360, 64800  (nrow, ncol, ncell)
  # resolution : 28564.46, 37331.58  (x, y)
  # extent     : -7976391, 2306816, 2493709, 9213392  (xmin, xmax, ymin, ymax)
  # crs        : +proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
  
  
  # Use that to reproject UTM to AEQD
  # ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), i.e. projectExtent(from, to, ...))
  newTemplate <- raster::projectExtent(xUTM, sp::proj4string(x.i))
  # change newTemplate to xAEQD####
  sp::proj4string(newTemplate) # "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster.
  # class      : RasterLayer 
  # dimensions : 97, 148, 14356  (nrow, ncol, ncell)
  # resolution : 96002.62, 115886.2  (x, y)
  # extent     : -7481255, 6727132, -3943493, 7297469  (xmin, xmax, ymin, ymax)
  # crs        : +proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
  
  # This changes xUTM res from 50x50 to newTemplate res 43.3 x 52.2 (x, y)####
  # Which breaks writeRaster later
  rm(x.i)
  
  # change xAEQD res so x & y match. Kills values
  # raster::res(newTemplate) <- rep(mean(raster::res(newTemplate)), 2)
  # Use user res instead of average of the two####
  raster::res(newTemplate) <- rep(rasterResolution, 2)
  
  # Give newTemplate some values. Make Rep equal to the ncell dimension
  ones <- rep(1, raster::ncell(newTemplate))
  xAEQD <- raster::setValues(newTemplate, ones)
  rm(newTemplate)
  rm(ones)
  xAEQD
  # class      : RasterLayer 
  # dimensions : 105, 135, 14175  (nrow, ncol, ncell)
  # resolution : 102405.1, 102405.1  (x, y)
  # extent     : -8727700, 5096989, -4485671, 6266865  (xmin, xmax, ymin, ymax)
  # crs        : +proj=aeqd +lat_0=44.99489997 +lon_0=-17.48575004 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
  # source     : memory
  # names      : layer 
  # values     : 1, 1  (min, max)
  
  # get resolution from raster in rasterlist, assign it object, squared
  rasterres <- (raster::res(xAEQD)[1]) ^ 2
  
  # Loop through all unique tags
  counter <- 0
  for (i in unique(data$ID)) { # i <- unique(data$ID)[1]
    
    # Print individual
    counter <- counter + 1
    print(paste0(
      "processing ", which(unique(data$ID) %in% i),
      " of ", length(unique(data$ID))
    ))
    
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
    sp::proj4string(move.i) # "+proj=longlat +datum=WGS84 +no_defs"
    
    # Incorporate uncertainty in the model by including a location error.
    # From communication with Rob, the gps location of a shark is estimated, from:
    # the boat coordinates, bearing and distance estimate, to be 1 m.
    
    if ("moveLocError" %in% colnames(data.i)) {
      move.i$LocationError <- data.i$moveLocError
    } else {
      if (exists("moveLocError")) {
        if (length(moveLocError) == 1) {
          move.i$LocationError <- moveLocError # single value, replicated down for each relocation
          # Robs: 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate
        } else { # else if multiple values
          if (length(moveLocError) == nrow(data.i)) { # should be same length as full dataset
            move.i$LocationError <- data.i %>% # take the full dataset,
              dplyr::bind_cols(moveLocError = moveLocError) %>% # cbind the full movelocerror
              dplyr::filter(ID == i) %>% # filter for just this ID
              dplyr::pull(moveLocError) # and pull just the movelocerror for this ID
          } else { # if not length 1 and not length of nrow(data)
            stop(print("moveLocError must be either length 1 or length(nrow(data))")) # if not stop and tell user
          } # close not length1 not same length as full dataset else
        } # close length1 or else
      } # close if exists movelocerror
    }
    
    # Convert projection to Azimuthal Equi-Distance projection (aeqd)
    # r.i <- spTransform(move.i, center = center)
    r.i <- sp::spTransform(move.i, CRSobj = rallCRS, center = center)
    
    
    # Make sure it changed correctly
    sp::proj4string(r.i) # "+proj=aeqd +lat_0=24.75136 +lon_0=-36.01593 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    # irish tuna: "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    # r.tmp <- projectExtent(r.i, proj4string(rall)) # kills the data in r.i
    # proj4string(r.tmp)
    
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
    
    
    # We need to make sure that the projection of the move object is in the same format as our raster. Check below.
    sp::proj4string(xAEQD) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    sp::proj4string(r.i) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    sp::proj4string(xAEQD) == sp::proj4string(r.i) # TRUE
    
    # Below we exclude the data points that are 'long' i.e. create a large time gap. move::burst ?
    bursted <- move::burst(
      r.i,
      c("normal", "long")[1 + (move::timeLag(r.i, units = timeDiffUnits) > timeDiffLong)]
    )
    rm(r.i)
    # There are 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command.
    
    if ("dbblocationerror" %in% colnames(data.i)) {
      dbblocationerror.i <- data.i$dbblocationerror
    } else {
      if (exists("dbblocationerror")) {
        if (length(dbblocationerror) == 1) {
          dbblocationerror.i <- dbblocationerror # single value, replicated down for each relocation
          # Robs: 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate
        } else { # else if multiple values
          if (length(dbblocationerror) == nrow(data)) { # should be same length as full dataset
            dbblocationerror.i <- data %>% # take the full dataset,
              dplyr::bind_cols(dbblocationerror = dbblocationerror) %>% # cbind the full movelocerror
              dplyr::filter(ID == i) %>% # filter for just this ID
              dplyr::pull(dbblocationerror) # and pull just the movelocerror for this ID
          } else { # if not length 1 and not length of nrow(data)
            stop(print("dbblocationerror must be either length 1 or length(nrow(data))")) # if not stop and tell user
          } # close not length1 not same length as full dataset else
        } # close length1 or else
      } # close if exists dbblocationerror
    }
    
    # location error needs to be a positive number. Replace zeroes with 0.00001
    dbblocationerror.i[which(dbblocationerror.i == 0)] <- 0.00001
    
    # Construct the model. The time.step should reflect the ping frequency of the tag (in minutes)
    bursted_dbbmm <- move::brownian.bridge.dyn(bursted,
                                               burstType = "normal",
                                               raster = xAEQD, # has to be AEQD for metre-based calculations, presuambly?
                                               # Error:The projection of the raster and the Move object are not equal.
                                               # Need bursted, r.i, to be the same projection as xAEQD
                                               location.error = dbblocationerror.i, # dbblocationerror.i
                                               ext = dbbext, # dbbext
                                               window.size = dbbwindowsize #  must be >=2*margin which is 11 so >=22, but odd so >=23
    )
    
    # data.i$NewEastingUTMmin <- data.i$NewEastingUTM - data.i$dbblocationerror
    # data.i$NewEastingUTMmax <- data.i$NewEastingUTM + data.i$dbblocationerror
    # data.i$NewNorthingUTMmin <- data.i$NewNorthingUTM - data.i$dbblocationerror
    # data.i$NewNorthingUTMmax <- data.i$NewNorthingUTM + data.i$dbblocationerror
    # 
    # data.i$EastRastExtentMin <- extent(xAEQD)[1]
    # data.i$EastRastExtentMax <- extent(xAEQD)[2]
    # data.i$NorthRastExtentMin <- extent(xAEQD)[3]
    # data.i$NorthRastExtentMax <- extent(xAEQD)[4]
    # 
    # TODO ext size error ####
    # Error in .local(object, raster, location.error = location.error, ext = ext,  : 
    # Higher x grid not large enough, consider extending the raster in that direction or enlarging the ext argument
    # Season==Winter Fish#6
    
    # > xAEQD
    # class      : RasterLayer
    # dimensions : 342, 402, 137484  (nrow, ncol, ncell)
    # resolution : 25705.69, 25705.69  (x, y)
    # extent     : -7587818, 2745872, -3540761, 5250587  (xmin, xmax, ymin, ymax)
    # crs        : +proj=aeqd +lat_0=42.43218337 +lon_0=-30.08680803 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
    # source     : memory
    # names      : layer
    # values     : 1, 1  (min, max)
    # 
    # > bursted
    # class       : MoveBurst
    # features    : 90
    # extent      : 1286895, 1910646, -1295349, 317435  (xmin, xmax, ymin, ymax)
    # crs         : +proj=aeqd +lat_0=42.43218337 +lon_0=-30.08680803 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
    # variables   : 9
    # names       :         Lat,          Lon,   Datetime,    MarineZone,     moveLocError, dbblocationerror,     NewEastingUTM,  NewNorthingUTM, TimeDiff
    # min values  : 29.16504378, -15.05470771, 1480550400, Bay_of_Biscay, 15963.1149120894, 15963.1149120894, -1675882.39631877, 3396669.2275379,        0
    # max values  : 44.09076622, -9.703392913, 1488240000,         Other, 317635.883473968, 317635.883473968, -1080176.75804222, 5479499.2258687,        1
    # timestamps  : 2016-11-30 16:00:00 ... 2017-02-27 16:00:00 Time difference of 89 days  (start ... end, duration)
    # sensors     : unknown
    # indiv. data : ID, Season
    # indiv. value: X5116032 Winter
    # bursts      : normal: 89
    # date created: 2022-01-03 16:40:03
    
    # But runs with location.error = 1 & ext = 0, even though some NewNorthingUTMs are outside xAEQD extents (xmin & ymax)
    # plot(x = data.i$NewEastingUTM,
    #      y = data.i$NewNorthingUTM)
    
    # plot(x = c(extent(xAEQD)[1], extent(xAEQD)[1], extent(xAEQD)[2], extent(xAEQD)[2]),
    #        y = c(extent(xAEQD)[4], extent(xAEQD)[3], extent(xAEQD)[4], extent(xAEQD)[3]),
    #        pch = 19, col = "red")
    # points(bursted@coords)
    # points(x = c(extent(bursted)[1], extent(bursted)[1], extent(bursted)[2], extent(bursted)[2]),
    #      y = c(extent(bursted)[4], extent(bursted)[3], extent(bursted)[4], extent(bursted)[3]),
    #      col = "green")
    # points(x = c(extent(bursted)[1] - max(data.i$dbblocationerror),
    #              extent(bursted)[1] - max(data.i$dbblocationerror),
    #              extent(bursted)[2] + max(data.i$dbblocationerror),
    #              extent(bursted)[2] + max(data.i$dbblocationerror)),
    #        y = c(extent(bursted)[4] + max(data.i$dbblocationerror),
    #              extent(bursted)[3] - max(data.i$dbblocationerror),
    #              extent(bursted)[4] + max(data.i$dbblocationerror),
    #              extent(bursted)[3] - max(data.i$dbblocationerror)),
    #        col = "blue")
    # # Bursted + max location error still isn't outside xAEQD extents
    
    # calculate the dynamic brownian motion variance of the gappy track
    # dbbv <- brownian.motion.variance.dyn(bursted,
    #                                      location.error=dbblocationerror.i,
    #                                      window.size=dbbwindowsize,
    #                                      margin=11
    #                                      )
    # 
    # bursted$Datetime
    # tl <- timeLag(bursted, units = "hours")
    # # posted online https://gitlab.com/bartk/move/-/issues/58 ####
    
    rm(bursted)
    
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
    # Note: Based on 'An introduciton to the 'move' package', it's the opposite of what you might expect: "A cell with a very high value in the UD raster will have a very low value in the contour raster"
    # area.50 <- sum(raster::values(move::getVolumeUD(bb) <= .50))
    # area.95 <- sum(raster::values(move::getVolumeUD(bb) <= .95))
    area.50 <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.5)), 4)
    area.95 <- round(sum(raster::values(bb) >= (max(bb@data@values) * 0.05)), 4)
    
    # Combine in single df
    area.ct <- data.frame(
      core.use = area.50,
      general.use = area.95
    )
    
    # Add ID id
    area.ct$ID <- i
    
    # Put in list
    bb.list[[counter]] <- area.ct
    bb.list
    
    # need to reproject bb back to UTM else the resolution is unequal in x & y
    # only an issue for ascii raster format
    
    # proj4string(bb) # "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    # bbUTM <- projectExtent(bb, xUTM)
    # proj4string(bbUTM) # "+proj=utm +zone=27 +datum=WGS84 +units=m +no_defs"
    # bbUTM
    # # class      : RasterLayer 
    # # dimensions : 97, 148, 14356  (nrow, ncol, ncell)
    # # resolution : 158649.6, 222066.9  (x, y)
    # # extent     : -16099541, 7380606, -5468203, 16072285  (xmin, xmax, ymin, ymax)
    # # crs        : +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs 
    # res(bb) # 96002.62 115886.21
    # res(xUTM) # 111000 111000
    # res(bbUTM) # 158649.6 222066.9
    # res(bbUTM) <- res(xUTM)
    
    names(bb) <- i # puts name in slot, means scaleraster can name its output properly with lapply
    bb
    # class      : .UD 
    # dimensions : 105, 135, 14175  (nrow, ncol, ncell)
    # resolution : 102405.1, 102405.1  (x, y)
    # extent     : -8727700, 5096989, -4485671, 6266865  (xmin, xmax, ymin, ymax)
    # crs        : +proj=aeqd +lat_0=44.99489997 +lon_0=-17.48575004 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
    # source     : memory
    # names      : X5103448 
    # values     : 0, 0.03363462  (min, max)
    
    # need to have all of these bb's in the same res, extent, dimensions.
    # plot(bb)
    # options(scipen = 5)
    
    # Export aggregated individual UDs to as ascii files for further exploration/spatial analysis in GIS software.
    # Needed for when the aim is to plot population-level and normalized UDs per species.
    raster::writeRaster(bb, # x has unequal horizontal and vertical resolutions. Such data cannot be stored in arc-ascii format
                        paste0(savedir, filename = i, writeRasterExtension),
                        format = writeRasterFormat, # "ascii"
                        datatype = writeRasterDatatype, # "FLT4S"
                        if (writeRasterFormat != "CDF") bylayer = TRUE, # bylayer kills ncdf4
                        overwrite = TRUE)
    # Error in .startAsciiWriting(x, filename, ...):x has unequal horizontal and vertical resolutions. Such data cannot be stored in arc-ascii format
    # writeRasterFormat = "ascii"
    # writeRasterExtension = ".asc"
    # writeRasterFormat = "raster"
    # writeRasterExtension = ".grd"
    # library(ncdf4)
    # writeRasterFormat = "CDF"
    # writeRasterExtension = ".nc"
    
    gc()
  } # close for i in unique data$ID
  
  # Put everything in a data.frame
  md <- dplyr::bind_rows(bb.list,
                         .id = "column_label"
  ) %>%
    dplyr::select(!column_label) # remove column_label column
  
  md$core.use <- (rasterres * md$core.use) / 1000000 # convert from cells/pixels to kilometres squared area
  md$general.use <- (rasterres * md$general.use) / 1000000

  write.csv(md,
            file = file.path(savedir, absVolumeAreaSaveName),
            row.names = FALSE
  )
  
  # # Scale and sum individual rasters, followed by a rescaling to obtain an aggregated UD (.asc)
  # print(Sys.info()["nodename"])
  # if (grepl("aurits", Sys.info()["nodename"])) {
  #   work.dir <- "~/Documents/Science/Projects/Rob Bullock - Bimini/dBBMM_HomeRange_output/"
  # }
  # if (grepl("nautilus", Sys.info()["nodename"]) | grepl("Poseidon", Sys.info()["nodename"]) | grepl("aquarius", Sys.info()["nodename"])) {
  #   work.dir <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/"
  # }
  
} # close function
