#' Automated Boosted Regression Tree modelling and mapping suite
#'
#' Automates delta log normal boosted regression trees abundance prediction.
#' Loops through all permutations of parameters provided (learning
#' rate, tree complexity, bag fraction), chooses the best, then simplifies it.
#' Generates line, dot and bar plots, and outputs these and the predictions
#' and a report of all variables used, statistics for tests, variable
#' interactions, predictors used and dropped, etc. If selected, generates
#' predicted abundance maps, and Unrepresentativeness surfaces.
#' See www.GitHub.com/SimonDedman/gbm.auto for issues, feedback, and development
#' suggestions. See SimonDedman.com for links to walkthrough paper, and papers
#' and thesis published using this package.
#'
#' @param grids Explantory data to predict to. Import with (e.g.) read.csv and
#' specify object name. Defaults to NULL (won't predict to grids).
#' 
#' @return Line, dot and bar plots, a report of all variables used, statistics
#' for tests, variable interactions, predictors used and dropped, etc. If
#' selected generates predicted abundance maps, and Unrepresentativeness surface
#'
#' @details Errors and their origins:
#' @examples
#' \donttest{
#' # Not run. Note: grids file was heavily cropped for CRAN upload so output map
#' # predictions only cover patchy chunks of the Irish Sea, not the whole area.
#' # Full versions of these files:
#' # https://drive.google.com/file/d/1WHYpftP3roozVKwi_R_IpW7tlZIhZA7r
#' # /view?usp=sharing
#' library(gbm.auto)
#' data(grids)
#' data(samples)
#' # Set your working directory
#' gbm.auto(grids = grids, samples = samples, expvar = c(4:8, 10), resvar = 11,
#' tc = c(2,7), lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)}
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export

#' @import dplyr
#' @import devtools
# install_git('https://gitlab.com/bartk/move.git') #Installs 'move' development version
#' @import move
#' @import ggplot2
#' @import maptools
#' @import circular
#' @import ggmap
#' @import mapproj
#' @import knitr
#' @import magrittr %<>%
#' @import tidyr drop_na
#' @import tidylog verbose functions
#' @import raster
#' @import sf
#' @import sp SpatialPoints and maybe others
#' @import knitr
#' @import kableExtra LaTeX format tables via pipe structure
#' @import lubridate today() and other datetime operations as needed.
#' @import stars ggplot rasters in plot sections.
#' @import starsExtra ggplot rasters in plot sections; trim2
#' @import devtools
#install_github("SimonDedman/gbm.auto")
#' @import gbm.auto

# importFrom example
#' @importFrom stringi stri_split_fixed


dBBMM_HomeRange <- function(
    data = NULL, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
    ID = NULL, # column name of IDs of individuals.
    Datetime = NULL, # name of Datetime column. Must be in POSIXct format.
    Lat = NULL, # name of Lat & Lon columns in data.
    Lon = NULL,
    Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
    dat.TZ = "US/Eastern", # timezone for as.POSIXct.
    proj = CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
    projectedCRS = "+init=epsg:32617", # EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
    sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
    moveLocError = 1, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
    timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
    timeDiffUnits = "hours", # units for time difference for move function.
    center = TRUE, # center move object within extent? See spTransform.
    buffpct = 0.3, # buffer extent for raster creation, proportion of 1.
    rasterCRS = CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation.
    rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres?
    bbdlocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
    bbdext = 3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
    bbdwindowsize = 21, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
    writeRasterFormat = "ascii",
    writeRasterExtension = ".asc",
    writeRasterDatatype = "FLT4S",
    absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
    savedir = tempdir(),  # save outputs to a temporary directory (default) else.
    # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
    alerts = TRUE # audio warning for failures
) {
  # ToDo:
  # option to use motion variance as the dependent variable, not the UD
  # clean up all notes, package elements, authorship, dependencies etc.
  # Add examples
  
  # Generalised Boosting Model / Boosted Regression Tree process chain automater.
  # Simon Dedman, 2012-6 simondedman@gmail.com GitHub.com/SimonDedman/gbm.auto
  
  # source("scaleraster.R") # might not be needed if included as function in same package
  
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts) options(error = function() {
    beep(9)# give warning noise if it fails
    graphics.off()# kill all graphics devices
    setwd(oldwd) # reinstate original working directory. Probably redundant given on.exit
  } # close options subcurly
  ) # close options
  
  library(sp)
  library(move)
  if (writeRasterFormat == "CDF") library(ncdf4)
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
  rename(Datetime = .data[[Datetime]],
         ID = .data[[ID]],
         Lat = .data[[Lat]],
         Lon = .data[[Lon]]) %>%
    mutate(ID = make.names(ID))
  
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
  cord.dec <- SpatialPoints(cbind(data$Lon, data$Lat),
                            proj4string = CRS("+proj=longlat")
  )
  
  # Transform to UTM by setting the EPSG to 32617 for WGS 84, UTM zone 17, northern hemisphere. This is where Bimini is located.
  cord.UTM <- as.data.frame(spTransform(cord.dec, CRS(projectedCRS)))
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
    group_by(ID) %>%
    summarise(relocations = length(Datetime))
  check2 <- filter(check1, relocations >= bbdwindowsize) # filter: removed 2 rows (14%), 12 rows remaining
  
  if (length(check1$ID) != length(check2$ID)) {
    data <- semi_join(data, check2) # Joining, by = "ID". semi_join: added no columns
    check1 <- data %>%
      group_by(ID) %>%
      summarise(relocations = length(Datetime))
    check2 <- filter(check1, relocations >= bbdwindowsize) # filter: no rows removed
    length(check1$ID) == length(check2$ID)
  } # data: 1253 x 7
  # ToDo: improve this####
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
    group_by(ID) %>%
    distinct(Datetime, .keep_all = TRUE) %>% # distinct (grouped): removed one row (<1%), 1,286 rows remaining
    # prevents duplicate Datetime crash in move() later
    ungroup() # 1253 x 7 after removing as.numeric above
  
  
  
  
  
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
  xUTM <- raster(
    xmn = min(data$NewEastingUTM, na.rm = TRUE) - ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
    xmx = max(data$NewEastingUTM, na.rm = TRUE) + ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
    ymn = min(data$NewNorthingUTM, na.rm = TRUE) - ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
    ymx = max(data$NewNorthingUTM, na.rm = TRUE) + ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
    crs = rasterCRS, # +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs
    resolution = rasterResolution # 50
  ) 
  proj4string(xUTM) # "+proj=utm +zone=27 +datum=WGS84 +units=m +no_defs"
  xUTM
  # class      : RasterLayer 
  # dimensions : 97, 148, 14356  (nrow, ncol, ncell)
  # resolution : 111000, 111000  (x, y)
  # extent     : -11061354, 5366646, 462297.3, 11229297  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=27 +datum=WGS84 +units=m +no_defs 
  
  
  
  # run move on all data together to generate global CRS for raster
  alldata <- data %>%
    arrange(Datetime) %>%
    group_by(Datetime) %>% # remove duplicates
    summarise(across(where(is.numeric), mean, na.rm = TRUE),
              across(where(~ is.character(.) | is.POSIXt(.)), first)) %>% # library(lubridate)
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)), #convert NaN to NA. POSIX needs lubridate
           across(where(~ is.character(.)), ~ ifelse(is.nan(.), NA, .))) %>% # https://community.rstudio.com/t/why-does-tidyrs-fill-work-with-nas-but-not-nans/25506/5
    ungroup # Joining, by = c("toppid", "Date")
  
  alldata <- as.data.frame(alldata)
  moveall <- move(
    x = alldata$Lon,
    y = alldata$Lat,
    time = alldata$Datetime,
    proj = proj,
    data = alldata,
    # animal = data$ID,
    sensor = sensor
  )
  # Convert projection to Azimuthal Equi-Distance projection (aeqd)
  rall <- spTransform(moveall, center = center)
  proj4string(rall) # "+proj=aeqd +lat_0=44.99489997 +lon_0=-17.48575004 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  rallCRS <- CRS(proj4string(rall))
  class(rallCRS) # CRS
  write.csv(proj4string(rall), paste0(savedir, "CRS.csv"), row.names = FALSE)
  saveRDS(rallCRS, file = paste0(savedir, "CRS.Rds"))
  rm(moveall)
  rm(alldata)
  
  
  # We now need to reproject this into AEQD. Make a dummy object to get the correct projection.
  x.i <- raster(
    xmn = min(data$NewEastingUTM),
    xmx = max(data$NewEastingUTM),
    ymn = min(data$NewNorthingUTM),
    ymx = max(data$NewNorthingUTM),
    crs = proj4string(rall) # if rasterCRS then same as above! trying proj4string(r.i) which is AEQD from spTransform
    # could add here: resolution = rasterResolution # 50 ####
  )
  
  # before loop, r.i doesn't exist
  # in loop, r.i is aeqd'd from move.i which creates different lon & lat centres each time, hence the rasters don't stack
  # need to create aeqd from data####
  
  # +proj=aeqd +lat_0=24.75136 +lon_0=-36.01593 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  proj4string(x.i) # "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  x.i
  # class      : RasterLayer 
  # dimensions : 180, 360, 64800  (nrow, ncol, ncell)
  # resolution : 28564.46, 37331.58  (x, y)
  # extent     : -7976391, 2306816, 2493709, 9213392  (xmin, xmax, ymin, ymax)
  # crs        : +proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
  
  
  # Use that to reproject UTM to AEQD
  # ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), i.e. projectExtent(from, to, ...))
  newTemplate <- projectExtent(xUTM, proj4string(x.i))
  # change newTemplate to xAEQD####
  proj4string(newTemplate) # "+proj=aeqd +lat_0=39.53798259 +lon_0=-39.4894484505 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
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
  res(newTemplate) <- rep(mean(res(newTemplate)), 2)
  
  # Give newTemplate some values. Make Rep equal to the ncell dimension
  ones <- rep(1, ncell(newTemplate))
  xAEQD <- setValues(newTemplate, ones)
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
    move.i <- move(
      x = data.i$Lon,
      y = data.i$Lat,
      time = data.i$Datetime,
      proj = proj,
      data = data.i,
      animal = data.i$ID,
      sensor = sensor
    )
    
    # Check the current projection
    proj4string(move.i) # "+proj=longlat +datum=WGS84 +no_defs"
    
    # Incorporate uncertainty in the model by including a location error.
    # From communication with Rob, the gps location of a shark is estimated, from:
    # the boat coordinates, bearing and distance estimate, to be 1 m.
    if (exists("moveLocError")) {
      if (length(moveLocError) == 1) {
        move.i$LocationError <- moveLocError # single value, replicated down for each relocation
        # Robs: 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate
      } else { # else if multiple values
        if (length(moveLocError) == nrow(data)) { # should be same length as full dataset
          move.i$LocationError <- data %>% # take the full dataset,
            bind_cols(moveLocError = moveLocError) %>% # cbind the full movelocerror
            filter(ID == i) %>% # filter for just this ID
            pull(moveLocError) # and pull just the movelocerror for this ID
        } else { # if not length 1 and not length of nrow(data)
          stop(print("moveLocError must be either length 1 or length(nrow(data))")) # if not stop and tell user
        } # close not length1 not same length as full dataset else
      } # close length1 or else
    } # close if exists movelocerror
    
    # Convert projection to Azimuthal Equi-Distance projection (aeqd)
    # r.i <- spTransform(move.i, center = center)
    r.i <- spTransform(move.i, CRSobj = rallCRS, center = center)
    
    
    # Make sure it changed correctly
    proj4string(r.i) # "+proj=aeqd +lat_0=24.75136 +lon_0=-36.01593 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
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
    TimeDiff <- timeLag(r.i, units = timeDiffUnits)
    r.i$TimeDiff <- append(0, TimeDiff)
    long <- which(r.i$TimeDiff > timeDiffLong)
    # if no timediffs are longer than timedifflong, long is integer(0), a potentially dangerous object.
    
    
    # We need to make sure that the projection of the move object is in the same format as our raster. Check below.
    proj4string(xAEQD) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    proj4string(r.i) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    proj4string(xAEQD) == proj4string(r.i) # TRUE
    
    # Below we exclude the data points that are 'long' i.e. create a large time gap. move::burst ?
    bursted <- burst(
      r.i,
      c("normal", "long")[1 + (timeLag(r.i, units = timeDiffUnits) > timeDiffLong)]
    )
    rm(r.i)
    # There are 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command.
    
    if (exists("bbdlocationerror")) {
      if (length(bbdlocationerror) == 1) {
        bbdlocationerror.i <- bbdlocationerror # single value, replicated down for each relocation
        # Robs: 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate
      } else { # else if multiple values
        if (length(bbdlocationerror) == nrow(data)) { # should be same length as full dataset
          bbdlocationerror.i <- data %>% # take the full dataset,
            bind_cols(bbdlocationerror = bbdlocationerror) %>% # cbind the full movelocerror
            filter(ID == i) %>% # filter for just this ID
            pull(bbdlocationerror) # and pull just the movelocerror for this ID
        } else { # if not length 1 and not length of nrow(data)
          stop(print("bbdlocationerror must be either length 1 or length(nrow(data))")) # if not stop and tell user
        } # close not length1 not same length as full dataset else
      } # close length1 or else
    } # close if exists movelocerror
    
    # location error needs to be a postive number. Replace zeroes with 0.00001
    bbdlocationerror.i[which(bbdlocationerror.i == 0)] <- 0.00001
    
    # Construct the model. The time.step should reflect the ping frequency of the tag (in minutes)
    bursted_dbbmm <- brownian.bridge.dyn(bursted,
                                         burstType = "normal",
                                         raster = xAEQD, # has to be AEQD for metre-based calculations, presuambly?
                                         # Error:The projection of the raster and the Move object are not equal.
                                         # Need bursted, r.i, to be the same projection as xAEQD
                                         location.error = bbdlocationerror.i,
                                         ext = bbdext,
                                         window.size = bbdwindowsize
    )
    rm(bursted)
    
    # Re-standardize (Dr. Kranstauber's trouble shooting solution).
    # Occurring errors are "due to limits of accuracy during calculations.
    # Given the error is really small (<.000001 %) I would not worry to much and re-standardize".
    # Then aggregate UD segments.
    tmp <- calc(bursted_dbbmm, sum) # raster::calc :Calculate values for a new Raster* object from another Raster* object
    # using a formula. returns a RasterLayer if fun returns a single value (e.g. sum)
    rm(bursted_dbbmm)
    bb <- new(".UD", tmp / sum(values(tmp))) # new() creates object from Class ("UD.")
    rm(tmp)
    
    # Calculate volume area (m^2) within 50% (core) and 95% (general use) contours. Note: absolute scale
    area.50 <- sum(values(getVolumeUD(bb) <= .50))
    area.95 <- sum(values(getVolumeUD(bb) <= .95))
    
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
    writeRaster(bb, # x has unequal horizontal and vertical resolutions. Such data cannot be stored in arc-ascii format
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
  md <- bind_rows(bb.list,
                  .id = "column_label"
  ) %>%
    dplyr::select(!column_label) # remove column_label column
  
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