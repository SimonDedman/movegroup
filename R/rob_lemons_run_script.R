# rob_lemons_run_script
# 2021-12-22 Simon Dedman & Mo VZB

# library(tidyverse)
# library(magrittr)
library(lubridate)
library(sp)
library(dplyr)
# library(tidylog)
# library(remotes)
# remotes::install_github("SimonDedman/dBBMMhomeRange", force = TRUE)
# library(dBBMMhomeRange)
# source("./R/dBBMM.build.R")
# library(move)
# source("./R/scaleraster.R")
# source("./R/dBBMM.plot.R")

utils::globalVariables("where") # https://github.com/r-lib/tidyselect/issues/201
# https://stackoverflow.com/questions/40251801/how-to-use-utilsglobalvariables
# https://github.com/r-lib/tidyselect/issues/248


saveloc <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/" # si
# DET <- read.csv("./data/TRACKS1.csv")
# 
# DET %<>%
#   mutate(Datetime = as.POSIXct(Datetime,
#                                format = "%m/%d/%y %H:%M",
#                                tz = "US/Eastern"
#   ),
#   Shark = make.names(Shark) # prefixes "X" to numerical-named sharks to avoid issues later
#   ) %>%
#   rename(
#     Lat = N,
#     Lon = W,
#     T.Ph = Tidal.Phase
#   ) %>%
#   select(Datetime, Shark, T.Ph, Lat, Lon) %>% # dropped 2 variables (Date, Time)
#   arrange(Shark, Datetime) # df size: 1308 x 5
# 
# write.csv(x = DET, file = paste0(saveloc, "TracksCleaned.csv"), row.names = FALSE) # Could load this directly here

DET <- read.csv(paste0(saveloc, "TracksCleaned.csv"))

dBBMMhomeRange(
  data = DET, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
  ID = "Shark", # column name of IDs of individuals.
  Datetime = "Datetime", # name of Datetime column. Must be in POSIXct format.
  Lat = "Lat", # name of Lat & Lon columns in data.
  Lon = "Lon",
  Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  dat.TZ = "US/Eastern", # timezone for as.POSIXct.
  proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
  projectedCRS = "+init=epsg:32617", # 32617 EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
  sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
  moveLocError = 1, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
  timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  timeDiffUnits = "hours", # units for time difference for move function.
  center = TRUE, # center move object within extent? See spTransform.
  buffpct = 3, # buffer extent for raster creation, proportion of 1.
  rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
  rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
  bbdlocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
  bbdext = 0.3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
  bbdwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = paste0(saveloc, "dBBMM ASCII/"),  # save outputs to a temporary directory (default) else.
  alerts = TRUE # audio warning for failures
) 

scaleraster(path = paste0(saveloc, "dBBMM ASCII/"), # Location of files created by dBBMM.build. No terminal slash.
            pattern = ".asc",
            format = "ascii",
            datatype = "FLT4S",
            bylayer = TRUE,
            overwrite = TRUE,
            scalefolder = "Scaled",
            summedname = "All_Rasters_Summed",
            scaledname = "All_Rasters_Scaled",
            crsloc = paste0(saveloc, "dBBMM ASCII/"), # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
            returnObj = FALSE)

dBBMM_plot(
  x = paste0(saveloc, "dBBMM ASCII/Scaled/All_Rasters_Scaled_LatLon.asc"), # path to scaled data
  # dataCRS = 2958, # one of (i) character: a string accepted by GDAL, (ii) integer, a valid EPSG value (numeric), or (iii) an object of class crs.
  myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
  # Default NULL, extents autocreated from data.
  # c(-79.3, 25.68331, -79.24, 25.78)
  googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion # FALSE
  # factors since Google Maps tiling zoom setup doesn't align to myLocation
  # extents.
  gmapsAPI = NULL, # enter your google maps API here, quoted character string
  expandfactor = 1.6, # extents expansion factor for basemap.
  # 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
  # 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
  mapzoom = 13, # 3 (continent) - 21 (building)
  mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present. # stamen
  maptype = "satellite", # Type of map for ggmap::get_map. # terrain
  contour1colour = "red", # colour for contour 1, typically 95%.
  contour2colour = "orange", # colour for contour 2, typically 50%.
  plottitle = "Aggregated 95% and 50% UD contours",
  # Can use the term 'home range' when an animal can be detected wherever it goes
  # i.e. using GPS, satellite or acoustic telemetry whereby it is known that acoustic
  # receivers cover the entire home range of the study species. 
  # This term is problematic when applied to a passive acoustic telemetry setting
  # where an array of non-overlapping receivers are used to assess local space use patterns
  # i.e. the home range is bigger than the coverage by the acoustic array; put in Details
  plotsubtitle = "Scaled contours. n = 13", # data %>% distinct(ID) %>% nrow() # 13
  legendtitle = "Percent UD Contours",
  plotcaption = paste0("dBBMMhomeRange, ", lubridate::today()),
  axisxlabel = "Longitude",
  axisylabel = "Latitude",
  legend.position = c(0.16, 0.92), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
  filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
  savedir = paste0(saveloc, "dBBMM ASCII/Plot/") # file.path(work.dir, out.dir, "Scaled")
)

# Tidal Cycle Loop ####
for (thistide in unique(DET$T.Ph)) {
  
  saveloc <- paste0("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/", thistide, "/") # si
  
  dBBMMhomeRange(
    data = DET %>% dplyr::filter(T.Ph == thistide), # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
    ID = "Shark", # column name of IDs of individuals.
    Datetime = "Datetime", # name of Datetime column. Must be in POSIXct format.
    Lat = "Lat", # name of Lat & Lon columns in data.
    Lon = "Lon",
    Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
    dat.TZ = "US/Eastern", # timezone for as.POSIXct.
    proj = CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
    projectedCRS = "+init=epsg:32617", # 32617 EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
    sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
    moveLocError = 1, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
    timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
    timeDiffUnits = "hours", # units for time difference for move function.
    center = TRUE, # center move object within extent? See spTransform.
    buffpct = 3, # buffer extent for raster creation, proportion of 1.
    rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
    rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
    bbdlocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
    bbdext = 0.3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
    bbdwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
    writeRasterFormat = "ascii",
    writeRasterExtension = ".asc",
    writeRasterDatatype = "FLT4S",
    absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
    savedir = paste0(saveloc, "dBBMM ASCII/"),  # save outputs to a temporary directory (default) else.
    alerts = TRUE # audio warning for failures
  ) 
  
  scaleraster(path = paste0(saveloc, "dBBMM ASCII/"), # Location of files created by dBBMM.build. No terminal slash.
              pattern = ".asc",
              format = "ascii",
              datatype = "FLT4S",
              bylayer = TRUE,
              overwrite = TRUE,
              scalefolder = "Scaled",
              summedname = "All_Rasters_Summed",
              scaledname = "All_Rasters_Scaled",
              crsloc = paste0(saveloc, "dBBMM ASCII/"), # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
              returnObj = FALSE)
  
  dBBMM_plot(
    x = paste0(saveloc, "dBBMM ASCII/Scaled/All_Rasters_Scaled_LatLon.asc"), # path to scaled data
    myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
    googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion # FALSE
    gmapsAPI = NULL, # enter your google maps API here, quoted character string
    expandfactor = 1.6, # extents expansion factor for basemap.
    mapzoom = 13, # 3 (continent) - 21 (building)
    mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present. # stamen
    maptype = "satellite", # Type of map for ggmap::get_map. # terrain
    contour1colour = "red", # colour for contour 1, typically 95%.
    contour2colour = "orange", # colour for contour 2, typically 50%.
    plottitle = "Aggregated 95% and 50% UD contours",
    plotsubtitle = "Scaled contours. n = 13", # data %>% distinct(ID) %>% nrow() # 13
    legendtitle = "Percent UD Contours",
    plotcaption = paste0("dBBMMhomeRange, ", lubridate::today()),
    axisxlabel = "Longitude",
    axisylabel = "Latitude",
    legend.position = c(0.16, 0.92), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
    filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
    savedir = paste0(saveloc, "dBBMM ASCII/Plot/") # file.path(work.dir, out.dir, "Scaled")
  )
} # close thistide loop