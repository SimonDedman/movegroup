# rob_lemons_run_script
# 2021-12-22 Simon Dedman & Mo VZB

# test

library(magrittr)
library(lubridate)
library(sp)
library(dplyr)
library(tidylog)
library(remotes)
# remotes::install_github("SimonDedman/dBBMMhomeRange", force = TRUE)
library(dBBMMhomeRange)

# saveloc <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/" # si
saveloc <- "~/github/A_day_in_the_life_analysis/dBBMMhomeRange/" # Mo

# DET <- read.csv("./data/TRACKS1.csv")
# DET <- read.csv("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMMhomeRange/Misc/TRACKS1.csv")
# 
# DET %<>%
#   dplyr::mutate(Datetime = as.POSIXct(Datetime,
#                                format = "%m/%d/%y %H:%M",
#                                tz = "US/Eastern"
#   ),
#   Shark = make.names(Shark) # prefixes "X" to numerical-named sharks to avoid issues later
#   ) %>%
#   dplyr::rename(
#     Lat = N,
#     Lon = W,
#     T.Ph = Tidal.Phase
#   ) %>%
#   dplyr::select(Datetime, Shark, T.Ph, Lat, Lon) %>% # dropped 2 variables (Date, Time)
#   dplyr::arrange(Shark, Datetime) # df size: 1308 x 5
# 
# write.csv(x = DET, file = paste0(saveloc, "TracksCleaned.csv"), row.names = FALSE) # Could load this directly here
# saveRDS(object = DET,
#         file = "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMMhomeRange/Misc/TracksCleaned.Rds")
# 
# DET <- read.csv(paste0(saveloc, "TracksCleaned.csv"))
DET <- get(load(file = "~/github/A_day_in_the_life_analysis/dBBMMhomeRange/data/TracksCleaned.Rdata")) # Mo
# DET <- filter(DET, T.Ph == "M")

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
  rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
  # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
  rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
  rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
  dbblocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
  dbbext = 0.3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
  dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = paste0(saveloc, "dBBMM ASCII/"),  # save outputs to a temporary directory (default) else.
  alerts = TRUE # audio warning for failures
) 

scaleraster(path = paste0(saveloc, "dBBMM ASCII/"), # Location of files created by dBBMM.build. No terminal slash.
            pathsubsets = paste0(saveloc, "dBBMM ASCII/"),
            pattern = ".asc",
            weighting = 1, # weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
            format = "ascii",
            datatype = "FLT4S",
            bylayer = TRUE,
            overwrite = TRUE,
            scalefolder = "Scaled",
            weightedsummedname = "All_Rasters_Weighted_Summed",
            scaledweightedname = "All_Rasters_Scaled_Weighted",
            crsloc = paste0(saveloc, "dBBMM ASCII/"), # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
            returnObj = FALSE)

dBBMMplot(
  x = paste0(saveloc, "dBBMM ASCII/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), # path to scaled data
  crsloc = paste0(saveloc, "dBBMM ASCII/"),
  trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE
  myLocation = c(-79.27, 25.72, -79.24, 25.75), # location for extents, format c(xmin, ymin, xmax, ymax).
  # Default NULL, extents autocreated from data.
  # c(-79.3, 25.68331, -79.24, 25.78)
  googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion
  # factors since Google Maps tiling zoom setup doesn't align to myLocation
  # extents.
  gmapsAPI = NULL, # enter your Google maps API here, quoted character string
  expandfactor = 1.6, # extents expansion factor for basemap.
  # 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
  # 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
  mapzoom = 12, # google: 3 (continent) - 21 (building). stamen: 0-18
  mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present.
  maptype = "satellite", # Type of map for ggmap::get_map.
  contour1colour = "orange", # colour for contour 1, typically 95%.
  contour2colour = "red", # colour for contour 2, typically 50%.
  plottitle = "Aggregated 95% and 50% UD contours",
  # Can use the term 'home range' when an animal can be detected wherever it goes
  # i.e. using GPS, satellite or acoustic telemetry whereby it is known that acoustic
  # receivers cover the entire home range of the study species. 
  # This term is problematic when applied to a passive acoustic telemetry setting
  # where an array of non-overlapping receivers are used to assess local space use patterns
  # i.e. the home range is bigger than the coverage by the acoustic array; put in Details
  plotsubtitle = "Scaled contours. n = 13", # data %>% distinct(ID) %>% nrow() # 13
  legendtitle = "Percent UD Contours",
  plotcaption = paste0("dBBMM_HomeRange, ", lubridate::today()),
  axisxlabel = "Longitude",
  axisylabel = "Latitude",
  legendposition = c(0.16, 0.78), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
  filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
  savedir = paste0(saveloc, "dBBMM ASCII/Plot/"), # file.path(work.dir, out.dir, "Scaled")
  receiverlats = NULL, # vector of latitudes for receivers to be plotted
  receiverlons = NULL, # vector of longitudes for receivers to be plotted
  receivernames = NULL, # vector of names for receivers to be plotted
  receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
  recpointscol = NULL, # Colour of receiver centrepoint outlines.
  recpointsfill = NULL, # Colour of receiver centrepoint fills.
  recpointsalpha = NULL, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
  recpointssize = NULL, # Size of receiver points.
  recpointsshape = NULL, # Shape of receiver points, default 21, circle with outline and fill.
  recbufcol = NULL, # Colour of the receiver buffer circle outlines.
  recbuffill = NULL, # Colour of the receiver buffer circle fills.
  recbufalpha = NULL,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
  reclabcol = NULL, # Receiver label text colour.
  reclabfill = NA, # Receiver label fill colour, NA for no fill.
  reclabnudgex = NULL, # Receiver label offset nudge in X dimension.
  reclabnudgey = NULL, # Receiver label offset nudge in Y dimension.
  reclabpad = NULL, # Receiver label padding in lines.
  reclabrad = NULL, # Receiver label radius in lines.
  reclabbord = NULL, # Receiver label border in mm.
  surface = TRUE
)

# Tidal Cycle Loop ####

## MO COMMENT 2022-10-11: the below loop within a loop no longer works as intended. first dbbmmhomerange() needs to be run for all tidal phases, then scaleraster()

for (thistide in unique(DET$T.Ph)) { # thistide <- unique(DET$T.Ph)[1]
  
  # saveloc <- paste0("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/", thistide, "/") # si
  saveloc <- paste0("~/github/A_day_in_the_life_analysis/dBBMMhomeRange/dBBMM ASCII/Tides/", thistide, "/") # si
  
  dir.create(saveloc)
  dir.create(paste0(saveloc, "Plot"))
  
  dBBMMhomeRange(
    data = DET %>% dplyr::filter(T.Ph == thistide), # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
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
    rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
    # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
    rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
    rasterResolution = 50, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
    dbblocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
    dbbext = 0.3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
    dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
    writeRasterFormat = "ascii",
    writeRasterExtension = ".asc",
    writeRasterDatatype = "FLT4S",
    absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
    savedir = paste0(saveloc),  # save outputs to a temporary directory (default) else.
    alerts = TRUE # audio warning for failures
  ) 
  
  scaleraster(path = paste0(saveloc), # Location of files created by dBBMM.build. No terminal slash.
              pathsubsets = "~/github/A_day_in_the_life_analysis/dBBMMhomeRange/dBBMM ASCII/Tides/",
              pattern = ".asc",
              weighting = 1, # weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
              format = "ascii",
              datatype = "FLT4S",
              bylayer = TRUE,
              overwrite = TRUE,
              scalefolder = "Scaled",
              weightedsummedname = "All_Rasters_Weighted_Summed",
              scaledweightedname = "All_Rasters_Scaled_Weighted",
              crsloc = paste0(saveloc), # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
              returnObj = FALSE)
  
  dBBMMplot(
    x = paste0(saveloc, "Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), # path to scaled data
    crsloc = saveloc,
    trim = TRUE,
    myLocation = c(-79.27, 25.72, -79.24, 25.75),
    googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion # FALSE
    gmapsAPI = NULL, # enter your google maps API here, quoted character string
    expandfactor = 1.6, # extents expansion factor for basemap.
    mapzoom = 12, # 3 (continent) - 21 (building)
    mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. # stamen
    maptype = "satellite", # Type of map for ggmap::get_map. # terrain
    contour1colour = "orange", # colour for contour 1, typically 95%.
    contour2colour = "red", # colour for contour 2, typically 50%.
    plottitle = "Aggregated 95% and 50% UD contours",
    plotsubtitle = paste0("Scaled contours. n = ",
                          length(list.files(path = saveloc, pattern = ".asc")),
                          ", Tide = ",
                          switch(thistide,
                                 L = "Low",
                                 M = "Mid",
                                 H = "High")), # data %>% distinct(ID) %>% nrow() # 13
    legendtitle = "Percent UD Contours",
    plotcaption = paste0("dBBMMhomeRange, ", lubridate::today()),
    axisxlabel = "Longitude",
    axisylabel = "Latitude",
    legendposition = c(0.16, 0.78), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
    filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
    savedir = paste0(saveloc, "Plot"), # file.path(work.dir, out.dir, "Scaled")
    receiverlats = NULL, # vector of latitudes for receivers to be plotted
    receiverlons = NULL, # vector of longitudes for receivers to be plotted
    receivernames = NULL, # vector of names for receivers to be plotted
    receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
    recpointscol = "black", # Colour of receiver centrepoint outlines.
    recpointsfill = "white", # Colour of receiver centrepoint fills.
    recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
    recpointssize = 1, # Size of receiver points.
    recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
    recbufcol = "black", # Colour of the receiver buffer circle outlines.
    recbuffill = "red", # Colour of the receiver buffer circle fills.
    recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
    reclabcol = "black", # Receiver label text colour.
    reclabfill = NA, # Receiver label fill colour, NA for no fill.
    reclabnudgex = 0, # Receiver label offset nudge in X dimension.
    reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
    reclabpad = 0, # Receiver label padding in lines.
    reclabrad = 0.15, # Receiver label radius in lines.
    reclabbord = 0, # Receiver label border in mm.
    surface = TRUE
  )
} # close thistide loop

