# rob_lemons_run_script
# 2021-12-22 Simon Dedman & Mo VZB

library(tidyverse)
library(magrittr)
library(lubridate)
library(tidylog)
source("./R/dBBMM.build.R")
source("./R/scaleraster.R")
source("./R/dBBMM.plot.R")
saveloc <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/" # si
DET <- read.csv("./data/TRACKS1.csv")

DET %<>%
  mutate(Datetime = as.POSIXct(Datetime,
                               format = "%m/%d/%y %H:%M",
                               tz = "US/Eastern"
  ),
  Shark = make.names(Shark) # prefixes "X" to numerical-named sharks to avoid issues later
  ) %>%
  rename(
    Lat = N,
    Lon = W,
    T.Ph = Tidal.Phase
  ) %>%
  select(Datetime, Shark, T.Ph, Lat, Lon) %>% # dropped 2 variables (Date, Time)
  arrange(Shark, Datetime) # df size: 1308 x 5

write.csv(x = DET, file = paste0(saveloc, "TracksCleaned.csv"), row.names = FALSE) # Could load this directly here

dBBMM_HomeRange(
  data = DET, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
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
  rasterCRS = CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
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
  plotcaption = paste0("dBBMM_HomeRange, ", today()),
  axisxlabel = "Longitude",
  axisylabel = "Latitude",
  legend.position = c(0.16, 0.92), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
  filesavename = paste0(today(), "_dBBMM-contours.png"),
  savedir = paste0(saveloc, "dBBMM ASCII/Plot/") # file.path(work.dir, out.dir, "Scaled")
)  

myLocation # -125.6578 -731.6077  635.0810  789.8699
# locations are wrong. all the way back to the initial asc files.
# build L500 r.i <- spTransform(move.i, CRSobj = rallCRS, center = center)
# > move.i
# class       : Move 
# features    : 39 
# extent      : -79.25101, -79.2461, 25.7315, 25.73574  (xmin, xmax, ymin, ymax)
# crs         : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
# variables   : 7
# names       :   Datetime, T.Ph,      Lat,       Lon,     NewEastingUTM,   NewNorthingUTM, LocationError 
# min values  : 1349195880,    L,  25.7315, -79.25101, -8822182.07805263, 2965863.89366829,             1 
# max values  : 1349299080,    M, 25.73574,  -79.2461, -8821635.49935284, 2966387.85337344,             1 
# timestamps  : 2012-10-02 12:38:00 ... 2012-10-03 17:18:00 Time difference of 1 days  (start ... end, duration) 
# sensors     : VR2W 
# indiv. data : ID 
# indiv. value: C69A 
# date created: 2021-12-18 00:59:02

# r.i
# class       : Move 
# features    : 39 
# extent      : -20.06822, 472.6107, -526.2208, -56.49803  (xmin, xmax, ymin, ymax)
# crs         : +proj=aeqd +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +lon_0=-79.25081 +lat_0=25.73625 
# variables   : 8
# names       :   Datetime, T.Ph,      Lat,       Lon,     NewEastingUTM,   NewNorthingUTM, LocationError,         TimeDiff 
# min values  : 1349195880,    L,  25.7315, -79.25101, -8822182.07805263, 2965863.89366829,             1,                0 
# max values  : 1349299080,    M, 25.73574,  -79.2461, -8821635.49935284, 2966387.85337344,             1, 24.5833333333333 
# timestamps  : 2012-10-02 12:38:00 ... 2012-10-03 17:18:00 Time difference of 1 days  (start ... end, duration) 
# sensors     : VR2W 
# indiv. data : ID 
# indiv. value: C69A 
# date created: 2021-12-18 00:59:02 

# r.i extent is out of "bounds"
# can only be the CRSobj projection?

# Switched to original CRS options: 
projectedCRS = "+init=epsg:32617"
rasterCRS = CRS("+proj=utm +zone=17 +datum=WGS84") # CRS for raster creation. 17
# but results in:
myLocation # -16.45209 -731.02838  592.14997  790.47677
# presumably should be longlat, -180:180, -90:90 ?
# try
x <- read_stars(x) %>% st_set_crs(2958) # 4326
myLocation # -16.45209 -731.02838  592.14997  790.47677
# made no difference.



x = paste0(saveloc, "dBBMM ASCII/Scaled/All_Rasters_Scaled_LatLon.asc")
y <- read_stars(x) %>% st_set_crs(2958) # 4326
y
st_crs(y) # 4326
st_bbox(y)
z <- raster(x)
extent(z)
# extents the same, correct, -79 25, lat lon.
# CRS is wrong, so x <- read_stars(x) %>% st_set_crs(4326) causes problems?
# st_set_crs(4326) is wrong so causes problems??

dataCRS <- readRDS(paste0(crsloc, "CRS.Rds"))
crs(z) # NA
crs(z) <- dataCRS
crs(z) # 
# crs = CRS("+proj=longlat")
y %<>% st_crs("+proj=longlat")
y$wkt
class(dataCRS)