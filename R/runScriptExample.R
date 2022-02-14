source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/dBBMM.build.R")
dBBMM_HomeRange(
  data = NULL, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
  ID = NULL, # column name of IDs of individuals.
  Datetime = NULL, # name of Datetime column. Must be in POSIXct format.
  Lat = NULL, # name of Lat & Lon columns in data.
  Lon = NULL,
  Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  dat.TZ = "US/Eastern", # timezone for as.POSIXct.
  proj = CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
  projectedCRS = "+init=epsg:32617", # EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
  sensor = "VR2W", # sensor for move function.
  moveLocError = 1, # location error in metres for move function.
  timeDiffUnits = "hours", # units for time difference for move function.
  timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  buffpct = 0.3, # buffer extent for raster creation, proportion of 1.
  rasterCRS = CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation.
  rasterResolution = 50,
  bbdlocationerror = "LocationError", # location.error param in brownian.bridge.dyn.
  bbdext = 3, # ext param in brownian.bridge.dyn.
  bbdwindowsize = 23, # window.size param in brownian.bridge.dyn.
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = tempdir()  # save outputs to a temporary directory (default) else.
  # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
)

source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/scaleraster.R")
# Warning in file(filename, "r", encoding = encoding) :  cannot open file './R/scaleraster.R': No such file or directory
# Error in file(filename, "r", encoding = encoding) : cannot open the connection
# Again, works fine in console
scaleraster(path = file.path(data.dir, out.dir)) # should be savedir ?


source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/dBBMM.plot.R")
dBBMM_plot(
  x = paste0(data.dir, "/dBBMM ASCII/Scaled/All_Rasters_Scaled.asc"), # path to scaled data
  myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
  # Default NULL, extents autocreated from data.
  # c(-79.3, 25.68331, -79.24, 25.78)
  googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion
  # factors since Google Maps timing zoom setup doesn't align to myLocation
  # extents.
  expandfactor = 1.6, # extents expansion factor for Google Map basemap.
  # 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
  # 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
  mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present.
  maptype = "satellite", # Type of map for ggmap::get_map.
  contour1colour = "red", # colour for contour 1, typically 95%.
  contour2colour = "orange", # colour for contour 2, typically 50%.
  plottitle = "Aggregated 95% and 50% contours, lemon sharks, Bimini",
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
  savedir = tempdir() # file.path(work.dir, out.dir, "Scaled")
)
