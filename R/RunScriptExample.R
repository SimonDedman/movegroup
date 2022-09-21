# dBBMM ####
# by month? season? region?

saveloc <- "/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/"

# moveLocError should be a vector of metres of error for each point. Need to create that from each polygon.
reproject <- function(x, coorda, coordb, latloncrs, projectedcrs) {
  x <- sf::st_as_sf(x, coords = c(coorda, coordb)) |>
    sf::st_set_crs(latloncrs) |> # latlon degrees sf object
    st_transform(projectedcrs) |> # eastings northings units metres
    dplyr::select(-everything()) # remove all columns. Geometry is protected and retained
  return(x)
}
centre <- reproject(x = AllDailies,
                    coorda = "lon",
                    coordb = "lat",
                    latloncrs = 4326,
                    projectedcrs = 3857)
right <- reproject(x = AllDailies,
                   coorda = "lon975",
                   coordb = "lat",
                   latloncrs = 4326,
                   projectedcrs = 3857)
down <- reproject(x = AllDailies,
                  coorda = "lon",
                  coordb = "lat025",
                  latloncrs = 4326,
                  projectedcrs = 3857)
left <- reproject(x = AllDailies,
                  coorda = "lon025",
                  coordb = "lat",
                  latloncrs = 4326,
                  projectedcrs = 3857)
up <- reproject(x = AllDailies,
                coorda = "lon",
                coordb = "lat975",
                latloncrs = 4326,
                projectedcrs = 3857)
rightdist <- st_distance(x = centre,
                         y = right,
                         by_element = TRUE)
downdist <- st_distance(x = centre,
                        y = down,
                        by_element = TRUE)
leftdist <- st_distance(x = centre,
                        y = left,
                        by_element = TRUE)
updist <- st_distance(x = centre,
                      y = up,
                      by_element = TRUE)
meandist <- cbind(rightdist, downdist, leftdist, updist)
AllDailies$meandist <- rowMeans(meandist, na.rm = TRUE)

# Datetime must be POSIXct
AllDailies$Datetime <- as.POSIXct.Date(AllDailies$Date)

# set rasterCRS
mean(AllDailies$lon, na.rm = TRUE) # -19.43032
mean(AllDailies$lat, na.rm = TRUE) # 44.2838
# https://tmackinnon.com/2005/images/utmworld.gif
# Mean point is zone 27
# rasterCRS = CRS("+proj=utm +zone=27 +datum=WGS84")

# bbdwindowsize
# A dBBMM is not run if total detections of individual < window size
AllDailies %>%
  dplyr::select(toppid, Date) %>%
  group_by(toppid) %>%
  summarise(n = n()) %>%
  arrange(n) %>% # min 29
  summarise(windowsize = first(n)) %>%
  pull()
# 29


source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/dBBMM.build.R")
source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/scaleraster.R")
source("~/Dropbox/Galway/Analysis/R/dBBMM_HomeRange/R/dBBMM.plot.R")

dBBMM_HomeRange(
  data = AllDailies %>% dplyr::select(toppid, lat, lon, Datetime, MarineZone), # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
  ID = "toppid", # column name of IDs of individuals.
  Datetime = "Datetime", # name of Datetime column. Must be in POSIXct format.
  Lat = "lat", # name of Lat & Lon columns in data.
  Lon = "lon",
  Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  dat.TZ = "US/Eastern", # timezone for as.POSIXct.
  proj = CRS("+proj=longlat +datum=WGS84"),             # CRS for move function.
  projectedCRS = "+init=epsg:3857", # EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
  sensor = "unknown", # sensor for move function.
  moveLocError = AllDailies$meandist, # location error in metres for move function.
  timeDiffLong = 20, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  timeDiffUnits = "days", # units for time difference for move function. 20 days = TOPP Nature paper
  center = TRUE,
  buffpct = 0.3, # buffer extent for raster creation, proportion of 1.
  rasterCRS = CRS("+proj=utm +zone=27 +datum=WGS84"),   # CRS for raster creation.
  rasterResolution = 111000 / 4, # 1 degree lat = 111km. Could maybe make a tad smaller.
  bbdlocationerror = AllDailies$meandist, # location.error param in brownian.bridge.dyn.
  bbdext = 0.3, # ext param in brownian.bridge.dyn.
  bbdwindowsize = 29, # window.size param in brownian.bridge.dyn.
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = "/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/",  # save outputs to a temporary directory (default) else.
  alerts = TRUE
)

scaleraster(path = "/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/",
            crsloc = "/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/")

rallCRS <- read.csv("/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/CRS.csv")
rallCRS <- CRS(rallCRS$x)

dBBMM_plot(
  x = paste0("/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/", "Scaled/All_Rasters_Scaled_LatLon.asc"), # path to scaled data
  myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
  googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion
  expandfactor = 1.1, # extents expansion factor for Google Map basemap.
  mapzoom = 5, # 3 (continent) - 21 (building)
  mapsource = "stamen", # Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present.
  maptype = "terrain", # Type of map for ggmap::get_map.
  contour1colour = "red", # colour for contour 1, typically 95%.
  contour2colour = "orange", # colour for contour 2, typically 50%.
  plottitle = "Aggregated 95% and 50% UD contours",
  plotsubtitle = "Scaled contours. Irish-tagged Atlantic bluefin tuna. n = 51", # data %>% distinct(ID) %>% nrow() # 13
  legendtitle = "Percent UD Contours",
  plotcaption = paste0("dBBMM_HomeRange, ", today()),
  axisxlabel = "Longitude",
  axisylabel = "Latitude",
  legend.position = c(0.11, 0.9), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
  filesavename = paste0(today(), "_dBBMM-contours.png"),
  savedir = "/home/simon/Dropbox/Blocklab Monterey/Data/IrishTags/dBBMM/Plots" # file.path(work.dir, out.dir, "Scaled")
)