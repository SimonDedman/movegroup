# SD problem solving scratchpad 2022-09-22

# load all variables & library quickly ####
library(ggplot2)
library(ggmap)
library(magrittr)
library(stars)
library(lubridate)
library(sf)
library(starsExtra)
library(fishtrack3d)

x = "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/All_Rasters_Scaled_Weighted_UDScaled.asc" # path to scaled data
# dataCRS = 2958 # one of (i) character: a string accepted by GDAL (ii) integer a valid EPSG value (numeric) or (iii) an object of class crs.
crsloc = NULL # Location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
trim = TRUE # remove NA & 0 values and crop to remaining date extents? Default TRUE
myLocation = NULL # location for extents format c(xmin ymin xmax ymax).
# Default NULL extents autocreated from data.
# c(-79.3 25.68331 -79.24 25.78)
googlemap = FALSE # If pulling basemap from Google maps this sets expansion
# factors since Google Maps tiling zoom setup doesn't align to myLocation
# extents.
gmapsAPI = NULL # enter your Google maps API here quoted character string
# expandfactor = 1.6 # extents expansion factor for basemap.
expandfactor = 0 # extents expansion factor for basemap.
# 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
# 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
mapzoom = 12 # google: 3 (continent) - 21 (building). stamen: 0-18
mapsource = "google" # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present.
maptype = "satellite" # Type of map for ggmap::get_map.
contour1colour = "red" # colour for contour 1 typically 95%.
contour2colour = "orange" # colour for contour 2 typically 50%.
plottitle = "Aggregated 95% and 50% UD contours"
# Can use the term 'home range' when an animal can be detected wherever it goes
# i.e. using GPS satellite or acoustic telemetry whereby it is known that acoustic
# receivers cover the entire home range of the study species. 
# This term is problematic when applied to a passive acoustic telemetry setting
# where an array of non-overlapping receivers are used to assess local space use patterns
# i.e. the home range is bigger than the coverage by the acoustic array; put in Details
plotsubtitle = "Scaled contours. n = 13" # data %>% distinct(ID) %>% nrow() # 13
legendtitle = "Percent UD Contours"
plotcaption = paste0("dBBMM_HomeRange ", lubridate::today())
axisxlabel = "Longitude"
axisylabel = "Latitude"
legendposition = c(0.11, 0.85) # Percent distance (of middle? of legend box) from L to R percent distance from Bottom to Top.
fontsize = 12
fontfamily = "Times New Roman"
filesavename = paste0(lubridate::today(), "_dBBMM-contours.png")
savedir = tempdir() # file.path(work.dir out.dir "Scaled")
receiverlats = NULL # vector of latitudes for receivers to be plotted
receiverlons = NULL # vector of longitudes for receivers to be plotted
receivernames = NULL # vector of names for receivers to be plotted
receiverrange = NULL # single (will be recycled) or vector of detection ranges in metres for receivers to be plotted
recpointscol = "black" # Colour of receiver centrepoint outlines.
recpointsfill = "white" # Colour of receiver centrepoint fills.
recpointsalpha = 0.5 # Alpha value of receiver centrepoint fills 0 (invisible) to 1 (fully visible).
recpointssize = 1 # Size of receiver points.
recpointsshape = 21 # Shape of receiver points default 21 circle with outline and fill.
recbufcol = "black" # Colour of the receiver buffer circle outlines.
recbuffill = "red" # Colour of the receiver buffer circle fills.
recbufalpha = 0.5  # Alpha value of receiver buffer fills 0 (invisible) to 1 (fully visible).
reclabcol = "black" # Receiver label text colour.
reclabfill = NA # Receiver label fill colour NA for no fill.
reclabnudgex = 0 # Receiver label offset nudge in X dimension.
reclabnudgey = -200 # Receiver label offset nudge in Y dimension.
reclabpad = 0 # Receiver label padding in lines.
reclabrad = 0.15 # Receiver label radius in lines.
reclabbord = 0 # Receiver label border in mm.
surface = TRUE # Plot complete UD surface as well as contours


# install fishtrack3d from github####
remotes::install_github("aspillaga/fishtrack3d")
library(fishtrack3d)


# dBBMM.plot: Stars needs CRS to plot in ggMap ####
library(stars)
xstars <- stars::read_stars("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/All_Rasters_Scaled_Weighted_UDScaled.asc")
plot(xstars, breaks = "equal") # works
dataCRS <- readRDS("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/CRS.Rds")
st_crs(xstars) # Coordinate Reference System: NA
class(dataCRS) # CRS, sp package
# ?crs : “For compatibility with sp you can use proj4string instead of crs”
proj4string(xstars) <- dataCRS # Error in (function (classes, fdef, mtable): unable to find an inherited method for function ‘proj4string<-’ for signature ‘"stars", "CRS"’
proj4string(dataCRS) # "+proj=aeqd +lat_0=25.6871 +lon_0=-79.29617 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
st_set_crs(xstars, proj4string(dataCRS))
# stars object with 2 dimensions and 1 attribute
# attribute(s):
#   Min. 1st Qu. Median         Mean 3rd Qu.        Max.
# All_Rasters_Scaled_Weighted_UD...     0       0      0 1.264063e-05       0 0.008182338
# dimension(s):
#   from  to   offset delta                       refsys point values x/y
# x    1 293 -29238.7   200 +proj=aeqd +lat_0=25.6871...    NA   NULL [x]
# y    1 270  27074.4  -200 +proj=aeqd +lat_0=25.6871...    NA   NULL [y]
st_crs(xstars) # Coordinate Reference System: NA
crs(xstars) # NA
# st_set_crs did nothing?
st_crs(xstars) <- dataCRS # Error in `st_crs<-.dimensions`(`*tmp*`, value = value): crs of class CRS not recognized
st_crs(xstars) <- proj4string(dataCRS) # no error
st_crs(xstars) #  above line works!

x <- stars::read_stars("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/All_Rasters_Scaled_Weighted_UDScaled.asc")

ggplot() + stars::geom_stars(data = surfaceUD %>% sf::st_transform(3857)) # works, centred on 0,0, in latlon
ggplot() + stars::geom_stars(data = surfaceUD) # works, centred on 0,0, stretched long/high, in metres
ggplot() + stars::geom_stars(data = x %>% sf::st_transform(3857)) # works, correct latlon position
ggplot() + stars::geom_stars(data = x) # works, centred on 0,0, stretched long/high, in metres
# but
ggmap::ggmap(myMap) + stars::geom_stars(data = x %>% sf::st_transform(3857), inherit.aes = FALSE)
# Coordinate system already present. Adding new coordinate system, which will replace the existing one.
# Error in FUN(X[[i]], ...) : object 'lon' not found
st_crs(myMap) # Coordinate Reference System: NA #HOW?!
crs(myMap) # NA
st_bbox(myMap)


xf <- raster::raster("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/All_Rasters_Scaled_Weighted_UDScaled.asc") # raster not in imports####
raster::crs(xf) <- dataCRS
xf <- fishtrack3d::volumeUD(xf)
spdf.50 <- fishtrack3d::contourPoly(xf, levels = c(0.5))
sf_50 <- sf::st_as_sf(spdf.50) %>% sf::st_transform(3857)
spdf.95 <- fishtrack3d::contourPoly(xf, levels = c(0.95))
sf_95 <- sf::st_as_sf(spdf.95) %>% sf::st_transform(3857)


crs(sf_95) # lots of stuff
st_crs(sf_95) # 3857, wgs84
st_crs(sf_95 %>% sf::st_transform(3857))

ggmap::ggmap(myMap) + # basemap CRS = 3857
  
  ### NOTE: Removed 79110 rows containing missing values (geom_raster): UD SURFACE DID NOT PLOT!
  # stars::geom_stars(data = x %>% sf::st_transform(3857)) +
  
  ggplot2::geom_sf(data = sf_95 %>%
                     sf::st_transform(3857), # Vector transform after st_contour
                   # already 3857 above so converting twice but it ain't broke
                   fill = NA, inherit.aes = FALSE,
                   ggplot2::aes(colour = "95% UD")) + # https://github.com/dkahle/ggmap/issues/160#issuecomment-966812818
  
  ggplot2::geom_sf(data = sf_50 %>%
                     sf::st_transform(3857),
                   # already 3857 above so converting twice but it ain't broke
                   fill = NA, inherit.aes = FALSE, 
                   ggplot2::aes(colour = "50% UD")) +
  
  ggplot2::scale_colour_manual(name = legendtitle, values = c("50% UD" = contour2colour, "95% UD" = contour1colour)) +
  # https://stackoverflow.com/questions/64425970/ggmap-in-r-keep-google-copyright-information-on-cropped-map
  # scale_x_continuous(limits = c(myLocation[1], myLocation[3]), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(myLocation[2], myLocation[4]), expand = c(0, 0)) +
  #   Coordinate system already present. Adding new coordinate system, which will replace the existing one.
  # Scale for 'x' is already present. Adding another scale for 'x', which will replace the existing scale.
  # Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
  # Warning message: Removed 1 rows containing missing values (geom_rect). 
  
  ggplot2::ggtitle(plottitle, subtitle = plotsubtitle) +
  ggplot2::labs(x = axisxlabel, y = axisylabel, caption = plotcaption) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = legendposition, #%dist (of middle? of legend box) from L to R, %dist from Bot to Top
    legend.spacing.x = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
    legend.spacing.y = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 8),
    legend.background = ggplot2::element_rect(fill = "white", colour = NA), # element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"), # white background
    plot.background = ggplot2::element_rect(fill = "white", colour = "grey50"), # white background
    legend.key = ggplot2::element_blank(), 
    text = ggplot2::element_text(size = fontsize,  family = fontfamily)
  ) # removed whitespace buffer around legend boxes which is nice


# UD% viridis scale values are 0.002 to 0.008, should be 0:100####
y[[1]] # same as y$All_Rasters_Scaled_Weighted_UDScaled.asc = values
range(y[[1]], na.rm = TRUE)
max(y[[1]], na.rm = TRUE)
(1 / max(y[[1]], na.rm = TRUE)) * 100
length(y[[1]]) # 3621
y[[1]] <- (y[[1]] / max(y[[1]], na.rm = TRUE)) * 100
length(y[[1]]) # 3621

# stars contours compare to base contour####
range(y[[1]], na.rm = TRUE) # 5 to 100
y[[1]][which(y[[1]] <= 50)] <- 0
y[[1]][which(y[[1]] >= 50)] <- 100
range(y[[1]], na.rm = TRUE) # 0 100
# compare graphics. stars is true, fishtrack3d is shit.