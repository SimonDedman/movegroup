# From SSM Exploration, Blocklab
# Simon Dedman, simondedman:gmail.com, 2023-10-12
# Plot all tracks
#https://jmlondon.github.io/crwexampleakbs/analysis.html

# library(sf)
# library(ggspatial)
# library(devtools)
# devtools::install_github('jmlondon/ptolemy')
# after install, load and follow prompts to download source data
# library(ptolemy) #not sure if I need this for my code or just their example

plotsTracksByID <- function(x = NULL, # data object, e.g. as read in by read.csv
                            Latitude = NULL, # name or number of latitude column in data, guessed if left blank,
                            Longitude = NULL, # name or number of longitude column in data, guessed if left blank,
                            DateTime = NULL, # name or number of DateTime column in data, guessed if left blank, to be parsed with as.POSIXct,
                            Tz = "America/NewYork", # name of timezone DateTimes were recorded in. Default is US east coast, "America/NewYork". See TIMEZONES URL
                            ID = NULL # name or number of ID column in data, guessed if left blank
                            # ADD PLOTRASTER GOOGLEMAP PARAMS ####
) {
# guess Lat Lon DateTime, ID from colnames
  latvariants <- c("Latitude", "latitude", "Lat", "lat")
  lonvariants <- c("Longitude", "longitude", "Lon", "lon")
  DateTimeVariants <- c("DateTime", "dateTime", "Datetime", "datetime", "Date", "date", "Time", "time") # add code to test for >1 match, if so, use the first (i.e. DateTimes, then dates, then times)
  IDvariants <- c("ID", "id", "iD", "Id", "eventid", "shark")
  # NEED CODE####
  # if any of the terms match, the match value is the column number to rename with Latitude, Longitude, DateTime, ID
  
  
  # group by ID
x %<>% dplyr::group_by(ID)
  # Necessary?####


# create multilinestring tracks from points. Lat & lon become geometry
sf_locs <- sf::st_as_sf(x, coords = c("longitude","latitude")) %>%
  sf::st_set_crs(4326) #using jlondon example, may be incorrect? Is N.Pac

sf_lines <- sf_locs %>%
  dplyr::arrange(ID, DateTime) %>% #already arranged by ID & date?
  dplyr::group_by(ID) %>% #already grouped above? NO. x grouped, grouping lost in sf_locs conversion
  dplyr::summarise(do_union = FALSE) %>%
  sf::st_cast("MULTILINESTRING")

# ptolemy/jmlondon reprojected from wgs84/4326 to projected (Bering) but mine works as wgs84 so keep
# sf_locs <- sf_locs %>% sf::st_transform(6931)
# sf_lines <- sf_lines %>% sf::st_transform(6931)


# INSERT GOOGLE MAP SECTION HERE ####
#
#
#
#
#


#library(gbm.auto)
# source('C:/Users/simon/Dropbox/Galway/Analysis/R/gbm.auto/R/gbm.basemap.R')
# natlantic <- gbm.basemap(grids = x,
#                          gridslat = 6,
#                          gridslon = 7,
#                          extrabounds = TRUE,
#                          getzip = "./GSHHS_shp/")
# setwd("C:/Users/simon/Documents/BL/iccat_SSM_data/outputs")
# setwd("/media/Seagate/Work/Blocklab/iccat_SSM_data/outputs")
# natlantic <- read_sf("./CroppedMap/Crop_Map.shp") #don't need since gbm.basemap update
#class(natlantic) # "sf"         "tbl_df"     "tbl"        "data.frame"


# Replace with plotraster google map elements ####
ggplot() +
  annotation_spatial(natlantic, fill = "grey", lwd = 0) +
  layer_spatial(sf_lines, size = 0.75, aes(color = ID)) +
  theme(legend.position = "none") +
  #scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(-180, 180, by = 5)) +
  ggtitle("Atlantic bluefin tuna tracks",
          subtitle = "1996:2019 satellite & archival tags")
# why does this add a huge top & bottom white buffer region around the plot
# when it doesn't for my data with bad basemap? Function of rstudio plot window
# export, unlink x&y, stretch, then save.
}