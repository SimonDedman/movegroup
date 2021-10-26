# Mo note 2021-10-09:
# I reviewed the ATT paper by Udyawer, and one of my criticisms was the lack of up-to-date utility of the package as well as moving beyond individual UDs.
# However, if you want to make strong inter-individual comparisons, it’s important to normalise among UDs.
# Otherwise you get a confounding effect by different durations of tag deployment.
# This is not addressed at all in their paper, nor by anyone who uses his package.
# It’s important we emphasise this in the methods.

# https://link.springer.com/article/10.1186/s40317-018-0162-2;
# A standardised framework for analysing animal detections from automated tracking arrays
# Vinay Udyawer, Ross G. Dwyer, Xavier Hoenner, Russell C. Babcock, Stephanie Brodie, Hamish A. Campbell, Robert G. Harcourt, Charlie Huveneers, Fabrice R. A. Jaine, Colin A. Simpfendorfer, Matthew D. Taylor & Michelle R. Heupel
# Animal Biotelemetry volume 6, Article number: 17 (2018)

# dBBMM text from Mo: /home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/Quantifying nursery habitat use in juvenile sharks ms draft_May 21_latest_MvZB.docx
# no dBBMM section in it though

# SD: I suspect we can turn the outputs into a map figure in R, without having to use GIS
# MvZB: when I tried I got to a dead end. For some reason there was always one UD for which I couldn’t normalise the scale. See bottom.

# Clear memory
rm(list = ls())

# load packages

require(move)
require(ggmap)
require(mapproj)
require(dplyr)
require(magrittr)
# require(plyr) # messes with dplyr
require(maptools)
require(circular)

### dBBMM analysis using the "move" package
# SOME NOTES REGARDING THE "MOVE" PACKAGE
# Several arguments need to be used to run the model: window size, margin size, extent and raster.
# 1) window size corresponds with number of locations and moves along a given trajectory to estimate the MA parameter within defined subsections of the path.
# This increases the ability to detect breakpoints where changes in behaviour occur. The window size should relate to what kind of behaviours the model is desired to identify
# e.g., a window size of 19 means the sliding window is moved every 19 locations or every 19 hours (has to do with sampling interval)

# 2) Margin: Motion variance based on only the middle section of the trajector; the ends of the movement trajectory where no changes are allowed because at some stage you
# want to have a few locations to base your estimation of the variance on and how many locations in either side of the window we use for this, is called the margin.
# Smaller values for window size and margin is expected to give a higher frequency of behavioural changes; make these large for looking at migrations.

# 3) The raster dictates the grid cell size for the UD to be calculated per grid cell per individual. Create the raster of certain size that matches with coordinates used to
# make the move object

# 4) The extent argument is incorporated if there are animal locations that border the edges of the raster.


# Preparing the dataset ####
setwd("~/Documents/Science/PhD/Data/Mariana Fuentes/dBBMM/tiger")
setwd("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/Mo_example_data/")
# import receiver location coordinates
vemcoloc <- read.table("vemcoreceiverlocations.txt",
                       sep = ",",
                       header = T,
                       na.strings = c("", " ", NA)) %>%
  select(GPS.N, GPS.W, location, Agency) %>%
  filter(Agency == "BBFSF.Bim")

df <- read.table("data.txt",
                 sep = ",",
                 header = T,
                 na.strings = c("", " ", NA)) %>%
  left_join(vemcoloc) %>% # link location coordinates to the object
  filter(Agency == "BBFSF.Bim", # make sure only data from Bimini is used
         # remove locations too far away from island i.e. most OTN locations
         !location %in% c("CAT", "Cat West", "Great Isaac's South East", "Great Isaac's South West",
                          "Great Iscas NE", "Great Isacs NW",  "GUN", "GUN WEST", "hesperus")) %>%
  select(time, station, elasmo, Species, location, GPS.N, GPS.W) %>% # get rid of unimportant columns
  # change column names i.e. timestamp, location.long, location.lat # colnames(df)[1:7]
  rename(Datetime = time, Lat = GPS.N, Lon = GPS.W) %>%
  mutate(Datetime = as.POSIXct(Datetime, format = "%Y-%m-%d %H:%M", tz = "US/Eastern"), # proper time format
         elasmo = sub("^", "e", elasmo)) %>%  # add letter to the front of acoustic ID, otherwise drama with script
  arrange(elasmo, Datetime) # sort by timestamp

# Converting coordinates to UTM
coords <- SpatialPoints(cbind(df$Lon, df$Lat), proj4string = CRS("+proj=longlat")) %>%
  spTransform(CRS("+init=epsg:32617")) %>% # # transform CRS to UTM. 32617 = WGS 84, UTM zone 17, northern hemisphere. Bimini location.
  as.data.frame() %>%
  rename(NewEastingUTM = coords.x1, NewNorthingUTM = coords.x2)

df %<>% bind_cols(coords)

# add residency events RE####
# add residency events which could be used to break up the dataset into subsets.
# This deals with very large motion variances that would result from long absences from the array.
# Current threshold set to: 24 h
RE <- numeric(dim(df)[1])
RE[1] <- 1 # why this?

for (i in 2:length(df$elasmo)) {
  if (df$elasmo[i] == df$elasmo[i - 1] & df$Datetime[i] <= df$Datetime[i - 1] + (3600 * 6)) {
    # df$elasmo[i] == df$elasmo[i - 1] means this shark same as last shark
    RE[i] <- RE[i - 1] # subtract 1. Removes 1 added above?
    # print(i/dim(df)[1]*100) #progress bar
  } else {
    if (df$elasmo[i] == df$elasmo[i - 1] & df$Datetime[i] > df$Datetime[i - 1] + (3600 * 6)) {
      RE[i] <- RE[i - 1] + 1
      # print(i/dim(df)[1]*100) #progress bar
    }
  }
}
df$RE <- RE

nbr.re.events <- df %>%
  group_by(elasmo, RE) %>%
  tally()

# What's the point of this section?####
# count.det not used
count.det <- df %>%
  group_by(elasmo) %>%
  summarise(det.count = length(Datetime),
            loc.count = length(unique(location))) %>%
  filter(loc.count > 1) # filter out IDs that were detected in 1 location



# a dBBMM is not run if the # locations < window size (default value, 31),
# therefore need to filter out IDs with insufficient locations
# 2021-10-19 only the first elasmo fails this
df %<>%
  group_by(elasmo) %>%
  mutate(n = n()) %>%
  filter(n > 31) %>%
  select(-n)

# START: STACK dBBMM - Multiple ID ####
for (uniqueelasmo in unique(df$elasmo)) { # uniqueelasmo <- "e16961"
  # Expecting a for loop here?####
  df.new <- df[which(df$elasmo == uniqueelasmo), ]

  # next 2 lines do the same as the 1 above, second one does nothing.
  # new.id <- data.frame(elasmo = count.det$elasmo)
  # df.new %<>% inner_join(new.id) # I presume this reduces df.new by removing rows filtered above
  # But can't test with these data as none filtered above. Semi_join works below.
  # df.new %<>% inner_join(count.det$elasmo) # could do this directly?

  move.sp <- move( # move::move create a move object
    x = df.new$Lon,
    y = df.new$Lat,
    time = df.new$Datetime,
    proj = CRS("+proj=longlat +datum=WGS84"),
    data = df.new, # extra data associated with the relocations
    animal = df.new$elasmo, # animal ID(s) or name(s)
    sensor = "VR2W" # Sensor name(s)
  )

  class(move.sp) # check class, = "move"
  move.sp
  # class       : Move
  # features    : 12
  # extent      : -79.3227, -79.31871, 25.65654, 25.67408  (xmin, xmax, ymin, ymax)
  # crs         : +proj=longlat +datum=WGS84 +no_defs
  # variables   : 8
  # names       :   Datetime, station,                location,      Lat,       Lon,    NewEastingUTM,   NewNorthingUTM, RE
  # min values  : 1489052580,  109655, south west south turtle, 25.65654,  -79.3227, 668360.513631175, 2838715.91143071,  1
  # max values  : 1490010360,  126409,        west west turtle, 25.67408, -79.31871, 668736.385947233, 2840663.87712246,  5
  # timestamps  : 2017-03-09 04:43:00 ... 2017-03-20 07:46:00 Time difference of 11 days  (start ... end, duration)
  # sensors     : VR2W
  # indiv. data : elasmo, Species
  # indiv. value: e16736 Galeocerdo cuvier
  # date created: 2021-08-25 01:42:43

  # check the currect projection
  proj4string(move.sp) # "+proj=longlat +datum=WGS84 +no_defs"
  # based on effective detection range of receivers (see Kessel et al. 2014). Taking most conservative estimate.
  move.sp$LocationError <- 225

  # Convert projection to Azimuthal Equi-Distance projection (AEGD)
  r.sp <- spTransform(move.sp, center = TRUE)

  # Make sure it changed correctly
  proj4string(r.sp) # "+proj=aeqd +lat_0=25.73195 +lon_0=-79.265615 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

  # Calculate time lag between consecutive detections. Obviously there is no time lag for the first detection, so we need to tell R this:
  # knowing time lags in acoustic telemetry is important, as the motion variance is based on the time difference b/w detections at consecutive locations, i.e. larger time lag creates
  # a larger gap where the animal could've been and therefore inflates the motion variance. Need to deal with this issue!
  TimeDiff <- timeLag(r.sp, units = "hours")

  # need to add a "0" to start of each sublist
  # for (i in 1:length(TimeDiff)) {
  #   TimeDiff[[i]] <- append(0, TimeDiff[[i]])
  # }
  # r.sp$TimeDiff <- unlist(TimeDiff) # unlist the list and add to move object. Done!

  # TimeDiff isn't a list it's a vector, just need to add a 0 with c(). Possibly because it's only 1 individual.
  r.sp$TimeDiff <- c(0, TimeDiff)

  long <- which(r.sp$TimeDiff > 6)
  length(long) # 4
  # View(r.sp[long,12])

  # Make a raster for the UD to plot into. Start with UTM.
  # These coordinates need to be big enough to cover your data.
  # May need to expand x and y ranges if encountering errors when making the DBBMM in the next chunk
  # NOTE: "xUTM" grid should be larger than the "x" i.e. lower minimums and higher
  # The variable 'e' influences the extent.
  # The resolution determines the grid size.
  # Finer resolution/grid size (in m) i.e. lower value means longer computation time. May need/want to play with this value too.
  buffpct = 30 # buffer extent as a %
  buffpct <- buffpct / 100
  xUTM.sp <- raster(xmn = min(df.new$NewEastingUTM, na.rm = TRUE) - ((max(df.new$NewEastingUTM, na.rm = TRUE) - min(df.new$NewEastingUTM, na.rm = TRUE)) * buffpct),
                   xmx = max(df.new$NewEastingUTM, na.rm = TRUE) + ((max(df.new$NewEastingUTM, na.rm = TRUE) - min(df.new$NewEastingUTM, na.rm = TRUE)) * buffpct),
                   ymn = min(df.new$NewNorthingUTM, na.rm = TRUE) - ((max(df.new$NewNorthingUTM, na.rm = TRUE) - min(df.new$NewNorthingUTM, na.rm = TRUE)) * buffpct),
                   ymx = max(df.new$NewNorthingUTM, na.rm = TRUE) + ((max(df.new$NewNorthingUTM, na.rm = TRUE) - min(df.new$NewNorthingUTM, na.rm = TRUE)) * buffpct),
                   crs = CRS("+proj=utm +zone=17 +datum=WGS84"),
                   resolution = 50) # resolution determines grid size i.e. 50 x 50 m
  # xUTM.sp
  # class      : RasterLayer
  # dimensions : 62, 12, 744  (nrow, ncol, ncell)
  # resolution : 50, 50  (x, y)
  # extent     : 668247.8, 668847.8, 2838148, 2841248  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs

  # We now need to reproject this into AEQD. Make a dummy object to get the correct projection
  # (Use the CRS string from the r object above)
  x.sp <- raster(
    xmn = min(df.new$NewEastingUTM),
    xmx = max(df.new$NewEastingUTM),
    ymn = min(df.new$NewNorthingUTM),
    ymx = max(df.new$NewNorthingUTM),
    crs = CRS("+proj=utm +zone=17 +datum=WGS84")
  )
  proj4string(x.sp) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

  # use that to reproject
  # ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), (projectRaster(from, to, ...)) > project UTM to AEQD
  newTemplate <- projectExtent(xUTM.sp, proj4string(x.sp))
  newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster
  # class      : RasterLayer
  # dimensions : 848, 707, 599536  (nrow, ncol, ncell)
  # resolution : 50, 50  (x, y)
  # extent     : 656360.5, 691710.5, 2825893, 2868293  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs

  # Give newTemplate some values. Make Rep equal to the ncell dimension
  ones <- rep(1, (newTemplate@ncols * newTemplate@nrows))
  xAEQD.sp <- setValues(newTemplate, ones)
  xAEQD.sp
  # class      : RasterLayer
  # dimensions : 62, 12, 744  (nrow, ncol, ncell)
  # resolution : 50, 50  (x, y)
  # extent     : 668247.8, 668847.8, 2838148, 2841248  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
  # source     : memory
  # names      : layer
  # values     : 1, 1  (min, max)

  # Reproject the move object r into AEQD
  rNew.sp <- spTransform(r.sp, proj4string(xAEQD.sp))
  rNew.sp
  # class      : RasterLayer
  # dimensions : 848, 707, 599536  (nrow, ncol, ncell)
  # resolution : 50, 50  (x, y)
  # extent     : 656360.5, 691710.5, 2825893, 2868293  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
  # source     : memory
  # names      : layer
  # values     : 1, 1  (min, max)

  # we need to make sure that the projection of the move object is in the same format as our raster
  proj4string(xAEQD.sp) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
  proj4string(rNew.sp)  # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
  proj4string(xAEQD.sp) == proj4string(rNew.sp) # TRUE

  # Plot xAEQD and r
  plot(xAEQD.sp) # rectangle at centre of plot. x lims too wide?####
  points(rNew.sp) # points added

  # a dBBMM is not run if the # locations < window size (default value, 31),
  # therefore need to filter out IDs with insufficient locations. Use dplyr:
  # check1 <- ddply(df.new, c("elasmo"), summarise, det.count = length(Datetime))
  # count.det <- df.new %>% # count.det used again, presumably ok?
  #   # group_by(elasmo) %>%
  #   summarise(det.count = length(Datetime)) %>%
  #   filter(det.count > 31)
  # IF FALSE: rerun shortened df?
  # df.new %<>% semi_join(count.det) # go back to create a new "move.sp" object
  # IDs <- levels(trackId(rNew.sp)) # rNew.sp is a ‘moveStack' object
  # IDs includes elasmo's filtered out by count.det/semi_join above i.e. e16736####

  # FromHere 2021-10-13: should something be populated into dbbmmlist?####
  # Can't diagnose without anything here

  # xAEQD.sp = raster template extents
  # rNew.sp = move object points as raster


  rNew.sp@data@max # Error: trying to get slot "max" from an object (class "data.frame") that is not an S4 object
  xAEQD.sp@data@max # 1 (looks like all values are 1)

  # Realistically we just need to pull all the rasters from /Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/
  # into a list
  # filter out those with <31 locations
  # extract the maxes
  # normalise
  # resave (lapply)

  dbbmmlist <- list()
  i <- 4
  dbbmmlist[[i]]@data@max # extracts the max cell value of the UD
  # > dbbmmlist[[i]]@data@values
  # logical(0) ?????????
  ## seems like for every list, the last ID doesn't seem to have values stored in the raster... how to fix this..

  # Normalise raster attempts####
  dbbmmlist[[i]]@data@values <- dbbmmlist[[i]]@data@values / max(dbbmmlist[[i]]@data@values)
  # normalize the raster by dividing by the largest value. It stretches the values from 0 to 1.
  plot(dbbmmlist[[1]])

} # close for uniqueelasmo
