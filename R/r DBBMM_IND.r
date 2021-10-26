# Mo 2021-10-22: I’d choose e18950 as shark individual


# Clear memory
rm(list = ls())

# load packages

require(move)
require(ggmap)
require(mapproj)
require(dplyr)
# require(plyr)
require(maptools)
require(circular)
require(ape)
library(rgdal)
library(sp)

###
### dBBMM analysis using the "move" package
###


# ====================================================================================================================================
# ====================================================================================================================================

### SOME NOTES REGARDING THE "MOVE" PACKAGE
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

# ====================================================================================================================================
# ====================================================================================================================================

#################################
#################################
##### Preparing the dataset #####
#################################
#################################

setwd("~/Documents/Science/PhD/Data/BiminiSharks/home range - dynamic Brownian bridge")
setwd("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09")

df <- read.table("ray.txt", sep = ",", dec = ".", header = T, na.strings = c("", " ", NA))
df <- read.table("Mo_example_data/hammer.txt", sep = ",", dec = ".", header = T, na.strings = c("", " ", NA))

head(df)

# import receiver location coordinates
vemcoloc <- read.table("vemcoreceiverlocations.txt", sep = ",", dec = ".", header = T, na.strings = c("", " ", NA))
vemcoloc <- read.table("Mo_example_data/vemcoreceiverlocations.txt", sep = ",", dec = ".", header = T, na.strings = c("", " ", NA))
vemcoloc <- dplyr::select(vemcoloc, GPS.N, GPS.W, location, Agency)
vemcoloc <- filter(vemcoloc, Agency == "BBFSF.Bim")

# link location coordinates to the object
df <- df %>%
  left_join(., vemcoloc)

# make sure only data from Bimini is used
df <- filter(df, Agency == "BBFSF.Bim")

head(df)

# remove locations too far away from island i.e. most OTN locations:
df <- filter(
  df, location != "CAT", location != "Cat West", location != "Great Isaac's South East", location != "Great Isaac's South West", location != "Great Iscas NE", location != "Great Isacs NW",
  location != "GUN", location != "GUN WEST", location != "hesperus"
)

# get rid of unimportant columns
df <- dplyr::select(df, time, station, elasmo, location, GPS.N, GPS.W, Species)

# change column names i.e. timestamp, location.long, location.lat
colnames(df)[1:7] <- c("Datetime", "station", "elasmo", "location", "Lat", "Lon", "Species")

df$elasmo <- sub("^", "e", df$elasmo) # add letter to the front of acoustic ID, otherwise drama with script

# sort by timestamp
df$Datetime <- as.POSIXct(df$Datetime, format = "%Y-%m-%d %H:%M", tz = "US/Eastern")

str(df)

df <- arrange(df, elasmo, Datetime)
head(df)
tail(df)

### Converting coordinates to UTM ###

cord.dec <- SpatialPoints(cbind(df$Lon, df$Lat), proj4string = CRS("+proj=longlat"))

# Now transform the CRS (coordinate reference system) to UTM by setting the EPSG to 32617 for WGS 84, UTM zone 17, northern hemisphere. This is where Bimini is located.
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32617"))
cord.UTM

# physical check:
# par(mfrow = c(1, 2))
# plot(cord.dec, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
# plot(cord.UTM, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)

test <- as.data.frame(cord.UTM)
colnames(test) <- c("NewEastingUTM", "NewNorthingUTM")

df <- cbind(df, test)

df <- dplyr::select(df, Datetime, station, elasmo, Species, location, Lat, Lon, NewEastingUTM, NewNorthingUTM)











# ====================================================================================================================================










#############################################
#############################################
### START: BURSTED_dBBMM - FULL TRAJ ########
#############################################
#############################################

df1 <- df

ind <- filter(df1, elasmo == "e20182")

# create move object
move.ind <- move(
  x = ind$Lon, y = ind$Lat, time = ind$Datetime,
  proj = CRS("+proj=longlat +datum=WGS84"),
  data = ind, animal = ind$elasmo, sensor = "VR2W"
)

# check class
class(move.ind)

move.ind

# check the currect projection
proj4string(move.ind)
head(move.ind)
move.ind$LocationError <- 211 # based on effective detection range of receivers (see Kessel et al. 2014). Taking most conservative estimate.
head(move.ind)


# Convert projection to Azimuthal Equi-Distance projection (aeqd)
r.ind <- spTransform(move.ind, center = TRUE)

# Make sure it changed correctly
proj4string(r.ind)

# Calculate time lag between consecutive detections. Obviously there is no time lag for the first detection, so we need to tell R this:
# knowing time lags in acoustic telemetry is important, as the motion variance is based on the time difference b/w detections at consecutive locations, i.e. larger time lag creates
# a larger gap where the animal could've been and therefore inflates the motion variance. Need to deal with this issue!
TimeDiff <- timeLag(r.ind, units = "hours")

r.ind$TimeDiff <- append(0, TimeDiff)

long <- which(r.ind$TimeDiff > 4)
length(long)

# Make a raster for the UD to plot into
# Start with UTM. These coordinates need to be big enough to cover your data.

# May need to expand x and y ranges if encountering errors when making the DBBMM in the next chunk
# NOTE: "xUTM" grid should be larger than the "x" i.e. lower minimums and higher
e <- 30000
xUTM.ind <- raster(
  xmn = min(ind$NewEastingUTM) - e, xmx = max(ind$NewEastingUTM) + e, ymn = min(ind$NewNorthingUTM) - e, ymx = max(ind$NewNorthingUTM) + e,
  crs = CRS("+proj=utm +zone=17 +datum=WGS84"), resolution = 1000
) # resolution determines grid size

# We now need to reproject this into aeqd
# Make a dummy object to get the correct projection
# (Use the CRS string from the r object above)
x.ind <- raster(
  xmn = min(df1$NewEastingUTM), xmx = max(df1$NewEastingUTM), ymn = min(df1$NewNorthingUTM), ymx = max(df1$NewNorthingUTM),
  crs = CRS("+proj=utm +zone=17 +datum=WGS84")
)
proj4string(x.ind)

# Ok now use that to reproject
# ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), (projectRaster(from, to, ...)) > project UTM to AEQD
newTemplate <- projectExtent(xUTM.ind, proj4string(x.ind))
newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster

# Give newTemplate some values. Make Rep equal to the ncell dimension
ones <- rep(1, (newTemplate@ncols * newTemplate@nrows))
xAEQD.ind <- setValues(newTemplate, ones)
xAEQD.ind

class(xAEQD.ind)

# Reproject the move object r into AEQD
rNew.ind <- spTransform(r.ind, proj4string(xAEQD.ind))
rNew.ind

# we need to make sure that the projection of the move object is in the same format as our raster
proj4string(xAEQD.ind)
proj4string(rNew.ind)

# separate short vs long time intervals b/w successive detections at consecutive locations. Cut-off is arbitrarily set to 5 hours.
bursted <- burst(rNew.ind, c("normal", "long")[1 + (timeLag(rNew.ind, units = "hours") > 4)]) # 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command

# run the dBBMM
bursted_dbbmm <- brownian.bridge.dyn(bursted, burstType = "normal", raster = xAEQD.ind, location.error = "LocationError", time.step = 2, ext = 3)
# note that the units for time step are MINUTES. 2 reflects the nominal delay of transmissions of acoustic tags

# re-standardize (Dr. Kranstauber's trouble shooting solution). Occurring errors are "due to limits of accuracy when calculating. Given the error is
# really small (< .000001 %) I would not worry to much and re standardize". then aggregate UD segments per ID
tmp <- calc(bursted_dbbmm, sum)
dbbmmlist <- new(".UD", tmp / sum(values(tmp)))

# export aggregated individual UDs to as ascii files for further exploration/spatial analysis in GIS software. Needed for when the aim is to plot population-level UDs per species.
# saveloc="~/Documents/Science/PhD/Data/BiminiSharks/Multi-species MPA paper/dBBMM/ascii files/nurse/"
asc <- writeRaster(dbbmmlist, paste(saveloc, filename = "e20182.asc", sep = ""), format = "ascii", datatype = "FLT4S", bylayer = T, overwrite = T)

check2

######################################################
########## START: BURSTED_dBBMM - PLOT DATA ##########
######################################################

###
### plot receiver locations and connections b/w them
###
# plot(move.ind) # plot a Move object
plot(move.ind, type = "o", xlab = "location_east", ylab = "location_north")
###
### Plot the DBBMM (no geographic reference) ###
###

# plot the bursted model
bursted_dbbmm <- new(".UD", calc(bursted_dbbmm, sum))
class(bursted_dbbmm)

tiff("e23156-UD.tiff", width = 140, height = 180, units = "mm", res = 500)

plot(sqrt(bursted_dbbmm))
points(rNew.ind, cex = 0.1)
contour(bursted_dbbmm, levels = c(.5, .95), col = c(6, 1), add = TRUE, lwd = .5)

dev.off()

###
### plot the DBBMM (with spatial reference) ###
###

# ?rasterToPolygons
DBBMM_poly.ind <- rasterToPolygons(trim(bursted_dbbmm, values = 0))
# 'trim' shrinks a raster object by removing outer rows and columns that have the same value (in this case '0's)

DBBMM_poly_lonlat.ind <- spTransform(DBBMM_poly.ind, CRS("+proj=longlat"))
# transforms coordinates back to longlat. The values in the column 'layer' are the UDs (in total they sum up to 1)

# Get map function makes a call to online databases and retrieves a map which is defined by the bounding box
map_elasmo.ind <- get_map(bbox(extent(DBBMM_poly_lonlat.ind) * .9), zoom = 12) ### setting maptype = "toner" creates a black background with white island)

# The fortify converts the spatialPolygonsDataFrame in a data frame that can be used in ggplot (this function is definded in that package),
# most extra columns can for example be used to indicate holes and are not very relevant for us except for the id that is a identifier for
# each cell. This function only converts spatial data and not the associated data.
DBBMM_data.ind <- fortify(DBBMM_poly_lonlat.ind)

# Here we associate the UD data with the dataframe we're using to plot
DBBMM_data.ind$col <- as.data.frame(DBBMM_poly_lonlat.ind)[DBBMM_data.ind$id, "layer"]

tiff("e23156-sr.tiff", width = 140, height = 120, units = "mm", res = 500)

ggmap(map_elasmo.ind) +
  coord_map() +
  geom_polygon(
    data = DBBMM_data.ind, aes(x = long, y = lat, group = id, alpha = sqrt(col)),
    fill = "purple"
  ) +
  scale_alpha(range = c(0, 1), "Utilization\nDistribution")

dev.off()

####################################################
########## END: BURSTED_dBBMM - PLOT DATA ##########
####################################################



#############################################################
########## START: BURSTED_dBBMM - CALC VOLUME AREA ##########
#############################################################

cont.ind <- getVolumeUD(bursted_dbbmm)
cont.ind

# Area of 95% Utilization distribution
cont95.ind <- cont.ind <= .95
cont95.ind

pixels95.ind <- sum(values(cont95.ind))
pixels95.ind

# Resolution of the raster map is given above if you run the cont object
# resolution  :   (x, y)
area95.ind <- 50.98023 * 50.42364 * pixels95.ind
area95.ind
#  m2

area95km.ind <- area95.ind / 1000000
area95km.ind
#  km2

# Area of 50% UD
cont50.ind <- cont.ind <= .5
pixels50.ind <- sum(values(cont50.ind))

# Resolution of the raster map is given above if you run the cont object
# resolution  :   (x, y)
area50.ind <- 50.98023 * 50.42364 * pixels50.ind
area50.ind
#  m2

area50km.ind <- area50.ind / 1000000
area50km.ind
#  km2

###########################################################
########## END: BURSTED_dBBMM - CALC VOLUME AREA ##########
###########################################################

###########################################
###########################################
### END: BURSTED_dBBMM - FULL TRAJ ########
###########################################
###########################################





######################################################################
######################################################################
########## START: BURSTED_dBBMM - DISCERN MOVEMENT SEGMENTS ##########
######################################################################
######################################################################

# I.E. SEGMENTS OF FULL TRAJECTORY (here: day vs night, depending on sign of value in the SunElev column; + is above, - is below horizon)
rNew.ind$SunElev <- solarpos(
  spTransform(rNew.ind, "+proj=longlat"), # solarpos calculates the position/angle of the sun, apply to time column
  timestamps(rNew.ind)
)[, 2] # a negative value in SunElev column indates sun is below the horizon and vice versa

dataBurst.ind <- burst(
  rNew.ind,
  factor(rNew.ind$SunElev[-n.locs(rNew.ind)]
  < 0)
)

burstUD.ind <- brownian.bridge.dyn(dataBurst.ind, raster = xAEQD.ind, location.error = "LocationError", time.step = 2, ext = 1.4)
# in this case it will estimate the variance for the whole trajectory at once, but will then calculate the Utilization density per part of this trajectory
### sum traj by day and night (true = night, false = day)

# see below help files:
# ?stackApply
# ?`RasterBrick-class`
# ?getZ # Initial functions for a somewhat more formal approach to get or set z values (e.g. time) associated with layers of Raster* objects. In development.

###################################################################################
########## START: BURSTED_dBBMM - day vs night (PLOT & VOLUME AREA CALC) ##########
###################################################################################

stackApply(x = burstUD.ind, indices = factor(names(getZ(burstUD.ind))), fun = sum)
plot(sqrt(stackApply(burstUD.ind, factor(names(getZ(burstUD.ind))), sum)))

###
### Plot individual levels and calculate HR contours of burst objects
###

# to do this, we need to:
# create a Burst object with factor levels,,
# run a dBBMM for each factor level independently,
# transform the DBBMMBurstStack class to .UD (a contour cannot be drawn on a DBBMMBurstStack object) and sum segments of UDs of same factor level
# plot contour lines
# calculate UD

dataBurst.ind <- burst(rNew.ind, c("night", "day")[(factor(rNew.ind$SunElev[-n.locs(rNew.ind)] > 0))])

### START: Nighttime UD ###

dataBurst.ind_dbbmm.night <- brownian.bridge.dyn(dataBurst.ind, burstType = "night", raster = xAEQD.ind, location.error = "LocationError", time.step = 120, ext = .3)
dataBurst.ind_dbbmm.night <- new(".UD", calc(dataBurst.ind_dbbmm.night, sum))

tiff("e23156 - UD (night).tiff", width = 140, height = 180, units = "mm", res = 500)
plot(sqrt(dataBurst.ind_dbbmm.night), main = "Bull 23786 (night)", xlab = "Easting (m)", ylab = "Northing (m)")
contour(dataBurst.ind_dbbmm.night, levels = c(.5, .95), col = c(6, 1), add = TRUE, lwd = .5)
points(rNew.ind, cex = 0.1)
dev.off()

# calculate volume areas

cont.ind.night <- getVolumeUD(dataBurst.ind_dbbmm.night)
cont.ind.night

# Area of 95% Utilization distribution
cont95.ind.night <- cont.ind.night <= .95
cont95.ind.night

pixels95.ind.night <- sum(values(cont95.ind.night))
pixels95.ind.night

# Resolution of the raster map is given above if you run the cont object
# resolution  : ,   (x, y)
area95.ind.night <- 51.51273 * 50.27405 * pixels95.ind.night
area95.ind.night
#  m2

area95km.ind.night <- area95.ind.night / 1000000
area95km.ind.night
#  km2

# Area of 50% UD
cont50.ind.night <- cont.ind.night <= .5
pixels50.ind.night <- sum(values(cont50.ind.night))

# Resolution of the raster map is given above if you run the cont object
# resolution  : ,   (x, y)
area50.ind.night <- 51.51273 * 50.27405 * pixels50.ind.night
area50.ind.night
#  m2

area50km.ind.night <- area50.ind.night / 1000000
area50km.ind.night
#  km2

### END: Nighttime UD ###


### START: Daytime UD ###

dataBurst.ind_dbbmm.day <- brownian.bridge.dyn(dataBurst.ind, burstType = "day", raster = xAEQD.ind, location.error = "LocationError", time.step = 2, ext = 1.4)
dataBurst.ind_dbbmm.day <- new(".UD", calc(dataBurst.ind_dbbmm.day, sum))

tiff("e23156 - UD (day).tiff", width = 140, height = 180, units = "mm", res = 500)
plot(sqrt(dataBurst.ind_dbbmm.day), main = "Bull 23786 (night)", xlab = "Easting (m)", ylab = "Northing (m)")
contour(dataBurst.ind_dbbmm.day, levels = c(.5, .95), col = c(6, 1), add = TRUE, lwd = .5)
points(rNew.ind, cex = 0.1)
dev.off()

# calculate volume areas

cont.ind.day <- getVolumeUD(dataBurst.ind_dbbmm.day)
cont.ind.day

# Area of 95% Utilization distribution
cont95.ind.day <- cont.ind.day <= .95
cont95.ind.day

pixels95.ind.day <- sum(values(cont95.ind.day))
pixels95.ind.day

# Resolution of the raster map is given above if you run the cont object
# resolution  : ,   (x, y)
area95.ind.day <- 51.51273 * 50.27405 * pixels95.ind.day
area95.ind.day
#  m2

area95km.ind.day <- area95.ind.day / 1000000
area95km.ind.day
#  km2

# Area of 50% UD
cont50.ind.day <- cont.ind.day <= .5
pixels50.ind.day <- sum(values(cont50.ind.day))

# Resolution of the raster map is given above if you run the cont object
# resolution  : ,   (x, y)
area50.ind.day <- 51.51273 * 50.27405 * pixels50.ind.day
area50.ind.day
#  m2

area50km.ind.day <- area50.ind.day / 1000000
area50km.ind.day
#  km2

### END: Daytime UD ###


###########################################################
########## START: BURSTED_dBBMM - Dendogram plot ##########
###########################################################

# compare of the space use of an individual over time bw and  among all days and nights (differences and similarities)

ind.utc <- filter(df, elasmo == "e20185")
# Mo 2021-10-22: I’d choose e18950 as shark individual (hammer dbase)

# convert datetime to UTC, otherwise the matching of day doesnt work correctly
rec.TZ <- "UTC"
dat.TZ <- "US/Eastern"
ind.utc$Datetime <- as.POSIXct(ind.utc$Datetime, format = "%Y-%m-%d %H:%M", tz = dat.TZ)
try <- ind.utc$Datetime
ind.utc$Datetime <- format(try, tz = rec.TZ)
ind.utc$Datetime <- as.POSIXct(ind.utc$Datetime, format = "%Y-%m-%d %H:%M", tz = rec.TZ)

m.ind.utc <- move(
  x = ind.utc$Lon, y = ind.utc$Lat, time = ind.utc$Datetime,
  proj = CRS("+proj=longlat +datum=WGS84"),
  data = ind.utc, animal = ind.utc$Species, sensor = "VR2W"
)

# check the currect projection
proj4string(m.ind.utc)
m.ind.utc$LocationError <- 225 # based on effective detection range of receivers (see Kessel et al. 2014). Taking most conservative estimate.

# Convert projection to Azimuthal Equi-Distance projection (aeqd)
r.ind.utc <- spTransform(m.ind.utc, center = TRUE)

# Make sure it changed correctly
proj4string(r.ind.utc)

# May need to expand x and y ranges if encountering errors when making the DBBMM in the next chunk
# NOTE: "xUTM" grid should be larger than the "x" i.e. lower minimums and higher
e <- 5000
xUTM.ind.utc <- raster(
  xmn = min(ind.utc$NewEastingUTM) - e, xmx = max(ind.utc$NewEastingUTM) + e, ymn = min(ind.utc$NewNorthingUTM) - e, ymx = max(ind.utc$NewNorthingUTM) + e,
  crs = CRS("+proj=utm +zone=17 +datum=WGS84"), resolution = 50
) # resolution determines grid size i.e. 50 x 50 m
xUTM.ind.utc
class(xUTM.ind.utc)

# We now need to reproject this into aeqd
# Make a dummy object to get the correct projection
# (Use the CRS string from the r object above)
x.ind.utc <- raster(
  xmn = min(ind.utc$NewEastingUTM), xmx = max(ind.utc$NewEastingUTM), ymn = min(ind.utc$NewNorthingUTM), ymx = max(ind.utc$NewNorthingUTM),
  crs = CRS(" +proj=aeqd +ellps=WGS84 +lon_0=-79.28215 +lat_0=25.676015")
)
proj4string(x.ind.utc)

# Ok now use that to reproject
# ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), (projectRaster(from, to, ...)) > project UTM to AEQD
newTemplate <- projectExtent(xUTM.ind.utc, proj4string(x.ind.utc))
newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster

# Give newTemplate some values. Make Rep equal to the ncell dimension
ones <- rep(1, (newTemplate@ncols * newTemplate@nrows))
xAEQD.ind.utc <- setValues(newTemplate, ones)
xAEQD.ind.utc

class(xAEQD.ind.utc)

# Reproject the move object r into AEQD
rNew.ind.utc <- spTransform(r.ind.utc, proj4string(xAEQD.ind.utc))
rNew.ind.utc

# we need to make sure that the projection of the move object is in the same format as our raster
proj4string(xAEQD.ind.utc)
proj4string(rNew.ind.utc)

# calculate sunrise and sunset times and use this to assign day and night
sunset <- sunriset(m.ind.utc, timestamps(m.ind.utc), POSIXct.out = T, direction = "sunset")$time
sunrise <- sunriset(m.ind.utc, timestamps(m.ind.utc), POSIXct.out = T, direction = "sunrise")$time
m.ind.utc$day <- ifelse(timestamps(m.ind.utc) > sunrise & timestamps(m.ind.utc) < sunset, "day", "night")

### define below your temporal resolution

# calculate day, locations before sunrise are counted to the previous day to not break up nighttime periods and keep the sections continuous
m.ind.utc$date <- as.Date(timestamps(m.ind.utc))
m.ind.utc$date <- m.ind.utc$date - (timestamps(m.ind.utc) < sunrise)

# ignore years (optional)
m.ind.utc$date <- sub("2014-", "", m.ind.utc$date)
m.ind.utc$date <- sub("2015-", "", m.ind.utc$date)
m.ind.utc$date <- sub("2016-", "", m.ind.utc$date)
m.ind.utc$date <- sub("2017-", "", m.ind.utc$date)

# create an combined identifier for day/night and date
id <- paste(m.ind.utc$day, m.ind.utc$date)

# burst m according to the identifier and assign any segment where the # identifier changes the label change
changeId <- ifelse(diff(as.numeric(factor(id))) == 0, id[-1], "Change")
mBursted <- burst(m.ind.utc, changeId)

# calculate the UD separated per behaviour
burstUD <- brownian.bridge.dyn(spTransform(mBursted, center = T), raster = xAEQD.ind.utc, location.error = "LocationError", ext = 1.4, burstType = validNames(unique(id))) # Standardize UDs
UDs <- UDStack(burstUD)

# Omit lowest probabilities of the UD that have negligible contribution to the UD
UDs[getVolumeUD(UDs) > 0.99] <- 0

# calculate earth movers distance using threshold of 50 m (distance above which the amount of work is considered equal to the threshold distance)
# if the EMD is divided by the threshold distance, the result can be standardized to values between 0 and 1
emdDists <- emd(UDs / cellStats(UDs, sum), threshold = 50)

# use ape to plot a circular diagram of similarities
tiff("Stingray - dendogram (per month).tiff", width = 300, height = 300, units = "mm", res = 500)
tree <- as.phylo(hclust(emdDists / 50))
plot(
  tree, "fan",
  no.margin = T, x.lim = c(-0.65, 0.65), cex = 0.7, tip.color = c("red", "blue")[grepl("day", tree$tip.label) + 1], label.offset = 0.01
)
dev.off()

#########################################################
########## END: BURSTED_dBBMM - Dendogram plot ##########
#########################################################

####################################################################
####################################################################
########## END: BURSTED_dBBMM - DISCERN MOVEMENT SEGMENTS ##########
####################################################################
####################################################################



###
### Calculate earth mover's distance to assess overlap (warning: will take a LONG while!)
###

emd(aggregate(st, 5, sum)) # overlap




# burst plot (quick and dirty exploration of UDs. Somehow does not plot more than 16 UDs in 1 pnale)
tiff("UD.burst-bull.tiff", width = 200, height = 200, units = "mm", res = 500)

X()
plot(sqrt(st))

plot(move.sp, type = "o", xlab = "location_east", ylab = "location_north")




###################################################
### EXPORT CONTOUR LINES TO USE IN GIS SOFTWARE ###
###################################################

writeOGR(conts, ".", "lines", "ESRI Shapefile")
