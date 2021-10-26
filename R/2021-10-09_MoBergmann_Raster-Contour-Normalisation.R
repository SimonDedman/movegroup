
# Clear memory
rm(list = ls())

# load packages

require(move)
require(ggmap)
require(mapproj)
require(dplyr)
require(plyr)
require(maptools)
require(circular)

### dBBMM analysis using the "move" package
###


#====================================================================================================================================
#====================================================================================================================================

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

#====================================================================================================================================
#====================================================================================================================================

#################################
#################################
##### Preparing the dataset #####
#################################
#################################

setwd("~/Documents/Science/PhD/Data/Mariana Fuentes/dBBMM/tiger")

df <- read.table("data.txt",sep=",",dec=".",header=T,na.strings=c(""," ",NA))

head(df)

# import receiver location coordinates
vemcoloc <- read.table("~/Documents/Science/PhD/Data/vemcoreceiverlocations.txt",sep=",",dec=".",header=T,na.strings=c(""," ",NA))
vemcoloc <- dplyr::select(vemcoloc, GPS.N,GPS.W,location,Agency)
vemcoloc <- filter(vemcoloc, Agency=="BBFSF.Bim")

# link location coordinates to the object
df <- df %>%
  left_join(.,vemcoloc)

# make sure only data from Bimini is used
df <- filter(df, Agency == "BBFSF.Bim")

head(df)

# remove locations too far away from island i.e. most OTN locations:
df <- filter(df,location!="CAT",location!="Cat West",location!="Great Isaac's South East",location!="Great Isaac's South West",location!="Great Iscas NE",location!="Great Isacs NW",
             location!="GUN",location!="GUN WEST",location!="hesperus")

# get rid of unimportant columns
df <- dplyr::select(df, time, station, elasmo, location, GPS.N, GPS.W,Species)

# change column names i.e. timestamp, location.long, location.lat
colnames(df)[1:7] <- c("Datetime","station","elasmo","location","Lat","Lon","Species")

df$elasmo <- sub("^", "e",df$elasmo) # add letter to the front of acoustic ID, otherwise drama with script

# sort by timestamp
df$Datetime <- as.POSIXct(df$Datetime,format="%Y-%m-%d %H:%M",tz="US/Eastern")

str(df)

df <- arrange(df, elasmo, Datetime)
head(df)
tail(df)

### Converting coordinates to UTM ###

cord.dec = SpatialPoints(cbind(df$Lon, df$Lat), proj4string=CRS("+proj=longlat"))

# Now transform the CRS (coordinate reference system) to UTM by setting the EPSG to 32617 for WGS 84, UTM zone 17, northern hemisphere. This is where Bimini is located.
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32617"))
cord.UTM

#physical check:
#par(mfrow = c(1, 2))
#plot(cord.dec, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
#plot(cord.UTM, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)

test <- as.data.frame(cord.UTM)
colnames(test) <- c("NewEastingUTM", "NewNorthingUTM")

df <- cbind(df,test)

df <- dplyr::select(df, Datetime, station, elasmo, Species, location, Lat, Lon, NewEastingUTM, NewNorthingUTM)

# add residency events which could be used to break up the dataset into subsets. This deals with very large motion variances that would result from long absenses from the array.
# Current threshold set to: 24 h
RE<-numeric(dim(df)[1])
RE[1]<-1

for (i in 2:length(df$elasmo)){
  if(df$elasmo[i]==df$elasmo[i-1] & df$Datetime[i]<=df$Datetime[i-1]+(3600*6)){
    RE[i] <- RE[i-1]
    # print(i/dim(df)[1]*100) #progress bar
  }
  else {
    if(df$elasmo[i]==df$elasmo[i-1] & df$Datetime[i]>df$Datetime[i-1]+(3600*6)){
      RE[i] <- RE[i-1]+1
      # print(i/dim(df)[1]*100) #progress bar
    }
  }
}

df$RE <- RE

nbr.re.events <- ddply(df, .(elasmo,RE), nrow)
# nbr.re.events

unique(df$elasmo)
count.det <- ddply(df, .(elasmo), summarise, det.count = length(Datetime))
count.det1 <- ddply(df, .(elasmo), summarise, count = length(unique(location)))
count.det$loc.count <- count.det1$count 
count.det










#====================================================================================================================================


















#############################################
#############################################
### START: STACK dBBMM - Multiple ID ########
#############################################
#############################################

df.new<-df

# filter out IDs that were detected in 1 location
r <- which(count.det$loc.count<2)
count.det <- count.det[-r,]
new.id<- data.frame(count.det$elasmo)
colnames(new.id)[1]<-"elasmo"

c.df.new <- inner_join(df.new,new.id)

df.new <- c.df.new

move.sp <- move(x=df.new$Lon, y=df.new$Lat, time=df.new$Datetime, 
                proj=CRS("+proj=longlat +datum=WGS84"),
                data=df.new, animal=df.new$elasmo, sensor="VR2W")

# check class
class(move.sp)

move.sp

# check the currect projection
proj4string(move.sp)
head(move.sp)
move.sp$LocationError = 225 # based on effective detection range of receivers (see Kessel et al. 2014). Taking most conservative estimate.
head(move.sp)

# Convert projection to Azimuthal Equi-Distance projection (aeqd)
r.sp = spTransform(move.sp, center=TRUE)

# Make sure it changed correctly
proj4string(r.sp)

# Calculate time lag between consecutive detections. Obviously there is no time lag for the first detection, so we need to tell R this:
# knowing time lags in acoustic telemetry is important, as the motion variance is based on the time difference b/w detections at consecutive locations, i.e. larger time lag creates 
# a larger gap where the animal could've been and therefore inflates the motion variance. Need to deal with this issue! 
TimeDiff <- timeLag(r.sp, units="hours")

# need to add a "0" to start of each sublist 
for (i in 1:length(TimeDiff)){
  TimeDiff[[i]] <- append(0,TimeDiff[[i]])
}

r.sp$TimeDiff <- unlist(TimeDiff) # unlist the list and add to move object. Done!

long <- which(r.sp$TimeDiff>6)
length(long)
#View(r.sp[long,12])

# Make a raster for the UD to plot into
# Start with UTM. These coordinates need to be big enough to cover your data.

# May need to expand x and y ranges if encountering errors when making the DBBMM in the next chunk
# NOTE: "xUTM" grid should be larger than the "x" i.e. lower minimums and higher
e<-12000
xUTM.sp = raster(xmn=min(df.new$NewEastingUTM)-e, xmx=max(df.new$NewEastingUTM)+e, ymn=min(df.new$NewNorthingUTM)-e, ymx=max(df.new$NewNorthingUTM)+e,
                 crs =CRS("+proj=utm +zone=17 +datum=WGS84"), resolution=50) # resolution determines grid size i.e. 50 x 50 m
xUTM.sp
class(xUTM.sp)

# We now need to reproject this into aeqd
# Make a dummy object to get the correct projection
# (Use the CRS string from the r object above)
x.sp = raster(xmn=min(df.new$NewEastingUTM), xmx=max(df.new$NewEastingUTM), ymn=min(df.new$NewNorthingUTM), ymx=max(df.new$NewNorthingUTM), 
              crs = CRS("+proj=utm +zone=17 +datum=WGS84"))
proj4string(x.sp)

# Ok now use that to reproject
# ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), (projectRaster(from, to, ...)) > project UTM to AEQD
newTemplate = projectExtent(xUTM.sp, proj4string(x.sp))
newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster

# Give newTemplate some values. Make Rep equal to the ncell dimension
ones = rep(1, (newTemplate@ncols*newTemplate@nrows))
xAEQD.sp = setValues(newTemplate, ones)
xAEQD.sp

class(xAEQD.sp) 

# Reproject the move object r into AEQD 
rNew.sp = spTransform(r.sp, proj4string(xAEQD.sp))
rNew.sp

#we need to make sure that the projection of the move object is in the same format as our raster
proj4string(xAEQD.sp)
proj4string(rNew.sp)

# Plot xAEQD and r
plot(xAEQD.sp)
points(rNew.sp)

# a dBBMM is not run if the # locations < window size (default value, 31), therfore need to filter out IDs with insufficient locations. Use dplyr:
check1 <- ddply(df.new, c("elasmo"), summarise, det.count = length(Datetime))
check1
count.det <- filter(check1,det.count>31)
count.det
length(check1$elasmo)==length(count.det$elasmo)

# IF FALSE: rerun shortened df
df.new <- semi_join(df.new,count.det) # go back to create a new "move.sp" object

IDs <- levels(trackId(rNew.sp)) #rNew.sp is a ‘moveStack' object
dbbmmlist<-list()




i=4
dbbmmlist[[i]]@data@max # extracts the max cell value of the UD
> dbbmmlist[[i]]@data@values
logical(0) ?????????

## seems like for every list, the last ID doesn't seem to have values stored in the raster... how to fix this..





dbbmmlist[[i]]@data@values <- dbbmmlist[[i]]@data@values/max(dbbmmlist[[i]]@data@values) # normalize the raster by dividing by the largest value. It stretches the values from 0 to 1.

plot(dbbmmlist[[1]])


