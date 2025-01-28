# Vital Heim, date unknown

## Quick comment for code:
## I used section [C] below to prefilter the dataframe before running it through the movegroup::movegroup()
## The data frame had the acoustic detections by individual (elasmo) sorted by timestamp as usual
## For each row I had two extra columns that contained the different grouping variables for 
## studyperiod (i.e. night, outside, within) and
## the group of the shark (provisioned (fed) vs. naive)
## but for the fishery overlap paper I used the grouping variables species, sex, and season.
## So section C can be adjusted for as many variables as the user would like.

## I then ran through section [C] and then ran movegroup::movegroup for the subsetted dataframe, ds.j
## This I repeated for each group combo, i.e. [C] to create the desired ds.j, then movegroup::mnovegroup()

## Once I ran all needed movegroup::movegroup() I then went back to re-define the variables for each group
## then run the movegroup::scaleraster()

## So it was a lot of scrolling up and down in the code and I haven't gotten around to automate it yet. And
## so far it had relatively few group combos, so it was doable manually.
## But I figured perhaps this idea could be used as a starting step for the grouping issue on Github.

### ....................................................................................................
### [C] Prefilter data for species, region and season ----
### ....................................................................................................

# C0: Selections - ONLY thing you need to change in section C ----

# CHANGE THIS ONLY ...........................................................................
i = 3 #Select for studyperiod: 1 = night, 2 = outside, 3 = within
#j = 1 #Select region: Central = 1, North = 2, South = 3
#k = 1 #Select day vs. night: day = 1, night = 2
#l = 2 #Select life stage: adult = 1, juv = 2
m = 2 #Select shark group: naive = 1, provisioned = 2, unknown = 3
# ............................................................................................

# C1.i: Studyperiod selection - DONT change this ----

spp <- sort(unique(DET$studyperiod))
#> spp
#[1] "night"   "outside" "within" 

# Select studyperiod
spp.f <- spp[i]

# Filter for studyperiod
DET.i <- filter(DET, studyperiod == spp.f)

# C2.l = Filter for shark group ----
#sort(unique(DET$Elasmo))

sgrp <- sort(unique(DET$group))
sgrp.f <- sgrp[m]

DET.i <- filter(DET.i, sgroup == sgrp.f)

# C3: summarise data and check before submitting to dBBMM calc ----

ds.j <- DET.i
table(ds.j$Elasmo)
length(unique(ds.j$Elasmo))
length(unique(ds.j$Station.Name))
#sort(unique(ds.j$Station.Name))

### Preparations are completed!
### Now, we can execute the movegroup package

### ....................................................................................................
### [D] Calcualte dBBMMs using movegroup package ----
### ....................................................................................................

# Pipeline
# 1. step: run dBBMM.build.R (change this filename to movegroup.R ?) for each species & region
# 2. step: run scaleraster on each array
# 3. step: Pull All_Rasters_Scaled.asc from each array/scaled subfolder
# 4. step: And weight each (of 3) scaleraster outputs by number of cells containing any receiver (a number the user calculates).
# 5. step: Save to a single folder.

# D1: Construct individual-level dBBMMs ----

## Remember to add daynight period to saveloc if you want to look at different dayperiods as well.

movegroup(
  data = ds.j, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
  ID = "Elasmo", # column name of IDs of individuals.
  Datetime = "Datetime", # name of Datetime column. Must be in POSIXct format.
  Lat = "Lat", # name of Lat & Lon columns in data.
  Lon = "Lon",
  # Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  dat.TZ = "US/Eastern", # timezone for as.POSIXct.
  proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function. 
  # projectedCRS = "+init=epsg:32617", # 32617 EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
  sensor = "VR2W", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
  moveLocError = 500, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
  timeDiffLong = 2, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  timeDiffUnits = "hours", # units for time difference for move function.
  center = TRUE, # center move object within extent? See spTransform.
  buffpct = 3, #3, # buffer extent for raster creation, proportion of 1.
  rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
  # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
  rasterCRS = sp::CRS("+proj=utm +zone=17 +datum=WGS84"), # CRS for raster creation. 17
  rasterResolution = 200, # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
  # dbblocationerror = "LocationError", # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
  dbbext = 0.3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
  dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = paste0(saveloc, sgrp.f, "/", spp.f, "/"),  # save outputs to a temporary directory (default) else.
  alerts = TRUE # audio warning for failures
) 
beep(4)

# D2: Scale the rasters from 0 to 1, and weigh based on regions ----

## WARNING: MAKE SURE ALL INDIVIDUAL-LEVEL UDs OF INTEREST HAVE BEEN CREATED FIRST BEFORE RUNNING SCALERASTER()
scaleraster(
  path = paste0(saveloc, sgrp.f, "/", spp.f, "/"), # same as savedir in dBBMM.build. Contains rasters to be scaled.
  pathsubsets = paste0(saveloc), # 'Parental folder': Location of ALL files created by dBBMM.build. No terminal slash.
  pattern = ".asc",
  weighting = 1,
  # weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
  format = "ascii",
  datatype = "FLT4S",
  bylayer = TRUE,
  overwrite = TRUE,
  scalefolder = "Scaled",
  # weightedsummedname = "All_Rasters_Non_Weighted_Summed",
  scaledweightedname = "All_Rasters_Scaled_Non_Weighted",
  crsloc = paste0(saveloc, sgrp.f, "/", spp.f, "/"),
  # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
  returnObj = FALSE
)
beep(4)