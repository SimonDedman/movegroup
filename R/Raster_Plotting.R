# 2021-10-12 Ggplot rasters created by brownian.bridge.dyn
# Simon Dedman, simondedman@gmail.com
library(raster)
setwd("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/")
filelist <- as.list(list.files(pattern = ".asc"))

# filelist[[1]]
# do as apply ####
# tmp <- raster(filelist[[1]])
rasterlist <- as.list(filelist)
names(rasterlist) <- filelist
rasterlist <- lapply(rasterlist, function(x) raster(x))





par("mar") # 5.1 4.1 4.1 2.1
par(mar = c(1,1,1,1))
plot(tmp)
tmp
# class      : RasterLayer
# dimensions : 234, 155, 36270  (nrow, ncol, ncell)
# resolution : 50, 50  (x, y)
# extent     : 671608.8, 679358.8, 2841800, 2853500  (xmin, xmax, ymin, ymax)
# crs        : NA
# source     : 267.asc
# names      : X267


# Need to set CRS so can change it later
# Convert to SF object
#
range(tmp[]) # 0.00000000 0.06514624
# numbers are all incredibly low, maybe not an issue, but plot is all blank i.e. 0
dim(tmp) # 234 155   1
hist(tmp[]) # all zeroes?
tmp2 <- tmp[][which(tmp[] != 0)]
hist(tmp2) # almost all < 0.005
# is the extent wrong? zoomed out to a huge area? The latter. Required for buffering to work in other script. Can edit later in ggplot/df



tmp <- raster(filelist[[2]])
plot(tmp) # same, blank. Max
range(tmp[]) # 0.000000000 0.007339591
tmp2 <- tmp[][which(tmp[] != 0)]
hist(tmp2) # almost all < 0.005 but more of a distribution trailing down to 0.003


# Plot rasters in ggplot ####
# https://rdrr.io/cran/ggspatial/man/layer_spatial.Raster.html 
# https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
# https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster
# also https://r-spatial.github.io/stars/articles/stars5.html
