# Raster Scaling function
# Simon Dedman simondedman@gmail.com 2021-10-19

#' Automated Boosted Regression Tree modelling and mapping suite
#'
#' Automates delta log normal boosted regression trees abundance prediction.
#' Loops through all permutations of parameters provided (learning
#' rate, tree complexity, bag fraction), chooses the best, then simplifies it.
#' Generates line, dot and bar plots, and outputs these and the predictions
#' and a report of all variables used, statistics for tests, variable
#' interactions, predictors used and dropped, etc. If selected, generates
#' predicted abundance maps, and Unrepresentativeness surfaces.
#' See www.GitHub.com/SimonDedman/gbm.auto for issues, feedback, and development
#' suggestions. See SimonDedman.com for links to walkthrough paper, and papers
#' and thesis published using this package.
#'
#' @param path No terminal slash.
#' @param pattern Default ".asc".
#' @param weighting Weighting to divide individual and summed-scaled rasters by, for unbalanced arrays. Individual, Scaled, and Scaled_Weighted rasters, and the volume areas csv, will have weightings applied, but NOT the summed raster.
#' @param format Default "ascii".
#' @param datatype Default "FLT4S".
#' @param bylayer Default TRUE.
#' @param overwrite Default TRUE.
#' @param scalefolder Default "Scaled".
#' @param summedname Default "All_Rasters_Summed".
#' @param scaledname Default "All_Rasters_Scaled".
#' @param crsloc Location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
#' @param returnObj Logical. Return the scaled object to the parent environment? Default FALSE.
#' 
#' @return Line, dot and bar plots, a report of all variables used, statistics
#' for tests, variable interactions, predictors used and dropped, etc. If
#' selected generates predicted abundance maps, and Unrepresentativeness surface
#'
#' @details Errors and their origins:
#' @examples
#' \donttest{
#' # Not run
#' }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export

#' @import magrittr
#' @importFrom raster raster setMinMax res maxValue writeRaster stack stackApply nlayers projectExtent crs projectRaster values
#' @importFrom stringr str_remove
#' @importFrom sp CRS

scaleraster <- function(path = NULL, # Location of files created by dBBMM.build. No terminal slash.
                        pattern = ".asc",
                        weighting = 1, # weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
                        format = "ascii",
                        datatype = "FLT4S",
                        bylayer = TRUE,
                        overwrite = TRUE,
                        scalefolder = "Scaled",
                        summedname = "All_Rasters_Summed",
                        scaledname = "All_Rasters_Scaled",
                        crsloc = NULL, # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
                        returnObj = FALSE) {
  # library(raster)
  # If path has a terminal slash, remove it, it's added later
  if (substr(x = path, start = nchar(path), stop = nchar(path)) == "/") path = substr(x = path, start = 1, stop = nchar(path) - 1)
  
  # Pull all raster names from path into a list
  filelist <- as.list(list.files(path = path, pattern = pattern))
  
  # Read in rasters and add to list
  rasterlist <-
    lapply(filelist, function(x) raster::raster(paste0(path, "/", x))) %>% # read in rasters
    lapply(function(x) raster::setMinMax(x)) # set minmax values
  names(rasterlist) <- stringr::str_remove(filelist, pattern = pattern) # Name the list object (raster); need to get rid of extension e.g. ".asc"
  
  # get resolution from first raster in rasterlist (they all have same res), assign it object, squared
  rasterres <- (raster::res(rasterlist[[1]])[1]) ^ 2
  
  # Get max of maxes
  scalemax <-
    lapply(rasterlist, function(x) raster::maxValue(x)) %>% # extract maxes
    unlist() %>% # to vector
    max(na.rm = TRUE) # get max of maxes
  
  # create new folder to save to
  dir.create(paste0(path, "/", scalefolder))
  
  rasterlist %<>%
    lapply(function(x) x / scalemax) %>% # scale to max of maxes
    lapply(function(x) x / weighting) %>% # divide by weighting value
    lapply(function(x) raster::writeRaster(x = x, # resave individual rasters
                                           filename = paste0(path, "/", scalefolder, "/", names(x)), # , pattern: removed ability to resave as different format
                                           # error: adds X to start of numerical named objects####
                                           format = format,
                                           datatype = datatype,
                                           if (format != "CDF") bylayer = bylayer,
                                           overwrite = overwrite))
  
  # extentlist <- lapply(rasterlist, function(x) data.frame(xmin = as.vector(extent(x))[1],
  #                                                         xmax = as.vector(extent(x))[2],
  #                                                         ymin = as.vector(extent(x))[3],
  #                                                         ymax = as.vector(extent(x))[4]))
  # extentdf <- do.call(rbind, extentlist)
  # fullextents <- c(min(extentdf$xmin, na.rm = TRUE),
  #                  max(extentdf$xmax, na.rm = TRUE),
  #                  min(extentdf$ymin, na.rm = TRUE),
  #                  max(extentdf$ymax, na.rm = TRUE))
  # rasterlist %<>% lapply(function(x) extent(x) <- fullextents)
  # 
  # extent(rasterlist[[2]]) <- fullextents
  # rasterlist[[2]]@extent <- fullextents
  
  # rasterstack <- do.call(mosaic, rasterlist)
  # tmp <- mosaic(x = rasterlist[[1]],
  #               y = rasterlist[[2]])
  # 
  # st_crs(rasterlist[[1]])
  # proj4string(rasterlist[[1]])
  # proj4string(rasterlist[[2]])
  
  # sum the normalised individual UDs
  rasterstack <- raster::stack(x = rasterlist)
  # for nc: Error in compareRaster(x) : different extent: still need to fix this in build.R####
  
  All_Rasters_Summed <- raster::stackApply(x = rasterstack, # Raster* object or list of
                                           indices = rep(1, raster::nlayers(rasterstack)), # Vector of length nlayers(x), performs the function (sum) PER UNIQUE index, i.e. 1:5 = 5 unique sums.
                                           fun = sum, # returns a single value, e.g. mean or min, and that takes a na.rm argument
                                           na.rm = TRUE, # If TRUE, NA cells are removed from calculations
                                           filename = paste0(path, "/", scalefolder, "/", summedname, pattern), # character. Optional output filename, causes file to be written
                                           format = format,
                                           datatype = datatype,
                                           bylayer = bylayer,
                                           overwrite = overwrite)
  
  # another rescaling from 0 to 1
  # Should result in a single aggregated or ‘group’ level UD
  All_Rasters_Summed %<>% raster::setMinMax()
  All_Rasters_Scaled <- All_Rasters_Summed / raster::maxValue(All_Rasters_Summed)
  raster::writeRaster(x = All_Rasters_Scaled, # resave individual rasters
                      filename = paste0(path, "/", scalefolder, "/", scaledname, pattern),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  All_Rasters_Scaled %<>% lapply(function(x) x / weighting) # divide by weighting value
  raster::writeRaster(x = All_Rasters_Scaled, # resave individual rasters
                      filename = paste0(path, "/", scalefolder, "/", scaledname, "_Weighted", pattern),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # change projection of All_Rasters_Scaled to latlon for proper plotting
  dataCRS <- readRDS(paste0(crsloc, "CRS.Rds"))
  raster::crs(All_Rasters_Scaled) <- dataCRS
  # proj = CRS("+proj=longlat +datum=WGS84")
  
  
  All_Rasters_Scaled_LatLon <- raster::projectExtent(object = All_Rasters_Scaled, crs = sp::CRS("+proj=longlat")) # crs = proj
  # returns RasterLayer with projected extent, but no values. Can be adjusted (e.g. by setting its 
  # resolution) and used as a template 'to' in projectRaster.
  
  # change res so x & y match. Kills values
  raster::res(All_Rasters_Scaled_LatLon) <- rep(mean(raster::res(All_Rasters_Scaled_LatLon)), 2)
  # TODO: increase res if blocky? do further up?####
  
  
  All_Rasters_Scaled_LatLon <- raster::projectRaster(from = All_Rasters_Scaled,
                                                     to = All_Rasters_Scaled_LatLon)
  
  # All_Rasters_Scaled_LatLon <- projectRaster(from = All_Rasters_Scaled, crs = proj) # 2958
  # with crs 2958:
  # class      : RasterLayer 
  # dimensions : 482, 284, 136888  (nrow, ncol, ncell)
  # resolution : 84500, 84400  (x, y)
  # extent     : -6893280, 17104720, -20266219, 20414581  (xmin, xmax, ymin, ymax)
  # crs        : +proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs 
  # source     : memory
  # names      : All_Rasters_Summed 
  # values     : 0, 0.9535157  (min, max)
  
  # with proj=longlat:
  # class      : RasterLayer 
  # dimensions : 140, 274, 38360  (nrow, ncol, ncell)
  # resolution : 1.35, 0.855  (x, y)
  # extent     : -186.3136, 183.5864, -24.71852, 94.98148  (xmin, xmax, ymin, ymax)
  # crs        : +proj=longlat +datum=WGS84 +no_defs 
  # source     : memory
  # names      : All_Rasters_Summed 
  # values     : 0, 0.8578524  (min, max)
  
  raster::writeRaster(x = All_Rasters_Scaled_LatLon, # resave individual rasters
                      filename = paste0(path, "/", scalefolder, "/", scaledname, "_Weighted_LatLon", pattern),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # Calculate individual scaled ("relative") UD areas
  area.50 <- rasterlist %>% sapply(function(x) length(which(raster::values(x) >= 0.5))) # 50%
  area.50.sd <- sd(area.50)
  area.50.sd <- round(area.50.sd * rasterres, 1)
  area.50 <- round(area.50 * rasterres, 1) # convert from cells to metres squared area
  area.95 <- rasterlist %>% sapply(function(x) length(which(raster::values(x) >= 0.05))) # 95%
  area.95.sd <- sd(area.95)
  area.95 <- round(area.95 * rasterres, 1)
  area.95.sd <- round(area.95.sd * rasterres, 1)
  
  area.ct <- data.frame(core.use = area.50, 
                        general.use = area.95
  ) # Combine in single df
  area.ct$ID <- row.names(area.ct) # create ID column from row.names
  row.names(area.ct) <- NULL # kill row.names, reverts to 1,2,3
  area.ct <- rbind(area.ct, 
                   c(round(length(which(raster::values(All_Rasters_Scaled) >= 0.5)) * rasterres, 1), # add a row for All_Rasters_Scaled
                     round(length(which(raster::values(All_Rasters_Scaled) >= 0.05)) * rasterres, 1),
                     "All_Rasters_Scaled_Sum"),
                   c(round(area.50.sd * rasterres, 1),
                     round(area.95.sd * rasterres, 1),
                     "All_Rasters_Scaled_SD")
  )
  
  write.csv(area.ct,
            file = paste0(path, "/", scalefolder, "/","VolumeAreas_ScaledAllFish.csv"),
            row.names = FALSE)
  
  if (returnObj) return(All_Rasters_Scaled)
}