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
#' @param grids Explanatory data to predict to. Import with (e.g.) read.csv and
#' specify object name. Defaults to NULL (won't predict to grids).
#' @param path no terminal slash.
#' @param pattern default ".asc".
#' @param format default "ascii".
#' @param datatype default "FLT4S".
#' @param bylayer default TRUE.
#' @param overwrite default TRUE.
#' @param scalefolder default "Scaled".
#' @param summedname default "All_Rasters_Summed".
#' @param scaledname default "All_Rasters_Scaled".
#' @param returnObj Logical. Return the scaled object to the parent environment? Default FALSE.
#' 
#' @return Line, dot and bar plots, a report of all variables used, statistics
#' for tests, variable interactions, predictors used and dropped, etc. If
#' selected generates predicted abundance maps, and Unrepresentativeness surface
#'
#' @details Errors and their origins:
#' @examples
#' \donttest{
#' # Not run. Note: grids file was heavily cropped for CRAN upload so output map
#' # predictions only cover patchy chunks of the Irish Sea, not the whole area.
#' # Full versions of these files:
#' # https://drive.google.com/file/d/1WHYpftP3roozVKwi_R_IpW7tlZIhZA7r
#' # /view?usp=sharing
#' library(gbm.auto)
#' data(grids)
#' data(samples)
#' # Set your working directory
#' gbm.auto(grids = grids, samples = samples, expvar = c(4:8, 10), resvar = 11,
#' tc = c(2,7), lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)}
#'
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export

#' @import raster
#' @import magrittr

scaleraster <- function(path = NULL, # Location of files created by dBBMM.build. No terminal slash.
                        pattern = ".asc",
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
    lapply(filelist, function(x) raster(paste0(path, "/", x))) %>% # read in rasters
    lapply(function(x) setMinMax(x)) # set minmax values
  names(rasterlist) <- str_remove(filelist, pattern = pattern) # Name the list object (raster); need to get rid of extension e.g. ".asc"

  # Get max of maxes
  scalemax <-
    lapply(rasterlist, function(x) maxValue(x)) %>% # extract maxes
    unlist() %>% # to vector
    max(na.rm = TRUE) # get max of maxes

  # create new folder to save to
  dir.create(paste0(path, "/", scalefolder))

  rasterlist %<>%
    lapply(function(x) x / scalemax) %>% # scale to max of maxes
    # lapply(function(x) x / cellStats(x, stat = 'max')) %>% # scale to individual max
    lapply(function(x) writeRaster(x = x, # resave individual rasters
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
  
  All_Rasters_Summed <- stackApply(x = rasterstack, # Raster* object or list of
                                   indices = rep(1, nlayers(rasterstack)), # Vector of length nlayers(x), performs the function (sum) PER UNIQUE index, i.e. 1:5 = 5 unique sums.
                                   fun = sum, # returns a single value, e.g. mean or min, and that takes a na.rm argument
                                   na.rm = TRUE, # If TRUE, NA cells are removed from calculations
                                   filename = paste0(path, "/", scalefolder, "/", summedname, pattern), # character. Optional output filename, causes file to be written
                                   format = format,
                                   datatype = datatype,
                                   bylayer = bylayer,
                                   overwrite = overwrite)

  # another rescaling from 0 to 1
  # Should result in a single aggregated or ‘group’ level UD
  All_Rasters_Summed %<>% setMinMax()
  All_Rasters_Scaled <- All_Rasters_Summed / maxValue(All_Rasters_Summed)
  writeRaster(x = All_Rasters_Scaled, # resave individual rasters
              filename = paste0(path, "/", scalefolder, "/", scaledname, pattern),
              format = format,
              datatype = datatype,
              bylayer = bylayer,
              overwrite = overwrite)
  
  # change projection of All_Rasters_Scaled to latlon for proper plotting
  dataCRS <- readRDS(paste0(crsloc, "CRS.Rds"))
  crs(All_Rasters_Scaled) <- dataCRS
  # proj = CRS("+proj=longlat +datum=WGS84")
  
  
  All_Rasters_Scaled_LatLon <- projectExtent(object = All_Rasters_Scaled, crs = CRS("+proj=longlat")) # crs = proj
  # returns RasterLayer with projected extent, but no values. Can be adjusted (e.g. by setting its 
  # resolution) and used as a template 'to' in projectRaster.
  
  # change res so x & y match. Kills values
  res(All_Rasters_Scaled_LatLon) <- rep(mean(res(All_Rasters_Scaled_LatLon)), 2)
  # TODO: increase res if blocky? do further up?####
  
  
  All_Rasters_Scaled_LatLon <- projectRaster(from = All_Rasters_Scaled,
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
  
  writeRaster(x = All_Rasters_Scaled_LatLon, # resave individual rasters
              filename = paste0(path, "/", scalefolder, "/", scaledname, "_LatLon", pattern),
              format = format,
              datatype = datatype,
              bylayer = bylayer,
              overwrite = overwrite)
  # Error in .startAsciiWriting(x, filename, ...) : x has unequal horizontal and vertical resolutions.
  # Such data cannot be stored in arc-ascii format
  
  # Calculate individual scaled ("relative") UD areas
  area.50 <- rasterlist %>% sapply(function(x) length(which(values(x) >= 0.5))) # 50%
  area.95 <- rasterlist %>% sapply(function(x) length(which(values(x) >= 0.05))) # 95%
  area.ct <- data.frame(core.use = area.50, general.use = area.95) # Combine in single df
  area.ct$ID <- row.names(area.ct) # create ID column from row.names
  row.names(area.ct) <- NULL # kill row.names, reverts to 1,2,3
  area.ct <- rbind(area.ct, c(length(which(values(All_Rasters_Scaled) >= 0.5)), # add a row for All_Rasters_Scaled
                              length(which(values(All_Rasters_Scaled) >= 0.05)),
                              "All_Rasters_Scaled"
                              )
                   )
  
  write.csv(area.ct,
            file = paste0(path, "/", scalefolder, "/","VolumeAreas_ScaledAllFish.csv"),
            row.names = FALSE)
  
  if (returnObj) return(All_Rasters_Scaled)
}
