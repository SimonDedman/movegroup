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
#' @param grids Explantory data to predict to. Import with (e.g.) read.csv and
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
#' @param returnObj Logical. Return the scaled objet to the parent environment? Default FALSE.
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

scaleraster <- function(path = NULL, # no terminal slash
                        pattern = ".asc",
                        format = "ascii",
                        datatype = "FLT4S",
                        bylayer = TRUE,
                        overwrite = TRUE,
                        scalefolder = "Scaled",
                        summedname = "All_Rasters_Summed",
                        scaledname = "All_Rasters_Scaled",
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
  names(rasterlist) <- filelist # need to get rid of extension e.g. ".asc"

  # Get max of maxes
  scalemax <-
    lapply(rasterlist, function(x) maxValue(x)) %>% # extract maxes
    unlist() %>% # to vector
    max(na.rm = TRUE) # get max of maxes

  # create new folder to save to
  dir.create(paste0(path, "/", scalefolder))

  # scale to max of maxes & write individual rasters
  rasterlist %<>%
    lapply(function(x) x / scalemax) %>% # scale
    lapply(function(x) writeRaster(x = x, # resave individual rasters
                                   filename = paste0(path, "/", scalefolder, "/", names(x)), # , pattern: removed ability to resave as different format
                                   # error: adds X to start of numerical named objects####
                                   format = format,
                                   datatype = datatype,
                                   bylayer = bylayer,
                                   overwrite = overwrite))

  # sum the normalised individual UDs
  rasterstack <- raster::stack(x = rasterlist)
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
  
  if (returnObj) return(All_Rasters_Scaled)
}