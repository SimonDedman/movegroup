#' Succinct title 8 words max
#'
#' Description paragraph: Automates delta log normal boosted regression trees abundance prediction.
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
#' @export
#' @param locations Character vector of locations of rasters to weight, e.g. locations of All_Rasters_Summed.asc. No terminal slash.
#' @param rasternames Same length character vector of names for rasters from locations.
#' @param pattern Names of the raster files in locations, e.g. "All_Rasters_Summed.asc".
#' @param weightings Same length numeric vector of weightings to apply to rasters. Rasters from locations are divided by their respective value from weightings e.g. raster 1 values divided by weightings 1 value.
#' @param saveloc Where to save output rasters. No terminal slash.
#' @param extension Extension format of output rasters.
#' @param format Raster format for output rasters.
#' @param datatype Data type for output rasters.
#' @param bylayer Raster saving detail.
#' @param overwrite Raster saving detail.
#' 
#' @importFrom raster raster setMinMax writeRaster
#' @import magrittr
#' 
# Mo receiver array doublescale####
# scale to each receiver array (loop of 3) then scale those together
# 1. run dBBMM.build.R (change this filename to dBBMMhomeRange.R ?) for each array
# 2. run scaleraster on each array
# new script: # 3. Pull All_Rasters_Scaled.asc from each array/scaled subfolder
# And weight each (of 3) scaleraster outputs by number of cells containing any receiver (a number the user calculates).
# Save to a single folder.
weightraster <- function(locations = c("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/H/Scaled",
                                       "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/L/Scaled",
                                       "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/M/Scaled"), # assumes they have subfolders called Scaled with files called All_Rasters_Scaled.asc
                         rasternames = c("H", "L", "M"), # names of rasters
                         pattern = "All_Rasters_Summed.asc", # 
                         weightings = c(1, 2, 3), # respective weightings per array
                         saveloc = "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/Weighted",
                         extension = ".asc",
                         format = "ascii",
                         datatype = "FLT4S",
                         bylayer = TRUE,
                         overwrite = TRUE)
{ # open function
  
  # user input checks
  if ((length(weightings) != length(locations)) & length(weightings) != 1) stop("length(weightings) must be 1 or length(locations)")
  if ((length(saveloc) != length(locations)) & length(saveloc) != 1) stop("length(saveloc) must be 1 or length(locations)")
  if ((length(rasternames) != length(locations)) & length(rasternames) != 1) stop("length(rasternames) must be 1 or length(locations)")
  
  # from here, edit these ####
  filelist <- as.list(list.files(path = locations, pattern = pattern, full.names = TRUE))
  
  # Read in rasters and add to list
  rasterlist <-
    lapply(filelist, function(x) raster::raster(x)) %>% # read in rasters
    lapply(function(x) raster::setMinMax(x)) # set minmax values
  
  rasterlist <- Map("/", rasterlist, weightings) # divide values by weightings
  # https://stackoverflow.com/questions/53319000/multiply-columns-of-a-data-table-by-a-vector
  
  rasterlist <- Map(setNames, rasterlist, nm = rasternames)
  # https://stackoverflow.com/questions/61655351/changing-column-names-of-multiple-raster-bricks-using-names-in-r
  
  rasterlist %<>%
    # lapply(function(x) x / cellStats(x, stat = 'max')) %>% # scale to individual max
    lapply(function(x) raster::writeRaster(x = x, # resave individual rasters
                                           filename = paste0(saveloc, "/", names(x), extension), # , pattern: removed ability to resave as different format
                                           # error: adds X to start of numerical named objects####
                                           format = format,
                                           datatype = datatype,
                                           if (format != "CDF") bylayer = bylayer,
                                           overwrite = overwrite))
  
  # Next steps in chain:
  # 4. run scaleraster on the 3 outputs (All_Rasters_Scaled.asc (which you weighted) * 3)
  # 5. run plot
  
} # close function

# weightraster(locations = c("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/H/Scaled",
#                            "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/L/Scaled",
#                            "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/M/Scaled"), # assumes they have subfolders called Scaled with files called All_Rasters_Scaled.asc
#              rasternames = c("H", "L", "M"), # names of rasters
#              pattern = "All_Rasters_Summed.asc", # 
#              weightings = c(1, 2, 3), # respective weightings per array
#              saveloc = "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/Weighted",
#              extension = ".asc",
#              format = "ascii",
#              datatype = "FLT4S",
#              bylayer = TRUE,
#              overwrite = TRUE)
