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
#' @param pathsubsets Location of parent folder that contains ALL files created by dBBMM.build. No terminal slash.
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

scaleraster <- function(path = NULL, # Location of files created by dBBMM.build within a subset. No terminal slash.
                        pathsubsets = NULL, # Location of parent folder that contains ALL files created by dBBMM.build. No terminal slash.
                        pattern = ".asc",
                        weighting = w, # Weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
                        format = "ascii",
                        datatype = "FLT4S",
                        bylayer = TRUE,
                        overwrite = TRUE,
                        scalefolder = "Scaled",
                        summedname = "All_Rasters_Summed",
                        scaledname = "All_Rasters_Scaled",
                        crsloc = NULL, # Location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
                        returnObj = FALSE) {
  
  # 1. Scale individual-level UD rasters and aggregate into one group-level UD raster ####
  
  # If path has a terminal slash, remove it, it's added later
  if (substr(x = path, start = nchar(path), stop = nchar(path)) == "/") path = substr(x = path, start = 1, stop = nchar(path) - 1)
  
  # Pull all raster names from path into a list
  filelist <- as.list(list.files(path = path, pattern = pattern)) 
  
  # Read in rasters and add to list
  rasterlist <-
    lapply(filelist, function(x) raster::raster(paste0(path, "/", x))) %>% # Read in rasters
    lapply(function(x) raster::setMinMax(x)) # Set minmax values
  names(rasterlist) <- stringr::str_remove(filelist, pattern = pattern) # Name the list object (raster); need to get rid of extension e.g. ".asc"
  
  # Get resolution from first raster in rasterlist (they all have same res), assign it object, squared
  rasterres <- (raster::res(rasterlist[[1]])[1]) ^ 2
  
  # Now extract the max of maxes across subsets. To do this, repeat the steps from above, but now import all individual-level rasters across subsets
  # Extract appropriate raster names i.e. do not import scaled rasters
  
  # NOTE: if pathsubsets = NULL, it will crash. can avoid by setting pathsubsets = path, but may be more elegant solution.
  if (substr(x = pathsubsets, start = nchar(pathsubsets), stop = nchar(pathsubsets)) == "/") pathsubsets = substr(x = pathsubsets, start = 1, stop = nchar(pathsubsets) - 1)
  
  filelist_subsets <- as.list(list.files(path = pathsubsets, pattern = pattern, recursive = TRUE))
  
  names(filelist_subsets) <-  stringr::str_remove(filelist_subsets, pattern = pattern) # Name elements of the list
  
  patt <- "Scaled" # Define pattern to identify which raster names to remove
  
  filelist_subsets <- filelist_subsets[!grepl(patt, (names(filelist_subsets)))] # Now remove list elements that have the pattern in them; return full list of no match found
  
  # Read in appropriate rasters and add to list
  rasterlist_subsets <- 
    lapply(filelist_subsets, function(x) raster::raster(paste0(pathsubsets, "/", x))) %>% # Read in only rasters that occur in the filtered list
    lapply(function(x) raster::setMinMax(x)) # Reintroduce values
  
  # Get max of maxes across subsets 
  scalemax <-
    lapply(rasterlist_subsets, function(x) raster::maxValue(x)) %>% # extract maxes
    unlist() %>% # to vector
    max(na.rm = TRUE) # extract max of maxes
  
  # Create new folder to save to
  dir.create(paste0(path, "/", scalefolder))
  
  # Scale all raster values to max of maxes (maximum value becomes 1)
  rasterlist %<>%
    lapply(function(x) x / scalemax) %>% # scaling occurs here
    lapply(function(x) raster::writeRaster(x = x, # save scaled individual rasters
                                           filename = paste0(path, "/", scalefolder, "/", names(x)),
                                           format = format,
                                           datatype = datatype,
                                           if (format != "CDF") bylayer = bylayer,
                                           overwrite = overwrite))
  
  # Create a raster stack so that rasters can be summed in the step below
  rasterstack <- raster::stack(x = rasterlist)
  
  # Sum the scaled individual UDs, which should result in a single aggregated or ‘group-level' UD
  All_Rasters_Summed <- raster::stackApply(x = rasterstack, # Raster* object or list of
                                           indices = rep(1, raster::nlayers(rasterstack)), # Vector of length nlayers(x), performs the function (sum) PER UNIQUE index, i.e. 1:5 = 5 unique sums.
                                           fun = sum, # returns a single value, e.g. mean or min, and that takes a na.rm argument
                                           na.rm = TRUE, # If TRUE, NA cells are removed from calculations
                                           # filename = paste0(path, "/", scalefolder, "/", summedname, pattern), # character. Optional output filename, causes file to be written
                                           format = format,
                                           datatype = datatype,
                                           bylayer = bylayer,
                                           overwrite = overwrite)
  
  # Put the values back in the raster object
  All_Rasters_Summed %<>% raster::setMinMax()
  
  # Now scale the group-level UD
  All_Rasters_Scaled <- All_Rasters_Summed / raster::maxValue(All_Rasters_Summed)
  
  # Save this raster
  # raster::writeRaster(x = All_Rasters_Scaled,
  #                     filename = paste0(path, "/", scalefolder, "/", scaledname, pattern),
  #                     format = format,
  #                     datatype = datatype,
  #                     bylayer = bylayer,
  #                     overwrite = overwrite)
  
  # Now weight the group-level UD raster
  All_Rasters_Scaled_Weighted <- All_Rasters_Scaled / weighting
  rm(All_Rasters_Scaled) # Remove, as not used again
  
  # Save this raster too
  raster::writeRaster(x = All_Rasters_Scaled_Weighted, # resave individual rasters
                      filename = paste0(path, "/", scalefolder, "/", scaledname, "_Weighted", pattern),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # Change projection of All_Rasters_Scaled_Weighted to latlon for proper plotting
  dataCRS <- readRDS(paste0(crsloc, "CRS.Rds"))
  raster::crs(All_Rasters_Scaled_Weighted) <- dataCRS
  
  # If any NA present, replace by 0 (safety)
  All_Rasters_Scaled_Weighted@data@values[is.na(All_Rasters_Scaled_Weighted@data@values)] <- 0
  # any(is.na(All_Rasters_Scaled_Weighted@data@values))
  
  
  
  # 2. Now we will deal with the creation of a group-level UD raster for plotting purposes ####
  
  # Standardize so the values within the raster sum to 1 (required to run the getVolumeUD() below)
  All_Rasters_Scaled_Weighted <- All_Rasters_Scaled_Weighted / sum(raster::values(All_Rasters_Scaled_Weighted))    
  
  # Change the crs to LatLong for plotting and calculation purposes
  All_Rasters_Scaled_Weighted_LatLon <- raster::projectExtent(object = All_Rasters_Scaled_Weighted, crs = sp::CRS("+proj=longlat")) # crs = proj
  
  # Change res so x & y match (kills values)
  raster::res(All_Rasters_Scaled_Weighted_LatLon) <- rep(mean(raster::res(All_Rasters_Scaled_Weighted_LatLon)), 2)
  
  # Now project raster
  All_Rasters_Scaled_Weighted_LatLon <- raster::projectRaster(from = All_Rasters_Scaled_Weighted,
                                                              to = All_Rasters_Scaled_Weighted_LatLon)
  
  # Save the raster
  raster::writeRaster(x = All_Rasters_Scaled_Weighted_LatLon, 
                      filename = paste0(path, "/", scalefolder, "/", scaledname, "_Weighted_LatLon", pattern),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # Ensure again that no NAs exist in the raster; replace by 0
  All_Rasters_Scaled_Weighted_LatLon@data@values[is.na(All_Rasters_Scaled_Weighted_LatLon@data@values)] <- 0
  
  # Ensure again that cell values within a raster sum to 1
  All_Rasters_Scaled_Weighted_LatLon <- All_Rasters_Scaled_Weighted_LatLon / sum(raster::values(All_Rasters_Scaled_Weighted_LatLon))
  
  
  
  # 3. Below we will deal with the calculation of volume areas ####
  # Replace any occurring NAs with 0s. These may be introduced if region-specific UD areas do not have the same extent.
  replaceNA <- function(x, na.rm, ...){
    if (is.na(x[1]))
      return(0)
    else
      return(x)
  }
  UDlist <- rasterlist %>% sapply(function(x) raster::calc(x, fun = replaceNA))
  
  # Convert the scaled individual-level rasters within rasterlist to class ".UD". Also ensure that values within a raster sum to 1 so that they can be fed into getVolumeUD()
  UDlist <- UDlist %>% sapply(function(x) new(".UD", x / sum(raster::values(x))))
  
  # Below code calculates 50% and 95% volume areas per UD, the mean and stdev across UDs, and finally core and home range volume area sizes of the group-level UD
  # A. individual core and home range volume area sizes
  area.50 <- UDlist %>% sapply(function(x) sum(raster::values(move::getVolumeUD(x) <= .50))) # 50% volume area
  area.50 <- round((area.50 * rasterres) / 1000000, 2) # Convert from m^2 to km^2
  
  area.95 <- UDlist %>% sapply(function(x) sum(raster::values(move::getVolumeUD(x) <= .95))) # 95% volume area
  area.95 <- round((area.95 * rasterres) / 1000000, 2) # Convert from m^2 to km^2
  
  # B. Mean and SD
  area.50.mean <- round(mean(area.50), 2) # 50% volume area mean
  area.50.sd <- round(sd(area.50), 2) # 50% volume area SD
  
  area.95.mean <- round(mean(area.95), 2) # 95% volume area mean
  area.95.sd <- round(sd(area.95), 2) # 95% volume area SD
  
  # 3. Group-level core and home range volume areas
  UDScaled <- new(".UD", All_Rasters_Scaled_Weighted) # This uses the aeqd raster for calculations
  UDScaled <- UDScaled / sum(raster::values(UDScaled)) # Safety 
  UDScaled <- new(".UD", All_Rasters_Scaled_Weighted) # Convert back to .UD so getVolumeUD can work
  
  group_area.50 <- round((sum(raster::values(move::getVolumeUD(UDScaled) <= .50)) * rasterres) / 1000000, 2)
  group_area.50
  
  group_area.95 <- round((sum(raster::values(move::getVolumeUD(UDScaled) <= .95)) * rasterres) / 1000000, 2)
  group_area.95
  
  # Combine in a single df
  area.ct <- data.frame(core.use = area.50,
                        general.use = area.95
  )
  
  # Create ID column from row.names and kill row names
  area.ct$ID <- row.names(area.ct)
  row.names(area.ct) <- NULL
  
  # Add mean, sd and group-level UD values
  area.ct <- rbind(area.ct,
                   c(area.50.mean,
                     area.95.mean,
                     "mean_across_UDs"),
                   c(area.50.sd,
                     area.95.sd,
                     "sd_across_UDs"),
                   c(group_area.50,
                     group_area.95,
                     "Group-level_UD")
  )
  
  write.csv(area.ct,
            file = paste0(path, "/", scalefolder, "/","VolumeAreas_ScaledAllFish.csv"),
            row.names = FALSE)
  
  if (returnObj) return(All_Rasters_Scaled_Weighted)
  
}