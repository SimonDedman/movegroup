#' Scales individual utilization distribution rasters and volume area estimates
#' 
#' Scales individual-level utilization distribution (UD) rasters from 0 to 1 to facilitate interpretation as relative 
#' intensity of utilization (as opposed to absolute), making comparisons across individuals and interpretations at 
#' the group level more straightforward. Subsequently, scaled individual-level rasters are aggregated to create a 
#' single group-level UD raster. See www.GitHub.com/SimonDedman/movegroup for issues, feedback, and development 
#' suggestions. There is an option to account for bias in acoustic receiver array spatial representation (see Details). 
#' 
#' Step 1. Scale rasters.
#' Individual-level UD rasters are scaled from 0 to 1 by dividing each raster by the maximum probability density value
#' occurring within the raster set.
#' 
#' Step 2. Aggregate into a group-level raster.
#' Scaled individual-level rasters are summed to create a single group-level UD raster. 
#' 
#' Step 3. Re-scale to 0 to 1.
#' The group-level raster is divided by its own maximum value.
#' 
#' Step 4. Weight raster (optional).
#' The scaled group-level UD raster is divided by the specified weighting factor(s). Note that this is only useful if you 
#' want to account for an unbalanced listening station (e.g., acoustic receivers) array and have split up the study site and receivers in regions, and have 
#' run the movegroup() for each regional data set separately. See van Zinnicq Bergmann et al. 2022 
#' (https://doi.org/10.1016/j.biocon.2022.109469) for example. If not applicable, choose a value of "1".
#' 
#' Step 5. Standardize raster.
#' Standardize the potentially weighted and scaled group-level UD raster so that its values sum to 1.
#' 
#' Step 6. Export as latlon CRS.
#' Change crs to latlon for plotting and calculation purposes, save file, continue.
#' 
#' Step 7. Estimate 50 and 95pct contour volume areas.
#' For each scaled individual-level UD raster, estimate 50 and 95pct contour volume areas, as well as their mean and standard
#' deviation. Additionally, the 50 and 95pct volume area is estimated for the group-level UD raster.
#' 
#' Step 8. Export the projected-CRS group-level raster.
#' 
#' @param path Path to directory where the individual-level UDs are saved. Likely the same as 
#' savedir from movegroup. Default NULL.
#' @param pathsubsets Path to parent directory that contains all UDs across spatial groups or 
#' subsets, i.e. if you ran movegroup multiple times for different areas in a connected system, this
#' would be the parent folder within which all the movegroup savedir's are located. Default NULL.
#' @param pattern Extension pattern used to read in all UDs in directory and pathsubsets directory. 
#' Default ".asc".
#' @param weighting Addresses unbalanced receiver array design after receivers have first been 
#' partitioned into regions, and group-level UDs estimated per region. Numeric. Weights 
#' area-specific scaled group-level UD raster by value. This then means that estimated scaled 
#' individual-level volume areas also become weighted. Default is 1 for no weighting.
#' @param format Character. Output file type for raster::writeRaster param format. Default "ascii".
#' @param datatype Character. Data type for writing values to disk for raster::writeRaster param 
#' Datatype. Default "FLT4S". 
#' @param bylayer For raster::writeRaster param bylayer. Default TRUE.
#' @param overwrite For raster::writeRaster param overwrite. Default TRUE.
#' @param scalefolder Folder to save outputs to. Default "Scaled".
#' @param scaledweightedname Name of chunk for scaled and weighted output rasters. Default 
#' "All_Rasters_Scaled_Weighted".
#' @param crsloc Location of saved CRS Rds file from movegroup.R. Should be same as path. Default 
#' NULL.
#' @param returnObj Logical. Return the scaled object to the parent environment? Default FALSE.
#' 
#' @details Errors and their origins:
#' 
#' 1. Error in (function (cond): error in evaluating the argument 'x' in selecting a method for 
#' function 'res': subscript out of bounds. Probably path can't find any files of type=pattern:
#' check you used a terminal slash in savedir in movegroup, and that path has files of type=pattern.
#' 
#' 2. Error in if (substr(x = pathsubsets, start = nchar(pathsubsets), stop = nchar(pathsubsets))==:
#' argument is of length zero: pathsubsets is wrong. Try setting to same as path. NULL does this.
#' 
#' 3. Error in gzfile(file, "rb"): cannot open compressed file 'CRS.Rds', probable reason 'No such 
#' file or directory': crsloc is wrong. Try setting to same as path. NULL does this.
#' 
#' 4. In min/max: No non-missing arguments to min; returning Inf: likely not enough memory, increase
#'  rasterResolution value.
#'
#' @return Scaled and weighted individual-level and group-level utilization distributions saved as 
#' rasters. Scaled 50 and 95pct contour volume area estimates (in km2) for individuals and the group
#' , saved in .csv format. Latlon raster.
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' # Having run the movegroup function example:
#' scaleraster(path = mysavedir)
#' 
#' # Weighted by number of positions per ID, fewer locations = lower Weighting value = higher final 
#' # UD values after dividing by Weighting. This scales all IDs up to match the group max.
#' Weighting <- TracksCleaned |>
#'  dplyr::group_by(Shark) |>
#'  dplyr::summarise(N = n()) |> 
#'  dplyr::filter(N > 23) |> 
#'  dplyr::mutate(N = N / max(N, na.rm = TRUE)) |> 
#'  dplyr::pull(N)
#'  
#'  scaleraster(path = mysavedir, weighting = Weighting)
#' }
#'
#' @export scaleraster

#' @importFrom raster raster setMinMax res maxValue writeRaster stack stackApply nlayers projectExtent crs projectRaster values
#' @importFrom purrr map2
#' @importFrom stats sd
#' @importFrom stringr str_remove
#' @importFrom sp CRS

scaleraster <- function(path = NULL,
                        pathsubsets = NULL,
                        pattern = ".asc",
                        weighting = 1,
                        format = "ascii",
                        datatype = "FLT4S",
                        bylayer = TRUE,
                        overwrite = TRUE,
                        scalefolder = "Scaled",
                        scaledweightedname = "All_Rasters_Scaled_Weighted",
                        crsloc = NULL,
                        returnObj = FALSE) {
  
  # 1. Scale individual-level UD rasters and aggregate into one group-level UD raster ####
  
  # If path has a terminal slash, remove it, it's added later
  if (substr(x = path, start = nchar(path), stop = nchar(path)) == "/") path = substr(x = path, start = 1, stop = nchar(path) - 1)
  
  # Pull all raster names from path into a list
  filelist <- as.list(list.files(path = path, pattern = pattern)) 
  
  # Read in rasters and add to list
  rasterlist <-
    lapply(filelist, function(x) raster::raster(file.path(path, x)))  |>  # Read in rasters
    lapply(function(x) raster::setMinMax(x)) # Set minmax values
  names(rasterlist) <- stringr::str_remove(filelist, pattern = pattern) # Name the list object (raster); need to get rid of extension e.g. ".asc"
  
  # Get resolution from first raster in rasterlist (they all have same res), assign it object, squared
  rasterres <- (raster::res(rasterlist[[1]])[1]) ^ 2
  
  # Now extract the max of maxes across subsets. To do this, repeat the steps from above, but now import all individual-level rasters across subsets
  # Extract appropriate raster names i.e. do not import scaled rasters
  
  # NOTE: if pathsubsets = NULL, it will crash, ditto crsloc
  if (is.null(pathsubsets)) pathsubsets <- path
  # path can be left NULL by user if using pathsubsets. If using path and not pathsubsets, above line has already made pathsubsets = path anyway.
  if (is.null(crsloc)) crsloc <- pathsubsets
  # If crsloc & pathsubsets have a terminal slash, remove it, it's added later
  if (substr(x = crsloc, start = nchar(crsloc), stop = nchar(crsloc)) == "/") crsloc = substr(x = crsloc, start = 1, stop = nchar(crsloc) - 1)
  if (substr(x = pathsubsets, start = nchar(pathsubsets), stop = nchar(pathsubsets)) == "/") pathsubsets = substr(x = pathsubsets, start = 1, stop = nchar(pathsubsets) - 1)
  
  filelist_subsets <- as.list(list.files(path = pathsubsets, pattern = pattern, recursive = TRUE))
  
  names(filelist_subsets) <-  stringr::str_remove(filelist_subsets, pattern = pattern) # Name elements of the list
  
  patt <- "Scaled" # Define pattern to identify which raster names to remove
  
  filelist_subsets <- filelist_subsets[!grepl(patt, (names(filelist_subsets)))] # Now remove list elements that have the pattern in them; return full list of no match found
  
  # Read in appropriate rasters and add to list
  rasterlist_subsets <- 
    lapply(filelist_subsets, function(x) raster::raster(file.path(pathsubsets, x)))  |>  # Read in only rasters that occur in the filtered list
    lapply(function(x) raster::setMinMax(x)) # Reintroduce values
  
  # Get max of maxes across subsets 
  scalemax <-
    lapply(rasterlist_subsets, function(x) raster::maxValue(x))  |>  # extract maxes
    unlist() |>  # to vector
    max(na.rm = TRUE) # extract max of maxes
  
  # Create new folder to save to
  dir.create(file.path(path, scalefolder))
  
  # Scale all raster values to max of maxes (maximum value becomes 1)
  rasterlist  <- lapply(rasterlist, function(x) x / scalemax)  |>  # scaling occurs here
    # lapply(function(x) x / weighting)  |>  # Weighting occurs here
    purrr::map2(.y = weighting, .f = `/`)  |>  # Weighting occurs here. .y will be recycled if length 1.
    lapply(function(x) raster::writeRaster(x = x, # save scaled individual rasters
                                           filename = file.path(path, scalefolder, "/", names(x)),
                                           format = format,
                                           datatype = datatype,
                                           if (format != "CDF") bylayer = bylayer,
                                           overwrite = overwrite))
  
  # Create a raster stack so that rasters can be summed in the step below
  rasterstack <- raster::stack(x = rasterlist)
  
  # Sum the scaled individual UDs, which should result in a single aggregated or ‘group-level' UD
  All_Rasters_Weighted_Summed <- raster::stackApply(x = rasterstack, # Raster* object or list of
                                                    indices = rep(1, raster::nlayers(rasterstack)), # Vector of length nlayers(x), performs the function (sum) PER UNIQUE index, i.e. 1:5 = 5 unique sums.
                                                    fun = sum, # returns a single value, e.g. mean or min, and that takes a na.rm argument
                                                    na.rm = TRUE, # If TRUE, NA cells are removed from calculations
                                                    # filename = paste0(path, "/", scalefolder, "/", summedname, pattern), # character. Optional output filename, causes file to be written
                                                    format = format,
                                                    datatype = datatype,
                                                    bylayer = bylayer,
                                                    overwrite = overwrite)
  
  # Put the values back in the raster object
  All_Rasters_Weighted_Summed <- raster::setMinMax(All_Rasters_Weighted_Summed)
  
  # Scale the group-level UD
  All_Rasters_Scaled_Weighted <- All_Rasters_Weighted_Summed / raster::maxValue(All_Rasters_Weighted_Summed)
  
  # Save this raster
  raster::writeRaster(x = All_Rasters_Scaled_Weighted,
                      filename = file.path(path, scalefolder, paste0(scaledweightedname, pattern)),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # Now weight the group-level UD raster
  # Change projection of All_Rasters_Scaled_Weighted to latlon for proper plotting
  dataCRS <- readRDS(file.path(crsloc, "CRS.Rds"))
  raster::crs(All_Rasters_Scaled_Weighted) <- dataCRS
  
  # If any NA present, replace by 0 (safety)
  All_Rasters_Scaled_Weighted@data@values[is.na(All_Rasters_Scaled_Weighted@data@values)] <- 0
  
  
  # 2. Deal with the creation of a group-level UD raster for plotting purposes ####
  # Standardize so the values within the raster sum to 1 (required to run the getVolumeUD() below)
  All_Rasters_Scaled_Weighted <- All_Rasters_Scaled_Weighted / sum(raster::values(All_Rasters_Scaled_Weighted))    
  
  # Change the crs to LatLong for plotting and calculation purposes
  All_Rasters_Scaled_Weighted_LatLon <- raster::projectExtent(object = All_Rasters_Scaled_Weighted,
                                                              crs = sp::CRS("+proj=longlat")) # crs = proj
  
  # Change res so x & y match (kills values)
  raster::res(All_Rasters_Scaled_Weighted_LatLon) <- rep(mean(raster::res(All_Rasters_Scaled_Weighted_LatLon)), 2)
  
  # Project old values to new raster
  All_Rasters_Scaled_Weighted_LatLon <- raster::projectRaster(from = All_Rasters_Scaled_Weighted,
                                                              to = All_Rasters_Scaled_Weighted_LatLon)
  
  # Ensure again that no NAs exist in the raster; replace by 0
  All_Rasters_Scaled_Weighted_LatLon@data@values[is.na(All_Rasters_Scaled_Weighted_LatLon@data@values)] <- 0
  
  # Ensure again that cell values within a raster sum to 1
  All_Rasters_Scaled_Weighted_LatLon <- All_Rasters_Scaled_Weighted_LatLon / sum(raster::values(All_Rasters_Scaled_Weighted_LatLon))
  
  # Save the raster
  raster::writeRaster(x = All_Rasters_Scaled_Weighted_LatLon, 
                      filename = file.path(path, scalefolder, paste0(scaledweightedname, "_LatLon", pattern)),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  
  # 3. Calculate volume areas ####
  # Replace any occurring NAs with 0s. These may be introduced if region-specific UD areas do not have the same extent.
  replaceNA <- function(x, na.rm, ...){
    if (is.na(x[1]))
      return(0)
    else
      return(x)
  }
  UDlist <- rasterlist |> sapply(function(x) raster::calc(x, fun = replaceNA))
  
  # Convert the scaled individual-level rasters within rasterlist to class ".UD". Also ensure that values within a raster sum to 1 so that they can be fed into getVolumeUD()
  UDlist <- UDlist |> sapply(function(x) new(".UD", x / sum(raster::values(x))))
  
  
  
  # Calculate 50% and 95% volume areas per UD, the mean and stdev across UDs, and finally core and home range volume area sizes of the group-level UD
  # A. individual core and home range volume area sizes
  # area.50 <- UDlist %>% sapply(function(x) sum(raster::values(move::getVolumeUD(x) <= .50))) # 50% volume area
  area.50 <- UDlist |> sapply(function(x) sum(raster::values(x) >= (max(x@data@values) * 0.5))) # 50% volume area
  area.50 <- round((area.50 * rasterres) / 1000000, 4) # Convert from m^2 to km^2
  
  # area.95 <- UDlist %>% sapply(function(x) sum(raster::values(move::getVolumeUD(x) <= .95))) # 95% volume area
  area.95 <- UDlist |> sapply(function(x) sum(raster::values(x) >= (max(x@data@values) * 0.05))) # 95% volume area
  area.95 <- round((area.95 * rasterres) / 1000000, 4) # Convert from m^2 to km^2
  
  # B. Mean and SD
  area.50.mean <- round(mean(area.50), 4) # 50% volume area mean
  area.50.sd <- round(stats::sd(area.50), 4) # 50% volume area SD
  
  area.95.mean <- round(mean(area.95), 4) # 95% volume area mean
  area.95.sd <- round(stats::sd(area.95), 4) # 95% volume area SD
  
  # 3. Group-level core and home range volume areas
  UDScaled <- new(".UD", All_Rasters_Scaled_Weighted) # This uses the aeqd raster for calculations
  # UDScaled <- UDScaled / sum(raster::values(UDScaled)) # Safety: ensures raster values sum to 1
  # UDScaled <- new(".UD", UDScaled) # Convert back to .UD so getVolumeUD can work
  
  # Save the raster
  raster::writeRaster(x = UDScaled, 
                      filename = file.path(path, scalefolder, paste0(scaledweightedname, "_UDScaled", pattern)),
                      format = format,
                      datatype = datatype,
                      bylayer = bylayer,
                      overwrite = overwrite)
  
  # group_area.50 <- round((sum(raster::values(move::getVolumeUD(UDScaled) <= .50)) * rasterres) / 1000000, 4)
  # group_area.95 <- round((sum(raster::values(move::getVolumeUD(UDScaled) <= .95)) * rasterres) / 1000000, 4)
  
  group_area.50 <- round((sum(raster::values(UDScaled) >= (max(UDScaled@data@values) * 0.5)) * rasterres) / 1000000, 4)
  group_area.95 <- round((sum(raster::values(UDScaled) >= (max(UDScaled@data@values) * 0.05)) * rasterres) / 1000000, 4)
  # See movegroup.R L748. Could switch above lines to this:
  # area.50 <- round(sum(raster::values(move::getVolumeUD(bb) <= .50)), 4)
  # area.95 <- round(sum(raster::values(move::getVolumeUD(bb) <= .95)), 4)
  
  
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
                   c(round(area.50.sd, 1),
                     round(area.95.sd, 1),
                     "sd_across_UDs"),
                   c(group_area.50,
                     group_area.95,
                     "Group-level_UD")
  )
  
  # 2023-10-04 Vital memory bug warning
  if ((round(area.50.sd, 1) == 0) | (round(area.95.sd, 1) == 0)) print("No standard deviation: all individual UDs identical. Possibly due to insufficient memory for raster calculations. Check rasterResolution in movegroup")
  
  write.csv(area.ct,
            file = file.path(path, scalefolder, "VolumeAreas_ScaledAllFish.csv"),
            row.names = FALSE)
  
  if (returnObj) return(All_Rasters_Scaled_Weighted)
  
}