#' Combines region-specific group-level UD rasters into a single raster.
#'
#' Extends the spatial extent of each area-specific group-level raster to the spatial extent shared by all rasters. 
#' This will only be required if you have multiple individuals (e.g. different sharks) divided amongst a few discrete areas 
#' (e.g. around different islands) and the effort (e.g. receiver coverage) is different among islands. 
#' Not required for multiple individuals all within the same region or sampling regime.
#'
#' @param folderroots Character vector of locations of folder roots output by movegroup. Function 
#' expects CRS.Rds file and a subfolder with the scaled raster.
#' @param foldernames Character vector names of folders corresponding to files in folderroots, i.e. 
#' the names of the objects, arrays, regions, etc.
#' @param pattern For input rasters from scaleraster. Default ".asc".
#' @param scalefolder For input rasters from scaleraster. Default "Scaled".
#' @param scaledweightedname For input rasters from scaleraster. Default "All_Rasters_Scaled".
#' @param savefolder Single character entry of folder to save outputs, no trailing slash.
#' @param format Character. Output file type for raster::writeRaster param format. Default ascii, 
#' other options have generally not worked well in SD's experience.
#' @param datatype Character. Data type for writing values to disk for raster::writeRaster param 
#' datatype. Default FLT4S.
#' @param bylayer For raster::writeRaster param bylayer. Default TRUE.
#' @param overwrite For raster::writeRaster param overwrite. Default TRUE.
#' @param returnObj Logical. Return the scaled object to the parent environment to be assigned as an
#'  object? Default FALSE.
#' 
#' @return Region-specific group-level UD rasters that share the same spatial extent.
#'
#' @details When used in a movegroup pipeline, the order would be: movegroup.R, scaleraster.R, 
#' alignraster.R if required, plotraster.R.
#' 
#' @examples
#' \donttest{
#' # Not run
#' }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}

#' @export

#' @importFrom sp bbox
#' @importFrom raster crs setMinMax raster extend writeRaster
#' @importFrom purrr map2
#' @importFrom terra project



# read in rasters & add to list####
alignraster <- function(folderroots = c("/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/H", # character vector of locations of folder roots output by movegroup. Function expects CRS.Rds file and a subfolder with the scaled raster.
                                        "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/L",
                                        "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/M"),
                        foldernames = c("H", "L", "M"), # character vector names of folders corresponding to files in folderroots, i.e. the names of the objects, arrays, regions, etc.
                        pattern = ".asc", # for input rasters from scaleraster
                        scalefolder = "Scaled", # for input rasters from scaleraster
                        scaledweightedname = "All_Rasters_Scaled", # for input rasters from scaleraster
                        savefolder = "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/dBBMM ASCII/Aligned", # single character entry, no trailing slash
                        format = "ascii", # save format
                        datatype = "FLT4S", # save format
                        bylayer = TRUE, # save format
                        overwrite = TRUE, # save format
                        returnObj = FALSE # return rasterlist object?
) {
  if (length(folderroots) != length(foldernames)) stop("length of folderroots and foldernames must be equal")
  # If folderroots or savefolder have a terminal slash, remove it, it's added later
  for (folders in folderroots) {
    if (substr(x = folders, start = nchar(folders), stop = nchar(folders)) == "/") folderroots[which(folderroots %in% folders)] = substr(x = folders, start = 1, stop = nchar(folders) - 1)
  }
  
  if (substr(x = savefolder, start = nchar(savefolder), stop = nchar(savefolder)) == "/") savefolder = substr(x = savefolder, start = 1, stop = nchar(savefolder) - 1)
  
  foldernames <- as.list(foldernames)
  
  # Read in CRS files as list
  crslist <- as.list(paste0(folderroots, "/CRS.Rds"))  |> 
    lapply(function(x) readRDS(x))
  names(crslist) <- foldernames # unnecessary?
  
  rasterlist <- 
    as.list(paste0(folderroots, "/", scalefolder, "/", scaledweightedname, pattern))  |> # Pull all raster names from folderroots into a list
    lapply(function(x) raster::raster(x))  |>  # read in rasters
    lapply(function(x) raster::setMinMax(x))  |>  # set minmax values
    # https://stackoverflow.com/questions/72063819/use-an-arrow-assignment-function-as-r-purrr-map2
    purrr::map2(crslist, ~ {raster::crs(.x) <- .y;.x})  |> 
    purrr::map2(foldernames, ~ {names(.x) <- .y;.x})
  
  # calculate full shared extent
  sharedextent <- lapply(rasterlist, function(x) as.vector(sp::bbox(x))) # xmin #ymin #xmax #ymax
  sharedextent <- data.frame(t(sapply(sharedextent, c)))
  sharedextent <- c(min(sharedextent[1]), # xmin
                    max(sharedextent[3]), # xmax
                    min(sharedextent[2]), # ymin
                    max(sharedextent[4])) # ymax
  
  # align to same spatial extent
  rasterlist <- lapply(rasterlist, function(x) raster::extend(x, sharedextent))
  
  # Convert to SpatRaster format to be used by {terra}
  rasterlist <- lapply(rasterlist, function(x) as(x, "SpatRaster"))
  
  # Reproject all rasters simultaneously
  rasterlist <- lapply(rasterlist, function(x) project(x, y = rasterlist[[length(rasterlist)]]))
  
  # Convert back to RasterLayer to save CRS
  rasterlist <- lapply(rasterlist, function(x) raster::raster(x))
  
  # Save CRS
  rasterlistCRS <- sp::CRS(sp::proj4string(rasterlist[[1]]))
  class(rasterlistCRS) # CRS
  write.csv(sp::proj4string(rasterlistCRS), paste0(savefolder, "/", "CRS.csv"), row.names = FALSE)
  saveRDS(rasterlistCRS, file = paste0(savefolder, "/", "CRS.Rds"))
  
  rasterlist <- lapply(rasterlist, function(x) raster::writeRaster(x = x, # resave individual rasters
                                                         filename = paste0(savefolder, "/", names(x)), # , pattern: removed ability to resave as different format
                                                         # error: adds X to start of numerical named objects####
                                                         format = format,
                                                         datatype = datatype,
                                                         if (format != "CDF") bylayer = bylayer,
                                                         overwrite = overwrite))
  
  if (returnObj) return(rasterlist)
} # close function