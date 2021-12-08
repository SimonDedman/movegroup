# Raster Scaling function
# Simon Dedman simondedman@gmail.com 2021-10-19
scaleraster <- function(path = NULL, # no terminal slash
                        pattern = ".asc",
                        format = "ascii",
                        datatype = "FLT4S",
                        bylayer = TRUE,
                        overwrite = TRUE,
                        scalefolder = "Scaled") {
  library(raster)
  # If path has a terminal slash, remove it, it's added later
  if (substr(x = path, start = nchar(path), stop = nchar(path)) == "/") path = substr(x = path, start = 1, stop = nchar(path) - 1)

  # Pull all raster names from path into a list
  filelist <- as.list(list.files(path = path, pattern = pattern))

  # Filter out those with <23 locations [already dealt with in {r filter_data} in r_Rob_lemons_SD.R]

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

  # scale to max of maxes
  rasterlist %<>%
    lapply(function(x) x / scalemax) %>% # scale
    lapply(function(x) writeRaster(x = x, # resave individual rasters
                                   filename = paste0(path, "/", scalefolder, "/", names(x)), # , pattern removed ability to resave as different format
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
                                   filename = paste0(path, "/", scalefolder, "/", "All_Rasters_Summed", pattern), # character. Optional output filename, causes file to be written
                                   format = format,
                                   datatype = datatype,
                                   bylayer = bylayer,
                                   overwrite = overwrite)

  # another rescaling from 0 to 1
  # Should result in a single aggregated or ‘group’ level UD
  All_Rasters_Summed %<>% setMinMax()
  All_Rasters_Scaled <- All_Rasters_Summed / maxValue(All_Rasters_Summed)
  writeRaster(x = All_Rasters_Scaled, # resave individual rasters
              filename = paste0(path, "/", scalefolder, "/", "All_Rasters_Scaled", pattern),
              format = format,
              datatype = datatype,
              bylayer = bylayer,
              overwrite = overwrite)
}
