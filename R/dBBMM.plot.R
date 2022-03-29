#' @import ggplot2
#' @importFrom ggmap get_map

dBBMM_plot <- function(
    x = paste0(data.dir, "Scaled/All_Rasters_Scaled.asc"), # path to scaled data
    # dataCRS = 2958, # one of (i) character: a string accepted by GDAL, (ii) integer, a valid EPSG value (numeric), or (iii) an object of class crs.
    trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE
    myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
    # Default NULL, extents autocreated from data.
    # c(-79.3, 25.68331, -79.24, 25.78)
    googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion
    # factors since Google Maps tiling zoom setup doesn't align to myLocation
    # extents.
    gmapsAPI = NULL, # enter your google maps API here, quoted character string
    expandfactor = 1.6, # extents expansion factor for basemap.
    # 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
    # 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
    mapzoom = 5, # google: 3 (continent) - 21 (building). stamen: 0-18
    mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present.
    maptype = "satellite", # Type of map for ggmap::get_map.
    contour1colour = "red", # colour for contour 1, typically 95%.
    contour2colour = "orange", # colour for contour 2, typically 50%.
    plottitle = "Aggregated 95% and 50% UD contours",
    # Can use the term 'home range' when an animal can be detected wherever it goes
    # i.e. using GPS, satellite or acoustic telemetry whereby it is known that acoustic
    # receivers cover the entire home range of the study species. 
    # This term is problematic when applied to a passive acoustic telemetry setting
    # where an array of non-overlapping receivers are used to assess local space use patterns
    # i.e. the home range is bigger than the coverage by the acoustic array; put in Details
    plotsubtitle = "Scaled contours. n = 13", # data %>% distinct(ID) %>% nrow() # 13
    legendtitle = "Percent UD Contours",
    plotcaption = paste0("dBBMM_HomeRange, ", today()),
    axisxlabel = "Longitude",
    axisylabel = "Latitude",
    legend.position = c(0.16, 0.92), #%dist (of middle? of legend box) from L to R, %dist from Bot to Top.
    filesavename = paste0(today(), "_dBBMM-contours.png"),
    savedir = tempdir() # file.path(work.dir, out.dir, "Scaled")
) {
  # ToDo
  # expandfactor could default to NULL and have a formula to set it based on the size of extents
  # gbm.basemap
  # gMaps API tutorial
  # trim section optional, depends magrittr
  # 50 & 95% breaks could be editable as function params but will be a bit of work
  
  library(stars) # read_stars
  library(ggmap) # get_map
  library(magrittr) # %<>%
  
  # tmp <- raster(x)
  # setMinMax(tmp)
  # plot(tmp)
  # is.na(tmp)
  # dataCRS <- readRDS(paste0(crsloc, "CRS.Rds"))
  # crs(tmp) <- dataCRS
  
  
  # Import raster
  x <- read_stars(x) %>% st_set_crs(4326) # 4326 2958
  # read_stars doens't have most of the info that raster() has
  # class(dataCRS)
  if (stars:::is_curvilinear(x)) stop(print("x is curvilinear; first reproject to planar"))
  
  y <- x # make dupe object else removing all data < 0.05 means the 0.05 contour doesn't work in ggplot
  if (trim) { # trim raster extent to data?
    is.na(y[[1]]) <- y[[1]] == 0 # replace char pattern (0) in whole df/tbl with NA
    is.na(y[[1]]) <- y[[1]] < (max(y[[1]], na.rm = TRUE) * 0.05) # replace anything < 95% contour with NA since it won't be drawn
  }
    y %<>% starsExtra::trim2() # remove NA columns, which were all zero columns. This changes the bbox accordingly
  
  if (is.null(myLocation)) myLocation <- st_bbox(y) %>% as.vector() # st_bbox(y %>% st_transform(4326))
  
  # Create basemap with gbm.auto####
  # # Remove gbm.auto from dependency at top if not using
  # dir.create(paste0(out.dir, "basemap"))
  # bounds <- myLocation[c(1, 3, 2, 4)] # c(xmin,xmax,ymin,ymax)
  # # run to generate basemap:
  # crop_map <- gbm.basemap(bounds = bounds,
  #                         res = "f",
  #                         # getzip = paste0(work.dir, out.dir, "basemap/GSHHS_shp"), # comment out first time, uncomment subsequent
  #                         savedir = paste0(work.dir, out.dir, "basemap"),
  #                         returnsf = TRUE)
  # # run to use generated basemap later:
  # crop_map <- st_read(dsn = paste0(work.dir, out.dir, "basemap/CroppedMap/Crop_Map.shp"),
  #                     layer = paste0("Crop_Map"),
  #                     quiet = TRUE) # read in worldmap
  
  # if (is.null(gmapsAPI)) register_google(key = gmapsAPI, # an api key
  #                                        account_type = "standard",
  #                                        write = TRUE)
  
  if (mapsource != "google") googlemap <- FALSE # in case user forgot to set both
  
  if (expandfactor != 0) { # grow bounds extents if requested
    xmid <- mean(myLocation[c(1,3)])
    ymid <- mean(myLocation[c(2,4)])
    xmax <- ((myLocation[3] - xmid) * expandfactor) + xmid #updated for sf/st
    xmin <- xmid - ((xmid - myLocation[1]) * expandfactor)
    ymax <- ((myLocation[4] - ymid) * expandfactor) + ymid
    ymin <- ymid - ((ymid - myLocation[2]) * expandfactor)
    myLocation <- c(xmin, ymin, xmax, ymax)
    if (googlemap) myLocation <- c(mean(c(myLocation[1], myLocation[3])), mean(c(myLocation[2], myLocation[4]))) # googlemap needs a center lon lat
  }
  
  myMap <- get_map(
    location = myLocation, # -62.57564  28.64368  33.78889  63.68533 # stamen etc want a bounding box
    zoom = mapzoom, # 3 (continent) - 21 (building)
    # scale = "auto", # default "auto", 1, 2, 4 all the same
    source = mapsource, # "google" # using stamen as fallback
    maptype = maptype, # "satellite"
    messaging = TRUE,
    crop = TRUE # google maps crs = 4326
  ) 
  
  # class(myMap) # "ggmap"  "raster"
  # tmp <- raster::crop(x = myMap, y = myLocation)
  # Error in (function (classes, fdef, mtable): unable to find an inherited method for function ‘crop’ for signature ‘"ggmap"’
  
  # Define a function to fix the bbox to be in EPSG:3857
  # https://stackoverflow.com/a/50844502/1736291
  # Fixes "error no lon value" in ggmap below
  # Move this to separate script & source it####
  ggmap_bbox <- function(map) {
    if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
    # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
    # and set the names to what sf::st_bbox expects:
    map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
    # Convert the bbox to an sf polygon, transform it to 3857, 
    # and convert back to a bbox (convoluted, but it works)
    bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
    # Overwrite the bbox of the ggmap object with the transformed coordinates 
    attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
    attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
    attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
    attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
    map
  }
  myMap <- ggmap_bbox(myMap) # Use the function
  
  # Automate width * height adjustments for different map extent / ratio
  # 6 (manually chosen width, below), divided by width range times by height range
  # Maintains ratio by scales height to width(6). Then *1.2 because it still wasn't perfect.
  # attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]] # longitude, x, width, bind as 6
  # attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]] # latitude, y, height
  autoheight <- (6 / (attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]])) * (attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]]) * 1.2
  
  ggmap(myMap) +
    geom_sf(data = st_contour(x = x,
                              contour_lines = TRUE, # makes lines not polys regardless of T or F
                              breaks = max(x[[1]], na.rm = TRUE) * 0.05
    ) %>%
      # breaks could be function param, but only allows 2 breaks. Whatevs ####
    st_transform(3857), # Vector transform after st_contour()  4326
    fill = NA, inherit.aes = FALSE,
    aes(colour = "UD_95_pct")) + # https://github.com/dkahle/ggmap/issues/160#issuecomment-966812818
    geom_sf(data = st_contour(x = x,
                              contour_lines = TRUE,
                              breaks = max(x[[1]], na.rm = TRUE) * 0.5
    ) %>%
      st_transform(3857), fill = NA, inherit.aes = FALSE, aes(colour = "UD_50_pct")) +
    scale_colour_manual(name = legendtitle, values = c(UD_95_pct = contour1colour, UD_50_pct = contour2colour)) +
    # https://stackoverflow.com/questions/64425970/ggmap-in-r-keep-google-copyright-information-on-cropped-map
    # scale_x_continuous(limits = c(myLocation[1], myLocation[3]), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(myLocation[2], myLocation[4]), expand = c(0, 0)) +
    #   Coordinate system already present. Adding new coordinate system, which will replace the existing one.
    # Scale for 'x' is already present. Adding another scale for 'x', which will replace the existing scale.
    # Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
    # Warning message:
    #   Removed 1 rows containing missing values (geom_rect). 
    ggtitle(plottitle, subtitle = plotsubtitle) +
    labs(x = axisxlabel, y = axisylabel, caption = plotcaption) +
    theme_minimal() +
    theme(legend.position = legend.position, #%dist (of middle? of legend box) from L to R, %dist from Bot to Top
          legend.spacing.x = unit(0, 'cm'), #compress spacing between legend items, this is min
          legend.spacing.y = unit(0, 'cm'), #compress spacing between legend items, this is min
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.background = element_rect(fill = "white", colour = NA), # element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50"), # white background
          plot.background = element_rect(fill = "white", colour = "grey50"), # white background
          legend.key = element_blank()) # removed whitespace buffer around legend boxes which is nice
  ggsave(filename = filesavename, plot = last_plot(), device = "png", path = savedir, scale = 1,
         #changes how big lines & legend items & axes & titles are relative to basemap. Smaller number = bigger items
         width = 6, height = autoheight, units = "in", dpi = 600, limitsize = TRUE)
} # close function