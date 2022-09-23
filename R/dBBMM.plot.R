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
#' @import ggplot2
#' @import ggmap
#' @import magrittr
#' @importFrom stars read_stars st_raster_type st_contour
#' @importFrom lubridate today
#' @importFrom sf st_set_crs st_bbox st_transform st_as_sfc st_as_sf
#' @importFrom starsExtra trim2
#' @importFrom fishtrack3d volumeUD contourPoly
#' @export
#' 
#' @param x Path to scaled data.
#' @param crsloc Location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
#' @param trim Remove NA & 0 values and crop to remaining date extents? Default TRUE.
#' @param myLocation Location for extents, format c(xmin, ymin, xmax, ymax). Default NULL, extents autocreated from data.
#' @param googlemap If pulling basemap from Google maps, this sets expansion factors since Google Maps tiling zoom setup doesn't align to myLocation extents.
#' @param gmapsAPI Enter your google maps API here, quoted character string.
#' @param expandfactor Extents expansion factor for basemap. 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out). 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
#' @param mapzoom Google: 3 (continent) - 21 (building). stamen: 0-18.
#' @param mapsource Source for ggmap::get_map; uses Stamen as fallback if no Goole Maps API present.
#' @param maptype Type of map for ggmap::get_map.
#' @param contour1colour Colour for contour 1, typically 95pct.
#' @param contour2colour Colour for contour 2, typically 50pct.
#' @param plottitle Can use the term 'home range' when an animal can be detected wherever it goes
#' i.e. using GPS, satellite or acoustic telemetry whereby it is known that acoustic receivers cover 
#' the entire home range of the study species. This term is problematic when applied to a passive 
#' acoustic telemetry setting where an array of non-overlapping receivers are used to assess local 
#' space use patterns i.e. the home range is bigger than the coverage by the acoustic array; put in Details.
#' @param plotsubtitle Plot subtitle.
#' @param legendtitle Legend title.
#' @param plotcaption Plot caption.
#' @param axisxlabel Longitude.
#' @param axisylabel Latitude.
#' @param legendposition Vector of 2, format c(1,2), Percent distance (of middle? of legend box) from L to R, percent distance from Bottom to Top, values 0 to 1.
#' @param fontsize Font size, default 12.
#' @param fontfamily = Font family, default "Times New Roman".
#' @param filesavename File savename.
#' @param savedir Save outputs to a temporary directory (default) else. Change to current directory e.g. "/home/me/folder". Do not use getwd() here. No terminal slash.
#' @param receiverlats Vector of latitudes for receivers to be plotted.
#' @param receiverlons Vector of longitudes for receivers to be plotted.
#' @param receivernames Vector of names for receivers to be plotted.
#' @param receiverrange Single (will be recycled), or vector of detection ranges in metres for receivers to be plotted.
#' @param recpointscol Colour of receiver centrepoint outlines.
#' @param recpointsfill Colour of receiver centrepoint fills.
#' @param recpointsalpha Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
#' @param recpointssize Size of receiver points.
#' @param recpointsshape Shape of receiver points, default 21, circle with outline and fill.
#' @param recbufcol Colour of the receiver buffer circle outlines.
#' @param recbuffill Colour of the receiver buffer circle fills.
#' @param recbufalpha Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
#' @param reclabcol Receiver label text colour.
#' @param reclabfill Receiver label fill colour, NA for no fill.
#' @param reclabnudgex Receiver label offset nudge in X dimension.
#' @param reclabnudgey Receiver label offset nudge in Y dimension.
#' @param reclabpad Receiver label padding in lines.
#' @param reclabrad Receiver label radius in lines.
#' @param reclabbord Receiver label border in mm.
#' @param surface Plot complete UD surface as well as contours.

dBBMMplot <- function(
    x = paste0("Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), # path to scaled data
    # dataCRS = 2958, # one of (i) character: a string accepted by GDAL, (ii) integer, a valid EPSG value (numeric), or (iii) an object of class crs.
    crsloc = NULL, # Location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
    trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE
    myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
    # Default NULL, extents autocreated from data.
    # c(-79.3, 25.68331, -79.24, 25.78)
    googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion
    # factors since Google Maps tiling zoom setup doesn't align to myLocation
    # extents.
    gmapsAPI = NULL, # enter your Google maps API here, quoted character string
    expandfactor = 1.6, # extents expansion factor for basemap.
    # 1.3 to 1.5 are the same zoom as 1. 1.6 is a big leap up in zoom (out).
    # 1.9 & maybe 1.7 or 1.8 is another step out. Ignored if not using Google Maps.
    mapzoom = 5, # google: 3 (continent) - 21 (building). stamen: 0-18
    mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present.
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
    plotcaption = paste0("dBBMM_HomeRange, ", lubridate::today()),
    axisxlabel = "Longitude",
    axisylabel = "Latitude",
    legendposition = c(0.11, 0.85), # Percent distance (of middle? of legend box) from L to R, percent distance from Bottom to Top.
    fontsize = 12,
    fontfamily = "Times New Roman",
    filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
    savedir = tempdir(), # file.path(work.dir, out.dir, "Scaled")
    receiverlats = NULL, # vector of latitudes for receivers to be plotted
    receiverlons = NULL, # vector of longitudes for receivers to be plotted
    receivernames = NULL, # vector of names for receivers to be plotted
    receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
    recpointscol = "black", # Colour of receiver centrepoint outlines.
    recpointsfill = "white", # Colour of receiver centrepoint fills.
    recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
    recpointssize = 1, # Size of receiver points.
    recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
    recbufcol = "black", # Colour of the receiver buffer circle outlines.
    recbuffill = "red", # Colour of the receiver buffer circle fills.
    recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
    reclabcol = "black", # Receiver label text colour.
    reclabfill = NA, # Receiver label fill colour, NA for no fill.
    reclabnudgex = 0, # Receiver label offset nudge in X dimension.
    reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
    reclabpad = 0, # Receiver label padding in lines.
    reclabrad = 0.15, # Receiver label radius in lines.
    reclabbord = 0, # Receiver label border in mm.
    surface = TRUE # Plot complete UD surface as well as contours
) {
  
  # ToDo
  # expandfactor could default to NULL and have a formula to set it based on the size of extents
  # gbm.basemap
  # gMaps API tutorial
  # trim section optional, depends magrittr
  # 50 & 95% breaks could be editable as function params but will be a bit of work
  
  # check receiver inputs are the correct lengths, if present.
  if (!is.null(receiverlats) & !is.null(receiverlons)) if (length(receiverlats) != length(receiverlons)) stop("length of receiverlats must equal length of receiverlons")
  if (!is.null(receiverlats) & !is.null(receivernames)) if (length(receiverlats) != length(receivernames)) stop("length of receivernames must equal length of receiverlats/lons")
  if (!is.null(receiverlats) & !is.null(receiverrange)) if (length(receiverrange) != length(receiverlons)) if (length(receiverrange) != 1) stop("length of receiverrange must equal length of receiverlats/lons, or 1")
  
  # raster version####
  # x <- raster::raster(x) # raster not in imports####
  # # x min max +- 29000
  # dataCRS <- readRDS(paste0(crsloc, "CRS.Rds")) # should be 3857
  # # dataCRS@projargs: "+proj=aeqd +lat_0=25.6871 +lon_0=-79.29617 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  # raster::crs(x) <- dataCRS
  
  # fishtrack3d version####
  # # use fishtrack3d function to convert to volumeUD raster
  # x <- fishtrack3d::volumeUD(x)
  # # Warning in install.packages : package ‘fishtrack3d’ is not available for this version of R (4.2)####
  # # need this until Edzer fix. Needed for MyLocation & contours
  # 
  # # create 50% and 95% contour spdf
  # spdf.50 <- fishtrack3d::contourPoly(x, levels = c(0.5))
  # sf_50 <- sf::st_as_sf(spdf.50) %>%  # convert to sf
  #   sf::st_transform(3857) # set CRS.
  # # 3857 = WGS 84 / Pseudo-Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI, World between 85.06°S and 85.06°N.
  # # Cartesian 2D CS. Axes: easting, northing (X,Y). Orientations: east, north. UoM: m.
  # 
  # spdf.95 <- fishtrack3d::contourPoly(x, levels = c(0.95))
  # sf_95 <- sf::st_as_sf(spdf.95) %>%  # convert to sf
  #   sf::st_transform(3857) # set CRS
  
  # stars version ####
  # Import raster
  x <- stars::read_stars(x)
  # %>% sf::st_set_crs(4326) # 4326 2958
  st_crs(x) <- proj4string(dataCRS) # set CRS
  
  # for plotting the surface UD on the map:
  # surfaceUD <- x %>% sf::st_set_crs(3857) # 4326 = WGS84. Ellipsoidal 2D CS. Axes: latitude, longitude. Orientations: north, east. UoM: degree
  # %>% sf::st_transform(3857)
  
  
  # if (stars::st_raster_type(x) == "curvilinear") stop(print("x is curvilinear; first reproject to planar"))
  # Warning message: object ‘is_curvilinear’ is not exported by 'namespace:stars'
  # https://github.com/r-spatial/stars/issues/464
  
  
  # Trim data to plot surface
  y <- x # make dupe object else removing all data < 0.05 means the 0.05 contour doesn't work in ggplot
  if (trim) { # trim raster extent to data?
    is.na(y[[1]]) <- y[[1]] == 0 # replace char pattern (0) in whole df/tbl with NA
    is.na(y[[1]]) <- y[[1]] < (max(y[[1]], na.rm = TRUE) * 0.05) # replace anything < 95% contour with NA since it won't be drawn
  }
  y %<>% starsExtra::trim2() # remove NA columns, which were all zero columns. This changes the bbox accordingly
  y[[1]] <- (y[[1]] / max(y[[1]], na.rm = TRUE)) * 100 # convert from raw values to 0:100 scale so legend is 0:100%
  
  # if (is.null(myLocation)) myLocation <- y %>% sf::st_transform(4326) %>% sf::st_bbox() %>% as.vector() # trimmer raster extents are smaller than UD95 contours
  # if (is.null(myLocation)) myLocation <- sf_95 %>% sf::st_transform(4326) %>% sf::st_bbox() %>% as.vector() # not using sf_95 any more
  if (is.null(myLocation)) myLocation <- stars::st_contour(x = x, contour_lines = TRUE, breaks = max(x[[1]], na.rm = TRUE) * 0.05) %>%
    sf::st_transform(4326) %>% sf::st_bbox() %>% as.vector()
  
  
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
  
  myMap <- ggmap::get_map(
    location = myLocation, # -62.57564  28.64368  33.78889  63.68533 # stamen etc want a bounding box
    zoom = mapzoom, # 3 (continent) - 21 (building). Stamen: 0-18
    # scale = "auto", # default "auto", 1, 2, 4 all the same
    source = mapsource, # "google" # using stamen as fallback
    maptype = maptype, # "satellite"
    messaging = TRUE,
    crop = TRUE # google maps crs = 4326
  ) 
  
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
    bbox_3857 <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(map_bbox, crs = 4326)), 3857))
    # Overwrite the bbox of the ggmap object with the transformed coordinates 
    attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
    attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
    attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
    attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
    map
  }
  myMap <- ggmap_bbox(myMap) # Use the function. Resulting map is CRS 3857
  
  # Automate width * height adjustments for different map extent / ratio
  # 6 (manually chosen width, below), divided by width range times by height range
  # Maintains ratio by scales height to width(6). Then *1.2 because it still wasn't perfect.
  # attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]] # longitude, x, width, bind as 6
  # attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]] # latitude, y, height
  autoheight <- (6 / (attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]])) * (attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]]) * 1.2
  
  # Create receiver objects
  if (!is.null(receiverlats) & !is.null(receiverlons)) {
    receiver <- data.frame(lon = receiverlons,
                           lat = receiverlats)
    receiver <- sf::st_as_sf(receiver, coords = c("lon","lat")) %>%
      sf::st_set_crs(4326) %>%
      sf::st_transform(3857)
    if (!is.null(receivernames)) {
      receiver <- cbind(receiver, receivernames)
    }
    if (!is.null(receiverrange)) {
      receiver <- cbind(receiver, receiverrange)
    }
  }
  
  # plot map ####
  ggmap::ggmap(myMap) + # basemap CRS = 3857
    
    # UD surface 
    {if (surface)
      stars::geom_stars(data = y %>% sf::st_transform(3857), inherit.aes = FALSE)
    } +
    
    # receiver centrepoints
    {if (!is.null(receiverlats) & !is.null(receiverlons))
      ggplot2::geom_sf(data = receiver %>%
                         sf::st_transform(3857), # Vector transform after st_contour
                       # already 3857 above so converting twice but it ain't broke
                       colour = recpointscol,
                       fill = recpointsfill,
                       alpha = recpointsalpha,
                       size = recpointssize,
                       shape = recpointsshape,
                       inherit.aes = FALSE,
      )
    } +
    
    # receiver buffer circles
    {if (!is.null(receiverlats) & !is.null(receiverlons) & !is.null(receiverrange))
      ggplot2::geom_sf(data = sf::st_buffer(receiver, dist = receiverrange) %>%
                         sf::st_transform(3857), # Vector transform after st_contour
                       # already 3857 above so converting twice but it ain't broke
                       colour = recbufcol,
                       fill = recbuffill,
                       alpha = recbufalpha,
                       inherit.aes = FALSE
      )
    } +
    
    # receiver labels
    {if (!is.null(receiverlats) & !is.null(receiverlons) & !is.null(receivernames))
      ggplot2::geom_sf_label(data = receiver %>%
                               sf::st_transform(3857), # Vector transform after st_contour
                             # already 3857 above so converting twice but it ain't broke
                             colour = reclabcol,
                             fill = reclabfill,
                             inherit.aes = FALSE,
                             nudge_x = reclabnudgex,
                             nudge_y = reclabnudgey,
                             label.padding = unit(reclabpad, "lines"), # 0.25
                             label.r = unit(reclabrad, "lines"),
                             label.size = reclabbord, # 0.25
                             ggplot2::aes(label = receivernames)
      )
    } +
    
    # 95% UD
    # ggplot2::geom_sf(data = sf_95 %>%
    #                    sf::st_transform(3857), # Vector transform after st_contour
    #                  # already 3857 above so converting twice but it ain't broke
    #                  fill = NA, inherit.aes = FALSE,
    #                  ggplot2::aes(colour = "95% UD")) + # https://github.com/dkahle/ggmap/issues/160#issuecomment-966812818
    
    ggplot2::geom_sf(data = stars::st_contour(x = x,
                                              contour_lines = TRUE,
                                              breaks = max(x[[1]],
                                                           na.rm = TRUE) * 0.05),
                     fill = NA,
                     inherit.aes = FALSE,
                     ggplot2::aes(colour = "95% UD")) +
    
    # 50% UD
    # ggplot2::geom_sf(data = sf_50 %>%
    #                    sf::st_transform(3857),
    #                  # already 3857 above so converting twice but it ain't broke
    #                  fill = NA, inherit.aes = FALSE,
    #                  ggplot2::aes(colour = "50% UD")) +
    
    ggplot2::geom_sf(data = stars::st_contour(x = x,
                                              contour_lines = TRUE,
                                              breaks = max(x[[1]],
                                                           na.rm = TRUE) * 0.5),
                     fill = NA,
                     inherit.aes = FALSE,
                     ggplot2::aes(colour = "50% UD")) +
    
    # UD surface colours
  viridis::scale_fill_viridis(
    alpha = 1, # 0:1
    begin = 0, # hue
    end = 1, # hue
    direction = 1, # colour order, 1 or -1
    discrete = FALSE, # false = continuous
    option = "D", # A magma B inferno C plasma D viridis E cividis F rocket G mako H turbo
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    # name = waiver(),
    name = "UD%",
    # limits = NA,
    # position = "left"
    position = "right"
  ) +
  
  # UD contour colours
  ggplot2::scale_colour_manual(name = legendtitle, values = c("50% UD" = contour2colour, "95% UD" = contour1colour)) +
    # https://stackoverflow.com/questions/64425970/ggmap-in-r-keep-google-copyright-information-on-cropped-map
    # scale_x_continuous(limits = c(myLocation[1], myLocation[3]), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(myLocation[2], myLocation[4]), expand = c(0, 0)) +
    #   Coordinate system already present. Adding new coordinate system, which will replace the existing one.
    # Scale for 'x' is already present. Adding another scale for 'x', which will replace the existing scale.
    # Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
    # Warning message: Removed 1 rows containing missing values (geom_rect). 
    
    ggplot2::ggtitle(plottitle, subtitle = plotsubtitle) +
    ggplot2::labs(x = axisxlabel, y = axisylabel, caption = plotcaption) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = legendposition, #%dist (of middle? of legend box) from L to R, %dist from Bot to Top
      legend.spacing.x = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
      legend.spacing.y = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA), # element_blank(),
      panel.background = ggplot2::element_rect(fill = "#99b3cc", colour = "grey50"), # Stamen ocean blue background
      plot.background = ggplot2::element_rect(fill = "white", colour = "grey50"), # white background
      legend.key = ggplot2::element_blank(), 
      text = ggplot2::element_text(size = fontsize,  family = fontfamily)
    ) # removed whitespace buffer around legend boxes which is nice
  
  ggplot2::ggsave(filename = filesavename, plot = ggplot2::last_plot(), device = "png", path = savedir, scale = 1,
                  #changes how big lines & legend items & axes & titles are relative to basemap. Smaller number = bigger items
                  width = 6, height = autoheight, units = "in", dpi = 600, limitsize = TRUE)
} # close function