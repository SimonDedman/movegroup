# Now that we have constructed individual movement models for complete movement trajectories, below we will do the same, but now we will discern different tidal phases when constructing the models.

# ```{r dir}
# set path to general folders where output files will be written
if (!dir.exists(file.path(data.dir, out.dir, "Tide"))) dir.create(file.path(data.dir, "dBBMM ASCII/Tide"))
out.dir_h <- file.path(data.dir, "dBBMM ASCII/Tide/H")
if (!dir.exists(out.dir_h)) dir.create(out.dir_h)
out.dir_l <- file.path(data.dir, "dBBMM ASCII/Tide/L")
if (!dir.exists(out.dir_l)) dir.create(out.dir_l)
out.dir_m <- file.path(data.dir, "dBBMM ASCII/Tide/M")
if (!dir.exists(out.dir_m)) dir.create(out.dir_m)
# ```

# ```{r construct_tidal_dbbmm}
# Loop through the different tidal phases
for (i in unique(data$T.Ph)) { # "H" "M" "L"

  # Filter individual tide data
  data.t <- data[data$T.Ph == i, ]

  # remove sharks with insufficient movement data (identical to code L113-L128)
  check1 <- data.t %>%
    group_by(ID) %>%
    summarise(relocations = length(Datetime))
  check1
  check2 <- filter(check1, relocations > 23)
  check2

  if (length(check1$ID) != length(check2$ID)) {
    data.t <- semi_join(data.t, check2)
    check1 <- data.t %>%
      group_by(ID) %>%
      summarise(relocations = length(Datetime))
    check1
    check2 <- filter(check1, relocations > 23)
    length(check1$ID) == length(check2$ID)
  } # data.t: 468 x 7



  bb <- list()
  bb.list <- list()

  data.t %<>%
    tidyr::drop_na(ID) %>%
    group_by(ID) %>%
    distinct(Datetime, .keep_all = TRUE) %>% # distinct (grouped): removed one row (<1%), 1,286 rows remaining
    # prevents duplicate Datetime crash in move() later
    ungroup() # 1253 x 7 after removing as.numeric above

  # Loop through all unique tags
  counter <- 0
  for (j in unique(data.t$ID)) {
    # Print individual
    counter <- counter + 1
    print(paste0(
      "processing ", which(unique(data.t$ID) %in% j),
      " of ", length(unique(data.t$ID)),
      " tags, tide = ", i,
      " (", which(unique(data$T.Ph) %in% i),
      " of ", length(unique(data$T.Ph)), ")"
    ))

    # Filter individual data
    data.i <- data.t[data.t$ID == j, ]

    # create move object
    move.i <- move(
      x = data.i$Lon,
      y = data.i$Lat,
      time = data.i$Datetime,
      proj = CRS("+proj=longlat +datum=WGS84"),
      data = data.i,
      animal = data.i$ID,
      sensor = "VR2W"
    )
    # Warning: Unknown or uninitialised column: `citation`.
    # Warning: Setting row names on a tibble is deprecated.

    # Check the current projection
    proj4string(move.i) # "+proj=longlat +datum=WGS84 +no_defs"

    # Incorporate uncertainty in the model by including a location error.
    # From communication with Rob, the gps location of a shark is estimated,
    # from the boat coordinates, bearing and distance estimate, to be 1 m.
    move.i$LocationError <- 1 # 1 m: gps loc of shark inferred from boat coord + bearing + distance estimate

    # Convert projection to Azimuthal Equi-Distance projection (aeqd)
    r.i <- spTransform(move.i, center = TRUE)

    # Make sure it changed correctly
    proj4string(r.i) # "+proj=aeqd +lat_0=25.73054 +lon_0=-79.248155 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

    # Calculate time lag between consecutive detections. Obviously there is no time lag for the first detection, so we need to tell R this. Knowing time lags in acoustic telemetry is important, as the motion variance is based on the time difference b/w detections at consecutive locations i.e. larger time lag creates a wider bridge where the animal could've been and therefore potentially inflates the motion variance. Below code deals with this issue by identifying the location and number of detections with large time gaps. Below we will use an arbitrary value of 2 h.
    TimeDiff <- timeLag(r.i, units = "hours")
    r.i$TimeDiff <- append(0, TimeDiff)
    long <- which(r.i$TimeDiff > 2)

    # Make a raster for the UD to plot into. Start with UTM.
    # These coordinates need to be big enough to cover your data.
    # May need to expand x and y ranges (i.e. the "e" variable below) if you encounter errors when constructing the DBBMM in the next code chunk.
    # NOTE: "xUTM" grid should be larger than the "x" i.e. lower minimums and higher.
    # The variable 'e' influences the extent.
    # The resolution determines the grid size.
    # Finer resolution/grid size (in m) i.e. lower value means longer computation time. May need/want to play with this value too.
    # e <- 80 * 1000 # make e a function of range of extent as %

    # Error in .local(object, raster, location.error = location.error, ext = ext,  :
    # Lower x grid not large enough, consider extending the raster in that direction or enlarging the ext argument
    # make buffpct larger
    buffpct <- 30 # buffer extent as a %
    buffpct <- buffpct / 100
    xUTM.i <- raster(
      xmn = min(data$NewEastingUTM, na.rm = TRUE) - ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
      xmx = max(data$NewEastingUTM, na.rm = TRUE) + ((max(data$NewEastingUTM, na.rm = TRUE) - min(data$NewEastingUTM, na.rm = TRUE)) * buffpct),
      ymn = min(data$NewNorthingUTM, na.rm = TRUE) - ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
      ymx = max(data$NewNorthingUTM, na.rm = TRUE) + ((max(data$NewNorthingUTM, na.rm = TRUE) - min(data$NewNorthingUTM, na.rm = TRUE)) * buffpct),
      crs = CRS("+proj=utm +zone=17 +datum=WGS84"),
      resolution = 50
    ) # 50

    # We now need to reproject this into AEQD. Make a dummy object to get the correct projection.
    x.i <- raster(
      xmn = min(data$NewEastingUTM),
      xmx = max(data$NewEastingUTM),
      ymn = min(data$NewNorthingUTM),
      ymx = max(data$NewNorthingUTM),
      crs = CRS("+proj=utm +zone=17 +datum=WGS84")
    )
    proj4string(x.i) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

    # Use that to reproject UTM to AEQD
    # ?projectExtent # project the values of a Raster object to a new Raster object with another projection (CRS), i.e. projectExtent(from, to, ...))
    newTemplate <- projectExtent(xUTM.i, proj4string(x.i))
    rm(x.i)
    newTemplate # get the number of pixels in newTemplate and write that number into the next function. This is an empty raster.
    # class      : RasterLayer
    # dimensions : 142615, 187156, 26691252940  (nrow, ncol, ncell)
    # resolution : 50, 50  (x, y)
    # extent     : 593395.5, 9951195, 2545826, 9676576  (xmin, xmax, ymin, ymax)
    # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs

    # Give newTemplate some values. Make Rep equal to the ncell dimension
    ones <- rep(1, ncell(newTemplate))
    xAEQD.i <- setValues(newTemplate, ones)
    rm(newTemplate)
    rm(ones)
    xAEQD.i
    # class      : RasterLayer
    # dimensions : 3270, 3245, 10611150  (nrow, ncol, ncell)
    # resolution : 50, 50  (x, y)
    # extent     : 594830.7, 757080.7, 2763971, 2927471  (xmin, xmax, ymin, ymax)
    # crs        : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
    # source     : memory
    # names      : layer
    # values     : 1, 1  (min, max)

    # Reproject the move object r into AEQD
    rNew.i <- spTransform(r.i, proj4string(xAEQD.i))
    rNew.i
    # class       : Move
    # features    : 37
    # extent      : 674830.7, 677094.8, 2843987, 2847471  (xmin, xmax, ymin, ymax)
    # crs         : +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
    # variables   : 7
    # names       :   Datetime,      Lat,       Lon,    NewEastingUTM,   NewNorthingUTM, LocationError,         TimeDiff
    # min values  : 1409777400, 25.70323, -79.25719, 674830.719845252, 2843986.90740363,             1,                0
    # max values  : 1410148680, 25.73457, -79.23457, 677094.847481247, 2847470.56026671,             1, 31.1666666666667
    # timestamps  : 2014-09-03 16:50:00 ... 2014-09-07 23:58:00 Time difference of 4 days  (start ... end, duration)
    # sensors     : VR2W
    # indiv. data : ID, T.Ph
    # indiv. value: E07 L
    # date created: 2021-08-25 01:42:43

    # We need to make sure that the projection of the move object is in the same format as our raster. Check below.
    proj4string(xAEQD.i) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    proj4string(rNew.i) # "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
    proj4string(xAEQD.i) == proj4string(rNew.i) # TRUE

    # Below we exclude the data points that are 'long' i.e. create a large time gap. move::burst ?
    bursted <- burst(
      rNew.i,
      c("normal", "long")[1 + (timeLag(rNew.i, units = "hours") > 2)]
    )
    rm(rNew.i)
    # There are 2 types of burst: "normal" and "long". You can select for these in the dbbmm by selecting the factor level in the burstType command.


    # Construct the model. The time.step should reflect the ping frequency of the tag (in minutes)
    bursted_dbbmm <- brownian.bridge.dyn(bursted,
      burstType = "normal",
      raster = xAEQD.i,
      location.error = "LocationError",
      ext = 3,
      window.size = 23
    )
    rm(bursted)
    rm(xAEQD.i)

    # Warning: Some burst are omitted for variance calculation since there are not segements of interest

    # Re-standardize (Dr. Kranstauber's trouble shooting solution). Occurring errors are
    # "due to limits of accuracy during calculations. Given the error is really small (<.000001 %)
    # I would not worry to much and re-standardize".
    # Then aggregate UD segments.
    tmp <- calc(bursted_dbbmm, sum)
    rm(bursted_dbbmm)
    bb <- new(".UD", tmp / sum(values(tmp)))
    rm(tmp)

    # Export aggregated individual UDs to as ascii files for further exploration/spatial analysis in GIS software.
    # Needed for when the aim is to plot population-level and normalized UDs per species.
    asc <- writeRaster(bb,
      file.path(data.dir, out.dir, "Tide", i, filename = paste0(j, ".asc")),
      format = "ascii",
      datatype = "FLT4S",
      bylayer = T,
      overwrite = T
    )

    # Calculate volume area (m^2) within 50% (core) and 95% (general use) contours
    area.50 <- sum(values(getVolumeUD(bb) <= .50))
    area.95 <- sum(values(getVolumeUD(bb) <= .95))
    rm(bb)
    # Combine in single df
    area.ct <- data.frame(core.use = area.50, general.use = area.95) # why change names?

    # Add ID id
    area.ct$ID <- j

    # Put in list
    bb.list[[counter]] <- area.ct
    bb.list
    gc()
  }

  # Put everything in a data.frame
  md <- bind_rows(bb.list, .id = "column_label")

  write.csv(
    x = md,
    file = file.path(data.dir, out.dir, "Tide", paste0("Lemon_VolumeArea_", i, ".csv")),
    row.names = FALSE
  )

  # Scale and sum individual rasters, followed by a rescaling to obtain an aggregated UD (.asc)
  print(Sys.info()["nodename"])
  if (grepl("aurits", Sys.info()["nodename"])) {
    work.dir <- "~/Documents/Science/Projects/Rob Bullock - Bimini/dBBMM_HomeRange_output/"
  }
  if (grepl("nautilus", Sys.info()["nodename"]) | grepl("Poseidon", Sys.info()["nodename"]) | grepl("aquarius", Sys.info()["nodename"])) {
    work.dir <- "/home/simon/Dropbox/PostDoc Work/Rob Bullock accelerometer Lemons 2020.09/"
  }

  source("./R/scaleraster.R")
  # Warning in file(filename, "r", encoding = encoding) :  cannot open file './R/scaleraster.R': No such file or directory
  # Error in file(filename, "r", encoding = encoding) : cannot open the connection
  # Again, works fine in console
  scaleraster(path = file.path(data.dir, out.dir, "Tide", i))
}
# ```

# plotting tidal phases

# ```{r tidal phase plots}
# Move this to separate script & source it####
# Import raster
# x <- raster(paste0(work.dir, "/data/dBBMM ASCII/Scaled/All_Rasters_Scaled.asc")) # mo
# x <- raster(paste0(work.dir, "/dBBMM ASCII/Scaled/All_Rasters_Scaled.asc")) # si
# x
# plot(sqrt(x)) # sqrt inflates the values to make it easier to view
# plot(x) # fix this nor plotting why?? SD

# x <- read_stars(paste0(work.dir, "/dBBMM ASCII/Scaled/All_Rasters_Scaled.asc")) %>% # si directory
#   st_set_crs(2958) %>%
#   st_transform(4326)
# stars:::is_curvilinear(x) # TRUE

for (i in unique(data$T.Ph)) {
  x <- read_stars(file.path(data.dir, out.dir, "Tide", i, "Scaled/All_Rasters_Scaled.asc")) %>%
    # si directory
    st_set_crs(2958)
  # stars:::is_curvilinear(x) # FALSE

  # maybe don't do this trimming bit? ####
  # is.na(x$All_Rasters_Scaled.asc) <- x$All_Rasters_Scaled.asc == 0 # replace char pattern (0) in whole df/tbl with NA
  # x %<>% trim2() # remove NA columns, which were all zero columns. This changes the bbox accordingly

  # myLocation to be function param, NULL by default so extents autocreated from data, but gives user the option to set extents####
  myLocation <- NULL
  if (is.null(myLocation)) myLocation <- st_bbox(x %>% st_transform(4326)) %>% as.vector() # xmin      ymin      xmax      ymax
  myLocation <- c(-79.3, 25.68331, -79.24, 25.78)

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

  # set googlemap as top parameter when this is a function ####
  googlemap <- FALSE
  if (googlemap) { # grow bounds extents if requested
    expandfactor <- 1.6 # 1.3 - 1.5 same zoom as 1. 1.6 is a big leap up in zoom (out). 1.9 (& maybe 1.7 or 1.8 is another step out, done manually for #2, 200368, because some points were out of bounds.)
    # Still need to improve this part
    xmid <- mean(myLocation[c(1, 3)])
    ymid <- mean(myLocation[c(2, 4)])
    xmax <- ((myLocation[3] - xmid) * expandfactor) + xmid # updated for sf/st
    xmin <- xmid - ((xmid - myLocation[1]) * expandfactor)
    ymax <- ((myLocation[4] - ymid) * expandfactor) + ymid
    ymin <- ymid - ((ymid - myLocation[2]) * expandfactor)
    myLocation <- c(xmin, ymin, xmax, ymax)
  }

  myMap <- get_map(
    location = myLocation,
    source = "google", # using stamen as fallback
    maptype = "satellite",
    crop = FALSE
  ) # google maps crs = 4326

  # Define a function to fix the bbox to be in EPSG:3857
  # https://stackoverflow.com/a/50844502/1736291
  # Fixes "error no lon value" in ggmap below
  # Move this to separate script & source it####
  ggmap_bbox <- function(map) {
    if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
    # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector,
    # and set the names to what sf::st_bbox expects:
    map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
    # Coonvert the bbox to an sf polygon, transform it to 3857,
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


  # ggplot() +
  ggmap(myMap) +
    # geom_stars(data = x %>% st_transform(3857), inherit.aes = FALSE) + # 4326
    # scale_fill_viridis_c() + # for stars
    geom_sf(
      data = st_contour(
        x = x,
        contour_lines = TRUE, # makes lines not polys regardless of T or F
        breaks = c(0.05)
      ) %>%
        # breaks could be function param, but only allows 2 breaks. Whatevs ####
        st_transform(3857), # Vector transform after st_contour()  4326
      fill = NA, inherit.aes = FALSE,
      aes(colour = "UD_95_pct")
    ) + # https://github.com/dkahle/ggmap/issues/160#issuecomment-966812818
    geom_sf(data = st_contour(x = x, contour_lines = TRUE, breaks = c(0.5)) %>%
      st_transform(3857), fill = NA, inherit.aes = FALSE, aes(colour = "UD_50_pct")) +
    scale_colour_manual(name = "Percent UD Contours", values = c(UD_95_pct = "red", UD_50_pct = "orange")) +
    ggtitle("Aggregated 95% and 50% contours, lemon sharks, Bimini",
      subtitle = "Scaled contours. n = 13"
    ) +
    # make title & subtitle function params ####
    # make n a dynamic variable
    # Can use the term 'home range' when an animal can be detected wherever it goes i.e. using GPS, satellite or acoustic telemetry whereby it is known that acoustic receivers cover the entire home range of the study species. This term is problematic when applied to a passive acoustic telemetry setting where an array of non-overlapping receivers are used to assess local space use patterns i.e. the home range is bigger than the coverage by the acoustic array; put in Details ####
    # data %>% distinct(ID) %>% nrow() # 13
    labs(x = "Longitude", y = "Latitude", caption = paste0("dBBMM_HomeRange, ", today())) +
    theme_minimal() +
    theme(
      legend.position = c(0.16, 0.92), # %dist (of middle? of legend box) from L to R, %dist from Bot to Top
      # Make legend.position as function param ####
      legend.spacing.x = unit(0, "cm"), # compress spacing between legend items, this is min
      legend.spacing.y = unit(0, "cm"), # compress spacing between legend items, this is min
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = "white", colour = NA), # element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey50"), # white background
      plot.background = element_rect(fill = "white", colour = "grey50"), # white background
      legend.key = element_blank()
    ) # removed whitespace buffer around legend boxes which is nice

  ggsave(paste0(today(), "_dBBMM-contours.png"),
    plot = last_plot(), device = "png", path = file.path(work.dir, out.dir, "Tide", i, "Scaled"), scale = 1,
    # changes how big lines & legend items & axes & titles are relative to basemap. Smaller number = bigger items
    width = 6, height = autoheight, units = "in", dpi = 600, limitsize = TRUE
  )
} # close function
