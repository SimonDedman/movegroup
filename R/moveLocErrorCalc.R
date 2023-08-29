#' moveLocError calculator for ARGOS or state space models resulting in 95percent latlon confidence intervals
#'
#' Builds a dataframe of original locations plus rowmeans of mean distance of location extremities 
#' lon975, lat; lon025, lat; lon, lat975; lon, lat025 from the centre point lon, lat.
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' 
#' @importFrom dplyr select
#' @importFrom purrr map_df
#' @importFrom rlang set_names
#' @importFrom sf st_as_sf st_set_crs st_transform st_distance
#' @importFrom tidyselect everything
#' 
#' @export moveLocErrorCalc
#' 
#' @param x Data frame or tibble with lats and lons and their high and low confidence interval 
#' counterparts.
#' @param loncol Name of longitude column in x, character. Default "lon".
#' @param latcol Name of latitude column in x, character. Default "lat".
#' @param latloncrs CRS of x, default 4326, numeric.
#' @param projectedcrs CRS to project to, should match your region, default 32617, numeric. See 
#' movegroup projectedCRS.
#' @param lon025 Name of low 2.5% confidence interval longitude column in x, character. Default 
#' "lon025".
#' @param lon975 Name of high 2.5% confidence interval longitude column in x, character. Default 
#' "lon975".
#' @param lat025 Name of low 2.5% confidence interval latitude column in x, character. Default 
#' "lat025".
#' @param lat975 Name of high 2.5% confidence interval latitude column in x, character. Default 
#' "lat975".
#' 
#' @return Dataframe of original locations plus rowmeans of mean distance of location extremities, 
#' for use in movegroup::movegroup(moveLocError).
#'
#' @details Use on your data object from movegroup::movegroup(data).
#' 
#' @examples
#' \dontrun{
#' data(argosFiltered)
#' myMoveLocError <- moveLocErrorCalc(argosFiltered)
#' }
#'

# moveLocError calculator for ARGOS / state space models resulting in 95% latlon confidence intervals
# Simon Dedman, simondedman@gmail.com, simondedman.com, 2023-07-26
moveLocErrorCalc <- function(x,
                             loncol = "lon",
                             latcol = "lat",
                             latloncrs = 4326,
                             projectedcrs = 32617,
                             lon025 = "lon025",
                             lon975 = "lon975",
                             lat025 = "lat025",
                             lat975 = "lat975"
) { # open moveLocErrorCalc function
  
  # build reproject function for later use
  reproject <- function(x,
                        loncol = loncol,
                        latcol = latcol,
                        latloncrs = latloncrs,
                        projectedcrs = projectedcrs) {
    x <- sf::st_as_sf(x, coords = c(loncol, latcol)) |>
      sf::st_set_crs(latloncrs) |> # latlon degrees sf object
      sf::st_transform(projectedcrs) |> # eastings northings units metres
      dplyr::select(-tidyselect::everything()) # remove all columns. Geometry is protected and retained
    return(x)
  }
  
  tracksfmean <- reproject(x = x,
                           loncol = loncol,
                           latcol = latcol,
                           latloncrs = latloncrs,
                           projectedcrs = projectedcrs)
  
  meanMoveLocDist <- list(
    c(x[,loncol], x[,lat975]), # U # were originally c(loncol, lat975), check format is right, 
    c(x[,lon975], x[,latcol]), # R
    c(x[,loncol], x[,lat025]), # D
    c(x[,lon025], x[,latcol]) # L
  ) |>
    lapply(function(x) reproject(x = x,
                                 loncol = x[1],
                                 latcol = x[2],
                                 latloncrs = latloncrs,
                                 projectedcrs = projectedcrs
    )) |>
    rlang::set_names(c("U", "R", "D", "L")) |> # set names of list elements
    lapply(
      function(vertextrack) { # distance from vertices to centre
        sf::st_distance(
          x = tracksfmean,
          y = vertextrack,
          by_element = TRUE
        )
      }
    ) |>
    purrr::map_df(~.x) |> # collapse list to df of 4 columns
    rowMeans()# make row means
  
  rm(tracksfmean)
  return(meanMoveLocDist)
} # close moveLocErrorCalc function