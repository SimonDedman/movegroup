#' Data: Tracks of lemon sharks off Bimini, Bahamas
#'
#' Tracks of 17 lemon sharks (Negaprion brevirostris) tagged off Bimini, Bahamas, 2012:2014, by 
#' Bimini Biological Field Station employees and volunteers, with accompanying tidal phase.
#'
#' @format A data frame with 1308 rows and 5 variables:
#' \describe{
#'   \item{Datetime}{POSIXct datetime, format YYYY-MM-DD HH:MM:SS.}
#'   \item{Shark}{Individual shark ID code.}
#'   \item{T.Ph}{Tidal phase, H M L High Medium Low.}
#'   \item{Lat}{Decimal latitudes.}
#'   \item{Lon}{Decimal longitudes.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name TracksCleaned
#' @usage data(TracksCleaned)
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Maurits van Zinnicq Bergmann, \email{mauritsvzb@@gmail.com}
#' @source \url{https://www.biminisharklab.com}
"TracksCleaned"