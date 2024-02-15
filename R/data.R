#' Data: Tracks of lemon sharks off Bimini, Bahamas
#'
#' Tracks of 17 lemon sharks (Negaprion brevirostris) tagged off Bimini, Bahamas, 2012:2014, by
#'  Bimini Biological Field Station employees and volunteers, with accompanying tidal phase.
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

#' Data: Tracks of two great hammerhead sharks with position confidence intervals
#'
#' Tracks of 2 great hammerhead sharks tagged in Jupiter, and The Keys, Florida, USA, in 2022 and 
#' 2023 respectively, by Saving The Blue (savingtheblue.org), filtered by argosfilter::sdafilter and
#'  with state space model applied using aniMotum package, using scripts by Vital Heim, see 
#'  https://github.com/SimonDedman/SavingTheBlue/blob/main/R/06A_Filter_SPOT_data.R and 
#'  https://github.com/SimonDedman/SavingTheBlue/blob/main/R/06B_CTCRW_SPOT_data_usin_animotum.R .
#'
#' @format A data frame with 1492 rows and 8 variables:
#' \describe{
#'   \item{id}{Character, shark ID.}
#'   \item{date}{POSIXct datetime, format YYYY-MM-DD HH:MM:SS.}
#'   \item{lon}{Decimal longitudes.}
#'   \item{lon025}{Decimal longitudes, lower 95% confidence interval bound.}
#'   \item{lon975}{Decimal longitudes, upper 95% confidence interval bound.}
#'   \item{lat}{Decimal latitudes.}
#'   \item{lat025}{Decimal latitudes, lower 95% confidence interval bound.}
#'   \item{lat975}{Decimal latitudes, upper 95% confidence interval bound.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name argosFiltered
#' @usage data(argosFiltered)
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Vital Heim, \email{vital.heim@@gmail.com}
#' @source \url{https://www.savingtheblue.org}
"argosFiltered"
