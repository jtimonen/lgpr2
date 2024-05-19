#' Simulated longitudinal test data
#'
#'
#' @format
#' A data frame with 176 rows and 7 columns:
#' \describe{
#'   \item{id}{Subject id, factor}
#'   \item{arm}{Trial arm, factor}
#'   \item{sex}{Sex, factor}
#'   \item{time}{Measurement time, numeric}
#'   \item{weight}{Weight, numeric}
#'   \item{y}{Response variable, numeric}
#'   \item{f_true}{True signal}
#' }
#' @family built-in datasets
"testdata"

# Load weather data in correct format
load_weather_data <- function() {
  library(fda)
  temp <- CanadianWeather$dailyAv[, , 1]
  precip <- CanadianWeather$dailyAv[, , 3]
  n_days <- nrow(temp)
  n_stations <- ncol(temp)
  regions <- rep(as.factor(CanadianWeather$region), each = n_days)
  stations <- rep(as.factor(CanadianWeather$place), each = n_days)
  day <- rep(c(1:n_days), times = n_stations)
  dat <- data.frame(day, as.vector(temp), as.vector(precip), stations, regions)
  colnames(dat) <- c("day", "temperature", "precipitation", "station", "region")
  dat$day <- as.numeric(dat$day)
  return(dat)
}
