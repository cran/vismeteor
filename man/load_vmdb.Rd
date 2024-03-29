% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{load_vmdb}
\alias{load_vmdb}
\alias{load_vmdb_rates}
\alias{load_vmdb_magnitudes}
\title{Loading visual meteor observations from the data base}
\usage{
load_vmdb_rates(
  dbcon,
  shower = NULL,
  period = NULL,
  sl = NULL,
  lim.magn = NULL,
  sun.alt.max = NULL,
  moon.alt.max = NULL,
  session.id = NULL,
  rate.id = NULL,
  withSessions = FALSE,
  withMagnitudes = FALSE
)

load_vmdb_magnitudes(
  dbcon,
  shower = NULL,
  period = NULL,
  sl = NULL,
  lim.magn = NULL,
  session.id = NULL,
  magn.id = NULL,
  withSessions = FALSE,
  withMagnitudes = TRUE
)
}
\arguments{
\item{dbcon}{database connection.}

\item{shower}{character; selects by meteor shower codes.
\code{NA} loads sporadic meteors.}

\item{period}{time; selects a time range by minimum/maximum.}

\item{sl}{numeric; selects a range of solar longitudes by minimum/maximum.}

\item{lim.magn}{numeric; selects a range of limiting magnitudes by minimum/maximum.}

\item{sun.alt.max}{numeric; selects the maximum altitude of the sun.}

\item{moon.alt.max}{numeric; selects the maximum altitude of the moon.}

\item{session.id}{integer; selects by session ids.}

\item{rate.id}{integer; selects rate observations by ids.}

\item{withSessions}{logical; if \code{TRUE}, also load the corresponding session data.}

\item{withMagnitudes}{logical; if \code{TRUE}, also load the corresponding magnitude observations.}

\item{magn.id}{integer; selects magnitude observations by ids.}
}
\value{
Both functions return a list, with

\tabular{ll}{
\code{observations} \tab data frame, rate or magnitude observations,\cr
\code{sessions} \tab data frame; session data of observations,\cr
\code{magnitudes} \tab table; contingency table of meteor magnitude frequencies.
}

\code{observations} depends on the function call. \code{load_vmdb_rates} returns a data frame, with

\tabular{ll}{
\code{rate.id} \tab unique identifier of the rate observation,\cr
\code{shower.code} \tab IAU code of the shower. It is \code{NA} in case of sporadic meteors.\cr
\code{period.start} \tab start of observation,\cr
\code{period.end} \tab end of observation,\cr
\code{sl.start} \tab solarlong at start of observation,\cr
\code{sl.end} \tab solarlong at start of observation,\cr
\code{session.id} \tab reference to the session,\cr
\code{freq} \tab count of observed meteors,\cr
\code{lim.magn} \tab limiting magnitude,\cr
\code{t.eff} \tab net observed time in hours,\cr
\code{f} \tab correction factor of cloud cover,\cr
\code{time.sidereal} \tab sidereal time,\cr
\code{sun.alt} \tab altitude of the sun,\cr
\code{sun.az} \tab azimuth of the sun,\cr
\code{moon.alt} \tab altitude of the moon,\cr
\code{moon.az} \tab azimuth of the moon,\cr
\code{moon.illum} \tab illumination of the moon (\verb{0.0 .. 1.0}),\cr
\code{field.alt} \tab altitude of the field of view (optional),\cr
\code{field.az} \tab azimuth of the field of view (optional),\cr
\code{radiant.alt} \tab altitude of the radiant (optional). The zenith attraction is already applied.\cr
\code{radiant.az} \tab azimuth of the radiant (optional),\cr
\code{magn.id} \tab reference to the magnitude observations (optional).
}

\code{load_vmdb_magnitudes} returns a \code{observations} data frame, with

\tabular{ll}{
\code{magn.id} \tab unique identifier of the magnitude observation,\cr
\code{shower.code} \tab IAU code of the shower. It is \code{NA} in case of sporadic meteors.\cr
\code{period.start} \tab start of observation,\cr
\code{period.end} \tab end of observation,\cr
\code{sl.start} \tab solarlong at start of observation,\cr
\code{sl.end} \tab solarlong at start of observation,\cr
\code{session.id} \tab reference to the session,\cr
\code{freq} \tab count of observed meteors,\cr
\code{magn.mean} \tab mean of magnitudes,\cr
\code{lim.magn} \tab limiting magnitude (optional).
}

The \code{sessions} data frame contains

\tabular{ll}{
\code{session.id} \tab unique identifier of the session,\cr
\code{longitude} \tab location’s longitude,\cr
\code{latitude} \tab location’s latitude,\cr
\code{elevation} \tab height above mean sea level in km,\cr
\code{country} \tab country name,\cr
\code{location.name} \tab location name,\cr
\code{observer.id} \tab observer id (optional),\cr
\code{observer.name} \tab observer name (optional).
}

\code{magnitudes} is a contingency table of meteor magnitude frequencies.
The row names refer to the id of magnitude observations.
The column names refer to the magnitude.
}
\description{
Loads the data of visual meteor observations from a data base created with
\href{https://pypi.org/project/imo-vmdb/}{imo-vmdb}.
}
\details{
\code{sl}, \code{period} and \code{lim.magn} expect a vector with successive minimum and maximum values.
\code{sun.alt.max} and \code{moon.alt.max} are expected to be scalar values.
}
\note{
Angle values are expected and returned in degrees.
}
\examples{
\dontrun{
# create a connection to the data base
con <- dbConnect(
    PostgreSQL(),
    dbname = "vmdb",
    host = "localhost",
    user = "vmdb"
)

# load rate observations including
# session data and magnitude observations
data <- load_vmdb_rates(
    con,
    shower = 'PER',
    sl = c(135.5, 145.5),
    period = c('2015-08-01', '2015-08-31'),
    lim.magn = c(5.3, 6.7),
    withMagnitudes = TRUE,
    withSessions = TRUE
)

# load magnitude observations including
# session data and magnitude observations
data <- load_vmdb_magnitudes(
    con,
    shower = 'PER',
    sl = c(135.5, 145.5),
    period = c('2015-08-01', '2015-08-31'),
    lim.magn = c(5.3, 6.7),
    withMagnitudes = TRUE,
    withSessions = TRUE
)
}
}
\references{
\url{https://pypi.org/project/imo-vmdb/}
}
