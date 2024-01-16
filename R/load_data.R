#' @name load_vmdb
#' @aliases load_vmdb_rates
#' @aliases load_vmdb_magnitudes
#' @import DBI
#' @title Loading visual meteor observations from the data base
#' @description
#' Loads the data of visual meteor observations from a data base created with
#' \href{https://pypi.org/project/imo-vmdb/}{imo-vmdb}.
#' @note Angle values are expected and returned in degrees.
#' @param dbcon database connection.
#' @param shower character; selects by meteor shower codes.
#'   `NA` loads sporadic meteors.
#' @param period time; selects a time range by minimum/maximum.
#' @param sl numeric; selects a range of solar longitudes by minimum/maximum.
#' @param lim.magn numeric; selects a range of limiting magnitudes by minimum/maximum.
#' @param sun.alt.max numeric; selects the maximum altitude of the sun.
#' @param moon.alt.max numeric; selects the maximum altitude of the moon.
#' @param session.id integer; selects by session ids.
#' @param rate.id integer; selects rate observations by ids.
#' @param magn.id integer; selects magnitude observations by ids.
#' @param withSessions logical; if `TRUE`, also load the corresponding session data.
#' @param withMagnitudes logical; if `TRUE`, also load the corresponding magnitude observations.
#' @details
#' `sl`, `period` and `lim.magn` expect a vector with successive minimum and maximum values.
#' `sun.alt.max` and `moon.alt.max` are expected to be scalar values.
#' @return
#' Both functions return a list, with
#'
#' \tabular{ll}{
#'     `observations` \tab data frame, rate or magnitude observations,\cr
#'     `sessions` \tab data frame; session data of observations,\cr
#'     `magnitudes` \tab table; contingency table of meteor magnitude frequencies.
#' }
#'
#' `observations` depends on the function call. `load_vmdb_rates` returns a data frame, with
#'
#' \tabular{ll}{
#' `rate.id` \tab unique identifier of the rate observation,\cr
#' `shower.code` \tab IAU code of the shower. It is `NA` in case of sporadic meteors.\cr
#' `period.start` \tab start of observation,\cr
#' `period.end` \tab end of observation,\cr
#' `sl.start` \tab solarlong at start of observation,\cr
#' `sl.end` \tab solarlong at start of observation,\cr
#' `session.id` \tab reference to the session,\cr
#' `freq` \tab count of observed meteors,\cr
#' `lim.magn` \tab limiting magnitude,\cr
#' `t.eff` \tab net observed time in hours,\cr
#' `f` \tab correction factor of cloud cover,\cr
#' `time.sidereal` \tab sidereal time,\cr
#' `sun.alt` \tab altitude of the sun,\cr
#' `sun.az` \tab azimuth of the sun,\cr
#' `moon.alt` \tab altitude of the moon,\cr
#' `moon.az` \tab azimuth of the moon,\cr
#' `moon.illum` \tab illumination of the moon (`0.0 .. 1.0`),\cr
#' `field.alt` \tab altitude of the field of view (optional),\cr
#' `field.az` \tab azimuth of the field of view (optional),\cr
#' `radiant.alt` \tab altitude of the radiant (optional). The zenith attraction is already applied.\cr
#' `radiant.az` \tab azimuth of the radiant (optional),\cr
#' `magn.id` \tab reference to the magnitude observations (optional).
#' }
#'
#' `load_vmdb_magnitudes` returns a `observations` data frame, with
#'
#' \tabular{ll}{
#' `magn.id` \tab unique identifier of the magnitude observation,\cr
#' `shower.code` \tab IAU code of the shower. It is `NA` in case of sporadic meteors.\cr
#' `period.start` \tab start of observation,\cr
#' `period.end` \tab end of observation,\cr
#' `sl.start` \tab solarlong at start of observation,\cr
#' `sl.end` \tab solarlong at start of observation,\cr
#' `session.id` \tab reference to the session,\cr
#' `freq` \tab count of observed meteors,\cr
#' `magn.mean` \tab mean of magnitudes,\cr
#' `lim.magn` \tab limiting magnitude (optional).
#' }
#'
#' The `sessions` data frame contains
#'
#' \tabular{ll}{
#' `session.id` \tab unique identifier of the session,\cr
#' `longitude` \tab location’s longitude,\cr
#' `latitude` \tab location’s latitude,\cr
#' `elevation` \tab height above mean sea level in km,\cr
#' `country` \tab country name,\cr
#' `location.name` \tab location name,\cr
#' `observer.id` \tab observer id (optional),\cr
#' `observer.name` \tab observer name (optional).
#' }
#'
#' `magnitudes` is a contingency table of meteor magnitude frequencies.
#' The row names refer to the id of magnitude observations.
#' The column names refer to the magnitude.
#'
#' @references \url{https://pypi.org/project/imo-vmdb/}
#' @examples
#' \dontrun{
#' # create a connection to the data base
#' con <- dbConnect(
#'     PostgreSQL(),
#'     dbname = "vmdb",
#'     host = "localhost",
#'     user = "vmdb"
#' )
#'
#' # load rate observations including
#' # session data and magnitude observations
#' data <- load_vmdb_rates(
#'     con,
#'     shower = 'PER',
#'     sl = c(135.5, 145.5),
#'     period = c('2015-08-01', '2015-08-31'),
#'     lim.magn = c(5.3, 6.7),
#'     withMagnitudes = TRUE,
#'     withSessions = TRUE
#' )
#'
#' # load magnitude observations including
#' # session data and magnitude observations
#' data <- load_vmdb_magnitudes(
#'     con,
#'     shower = 'PER',
#'     sl = c(135.5, 145.5),
#'     period = c('2015-08-01', '2015-08-31'),
#'     lim.magn = c(5.3, 6.7),
#'     withMagnitudes = TRUE,
#'     withSessions = TRUE
#' )
#' }

#' @rdname load_vmdb
#' @export
load_vmdb_rates <- function(
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
) {
    shower.filter <- ''
    if (!is.null(shower)) {
        shower <- unique(shower)
        shower.is_na <- is.na(shower)
        with_spo <- length(shower[shower.is_na]) > 0
        shower <- shower[!shower.is_na]
        with_showers <- length(shower) > 0

        if (with_spo) {
            spo.filter <- "r.shower IS NULL"
        } else {
            spo.filter <- "FALSE"
        }

        if (with_showers) {
            shower.filter <- paste(shower, collapse="','")
            shower.filter <- paste0("r.shower IN ('", shower.filter , "')" )
        } else {
            shower.filter <- "FALSE"
        }

        if (with_showers | with_spo) {
            shower.filter <- paste(
                '(', shower.filter, 'OR', spo.filter , ') AND '
            )
        } else {
            shower.filter <- ''
        }
    }

    period.filter <- ''
    if (!is.null(period)) {
        period <- matrix(period, ncol=2)
        period.filter <- apply(period, 1, function(period){
            paste0(
                "r.period_start >= '", period[1] ,
                "' AND r.period_end <= '", period[2] , "'"
            )
        })
        period.filter <- paste('(', paste(period.filter, collapse=') OR ('), ') AND ')
    }

    sl.filter <- ''
    if (!is.null(sl)) {
        sl <- matrix(sl, ncol=2)
        sl.filter <- apply(sl, 1, function(sl){
            if (sl[1] > sl[2]) {
                paste0(
                    "r.sl_start BETWEEN ", sl[1], " AND 360.0 AND ",
                    "r.sl_end BETWEEN 0.0 AND ", sl[2]
                )
            } else {
                paste0("r.sl_start > ", sl[1], " AND r.sl_end <= ", sl[2])
            }
        })
        sl.filter <- paste('(', paste(sl.filter, collapse=') OR ('), ') AND ')
    }

    lim.magn.filter <- ''
    if (!is.null(lim.magn)) {
        lim.magn <- matrix(lim.magn, ncol=2)
        lim.magn.filter <- apply(lim.magn, 1, function(lim.magn){
            paste0("r.lim_mag BETWEEN ", lim.magn[1], " AND ", lim.magn[2])
        })
        lim.magn.filter <- paste('(', paste(lim.magn.filter, collapse=') OR ('), ') AND ')
    }

    sun_alt.filter <- ''
    if (!is.null(sun.alt.max)) {
        sun_alt.filter <- paste0("r.sun_alt <= " , sun.alt.max , " AND " )
    }

    moon_alt.filter <- ''
    if (!is.null(moon.alt.max)) {
        moon_alt.filter <- paste0("r.moon_alt <= " , moon.alt.max , " AND " )
    }

    session.id.filter <- ''
    if (!is.null(session.id)) {
        session.id.filter <- paste(session.id, collapse=",")
        session.id.filter <- paste0(
            "r.session_id IN (", session.id.filter , ") AND "
        )
    }

    rate.id.filter <- ''
    if (!is.null(rate.id)) {
        rate.id.filter <- paste(rate.id, collapse=",")
        rate.id.filter <- paste0("r.id IN (", rate.id.filter , ") AND " )
    }

    # query the data from SQL-Database
    with_query <- paste0("
WITH selection as (
    SELECT
        r.id as \"rate.id\",
        r.shower as \"shower.code\",
        r.period_start as \"period.start\",
        r.period_end as \"period.end\",
        r.sl_start as \"sl.start\",
        r.sl_end as \"sl.end\",
        r.session_id as \"session.id\",
        r.freq as \"freq\",
        r.lim_mag as \"lim.magn\",
        r.t_eff as \"t.eff\",
        r.f as \"f\",
        r.sidereal_time as \"time.sidereal\",
        r.sun_alt as \"sun.alt\",
        r.sun_az as \"sun.az\",
        r.moon_alt as \"moon.alt\",
        r.moon_az as \"moon.az\",
        r.moon_illum as \"moon.illum\",
        r.field_alt as \"field.alt\",
        r.field_az as \"field.az\",
        r.rad_alt as \"radiant.alt\",
        r.rad_az as \"radiant.az\",
        rm.magn_id as \"magn.id\"
    FROM rate as r
    LEFT JOIN rate_magnitude rm ON r.id = rm.rate_id
    WHERE
        ", period.filter, "
        ", sl.filter, "
        ", lim.magn.filter, "
        ", sun_alt.filter, "
        ", moon_alt.filter, "
        ", session.id.filter, "
        ", rate.id.filter, "
        ", shower.filter, " TRUE
)
    ")
    query <- paste0(with_query, "SELECT * FROM selection")
    observations <- DBI::dbGetQuery(dbcon, query)
    row.names(observations) <- observations$rate.id
    observations$shower.code <- factor(observations$shower.code)
    observations$session.id <- factor(observations$session.id)
    observations$magn.id <- factor(observations$magn.id)

    magnitudes <- NULL
    if (withMagnitudes) {
        query <- paste0(with_query, "
                SELECT
                    id as \"magn.id\",
                    magn as \"magn\",
                    freq as \"freq\"
                FROM magnitude_detail
                WHERE id IN (
                    SELECT \"magn.id\" FROM selection
                    WHERE \"magn.id\" IS NOT NULL
                )
            ")
        magnitudes <- DBI::dbGetQuery(dbcon, query)
        magnitudes$magn <- factor(
            magnitudes$magn,
            levels = sort(unique(magnitudes$magn), decreasing = TRUE),
            ordered = TRUE
        )
        magnitudes <- stats::xtabs(freq ~ magn.id + magn, data = magnitudes)
    }

    sessions <- NULL
    if (withSessions) {
        query <- paste0(with_query, "
                SELECT
                    id AS \"session.id\",
                    longitude,
                    latitude,
                    elevation,
                    country,
                    city as \"location.name\",
                    observer_id AS \"observer.id\",
                    observer_name AS \"observer.name\"
                FROM obs_session
                WHERE \"id\" IN (
                    SELECT \"session.id\" FROM selection
                    WHERE \"session.id\" IS NOT NULL
                )
            ")

        sessions <- DBI::dbGetQuery(dbcon, query)
        row.names(sessions) <- sessions$session.id
        sessions$country <- factor(sessions$country)
        sessions$location.name <- factor(sessions$location.name)
        sessions$observer.id <- factor(sessions$observer.id)
        sessions$observer.name <- factor(sessions$observer.name)
    }

    list(observations=observations, sessions=sessions, magnitudes=magnitudes)
}

#' @rdname load_vmdb
#' @export
load_vmdb_magnitudes <- function(
        dbcon,
        shower = NULL,
        period = NULL,
        sl = NULL,
        lim.magn = NULL,
        session.id = NULL,
        magn.id = NULL,
        withSessions = FALSE,
        withMagnitudes = TRUE
) {
    shower.filter <- ''
    if (!is.null(shower)) {
        shower <- unique(shower)
        shower.is_na <- is.na(shower)
        with_spo <- length(shower[shower.is_na]) > 0
        shower <- shower[!shower.is_na]
        with_showers <- length(shower) > 0

        if (with_spo) {
            spo.filter <- "shower IS NULL"
        } else {
            spo.filter <- "FALSE"
        }

        if (with_showers) {
            shower.filter <- paste(shower, collapse="','")
            shower.filter <- paste0("shower IN ('", shower.filter , "')" )
        } else {
            shower.filter <- "FALSE"
        }

        if (with_showers | with_spo) {
            shower.filter <- paste(
                '(', shower.filter, 'OR', spo.filter , ') AND '
            )
        } else {
            shower.filter <- ''
        }
    }

    period.filter <- ''
    if (!is.null(period)) {
        period <- matrix(period, ncol=2)
        period.filter <- apply(period, 1, function(period){
            paste0(
                "period_start >= '", period[1] ,
                "' AND period_end <= '", period[2] , "'"
            )
        })
        period.filter <- paste('(', paste(period.filter, collapse=') OR ('), ') AND ')
    }

    sl.filter <- ''
    if (!is.null(sl)) {
        sl <- matrix(sl, ncol=2)
        sl.filter <- apply(sl, 1, function(sl){
            if (sl[1] > sl[2]) {
                paste0(
                    "sl_start BETWEEN ", sl[1], " AND 360.0 AND ",
                    "sl_end BETWEEN 0.0 AND ", sl[2]
                )
            } else {
                paste0("sl_start > ", sl[1], " AND sl_end <= ", sl[2])
            }
        })
        sl.filter <- paste('(', paste(sl.filter, collapse=') OR ('), ') AND ')
    }

    lim.magn.filter <- ''
    if (!is.null(lim.magn)) {
        lim.magn <- matrix(lim.magn, ncol=2)
        lim.magn.filter <- apply(lim.magn, 1, function(lim.magn){
            paste0("lim_mag BETWEEN ", lim.magn[1], " AND ", lim.magn[2])
        })
        lim.magn.filter <- paste('(', paste(lim.magn.filter, collapse=') OR ('), ') AND ')
    }

    session.id.filter <- ''
    if (!is.null(session.id)) {
        session.id.filter <- paste(session.id, collapse=",")
        session.id.filter <- paste0(
            "\"session_id\" IN (", session.id.filter , ") AND "
        )
    }

    magn.id.filter <- ''
    if (!is.null(magn.id)) {
        magn.id.filter <- paste(magn.id, collapse=",")
        magn.id.filter <- paste0("id IN (", magn.id.filter , ") AND " )
    }

    # query the data from PostgreSQL
    with_query <- paste0("
WITH selection as (
    SELECT
        id as \"magn.id\",
        shower as \"shower.code\",
        period_start as \"period.start\",
        period_end as \"period.end\",
        sl_start as \"sl.start\",
        sl_end as \"sl.end\",
        session_id as \"session.id\",
        freq as \"freq\",
        mean as \"magn.mean\",
        lim_mag as \"lim.magn\"
    FROM magnitude
    WHERE
        ", period.filter, "
        ", sl.filter, "
        ", lim.magn.filter, "
        ", session.id.filter, "
        ", magn.id.filter, "
        ", shower.filter, " TRUE
)
    ")
    query <- paste0(with_query, "SELECT * FROM selection")
    observations <- DBI::dbGetQuery(dbcon, query)
    row.names(observations) <- observations$magn.id
    observations$shower.code <- factor(observations$shower.code)
    observations$session.id <- factor(observations$session.id)

    magnitudes <- NULL
    if (withMagnitudes) {
        query <- paste0(with_query, "
                SELECT
                    id as \"magn.id\",
                    magn as \"magn\",
                    freq as \"freq\"
                FROM magnitude_detail
                WHERE id IN (
                    SELECT \"magn.id\" FROM selection
                )
            ")
        magnitudes <- DBI::dbGetQuery(dbcon, query)
        magnitudes$magn <- factor(
            magnitudes$magn,
            levels = sort(unique(magnitudes$magn), decreasing = TRUE),
            ordered = TRUE
        )
        magnitudes <- stats::xtabs(freq ~ magn.id + magn, data = magnitudes)
    }

    sessions <- NULL
    if (withSessions) {
        query <- paste0(with_query, "
                SELECT
                    id AS \"session.id\",
                    longitude,
                    latitude,
                    elevation,
                    country,
                    city as \"location.name\",
                    observer_id AS \"observer.id\",
                    observer_name AS \"observer.name\"
                FROM obs_session
                WHERE \"id\" IN (
                    SELECT \"session.id\" FROM selection
                    WHERE \"session.id\" IS NOT NULL
                )
            ")
        sessions <- DBI::dbGetQuery(dbcon, query)
        row.names(sessions) <- sessions$session.id
        sessions$country <- factor(sessions$country)
        sessions$location.name <- factor(sessions$location.name)
        sessions$observer.id <- factor(sessions$observer.id)
        sessions$observer.name <- factor(sessions$observer.name)
    }

    list(observations=observations, sessions=sessions, magnitudes=magnitudes)
}
