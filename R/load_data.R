#' @name load_vmdb
#' @aliases load_vmdb_rates
#' @aliases load_vmdb_magnitudes
#' @import httr2
#' @title Loading visual meteor observations via the imo-vmdb REST API
#' @description
#' Loads visual meteor observations from an
#' \href{https://pypi.org/project/imo-vmdb/}{imo-vmdb} web server via its
#' REST API.
#' @note Angle values are expected and returned in degrees.
#' @param base_url character; base URL of the imo-vmdb API, e.g.
#'   \code{"http://localhost:8000/api/v1"}.
#' @param shower character; selects by meteor shower codes.
#'   \code{NA} loads sporadic meteors.
#' @param period time; selects a time range by minimum/maximum.
#' @param sl numeric; selects a range of solar longitudes by minimum/maximum.
#' @param lim.magn numeric; selects a range of limiting magnitudes by
#'   minimum/maximum.
#' @param sun.alt.max numeric; selects the maximum altitude of the sun
#'   (rates only).
#' @param moon.alt.max numeric; selects the maximum altitude of the moon
#'   (rates only).
#' @param session.id integer; selects by session ids.
#' @param rate.id integer; selects rate observations by ids.
#' @param magn.id integer; selects magnitude observations by ids.
#' @param withSessions logical; if \code{TRUE}, also load the corresponding
#'   session data.
#' @param withMagnitudes logical; if \code{TRUE}, also load the corresponding
#'   magnitude observations.
#' @details
#' \code{sl}, \code{period} and \code{lim.magn} expect a vector with
#' successive minimum and maximum values.
#' \code{sun.alt.max} and \code{moon.alt.max} are expected to be scalar values.
#'
#' \strong{Note:} Unlike the previous DBI-based version, only a single range
#' per filter parameter is supported.  If you previously passed a matrix with
#' multiple rows to \code{period}, \code{sl}, or \code{lim.magn}, flatten
#' them to a single min/max pair or issue multiple calls and combine with
#' \code{rbind()}.
#' @return
#' Both functions return a list, with
#'
#' \tabular{ll}{
#'     \code{observations} \tab data frame, rate or magnitude observations,\cr
#'     \code{sessions} \tab data frame; session data of observations,\cr
#'     \code{magnitudes} \tab table; contingency table of meteor magnitude
#'       frequencies.
#' }
#'
#' \code{observations} depends on the function call.
#' \code{load_vmdb_rates} returns a data frame with columns:
#'
#' \tabular{ll}{
#' \code{rate.id} \tab unique identifier of the rate observation,\cr
#' \code{shower.code} \tab IAU code of the shower. \code{NA} for sporadic.\cr
#' \code{period.start} \tab start of observation,\cr
#' \code{period.end} \tab end of observation,\cr
#' \code{sl.start} \tab solar longitude at start,\cr
#' \code{sl.end} \tab solar longitude at end,\cr
#' \code{session.id} \tab reference to the session,\cr
#' \code{freq} \tab count of observed meteors,\cr
#' \code{lim.magn} \tab limiting magnitude,\cr
#' \code{t.eff} \tab net observed time in hours,\cr
#' \code{f} \tab correction factor of cloud cover,\cr
#' \code{time.sidereal} \tab sidereal time,\cr
#' \code{sun.alt} \tab altitude of the sun,\cr
#' \code{sun.az} \tab azimuth of the sun,\cr
#' \code{moon.alt} \tab altitude of the moon,\cr
#' \code{moon.az} \tab azimuth of the moon,\cr
#' \code{moon.illum} \tab illumination of the moon (\code{0.0 .. 1.0}),\cr
#' \code{field.alt} \tab altitude of the field of view (optional),\cr
#' \code{field.az} \tab azimuth of the field of view (optional),\cr
#' \code{radiant.alt} \tab altitude of the radiant (optional),\cr
#' \code{radiant.az} \tab azimuth of the radiant (optional),\cr
#' \code{magn.id} \tab reference to the magnitude observations (optional).
#' }
#'
#' \code{load_vmdb_magnitudes} returns an \code{observations} data frame with:
#'
#' \tabular{ll}{
#' \code{magn.id} \tab unique identifier of the magnitude observation,\cr
#' \code{shower.code} \tab IAU code of the shower. \code{NA} for sporadic.\cr
#' \code{period.start} \tab start of observation,\cr
#' \code{period.end} \tab end of observation,\cr
#' \code{sl.start} \tab solar longitude at start,\cr
#' \code{sl.end} \tab solar longitude at end,\cr
#' \code{session.id} \tab reference to the session,\cr
#' \code{freq} \tab count of observed meteors,\cr
#' \code{magn.mean} \tab mean magnitude,\cr
#' \code{lim.magn} \tab limiting magnitude (optional).
#' }
#'
#' The \code{sessions} data frame contains
#'
#' \tabular{ll}{
#' \code{session.id} \tab unique identifier of the session,\cr
#' \code{longitude} \tab location's longitude,\cr
#' \code{latitude} \tab location's latitude,\cr
#' \code{elevation} \tab height above mean sea level in km,\cr
#' \code{country} \tab country name,\cr
#' \code{location.name} \tab location name,\cr
#' \code{observer.id} \tab observer id (optional),\cr
#' \code{observer.name} \tab observer name (optional).
#' }
#'
#' \code{magnitudes} is a contingency table of meteor magnitude frequencies.
#' Row names are magnitude observation IDs; column names are magnitude classes.
#'
#' @references \url{https://pypi.org/project/imo-vmdb/}
#' @examples
#' \dontrun{
#' # Load rate observations including session and magnitude data
#' data <- load_vmdb_rates(
#'     base_url       = "http://localhost:8000/api/v1",
#'     shower         = "PER",
#'     sl             = c(135.5, 145.5),
#'     period         = c("2015-08-01", "2015-08-31"),
#'     lim.magn       = c(5.3, 6.7),
#'     withMagnitudes = TRUE,
#'     withSessions   = TRUE
#' )
#'
#' # Load magnitude observations
#' data <- load_vmdb_magnitudes(
#'     base_url     = "http://localhost:8000/api/v1",
#'     shower       = "PER",
#'     sl           = c(135.5, 145.5),
#'     period       = c("2015-08-01", "2015-08-31"),
#'     lim.magn     = c(5.3, 6.7),
#'     withSessions = TRUE
#' )
#' }

#' @rdname load_vmdb
#' @export
load_vmdb_rates <- function(
    base_url,
    shower       = NULL,
    period       = NULL,
    sl           = NULL,
    lim.magn     = NULL,
    sun.alt.max  = NULL,
    moon.alt.max = NULL,
    session.id   = NULL,
    rate.id      = NULL,
    withSessions   = FALSE,
    withMagnitudes = FALSE
) {
    p <- .build_params(
        shower, period, sl, lim.magn,
        sun.alt.max  = sun.alt.max,
        moon.alt.max = moon.alt.max,
        session.id   = session.id,
        id.param     = "rate_id",
        id.values    = rate.id,
        withSessions   = withSessions,
        withMagnitudes = withMagnitudes
    )

    body <- .api_get(base_url, "rates", p$scalar, p$multi)

    if (length(body$observations) == 0) {
        observations <- data.frame()
    } else {
        observations <- .remap_cols(as.data.frame(body$observations), .rate_col_map)
        observations$shower.code <- factor(observations$shower.code)
        observations$session.id  <- factor(observations$session.id)
        observations$magn.id     <- factor(observations$magn.id)
        row.names(observations)  <- observations$rate.id
    }

    list(
        observations = observations,
        sessions     = .parse_sessions(body$sessions),
        magnitudes   = .parse_magnitudes(body$magnitudes)
    )
}

#' @rdname load_vmdb
#' @export
load_vmdb_magnitudes <- function(
    base_url,
    shower       = NULL,
    period       = NULL,
    sl           = NULL,
    lim.magn     = NULL,
    session.id   = NULL,
    magn.id      = NULL,
    withSessions   = FALSE,
    withMagnitudes = TRUE
) {
    p <- .build_params(
        shower, period, sl, lim.magn,
        session.id   = session.id,
        id.param     = "magn_id",
        id.values    = magn.id,
        withSessions   = withSessions,
        withMagnitudes = withMagnitudes
    )

    body <- .api_get(base_url, "magnitudes", p$scalar, p$multi)

    if (length(body$observations) == 0) {
        observations <- data.frame()
    } else {
        observations <- .remap_cols(as.data.frame(body$observations), .magn_col_map)
        observations$shower.code <- factor(observations$shower.code)
        observations$session.id  <- factor(observations$session.id)
        row.names(observations)  <- observations$magn.id
    }

    list(
        observations = observations,
        sessions     = .parse_sessions(body$sessions),
        magnitudes   = .parse_magnitudes(body$magnitudes)
    )
}

# Build a named list of scalar query parameters from filter arguments.
# Multi-value parameters (shower, session_id, rate_id, magn_id) are returned
# separately as a named list of vectors.
.build_params <- function(
    shower, period, sl, lim.magn,
    sun.alt.max = NULL, moon.alt.max = NULL,
    session.id = NULL, id.param = NULL, id.values = NULL,
    withSessions = FALSE, withMagnitudes = FALSE
) {
    params <- list()
    multi  <- list()

    if (!is.null(shower)) {
        multi$shower <- ifelse(is.na(shower), "SPO", shower)
    }

    if (!is.null(period)) {
        period <- matrix(period, ncol = 2)
        params$period_start <- min(period[, 1])
        params$period_end   <- max(period[, 2])
    }

    if (!is.null(sl)) {
        sl <- matrix(sl, ncol = 2)
        params$sl_min <- min(sl[, 1])
        params$sl_max <- max(sl[, 2])
    }

    if (!is.null(lim.magn)) {
        lim.magn <- matrix(lim.magn, ncol = 2)
        params$lim_magn_min <- min(lim.magn[, 1])
        params$lim_magn_max <- max(lim.magn[, 2])
    }

    if (!is.null(sun.alt.max))  params$sun_alt_max  <- sun.alt.max
    if (!is.null(moon.alt.max)) params$moon_alt_max <- moon.alt.max

    if (!is.null(session.id))   multi$session_id    <- as.integer(session.id)
    if (!is.null(id.values))    multi[[id.param]]   <- as.integer(id.values)

    include <- character(0)
    if (withSessions)   include <- c(include, "sessions")
    if (withMagnitudes) include <- c(include, "magnitudes")
    if (length(include) > 0) params$include <- paste(include, collapse = ",")

    list(scalar = params, multi = multi)
}


# Send a GET request to the API and return the parsed JSON body.
.api_get <- function(base_url, path, scalar_params, multi_params) {
    req <- httr2::request(base_url) |>
        httr2::req_url_path_append(path)

    if (length(scalar_params) > 0) {
        req <- do.call(httr2::req_url_query, c(list(req), scalar_params))
    }

    for (name in names(multi_params)) {
        vals <- multi_params[[name]]
        req <- do.call(
            httr2::req_url_query,
            c(list(req, .multi = "explode"), stats::setNames(list(vals), name))
        )
    }

    resp <- httr2::req_perform(req)
    httr2::resp_body_json(resp, simplifyVector = TRUE)
}


# Explicit column maps: API name (DB column) -> R output name.
.rate_col_map <- c(
    id            = "rate.id",
    shower        = "shower.code",
    period_start  = "period.start",
    period_end    = "period.end",
    sl_start      = "sl.start",
    sl_end        = "sl.end",
    session_id    = "session.id",
    freq          = "freq",
    lim_mag       = "lim.magn",
    t_eff         = "t.eff",
    f             = "f",
    sidereal_time = "time.sidereal",
    sun_alt       = "sun.alt",
    sun_az        = "sun.az",
    moon_alt      = "moon.alt",
    moon_az       = "moon.az",
    moon_illum    = "moon.illum",
    field_alt     = "field.alt",
    field_az      = "field.az",
    rad_alt       = "radiant.alt",
    rad_az        = "radiant.az",
    magn_id       = "magn.id"
)

.magn_col_map <- c(
    id           = "magn.id",
    shower       = "shower.code",
    period_start = "period.start",
    period_end   = "period.end",
    sl_start     = "sl.start",
    sl_end       = "sl.end",
    session_id   = "session.id",
    freq         = "freq",
    mean         = "magn.mean",
    lim_mag      = "lim.magn"
)

.session_col_map <- c(
    id            = "session.id",
    longitude     = "longitude",
    latitude      = "latitude",
    elevation     = "elevation",
    country       = "country",
    city          = "location.name",
    observer_id   = "observer.id",
    observer_name = "observer.name"
)

.remap_cols <- function(df, col_map) {
    names(df) <- col_map[names(df)]
    df
}


# Parse the sessions array from the API response into a data.frame.
.parse_sessions <- function(sessions_list) {
    if (is.null(sessions_list) || length(sessions_list) == 0) return(NULL)
    s <- .remap_cols(as.data.frame(sessions_list), .session_col_map)
    s$country       <- factor(s$country)
    s$location.name <- factor(s$location.name)
    s$observer.id   <- factor(s$observer.id)
    s$observer.name <- factor(s$observer.name)
    row.names(s) <- s$session.id
    s
}


# Parse the magnitudes array from the API response into an xtabs table.
.parse_magnitudes <- function(magnitudes_list) {
    if (is.null(magnitudes_list) || length(magnitudes_list) == 0) return(NULL)
    m <- as.data.frame(magnitudes_list)
    if (nrow(m) == 0) return(NULL)
    m$magn <- factor(
        m$magn,
        levels  = sort(unique(m$magn), decreasing = TRUE),
        ordered = TRUE
    )
    names(m)[names(m) == "id"] <- "magn.id"
    stats::xtabs(freq ~ magn.id + magn, data = m)
}


