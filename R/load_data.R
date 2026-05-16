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
#' @param lim_magn numeric; selects a range of limiting magnitudes by
#'   minimum/maximum.
#' @param sun_alt_max numeric; selects the maximum altitude of the sun
#'   (rates only).
#' @param moon_alt_max numeric; selects the maximum altitude of the moon
#'   (rates only).
#' @param session_id integer; selects by session ids.
#' @param rate_id integer; selects rate observations by ids.
#' @param magn_id integer; selects magnitude observations by ids.
#' @param with_sessions logical; if \code{TRUE}, also load the corresponding
#'   session data.
#' @param with_magnitudes logical; if \code{TRUE}, also load the corresponding
#'   magnitude observations.
#' @details
#' \code{sl}, \code{period} and \code{lim_magn} expect a vector with
#' successive minimum and maximum values.
#' \code{sun_alt_max} and \code{moon_alt_max} are expected to be scalar values.
#'
#' \strong{Note:} Unlike the previous DBI-based version, only a single range
#' per filter parameter is supported.  If you previously passed a matrix with
#' multiple rows to \code{period}, \code{sl}, or \code{lim_magn}, flatten
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
#' \code{rate_id} \tab unique identifier of the rate observation,\cr
#' \code{shower} \tab IAU code of the shower. \code{NA} for sporadic.\cr
#' \code{period_start} \tab start of observation,\cr
#' \code{period_end} \tab end of observation,\cr
#' \code{sl_start} \tab solar longitude at start,\cr
#' \code{sl_end} \tab solar longitude at end,\cr
#' \code{session_id} \tab reference to the session,\cr
#' \code{freq} \tab count of observed meteors,\cr
#' \code{lim_magn} \tab limiting magnitude,\cr
#' \code{t_eff} \tab net observed time in hours,\cr
#' \code{f} \tab correction factor of cloud cover,\cr
#' \code{sidereal_time} \tab sidereal time,\cr
#' \code{sun_alt} \tab altitude of the sun,\cr
#' \code{sun_az} \tab azimuth of the sun,\cr
#' \code{moon_alt} \tab altitude of the moon,\cr
#' \code{moon_az} \tab azimuth of the moon,\cr
#' \code{moon_illum} \tab illumination of the moon (\code{0.0 .. 1.0}),\cr
#' \code{field_alt} \tab altitude of the field of view (optional),\cr
#' \code{field_az} \tab azimuth of the field of view (optional),\cr
#' \code{rad_alt} \tab altitude of the radiant (optional),\cr
#' \code{rad_az} \tab azimuth of the radiant (optional),\cr
#' \code{magn_id} \tab reference to the magnitude observations (optional).
#' }
#'
#' \code{load_vmdb_magnitudes} returns an \code{observations} data frame with:
#'
#' \tabular{ll}{
#' \code{magn_id} \tab unique identifier of the magnitude observation,\cr
#' \code{shower} \tab IAU code of the shower. \code{NA} for sporadic.\cr
#' \code{period_start} \tab start of observation,\cr
#' \code{period_end} \tab end of observation,\cr
#' \code{sl_start} \tab solar longitude at start,\cr
#' \code{sl_end} \tab solar longitude at end,\cr
#' \code{session_id} \tab reference to the session,\cr
#' \code{freq} \tab count of observed meteors,\cr
#' \code{magn_mean} \tab mean magnitude,\cr
#' \code{lim_magn} \tab limiting magnitude (optional).
#' }
#'
#' The \code{sessions} data frame contains
#'
#' \tabular{ll}{
#' \code{session_id} \tab unique identifier of the session,\cr
#' \code{longitude} \tab location's longitude,\cr
#' \code{latitude} \tab location's latitude,\cr
#' \code{elevation} \tab height above mean sea level in km,\cr
#' \code{country} \tab country name,\cr
#' \code{location_name} \tab location name,\cr
#' \code{observer_id} \tab observer id (optional),\cr
#' \code{observer_name} \tab observer name (optional).
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
#'     lim_magn       = c(5.3, 6.7),
#'     with_magnitudes = TRUE,
#'     with_sessions   = TRUE
#' )
#'
#' # Load magnitude observations
#' data <- load_vmdb_magnitudes(
#'     base_url     = "http://localhost:8000/api/v1",
#'     shower       = "PER",
#'     sl           = c(135.5, 145.5),
#'     period       = c("2015-08-01", "2015-08-31"),
#'     lim_magn     = c(5.3, 6.7),
#'     with_sessions = TRUE
#' )
#' }

#' @rdname load_vmdb
#' @export
load_vmdb_rates <- function(
  base_url,
  shower = NULL,
  period = NULL,
  sl = NULL,
  lim_magn = NULL,
  sun_alt_max = NULL,
  moon_alt_max = NULL,
  session_id = NULL,
  rate_id = NULL,
  with_sessions = FALSE,
  with_magnitudes = FALSE
) {
    p <- .build_params(
        shower, period, sl, lim_magn,
        sun_alt_max = sun_alt_max,
        moon_alt_max = moon_alt_max,
        session_id = session_id,
        id_param = "rate_id",
        id_values = rate_id,
        with_sessions = with_sessions,
        with_magnitudes = with_magnitudes
    )

    body <- .api_get(base_url, "rates", p$scalar, p$multi)

    if (length(body$observations) == 0) {
        observations <- data.frame()
    } else {
        observations <- .remap_cols(as.data.frame(body$observations), .rate_col_map)
        observations$shower <- factor(observations$shower)
        observations$session_id <- factor(observations$session_id)
        observations$magn_id <- factor(observations$magn_id)
        row.names(observations) <- observations$rate_id
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
  shower = NULL,
  period = NULL,
  sl = NULL,
  lim_magn = NULL,
  session_id = NULL,
  magn_id = NULL,
  with_sessions = FALSE,
  with_magnitudes = TRUE
) {
    p <- .build_params(
        shower, period, sl, lim_magn,
        session_id = session_id,
        id_param = "magn_id",
        id_values = magn_id,
        with_sessions = with_sessions,
        with_magnitudes = with_magnitudes
    )

    body <- .api_get(base_url, "magnitudes", p$scalar, p$multi)

    if (length(body$observations) == 0) {
        observations <- data.frame()
    } else {
        observations <- .remap_cols(as.data.frame(body$observations), .magn_col_map)
        observations$shower <- factor(observations$shower)
        observations$session_id <- factor(observations$session_id)
        row.names(observations) <- observations$magn_id
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
  shower, period, sl, lim_magn,
  sun_alt_max = NULL, moon_alt_max = NULL,
  session_id = NULL, id_param = NULL, id_values = NULL,
  with_sessions = FALSE, with_magnitudes = FALSE
) {
    params <- list()
    multi <- list()

    if (!is.null(shower)) {
        multi$shower <- ifelse(is.na(shower), "SPO", shower)
    }

    if (!is.null(period)) {
        period <- matrix(period, ncol = 2)
        params$period_start <- min(period[, 1])
        params$period_end <- max(period[, 2])
    }

    if (!is.null(sl)) {
        sl <- matrix(sl, ncol = 2)
        params$sl_min <- min(sl[, 1])
        params$sl_max <- max(sl[, 2])
    }

    if (!is.null(lim_magn)) {
        lim_magn <- matrix(lim_magn, ncol = 2)
        params$lim_magn_min <- min(lim_magn[, 1])
        params$lim_magn_max <- max(lim_magn[, 2])
    }

    if (!is.null(sun_alt_max)) params$sun_alt_max <- sun_alt_max
    if (!is.null(moon_alt_max)) params$moon_alt_max <- moon_alt_max

    if (!is.null(session_id)) multi$session_id <- as.integer(session_id)
    if (!is.null(id_values)) multi[[id_param]] <- as.integer(id_values)

    include <- character(0)
    if (with_sessions) include <- c(include, "sessions")
    if (with_magnitudes) include <- c(include, "magnitudes")
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


# Semantic renames only: disambiguate the generic `id` column against
# foreign-key columns in the same data frame, and avoid shadowing base R's
# `mean()`. Columns not listed here are passed through from the imo-vmdb
# API unchanged (it already uses snake_case).
.rate_col_map <- c(
    id = "rate_id"
)

.magn_col_map <- c(
    id   = "magn_id",
    mean = "magn_mean"
)

.session_col_map <- c(
    id = "session_id"
)

.remap_cols <- function(df, col_map) {
    new_names <- col_map[names(df)]
    names(df) <- ifelse(is.na(new_names), names(df), new_names)
    df
}


# Parse the sessions array from the API response into a data.frame.
.parse_sessions <- function(sessions_list) {
    if (is.null(sessions_list) || length(sessions_list) == 0) {
        return(NULL)
    }
    s <- .remap_cols(as.data.frame(sessions_list), .session_col_map)
    s$country <- factor(s$country)
    s$location_name <- factor(s$location_name)
    s$observer_id <- factor(s$observer_id)
    s$observer_name <- factor(s$observer_name)
    row.names(s) <- s$session_id
    s
}


# Parse the magnitudes array from the API response into an xtabs table.
.parse_magnitudes <- function(magnitudes_list) {
    if (is.null(magnitudes_list) || length(magnitudes_list) == 0) {
        return(NULL)
    }
    m <- as.data.frame(magnitudes_list)
    if (nrow(m) == 0) {
        return(NULL)
    }
    m$magn <- factor(
        m$magn,
        levels  = sort(unique(m$magn), decreasing = TRUE),
        ordered = TRUE
    )
    names(m)[names(m) == "id"] <- "magn_id"
    stats::xtabs(freq ~ magn_id + magn, data = m)
}
