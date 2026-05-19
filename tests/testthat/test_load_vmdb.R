test_that(".build_params: shower encoding", {
    p <- vismeteor:::.build_params("PER", NULL, NULL, NULL)
    expect_equal(p$multi$shower, "PER")

    # NA → SPO (sporadic)
    p <- vismeteor:::.build_params(NA, NULL, NULL, NULL)
    expect_equal(p$multi$shower, "SPO")

    # mixed: NA becomes SPO, named entries kept
    p <- vismeteor:::.build_params(c("PER", NA), NULL, NULL, NULL)
    expect_equal(p$multi$shower, c("PER", "SPO"))

    # NULL → no shower param
    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL)
    expect_null(p$multi$shower)
})

test_that(".build_params: range params", {
    # Date-only character input is expanded to strict ISO-T form:
    # lower bound -> midnight, upper bound -> 23:59:59 UTC.
    p <- vismeteor:::.build_params(NULL, c("2015-08-01", "2015-08-31"), NULL, NULL)
    expect_equal(p$scalar$period_start, "2015-08-01T00:00:00")
    expect_equal(p$scalar$period_end, "2015-08-31T23:59:59")

    p <- vismeteor:::.build_params(NULL, NULL, c(135.5, 145.5), NULL)
    expect_equal(p$scalar$sl_min, 135.5)
    expect_equal(p$scalar$sl_max, 145.5)

    p <- vismeteor:::.build_params(NULL, NULL, NULL, c(5.5, 6.5))
    expect_equal(p$scalar$lim_magn_min, 5.5)
    expect_equal(p$scalar$lim_magn_max, 6.5)
})

test_that(".parse_dt: strict ISO-T parser", {
    v <- vismeteor:::.parse_dt("2015-08-12T22:00:00")
    expect_s3_class(v, "POSIXct")
    expect_equal(attr(v, "tzone"), "UTC")
    expect_equal(
        as.numeric(v),
        as.numeric(as.POSIXct("2015-08-12 22:00:00", tz = "UTC"))
    )

    # Non-canonical forms are NOT tolerated.
    expect_true(is.na(vismeteor:::.parse_dt("2015-08-12 22:00:00")))
    expect_true(is.na(vismeteor:::.parse_dt("2015-08-12")))
    expect_true(is.na(vismeteor:::.parse_dt("not a date")))
})

test_that(".fmt_period: strict ISO-T formatter", {
    # POSIXct in -> T-form out, UTC.
    pt <- as.POSIXct("2015-08-12 22:00:00", tz = "UTC")
    expect_equal(vismeteor:::.fmt_period(pt, "lower"), "2015-08-12T22:00:00")
    expect_equal(vismeteor:::.fmt_period(pt, "upper"), "2015-08-12T22:00:00")

    # Date -> midnight (lower) / 23:59:59 (upper).
    d <- as.Date("2015-08-12")
    expect_equal(vismeteor:::.fmt_period(d, "lower"), "2015-08-12T00:00:00")
    expect_equal(vismeteor:::.fmt_period(d, "upper"), "2015-08-12T23:59:59")

    # character T-form -> identity.
    expect_equal(
        vismeteor:::.fmt_period("2015-08-12T22:00:00", "lower"),
        "2015-08-12T22:00:00"
    )

    # character space-form -> reformatted with T.
    expect_equal(
        vismeteor:::.fmt_period("2015-08-12 22:00:00", "lower"),
        "2015-08-12T22:00:00"
    )

    # character date-only -> bound expansion.
    expect_equal(
        vismeteor:::.fmt_period("2015-08-12", "lower"),
        "2015-08-12T00:00:00"
    )
    expect_equal(
        vismeteor:::.fmt_period("2015-08-12", "upper"),
        "2015-08-12T23:59:59"
    )

    # Unparseable -> error.
    expect_error(vismeteor:::.fmt_period("not a date", "lower"))
    expect_error(vismeteor:::.fmt_period(list(), "lower"))
})

test_that(".build_params: scalar altitude / id params", {
    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, sun_alt_max = -10)
    expect_equal(p$scalar$sun_alt_max, -10)

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, moon_alt_max = 5)
    expect_equal(p$scalar$moon_alt_max, 5)

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, session_id = c(1L, 2L))
    expect_equal(p$multi$session_id, c(1L, 2L))

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL,
        id_param = "rate_id", id_values = 100L
    )
    expect_equal(p$multi$rate_id, 100L)
})

test_that(".build_params: include parameter", {
    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL)
    expect_null(p$scalar$include)

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, with_sessions = TRUE)
    expect_equal(p$scalar$include, "sessions")

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, with_magnitudes = TRUE)
    expect_equal(p$scalar$include, "magnitude_details")

    p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL,
        with_sessions = TRUE, with_magnitudes = TRUE
    )
    expect_equal(p$scalar$include, "sessions,magnitude_details")
})

test_that(".parse_sessions: correct data.frame with factors and row names", {
    df <- data.frame(
        id = 1L,
        longitude = 10.0,
        latitude = 50.0,
        elevation = 0.3,
        country = "DE",
        location_name = "Somewhere",
        observer_id = "XX",
        observer_name = "Doe, J.",
        stringsAsFactors = FALSE
    )
    s <- vismeteor:::.parse_sessions(df)
    expect_true(is.data.frame(s))
    expect_equal(names(s), c(
        "session_id", "longitude", "latitude", "elevation",
        "country", "location_name", "observer_id", "observer_name"
    ))
    expect_equal(s$session_id, 1L)
    expect_true(is.factor(s$country))
    expect_true(is.factor(s$location_name))
    expect_true(is.factor(s$observer_id))
    expect_true(is.factor(s$observer_name))
    expect_equal(row.names(s), "1")

    expect_null(vismeteor:::.parse_sessions(NULL))
    expect_null(vismeteor:::.parse_sessions(list()))
})

test_that(".parse_magnitudes: builds xtabs table with correct dims", {
    df <- data.frame(
        id = c(200L, 200L, 200L, 201L, 201L),
        magn = c(3L, 2L, 1L, 3L, 2L),
        freq = c(2.5, 4.0, 5.5, 3.0, 5.0),
        stringsAsFactors = FALSE
    )
    m <- vismeteor:::.parse_magnitudes(df)
    expect_true(methods::is(m, "table"))
    expect_setequal(row.names(m), c("200", "201"))
    expect_setequal(colnames(m), c("1", "2", "3"))
    expect_equal(m["200", "2"], 4.0)
    expect_equal(m["201", "3"], 3.0)

    expect_null(vismeteor:::.parse_magnitudes(NULL))
    expect_null(vismeteor:::.parse_magnitudes(as.data.frame(list())))
})

# Integration tests: full pipeline with mocked HTTP responses.
# Fixture files live under tests/testthat/fixtures/{scenario}/
# and map to http://example.com/api/v1/{endpoint}.json

test_that("load_vmdb_rates: empty observations", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/empty", {
        res <- load_vmdb_rates("http://example.com/api/v1")
        expect_type(res, "list")
        expect_true(is.data.frame(res$observations))
        expect_equal(nrow(res$observations), 0)
        expect_null(res$sessions)
        expect_null(res$magnitudes)
    })
})

test_that("load_vmdb_rates: parses observations, sessions, magnitudes", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/full", {
        res <- load_vmdb_rates("http://example.com/api/v1")
        expect_type(res, "list")

        obs <- res$observations
        expect_true(is.data.frame(obs))
        expect_equal(nrow(obs), 1)
        expect_equal(obs$rate_id, 100L)
        expect_true(is.factor(obs$shower))
        expect_equal(as.character(obs$shower), "PER")
        expect_true(is.factor(obs$session_id))
        expect_true(is.factor(obs$magn_id))
        expect_true(is.logical(obs$magn_solo))
        expect_true(obs$magn_solo)
        expect_s3_class(obs$period_start, "POSIXct")
        expect_s3_class(obs$period_end, "POSIXct")
        expect_equal(attr(obs$period_start, "tzone"), "UTC")
        expect_equal(attr(obs$period_end, "tzone"), "UTC")
        expect_equal(row.names(obs), "100")

        sess <- res$sessions
        expect_true(is.data.frame(sess))
        expect_equal(row.names(sess), "1")
        expect_true(is.factor(sess$country))

        magn <- res$magnitudes
        expect_true(methods::is(magn, "table"))
        expect_equal(row.names(magn), "200")
        expect_setequal(colnames(magn), c("1", "2", "3"))
    })
})

test_that("load_vmdb_rates: all-sporadic response maps to SPO factor", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/sporadic_only", {
        res <- load_vmdb_rates("http://example.com/api/v1")
        obs <- res$observations
        expect_equal(nrow(obs), 2)
        expect_true(is.factor(obs$shower))
        expect_true("SPO" %in% levels(obs$shower))
        expect_equal(as.character(obs$shower), c("SPO", "SPO"))
        expect_false(any(is.na(obs$shower)))
    })
})

test_that("load_vmdb_rates: mixed PER/sporadic response gets SPO mapping", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/sporadic_mixed", {
        res <- load_vmdb_rates("http://example.com/api/v1")
        obs <- res$observations
        expect_equal(nrow(obs), 2)
        expect_true(is.factor(obs$shower))
        expect_setequal(as.character(obs$shower), c("PER", "SPO"))
        expect_false(any(is.na(obs$shower)))
    })
})

test_that("load_vmdb_magnitudes: all-sporadic response maps to SPO factor", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/sporadic_only_magn", {
        res <- load_vmdb_magnitudes("http://example.com/api/v1",
                                    with_magnitudes = FALSE)
        obs <- res$observations
        expect_equal(nrow(obs), 2)
        expect_true(is.factor(obs$shower))
        expect_equal(as.character(obs$shower), c("SPO", "SPO"))
    })
})

test_that("load_vmdb_magnitudes: mixed PER/sporadic response gets SPO mapping", {
    testthat::skip_if_not_installed("httptest2")
    httptest2::with_mock_dir("fixtures/sporadic_mixed_magn", {
        res <- load_vmdb_magnitudes("http://example.com/api/v1",
                                    with_magnitudes = FALSE)
        obs <- res$observations
        expect_equal(nrow(obs), 2)
        expect_setequal(as.character(obs$shower), c("PER", "SPO"))
    })
})

test_that(".build_params: POSIXct period preserves class", {
    pt <- as.POSIXct(c("2015-08-12 00:00:00", "2015-08-13 23:59:59"),
                     tz = "UTC")
    p <- vismeteor:::.build_params(NULL, pt, NULL, NULL)
    expect_equal(p$scalar$period_start, "2015-08-12T00:00:00")
    expect_equal(p$scalar$period_end, "2015-08-13T23:59:59")
})

test_that("load_vmdb_magnitudes: parses observations, sessions, magnitudes", {
    testthat::skip_if_not_installed("httptest2")
    # with_magnitudes defaults to TRUE → suppress include param for clean fixture path
    httptest2::with_mock_dir("fixtures/full", {
        res <- load_vmdb_magnitudes("http://example.com/api/v1", with_magnitudes = FALSE)
        expect_type(res, "list")

        obs <- res$observations
        expect_true(is.data.frame(obs))
        expect_equal(nrow(obs), 1)
        expect_equal(obs$magn_id, 200L)
        expect_true(is.factor(obs$shower))
        expect_equal(as.character(obs$shower), "PER")
        expect_s3_class(obs$period_start, "POSIXct")
        expect_s3_class(obs$period_end, "POSIXct")
        expect_equal(attr(obs$period_start, "tzone"), "UTC")
        expect_equal(attr(obs$period_end, "tzone"), "UTC")
        expect_equal(row.names(obs), "200")

        expect_true(is.data.frame(res$sessions))
        expect_true(methods::is(res$magnitudes, "table"))
    })
})
