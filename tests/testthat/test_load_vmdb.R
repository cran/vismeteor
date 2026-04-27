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
  p <- vismeteor:::.build_params(NULL, c("2015-08-01", "2015-08-31"), NULL, NULL)
  expect_equal(p$scalar$period_start, "2015-08-01")
  expect_equal(p$scalar$period_end,   "2015-08-31")

  p <- vismeteor:::.build_params(NULL, NULL, c(135.5, 145.5), NULL)
  expect_equal(p$scalar$sl_min, 135.5)
  expect_equal(p$scalar$sl_max, 145.5)

  p <- vismeteor:::.build_params(NULL, NULL, NULL, c(5.5, 6.5))
  expect_equal(p$scalar$lim_magn_min, 5.5)
  expect_equal(p$scalar$lim_magn_max, 6.5)
})

test_that(".build_params: scalar altitude / id params", {
  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, sun.alt.max = -10)
  expect_equal(p$scalar$sun_alt_max, -10)

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, moon.alt.max = 5)
  expect_equal(p$scalar$moon_alt_max, 5)

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, session.id = c(1L, 2L))
  expect_equal(p$multi$session_id, c(1L, 2L))

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL,
                                  id.param = "rate_id", id.values = 100L)
  expect_equal(p$multi$rate_id, 100L)
})

test_that(".build_params: include parameter", {
  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL)
  expect_null(p$scalar$include)

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, withSessions = TRUE)
  expect_equal(p$scalar$include, "sessions")

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL, withMagnitudes = TRUE)
  expect_equal(p$scalar$include, "magnitudes")

  p <- vismeteor:::.build_params(NULL, NULL, NULL, NULL,
                                  withSessions = TRUE, withMagnitudes = TRUE)
  expect_equal(p$scalar$include, "sessions,magnitudes")
})

test_that(".parse_sessions: correct data.frame with factors and row names", {
  df <- data.frame(
    id            = 1L,
    longitude     = 10.0,
    latitude      = 50.0,
    elevation     = 0.3,
    country       = "DE",
    city          = "Somewhere",
    observer_id   = "XX",
    observer_name = "Doe, J.",
    stringsAsFactors = FALSE
  )
  s <- vismeteor:::.parse_sessions(df)
  expect_true(is.data.frame(s))
  expect_equal(names(s), c("session.id", "longitude", "latitude", "elevation",
                            "country", "location.name", "observer.id", "observer.name"))
  expect_equal(s$session.id, 1L)
  expect_true(is.factor(s$country))
  expect_true(is.factor(s$location.name))
  expect_true(is.factor(s$observer.id))
  expect_true(is.factor(s$observer.name))
  expect_equal(row.names(s), "1")

  expect_null(vismeteor:::.parse_sessions(NULL))
  expect_null(vismeteor:::.parse_sessions(list()))
})

test_that(".parse_magnitudes: builds xtabs table with correct dims", {
  df <- data.frame(
    id   = c(200L, 200L, 200L, 201L, 201L),
    magn = c(3L,   2L,   1L,   3L,   2L),
    freq = c(2.5,  4.0,  5.5,  3.0,  5.0),
    stringsAsFactors = FALSE
  )
  m <- vismeteor:::.parse_magnitudes(df)
  expect_true(methods::is(m, "table"))
  expect_setequal(row.names(m), c("200", "201"))
  expect_setequal(colnames(m),  c("1", "2", "3"))
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
    expect_equal(obs$rate.id,      100L)
    expect_true(is.factor(obs$shower.code))
    expect_equal(as.character(obs$shower.code), "PER")
    expect_true(is.factor(obs$session.id))
    expect_true(is.factor(obs$magn.id))
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

test_that("load_vmdb_magnitudes: parses observations, sessions, magnitudes", {
  testthat::skip_if_not_installed("httptest2")
  # withMagnitudes defaults to TRUE → suppress include param for clean fixture path
  httptest2::with_mock_dir("fixtures/full", {
    res <- load_vmdb_magnitudes("http://example.com/api/v1", withMagnitudes = FALSE)
    expect_type(res, "list")

    obs <- res$observations
    expect_true(is.data.frame(obs))
    expect_equal(nrow(obs), 1)
    expect_equal(obs$magn.id, 200L)
    expect_true(is.factor(obs$shower.code))
    expect_equal(as.character(obs$shower.code), "PER")
    expect_equal(row.names(obs), "200")

    expect_true(is.data.frame(res$sessions))
    expect_true(methods::is(res$magnitudes, "table"))
  })
})
