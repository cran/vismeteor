test_that("per_2015", {
    if (FALSE) {
        base_url <- "http://localhost:8000/api/v1"
        data_dir <- system.file("data", package = "vismeteor")

        PER_2015_rates <- load_vmdb_rates(
            base_url,
            shower = "PER",
            period = c("2015-01-01", "2015-12-31"),
            with_magnitudes = TRUE,
            with_sessions = TRUE
        )
        save(
            PER_2015_rates,
            file = normalizePath(paste0(data_dir, "/PER_2015_rates.RData")),
            compress = "xz"
        )

        PER_2015_magn <- load_vmdb_magnitudes(
            base_url,
            shower = "PER",
            period = c("2015-01-01", "2015-12-31"),
            with_magnitudes = TRUE,
            with_sessions = TRUE
        )
        save(
            PER_2015_magn,
            file = normalizePath(paste0(data_dir, "/PER_2015_magn.RData")),
            compress = "xz"
        )
    }

    expect_type(PER_2015_rates, "list")
    sessions <- PER_2015_rates$sessions
    expect_type(sessions, "list")
    expect_true(is(sessions, "data.frame"))
    observations <- PER_2015_rates$observations
    expect_type(observations, "list")
    expect_true(is(observations, "data.frame"))
    magnitudes <- PER_2015_rates$magnitudes
    expect_type(magnitudes, "double")
    expect_true(is(magnitudes, "table"))

    expect_type(PER_2015_magn, "list")
    sessions <- PER_2015_magn$sessions
    expect_type(sessions, "list")
    expect_true(is(sessions, "data.frame"))
    observations <- PER_2015_magn$observations
    expect_type(observations, "list")
    expect_true(is(observations, "data.frame"))
    magnitudes <- PER_2015_magn$magnitudes
    expect_type(magnitudes, "double")
    expect_true(is(magnitudes, "table"))
})
