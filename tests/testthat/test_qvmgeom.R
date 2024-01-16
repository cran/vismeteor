test_that("qvmgeom", {
    r <- 1.8

    p <- vismeteor::pvmgeom(6, 5.5, r, lower.tail = FALSE)
    q <- suppressWarnings(
        vismeteor::qvmgeom(c(-1.0, 0.0, p - 1e-06, p, p + 1e-06, 1.0), 5.5, r, lower.tail = FALSE)
    )
    expect_equal(q, c(NA, 5, NA, 5, 5, -Inf))

    p <- vismeteor::pvmgeom(6, 5.5, r, lower.tail = TRUE)
    q <- suppressWarnings(
        vismeteor::qvmgeom(c(-1.0, 1.0, p + 1e-06, p, p - 1e-06, 0.0), 5.5, r, lower.tail = TRUE)
    )
    expect_equal(q, c(NA, 5, NA, 5, 5, -Inf))

    lms <- seq(5.6, 6.5, 0.1)

    for (lm in lms) {
        q <- suppressWarnings(
            vismeteor::qvmgeom(c(-1.0, 0.0, 1.0), lm, r, lower.tail = FALSE)
        )
        expect_equal(q, c(NA, 6, -Inf), label = paste('lm =', lm))

        q <- suppressWarnings(
            vismeteor::qvmgeom(c(-1.0, 0.0, 1.0), lm, r, lower.tail = TRUE)
        )
        expect_equal(q, c(NA, -Inf, 6), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- vismeteor::pvmgeom(6, lm, r, lower.tail = FALSE)
        q <- vismeteor::qvmgeom(c(p - 1e-06, p, p + 1e-06), lm, r, lower.tail = FALSE)
        expect_equal(q, c(6, 6, 5), label = paste('lm =', lm))

        p <- vismeteor::pvmgeom(6, lm, r, lower.tail = TRUE)
        q <- vismeteor::qvmgeom(c(p + 1e-06, p -1e-06), lm, r, lower.tail = TRUE)
        expect_equal(q, c(6, 5), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- vismeteor::pvmgeom(5, lm, r, lower.tail = FALSE)
        q <- vismeteor::qvmgeom(c(p - 1e-06, p, p + 1e-06), lm, r, lower.tail = FALSE)
        expect_equal(q, c(5, 5, 4), label = paste('lm =', lm))

        p <- vismeteor::pvmgeom(5, lm, r, lower.tail = TRUE)
        q <- vismeteor::qvmgeom(c(p + 1e-06, p - 1e-06), lm, r, lower.tail = TRUE)
        expect_equal(q, c(5, 4), label = paste('lm =', lm))
    }

    # quantile of meteor magnitudes equals geometric distribution

    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    m <- as.integer(seq(0, 17, 1))
    p <- round(vismeteor::pvmgeom(6 - m, 6.5, r, lower.tail = FALSE, perception.fun = perception.const), 6)
    m <- vismeteor::qvmgeom(p, 6.5, r, lower.tail = FALSE, perception.fun = perception.const)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(stats::qgeom(p, 1 - 1/r, lower.tail = TRUE), 6 - m)

    m <- as.integer(seq(10, 30, 1))
    p <- round(vismeteor::pvmgeom(6 - m, 6.5, r, lower.tail = TRUE, perception.fun = perception.const), 6)
    m <- vismeteor::qvmgeom(p, 6.5, r, lower.tail = TRUE, perception.fun = perception.const)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(stats::qgeom(p, 1 - 1/r, lower.tail = FALSE), 6 - m)
})
