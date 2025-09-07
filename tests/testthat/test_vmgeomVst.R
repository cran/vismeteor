test_that("vmgeomVst", {
    r <- seq(1.5, 3.3, 0.1)
    limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
    m <- seq(-200, 6, 1)
    model <- expand.grid(r = r, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(r, limmag) {
            p <- dvmgeom(m, limmag, r)
            t <- vmgeomVstFromMagn(m, limmag)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                r = r,
                q = log(r),
                limmag = limmag,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$r, model$limmag , SIMPLIFY = FALSE)
    )

    #
    # test vmgeomVstFromMagn
    #
    expect_true(vmgeomVstFromMagn(6, 5.5) > 0.0)
    expect_true(all(abs(model$t.var - 1.0) < 0.018))

    #
    # test vmgeomVstToR
    #
    model$r.est <- vmgeomVstToR(model$t.mean)
    expect_true(all(abs(model$r - model$r.est) < 0.013))

    # test with non-exotic values
    model0 <- with(model, {
        subset(model, r < 2.7)
    })
    expect_true(all(abs(model0$r - model0$r.est) < 0.007))

    # test first derivative
    f <- function(x) {
        vmgeomVstToR(x, deriv.degree = 1L)
    }
    y <- vmgeomVstToR(5.5) - vmgeomVstToR(4.5)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)

    # test second derivative
    f <- function(x) {
        vmgeomVstToR(x, deriv.degree = 2L)
    }
    y <- vmgeomVstToR(5.5, deriv.degree = 1L) - vmgeomVstToR(4.5, deriv.degree = 1L)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)

    # log ...
    model$q.est <- vmgeomVstToR(model$t.mean, log = TRUE)
    expect_true(all(abs(model$q - model$q.est) < 0.004))

    # test with non-exotic values
    model0 <- with(model, {
        subset(model, r < 2.7)
    })
    expect_true(all(abs(model0$q - model0$q.est) < 0.004))

    # test log first derivative
    f <- function(x) {
        vmgeomVstToR(x, log = TRUE, deriv.degree = 1L)
    }
    y <- vmgeomVstToR(5.5, log = TRUE) - vmgeomVstToR(4.5, log = TRUE)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)

    # test log second derivative
    f <- function(x) {
        vmgeomVstToR(x, log = TRUE, deriv.degree = 2L)
    }
    y <- vmgeomVstToR(5.5, log = TRUE, deriv.degree = 1L) -
        vmgeomVstToR(4.5, log = TRUE, deriv.degree = 1L)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)
})
