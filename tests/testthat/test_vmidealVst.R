test_that("vmidealVst", {
    psi <- seq(-3.5, 6.5, 0.5)
    limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
    m <- seq(-200, 6, 1)

    model <- expand.grid(psi = psi, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(psi, limmag) {
            p <- dvmideal(m, limmag, psi)
            t <- vmidealVstFromMagn(m, limmag)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                psi = psi,
                limmag = limmag,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$psi, model$limmag , SIMPLIFY = FALSE)
    )

    #
    # test vmidealVstFromMagn
    #
    expect_true(all(abs(model$t.var - 1.0) < 0.0017))

    #
    # test vmidealVstToPsi
    #
    model$psi.est <- vmidealVstToPsi(model$t.mean, model$limmag)
    expect_true(mean(abs(model$psi - model$psi.est)) < 0.0072)
    expect_true(all(abs(model$psi - model$psi.est) < 0.03))

    # test with low psi values
    model0 <- with(model, {
        subset(model, psi < 5.0)
    })
    expect_true(mean(abs(model0$psi - model0$psi.est)) < 0.006)
    expect_true(all(abs(model0$psi - model0$psi.est) < 0.023))

    # test first derivative
    f <- function(x) {
        vmidealVstToPsi(x, 6.0, deriv.degree = 1L)
    }
    y <- vmidealVstToPsi(5.5, 6.0) - vmidealVstToPsi(4.5, 6.0)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)

    # test second derivative
    f <- function(x) {
        vmidealVstToPsi(x, 6.0, deriv.degree = 2L)
    }
    y <- vmidealVstToPsi(5.5, 6.0, deriv.degree = 1L) - vmidealVstToPsi(4.5, 6.0, deriv.degree = 1L)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)
})
