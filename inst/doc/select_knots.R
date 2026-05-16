## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(vismeteor)

## ----echo=TRUE----------------------------------------------------------------
data(PER_2015_rates)
obs <- PER_2015_rates$observations
obs <- subset(
    obs,
    rad_alt > 15 & moon_alt < 0 & sun_alt < -5 &
        sl_start >= 139.2 & sl_end <= 140.6
)

deg2rad <- \(x) x * pi / 180
obs$sl <- (obs$sl_start + obs$sl_end) / 2
obs$t <- with(obs, t_eff * sin(deg2rad(rad_alt)) / f)

rates <- data.frame(
    sl     = obs$sl,
    t      = obs$t,
    limmag = obs$lim_magn,
    freq   = obs$freq
)
rates <- rates[order(rates$sl), ]
nrow(rates)

## ----echo=TRUE----------------------------------------------------------------
fit_glm <- \(rates, knots) {
    f <- if (length(knots) == 0L) {
        freq ~ splines::ns(sl) + limmag + offset(log(t))
    } else {
        freq ~ splines::ns(sl, knots = knots) +
            limmag + offset(log(t))
    }
    stats::glm(f, data = rates, family = stats::poisson())
}
score_bic <- \(rates, knots) stats::BIC(fit_glm(rates, knots))

## ----echo=TRUE----------------------------------------------------------------
rates$q <- vismeteor::freq_quantile(rates$freq, 12) # >= 12 meteors per bin
rates$qsl <- with(
    rates,
    ave(freq * sl, q, FUN = sum) /
        ave(freq, q, FUN = sum)
)
cand <- unique(round(rates$qsl, 2))
length(cand)

## ----echo=TRUE----------------------------------------------------------------
res_fwd <- vismeteor::select_knots(rates, cand, score_bic)
sort(c(res_fwd$knots, res_fwd$fixed_knots))
res_fwd$score

## ----echo=TRUE, fig.alt="BIC trajectory of the forward search over the number of interior knots, with the score optimum at the endpoint"----
plot(res_fwd$history$n_knots, res_fwd$history$score,
    type = "b", pch = 19, col = "blue",
    xlab = "number of interior knots", ylab = "BIC",
    main = "Forward selection: BIC trajectory"
)
points(tail(res_fwd$history$n_knots, 1L), res_fwd$score,
    pch = 19, col = "red", cex = 1.5
)
abline(h = res_fwd$score, lty = 2, col = "grey50")

## ----echo=TRUE----------------------------------------------------------------
res_fwd_plus1 <- vismeteor::select_knots(rates, cand, score_bic,
    fixed_knots = c(
        res_fwd$knots,
        res_fwd$fixed_knots
    ),
    n_steps = 1L
)
res_fwd_plus1$score - res_fwd$score # always >= 0: how much worse one step makes it

## ----echo=TRUE----------------------------------------------------------------
res_bwd <- vismeteor::select_knots(rates, cand, score_bic, backward = TRUE)
sort(c(res_bwd$knots, res_bwd$fixed_knots))
res_bwd$score

## ----echo=TRUE, fig.alt="BIC trajectory of the backward search over the number of interior knots, with the score optimum at the endpoint"----
plot(res_bwd$history$n_knots, res_bwd$history$score,
    type = "b", pch = 19, col = "blue",
    xlab = "number of interior knots", ylab = "BIC",
    main = "Backward elimination: BIC trajectory"
)
points(tail(res_bwd$history$n_knots, 1L), res_bwd$score,
    pch = 19, col = "red", cex = 1.5
)
abline(h = res_bwd$score, lty = 2, col = "grey50")

## ----echo=TRUE----------------------------------------------------------------
best <- if (res_fwd$score <= res_bwd$score) {
    c(res_fwd$knots, res_fwd$fixed_knots)
} else {
    c(res_bwd$knots, res_bwd$fixed_knots)
}
fit <- fit_glm(rates, sort(best))

# Data-driven population index of the shower, with 95% Wald CI
# (log-scale CI on the limmag coefficient, exponentiated).
limmag_row <- summary(fit)$coefficients["limmag", ]
r_hat <- exp(limmag_row["Estimate"])
r_ci <- exp(limmag_row["Estimate"] +
    c(-1, 1) * 1.96 * limmag_row["Std. Error"])
c(
    estimate = unname(r_hat),
    lower = unname(r_ci[1]),
    upper = unname(r_ci[2])
)

## ----echo=TRUE, fig.alt="Fitted ZHR curve over solar longitude with a 95% confidence band and selected knots marked on the axis"----
sl_grid <- seq(139.2, 140.6, length.out = 400)
pred_df <- data.frame(sl = sl_grid, limmag = 6.5, t = 1)
link <- predict(fit, newdata = pred_df, type = "link", se.fit = TRUE)
zhr <- exp(link$fit)
zhr_lo <- exp(link$fit - 1.96 * link$se.fit)
zhr_hi <- exp(link$fit + 1.96 * link$se.fit)

plot(sl_grid, zhr,
    type = "n",
    ylim = range(c(zhr_lo, zhr_hi)),
    xlab = "solar longitude (deg)", ylab = "ZHR",
    main = "PER 2015 - fitted activity profile (95% CI)"
)
polygon(c(sl_grid, rev(sl_grid)),
    c(zhr_lo, rev(zhr_hi)),
    col = adjustcolor("blue", alpha.f = 0.2), border = NA
)
lines(sl_grid, zhr, col = "blue", lwd = 2)
rug(best, col = "red")

## ----echo=TRUE, fig.alt="Binned observed-over-expected ratios with one-sigma error bars over solar longitude"----
# Use a coarser binning here than for the candidate grid above
# (>= 500 meteors per bin) so that bin uncertainties are small enough
# to read systematic deviations from the diagnostic plot.
rates$qres <- vismeteor::freq_quantile(rates$freq, 500)
rates$pred <- predict(fit, type = "response")

bins <- with(rates, {
    obs_n <- tapply(freq, qres, sum)
    pred_n <- tapply(pred, qres, sum)
    data.frame(
        sl    = tapply(sl * freq, qres, sum) / obs_n,
        ratio = obs_n / pred_n,
        se    = sqrt(obs_n) / pred_n
    )
})
nrow(bins)

plot(bins$sl, bins$ratio,
    pch = 19, col = "blue",
    ylim = range(c(
        bins$ratio - 2 * bins$se,
        bins$ratio + 2 * bins$se
    )),
    xlab = "solar longitude (deg)",
    ylab = "observed / expected",
    main = "PER 2015 - binned residual ratios (+/- 1 sigma)"
)
arrows(bins$sl, bins$ratio - bins$se,
    bins$sl, bins$ratio + bins$se,
    angle = 90, code = 3, length = 0.03, col = "blue"
)
abline(h = 1, lty = 2, col = "grey50")

