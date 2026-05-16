## ----setup, include = FALSE---------------------------------------------------
library(vismeteor)
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----echo=TRUE, results='hide'------------------------------------------------
observations <- with(PER_2015_magn$observations, {
    idx <- !is.na(lim_magn) & sl_start > 135.81 & sl_end < 135.87
    data.frame(
        magn_id = magn_id[idx],
        lim_magn = lim_magn[idx]
    )
})
head(observations, 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(observations, 5))

## ----echo=TRUE, results='hide'------------------------------------------------
magnitudes <- with(new.env(), {
    magnitudes <- merge(
        observations,
        as.data.frame(PER_2015_magn$magnitudes),
        by = "magn_id"
    )
    magnitudes$magn <- as.integer(as.character(magnitudes$magn))
    subset(magnitudes, (magnitudes$lim_magn - magnitudes$magn) > -0.5)
})
head(magnitudes, 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(magnitudes[magnitudes$Freq > 0, ], 5))

## ----echo=TRUE, results='hide'------------------------------------------------
# maximum likelihood estimation (MLE) of r
result_ml <- with(magnitudes, {
    # log likelihood function
    ll <- function(r) -sum(Freq * dvmgeom(magn, lim_magn, r, log = TRUE))
    r_start <- 2.0 # starting value
    r_lower <- 1.2 # lowest expected value
    r_upper <- 4.0 # highest expected value
    # find minimum
    optim(r_start, ll, method = "Brent", lower = r_lower, upper = r_upper, hessian = TRUE)
})

## ----echo=TRUE----------------------------------------------------------------
print(result_ml$par) # mean of r
print(1 / result_ml$hessian[1][1]) # variance of r

## ----echo=TRUE, results='hide'------------------------------------------------
with(new.env(), {
    data_plot <- data.frame(r = seq(2.0, 2.8, 0.01))
    data_plot$ll <- mapply(function(r) {
        with(magnitudes, {
            # log likelihood function
            sum(Freq * dvmgeom(magn, lim_magn, r, log = TRUE))
        })
    }, data_plot$r)
    data_plot$l <- exp(data_plot$ll - max(data_plot$ll))
    data_plot$l <- data_plot$l / sum(data_plot$l)
    brks <- seq(min(data_plot$r) - 0.02, max(data_plot$r) + 0.02, by = 0.02)
    plot(data_plot$r, data_plot$l,
        # breaks = brks,
        type = "l",
        col = "blue",
        xlab = "r",
        xaxt = "n",
        ylab = "likelihood"
    )
    xlabels <- seq(min(round(data_plot$r, 1)) - 0.1, max(round(data_plot$r, 1)) + 0.1, by = 0.1)
    axis(
        side = 1,
        at = xlabels,
        labels = sprintf("%.1f", xlabels)
    )
    abline(v = result_ml$par, col = "red", lwd = 1)
})

## ----echo=TRUE, results='hide'------------------------------------------------
result_rate <- with(magnitudes, {
    N <- sum(Freq)
    a <- vmperception(lim_magn - magn - 1) / vmperception(lim_magn - magn)
    a_mean <- as.numeric(weighted.mean(a, w = Freq))
    a_var <- as.numeric(cov.wt(cbind(a), wt = Freq)$cov) / N
    # apply the delta method and return the result
    list(
        "mean" = 1 / a_mean - a_var / a_mean^3,
        "var" = a_var / a_mean^4
    )
})

## ----echo=TRUE----------------------------------------------------------------
print(result_rate$mean) # mean of r
print(result_rate$var) # variance of r

## ----echo=TRUE, results='hide'------------------------------------------------
# Bootstrapping Method
r_means <- with(magnitudes, {
    N <- sum(Freq)
    a <- vmperception(lim_magn - magn - 1) / vmperception(lim_magn - magn)
    replicate(50000, {
        s <- sample(a, size = N, replace = TRUE, prob = Freq)
        s_mean <- mean(s)
        s_var <- var(s) / N
        1 / s_mean - s_var / s_mean^3
    })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
    r_sd <- sqrt(result_rate$var)
    r_min <- result_rate$mean - 3 * r_sd
    r_max <- result_rate$mean + 3 * r_sd
    r <- subset(r_means, r_means > r_min & r_means < r_max)
    brks <- seq(min(r) - 0.02, max(r) + 0.02, by = 0.02)
    hist(r,
        breaks = brks,
        col = "skyblue",
        border = "black",
        main = "Histogram of mean r",
        xlab = "r",
        xaxt = "n",
        ylab = "count"
    )
    xlabels <- seq(min(round(r, 1)) - 0.1, max(round(r, 1)) + 0.1, by = 0.1)
    axis(
        side = 1,
        at = xlabels,
        labels = sprintf("%.1f", xlabels)
    )
    abline(v = result_rate$mean, col = "red", lwd = 1)
})

## ----echo=TRUE, results='show'------------------------------------------------
result_vs <- with(magnitudes, {
    N <- sum(Freq)
    tm <- vmgeom_vst_from_magn(magn, lim_magn)
    tm_mean <- sum(Freq * tm) / N
    tm_var <- sum(Freq * (tm - tm_mean)^2) / (N - 1)
    tm_mean_var <- tm_var / N

    # Delta method: variance and variance of r
    r_hat <- vmgeom_vst_to_r(tm_mean)
    dr_dtm <- vmgeom_vst_to_r(tm_mean, deriv_degree = 1L)
    d2r_dtm2 <- vmgeom_vst_to_r(tm_mean, deriv_degree = 2L)
    r_hat <- r_hat - 0.5 * d2r_dtm2 * tm_mean_var
    r_var <- dr_dtm^2 * tm_mean_var + 0.5 * d2r_dtm2^2 * tm_mean_var^2

    list("mean" = r_hat, "var" = r_var)
})

## ----echo=TRUE, results='show'------------------------------------------------
print(paste("r mean:", result_vs$mean))
print(paste("r var:", result_vs$var))

## ----echo=TRUE, results='hide'------------------------------------------------
# Bootstrapping Method
r_means <- with(magnitudes, {
    N <- sum(Freq)
    tm <- vmgeom_vst_from_magn(magn, lim_magn)
    replicate(50000, {
        s <- sample(tm, size = N, replace = TRUE, prob = Freq)
        s_mean <- mean(s)
        s_var <- var(s) / N

        r_hat <- vmgeom_vst_to_r(s_mean)
        d2r_ds2 <- vmgeom_vst_to_r(s_mean, deriv_degree = 2L)
        r_hat - 0.5 * d2r_ds2 * s_var
    })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
    r_sd <- sqrt(result_vs$var)
    r_min <- as.vector(result_vs$mean - 3 * r_sd)
    r_max <- as.vector(result_vs$mean + 3 * r_sd)
    r <- subset(r_means, r_means > r_min & r_means < r_max)
    brks <- seq(min(r) - 0.02, max(r) + 0.02, by = 0.02)
    hist(r,
        breaks = brks,
        col = "skyblue",
        border = "black",
        main = "Histogram of mean r",
        xlab = "r",
        xaxt = "n",
        ylab = "count"
    )
    xlabels <- seq(min(round(r, 1)) - 0.1, max(round(r, 1)) + 0.1, by = 0.1)
    axis(
        side = 1,
        at = xlabels,
        labels = sprintf("%.1f", xlabels)
    )
    abline(v = result_vs$mean, col = "red", lwd = 1)
})

## ----echo=TRUE, results='hide'------------------------------------------------
magnitudes$p <- with(magnitudes, dvmgeom(m = magn, lm = lim_magn, result_rate$mean))

## ----echo=TRUE, results='hide'------------------------------------------------
magn_min <- min(magnitudes$magn)

## ----echo=TRUE, results='asis'------------------------------------------------
idx <- magnitudes$magn == magn_min
magnitudes$p[idx] <- with(
    magnitudes[idx, ],
    pvmgeom(m = magn + 1L, lm = lim_magn, result_rate$mean, lower.tail = TRUE)
)

## ----echo=TRUE----------------------------------------------------------------
magnitutes_observed <- xtabs(Freq ~ magn_id + magn, data = magnitudes)
magnitutes_observed_mt <- margin.table(magnitutes_observed, margin = 2)
print(magnitutes_observed_mt)

## ----echo=TRUE----------------------------------------------------------------
magnitudes$magn[magnitudes$magn <= 0] <- "0-"
magnitudes$magn[magnitudes$magn >= 4] <- "4+"
magnitutes_observed <- xtabs(Freq ~ magn_id + magn, data = magnitudes)
print(margin.table(magnitutes_observed, margin = 2))

## ----echo=TRUE----------------------------------------------------------------
magnitutes_expected <- xtabs(p ~ magn_id + magn, data = magnitudes)
magnitutes_row_freq <- margin.table(magnitutes_observed, margin = 1)
magnitutes_expected <- sweep(magnitutes_expected, 1, magnitutes_row_freq, `*`)
magnitutes_expected <- magnitutes_expected / sum(magnitutes_expected)
print(sum(magnitudes$Freq) * margin.table(magnitutes_expected, margin = 2))

## ----echo=TRUE, results='asis'------------------------------------------------
chisq_test_result <- chisq.test(
    x = margin.table(magnitutes_observed, margin = 2),
    p = margin.table(magnitutes_expected, margin = 2)
)

## ----echo=TRUE----------------------------------------------------------------
chi2_df <- chisq_test_result$parameter - 1
chi2_pval <- pchisq(chisq_test_result$statistic, df = chi2_df, lower.tail = FALSE)
print(chi2_pval)

## ----fig.show='hold'----------------------------------------------------------
chisq_test_residuals <- with(new.env(), {
    chisq_test_residuals <- residuals(chisq_test_result)
    v <- as.vector(chisq_test_residuals)
    names(v) <- names(chisq_test_residuals)
    v
})
plot(
    chisq_test_residuals,
    main = "Residuals of the chi-square goodness-of-fit test",
    xlab = "m",
    ylab = "Residuals",
    ylim = c(-3, 3),
    xaxt = "n"
)
abline(h = 0.0, lwd = 2)
axis(1, at = seq_along(chisq_test_residuals), labels = names(chisq_test_residuals))

