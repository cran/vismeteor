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
head(magnitudes[magnitudes$Freq > 0, ], 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(magnitudes[magnitudes$Freq > 0, ], 5))

## ----echo=TRUE, results='hide'------------------------------------------------
# maximum likelihood estimation (MLE) of psi
result <- with(magnitudes, {
    # log likelihood function
    ll <- function(psi) -sum(Freq * dvmideal(magn, lim_magn, psi, log = TRUE))
    psi_start <- 6.0 # starting value
    psi_lower <- 4.0 # lowest expected value
    psi_upper <- 10.0 # highest expected value
    # find minimum
    optim(psi_start, ll, method = "Brent", lower = psi_lower, upper = psi_upper, hessian = TRUE)
})

## ----echo=TRUE----------------------------------------------------------------
psi_mean <- result$par # mean of psi
print(psi_mean)
psi_var <- 1 / result$hessian[1][1] # variance of psi
print(psi_var)

## ----echo=TRUE, results='hide'------------------------------------------------
with(new.env(), {
    data_plot <- data.frame(psi = seq(4.0, 11, 0.1))
    data_plot$ll <- mapply(function(psi) {
        with(magnitudes, {
            # log likelihood function
            sum(Freq * dvmideal(magn, lim_magn, psi, log = TRUE))
        })
    }, data_plot$psi)
    data_plot$l <- exp(data_plot$ll - max(data_plot$ll))
    data_plot$l <- data_plot$l / sum(data_plot$l)
    plot(data_plot$psi, data_plot$l,
        type = "l",
        col = "blue",
        xlab = "psi",
        ylab = "likelihood"
    )
    abline(v = result$par, col = "red", lwd = 1)
})

## ----echo=TRUE, results='show'------------------------------------------------
tm_mean <- with(magnitudes, {
    N <- sum(Freq)
    tm <- vmideal_vst_from_magn(magn, lim_magn)
    tm_mean <- sum(Freq * tm) / N
    tm_var <- sum(Freq * (tm - tm_mean)^2) / (N - 1)
    tm_mean_var <- tm_var / N
    list("val" = tm_mean, "var" = tm_mean_var, "sd" = sqrt(tm_mean_var))
})

## ----echo=TRUE, results='show'------------------------------------------------
print(paste("tm mean:", tm_mean$val))
print(paste("tm var:", tm_mean$var))

## ----echo=TRUE, results='show'------------------------------------------------
tm_means <- with(magnitudes, {
    N <- sum(Freq)
    tm <- vmideal_vst_from_magn(magn, lim_magn)
    replicate(50000, {
        mean(sample(tm, size = N, replace = TRUE, prob = Freq))
    })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
    tm_min <- tm_mean$val - 3 * tm_mean$sd
    tm_max <- tm_mean$val + 3 * tm_mean$sd
    tm_means <- subset(tm_means, tm_means > tm_min & tm_means < tm_max)
    brks <- seq(min(tm_means) - 0.02, max(tm_means) + 0.02, by = 0.02)
    hist(tm_means,
        breaks = brks,
        col = "skyblue",
        border = "black",
        main = "Histogram of mean tm",
        xlab = "tm",
        ylab = "count",
        xaxt = "n"
    )
    axis(1, at = seq(round(min(brks), 1), round(max(brks), 1) + 0.1, by = 0.1))
    abline(v = 0, col = "red", lwd = 1)
})

## ----echo=TRUE----------------------------------------------------------------
lim_magn_mean <- with(magnitudes, {
    N <- sum(Freq)
    sum(Freq * lim_magn) / N
})
print(paste("lim_magn_mean:", lim_magn_mean))
print(paste("mean psi:", vmideal_vst_to_psi(tm_mean$val, lim_magn_mean)))

## ----echo=TRUE----------------------------------------------------------------
print(vmideal_vst_to_psi(qnorm(0.90, tm_mean$val, tm_mean$sd), lim_magn_mean))

## ----echo=TRUE, results='asis'------------------------------------------------
psi_mean <- vmideal_vst_to_psi(tm_mean$val, lim_magn_mean)
magnitudes$p <- with(magnitudes, dvmideal(m = magn, lm = lim_magn, psi_mean))

## ----echo=TRUE, results='hide'------------------------------------------------
magn_min <- min(magnitudes$magn)

## ----echo=TRUE, results='asis'------------------------------------------------
idx <- magnitudes$magn == magn_min
magnitudes$p[idx] <- with(
    magnitudes[idx, ],
    pvmideal(m = magn + 1L, lm = lim_magn, psi_mean, lower.tail = TRUE)
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

