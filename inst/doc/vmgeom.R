## ----setup, include = FALSE---------------------------------------------------
library(vismeteor)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, results='hide'------------------------------------------------
observations <- with(PER_2015_magn$observations, {
    idx <- !is.na(lim.magn) & sl.start > 135.81 & sl.end < 135.87
    data.frame(
        magn.id = magn.id[idx],
        lim.magn = lim.magn[idx]
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
    by = 'magn.id'
  )
  magnitudes$magn <- as.integer(as.character(magnitudes$magn))
  subset(magnitudes, (magnitudes$lim.magn - magnitudes$magn) > -0.5)
})
head(magnitudes, 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(magnitudes[magnitudes$Freq>0,], 5))

## ----echo=TRUE, results='hide'------------------------------------------------
# maximum likelihood estimation (MLE) of r
result.ml <- with(magnitudes, {
    # log likelihood function
    ll <- function(r) -sum(Freq * dvmgeom(magn, lim.magn, r, log=TRUE))
    r.start <- 2.0 # starting value
    r.lower <- 1.2 # lowest expected value
    r.upper <- 4.0 # highest expected value
    # find minimum
    optim(r.start, ll, method='Brent', lower=r.lower, upper=r.upper, hessian=TRUE)
})

## ----echo=TRUE----------------------------------------------------------------
print(result.ml$par) # mean of r
print(1/result.ml$hessian[1][1]) # variance of r

## ----echo=TRUE, results='hide'------------------------------------------------
with(new.env(), {
  data.plot <- data.frame(r = seq(2.0, 2.8, 0.01))
  data.plot$ll <- mapply(function(r){
    with(magnitudes, {
      # log likelihood function
      sum(Freq * dvmgeom(magn, lim.magn, r, log = TRUE))
    })
  }, data.plot$r)
  data.plot$l <- exp(data.plot$ll - max(data.plot$ll))
  data.plot$l <- data.plot$l / sum(data.plot$l)
  brks <- seq(min(data.plot$r) - 0.02, max(data.plot$r) + 0.02, by = 0.02)
  plot(data.plot$r, data.plot$l,
    #breaks = brks,
    type = "l",
    col = "blue",
    xlab = "r",
    xaxt = "n",
    ylab = "likelihood"
  )
  xlabels = seq(min(round(data.plot$r, 1)) - 0.1, max(round(data.plot$r, 1)) + 0.1, by = 0.1)
  axis(
    side = 1,
    at = xlabels,
    labels = sprintf("%.1f", xlabels)
  )
  abline(v = result.ml$par, col = "red", lwd = 1)
})

## ----echo=TRUE, results='hide'------------------------------------------------
result.rate <- with(magnitudes, {
    N <- sum(Freq)
    a <- vmperception(lim.magn - magn - 1)/vmperception(lim.magn - magn)
    a.mean <- as.numeric(weighted.mean(a, w = Freq))
    a.var <- as.numeric(cov.wt(cbind(a), wt = Freq)$cov) / N
    # apply the delta method and return the result
    list(
        'mean' = 1/a.mean - a.var/a.mean^3,
        'var' = a.var/a.mean^4
    )
})

## ----echo=TRUE----------------------------------------------------------------
print(result.rate$mean) # mean of r
print(result.rate$var) # variance of r

## ----echo=TRUE, results='hide'------------------------------------------------
# Bootstrapping Method
r.means <- with(magnitudes, {
  N <- sum(Freq)
  a <- vmperception(lim.magn - magn - 1)/vmperception(lim.magn - magn)
  replicate(50000, {
    s <- sample(a, size = N, replace = TRUE, prob = Freq)
    s.mean <- mean(s)
    s.var <- var(s)/N
    1/s.mean - s.var/s.mean^3
  })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
  r.sd <- sqrt(result.rate$var)
  r.min <- result.rate$mean - 3 * r.sd
  r.max <- result.rate$mean + 3 * r.sd
  r <- subset(r.means, r.means > r.min & r.means < r.max)
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
  xlabels = seq(min(round(r, 1)) - 0.1, max(round(r, 1)) + 0.1, by = 0.1)
  axis(
    side = 1,
    at = xlabels,
    labels = sprintf("%.1f", xlabels)
  )
  abline(v = result.rate$mean, col = "red", lwd = 1)
})

## ----echo=TRUE, results='show'------------------------------------------------
result.vs <- with(magnitudes, {
  N <- sum(Freq)
  tm <- vmgeomVstFromMagn(magn, lim.magn)
  tm.mean <- sum(Freq * tm)/N
  tm.var <- sum(Freq * (tm - tm.mean)^2)/(N-1)
  tm.mean.var <- tm.var / N

  # Delta method: variance and variance of r
  r.hat <- vmgeomVstToR(tm.mean)
  dr_dtm <- vmgeomVstToR(tm.mean, deriv.degree = 1L)
  d2r_dtm2 <- vmgeomVstToR(tm.mean, deriv.degree = 2L)
  r.hat <- r.hat - 0.5 * d2r_dtm2 * tm.mean.var
  r.var <- dr_dtm^2 * tm.mean.var + 0.5 * d2r_dtm2^2 * tm.mean.var^2

  list('mean' = r.hat, 'var' = r.var)
})

## ----echo=TRUE, results='show'------------------------------------------------
print(paste('r mean:', result.vs$mean))
print(paste('r var:', result.vs$var))

## ----echo=TRUE, results='hide'------------------------------------------------
# Bootstrapping Method
r.means <- with(magnitudes, {
  N <- sum(Freq)
  tm <- vmgeomVstFromMagn(magn, lim.magn)
  replicate(50000, {
    s <- sample(tm, size = N, replace = TRUE, prob = Freq)
    s.mean <- mean(s)
    s.var <- var(s)/N

    r.hat <- vmgeomVstToR(s.mean)
    d2r_ds2 <- vmgeomVstToR(s.mean, deriv.degree = 2L)
    r.hat - 0.5 * d2r_ds2 * s.var
  })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
  r.sd <- sqrt(result.vs$var)
  r.min <- as.vector(result.vs$mean - 3 * r.sd)
  r.max <- as.vector(result.vs$mean + 3 * r.sd)
  r <- subset(r.means, r.means > r.min & r.means < r.max)
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
  xlabels = seq(min(round(r, 1)) - 0.1, max(round(r, 1)) + 0.1, by = 0.1)
  axis(
    side = 1,
    at = xlabels,
    labels = sprintf("%.1f", xlabels)
  )
  abline(v = result.vs$mean, col = "red", lwd = 1)
})

## ----echo=TRUE, results='hide'------------------------------------------------
magnitudes$p <- with(magnitudes, dvmgeom(m = magn, lm = lim.magn, result.rate$mean))

## ----echo=TRUE, results='hide'------------------------------------------------
magn.min <- min(magnitudes$magn)

## ----echo=TRUE, results='asis'------------------------------------------------
idx <- magnitudes$magn == magn.min
magnitudes$p[idx] <- with(
    magnitudes[idx,],
    pvmgeom(m = magn + 1L, lm = lim.magn, result.rate$mean, lower.tail = TRUE)
)

## ----echo=TRUE----------------------------------------------------------------
magnitutes.observed <- xtabs(Freq ~ magn.id + magn, data = magnitudes)
magnitutes.observed.mt <- margin.table(magnitutes.observed, margin = 2) 
print(magnitutes.observed.mt)

## ----echo=TRUE----------------------------------------------------------------
magnitudes$magn[magnitudes$magn <= 0] <- '0-'
magnitudes$magn[magnitudes$magn >= 4] <- '4+'
magnitutes.observed <- xtabs(Freq ~ magn.id + magn, data = magnitudes)
print(margin.table(magnitutes.observed, margin = 2))

## ----echo=TRUE----------------------------------------------------------------
magnitutes.expected <- xtabs(p ~ magn.id + magn, data = magnitudes)
magnitutes.row_freq <- margin.table(magnitutes.observed, margin = 1)
magnitutes.expected <- sweep(magnitutes.expected, 1, magnitutes.row_freq, `*`)
magnitutes.expected <- magnitutes.expected/sum(magnitutes.expected)
print(sum(magnitudes$Freq) * margin.table(magnitutes.expected, margin = 2))

## ----echo=TRUE, results='asis'------------------------------------------------
chisq.test.result <- chisq.test(
    x = margin.table(magnitutes.observed, margin = 2),
    p = margin.table(magnitutes.expected, margin = 2)
)

## ----echo=TRUE----------------------------------------------------------------
chi2.df <- chisq.test.result$parameter - 1
chi2.pval <- pchisq(chisq.test.result$statistic, df = chi2.df, lower.tail = FALSE)
print(chi2.pval)

## ----fig.show='hold'----------------------------------------------------------
chisq.test.residuals <- with(new.env(), {
    chisq.test.residuals <- residuals(chisq.test.result)
    v <- as.vector(chisq.test.residuals)
    names(v) <- names(chisq.test.residuals)
    v
})
plot(
    chisq.test.residuals,
    main="Residuals of the chi-square goodness-of-fit test",
    xlab="m",
    ylab="Residuals",
    ylim=c(-3, 3),
    xaxt = "n"
)
abline(h=0.0, lwd=2)
axis(1, at = seq_along(chisq.test.residuals), labels = names(chisq.test.residuals))

