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
head(magnitudes[magnitudes$Freq>0,], 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(magnitudes[magnitudes$Freq>0,], 5))

## ----echo=TRUE, results='hide'------------------------------------------------
# maximum likelihood estimation (MLE) of psi
result <- with(magnitudes, {
    # log likelihood function
    ll <- function(psi) -sum(Freq * dvmideal(magn, lim.magn, psi, log=TRUE))
    psi.start <- 6.0 # starting value
    psi.lower <- 4.0 # lowest expected value
    psi.upper <- 10.0 # highest expected value
    # find minimum
    optim(psi.start, ll, method='Brent', lower=psi.lower, upper=psi.upper, hessian=TRUE)
})

## ----echo=TRUE----------------------------------------------------------------
psi.mean <- result$par # mean of psi
print(psi.mean)
psi.var <- 1/result$hessian[1][1] # variance of psi
print(psi.var)

## ----echo=TRUE, results='hide'------------------------------------------------
with(new.env(), {
  data.plot <- data.frame(psi = seq(4.0, 11, 0.1))
  data.plot$ll <- mapply(function(psi){
    with(magnitudes, {
      # log likelihood function
      sum(Freq * dvmideal(magn, lim.magn, psi, log=TRUE))
    })
  }, data.plot$psi)
  data.plot$l <- exp(data.plot$ll - max(data.plot$ll))
  data.plot$l <- data.plot$l / sum(data.plot$l)
  plot(data.plot$psi, data.plot$l,
    type = "l",
    col = "blue",
    xlab = "psi",
    ylab = "likelihood"
  )
  abline(v = result$par, col = "red", lwd = 1)
})

## ----echo=TRUE, results='show'------------------------------------------------
tm.mean <- with(magnitudes, {
  N <- sum(Freq)
  tm <- vmidealVstFromMagn(magn, lim.magn)
  tm.mean <- sum(Freq * tm)/N
  tm.var <- sum(Freq * (tm - tm.mean)^2)/(N-1)
  tm.mean.var <- tm.var / N
  list('val' = tm.mean, 'var' = tm.mean.var, 'sd' = sqrt(tm.mean.var))
})

## ----echo=TRUE, results='show'------------------------------------------------
print(paste('tm mean:', tm.mean$val))
print(paste('tm var:', tm.mean$var))

## ----echo=TRUE, results='show'------------------------------------------------
tm.means <- with(magnitudes, {
    N <- sum(Freq)
    tm <- vmidealVstFromMagn(magn, lim.magn)
    replicate(50000, {
      mean(sample(tm, size = N, replace = TRUE, prob = Freq))
    })
})

## ----echo=TRUE, results='show'------------------------------------------------
with(new.env(), {
  tm.min <- tm.mean$val - 3 * tm.mean$sd
  tm.max <- tm.mean$val + 3 * tm.mean$sd
  tm.means <- subset(tm.means, tm.means > tm.min & tm.means < tm.max)
  brks <- seq(min(tm.means) - 0.02, max(tm.means) + 0.02, by = 0.02)
  hist(tm.means,
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
lim.magn.mean <- with(magnitudes, {
    N <- sum(Freq)
    sum(Freq * lim.magn)/N
})
print(paste('lim.magn.mean:', lim.magn.mean))
print(paste('mean psi:', vmidealVstToPsi(tm.mean$val, lim.magn.mean)))

## ----echo=TRUE----------------------------------------------------------------
print(vmidealVstToPsi(qnorm(0.90, tm.mean$val, tm.mean$sd), lim.magn.mean))

## ----echo=TRUE, results='asis'------------------------------------------------
psi.mean <- vmidealVstToPsi(tm.mean$val, lim.magn.mean)
magnitudes$p <- with(magnitudes, dvmideal(m = magn, lm = lim.magn, psi.mean))

## ----echo=TRUE, results='hide'------------------------------------------------
magn.min <- min(magnitudes$magn)

## ----echo=TRUE, results='asis'------------------------------------------------
idx <- magnitudes$magn == magn.min
magnitudes$p[idx] <- with(
    magnitudes[idx,],
    pvmideal(m = magn + 1L, lm = lim.magn, psi.mean, lower.tail = TRUE)
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
magnitutes.expected <- magnitutes.expected/nrow(magnitutes.expected)
print(sum(magnitudes$Freq) * margin.table(magnitutes.expected, margin = 2))

## ----echo=TRUE, results='asis'------------------------------------------------
chisq.test.result <- chisq.test(
    x = margin.table(magnitutes.observed, margin = 2),
    p = margin.table(magnitutes.expected, margin = 2)
)

## ----echo=TRUE----------------------------------------------------------------
print(chisq.test.result$p.value)

## ----fig.show='hold'----------------------------------------------------------
chisq.test.residuals <- with(new.env(), {
    chisq.test.residuals <- residuals(chisq.test.result)
    v <- as.vector(chisq.test.residuals)
    names(v) <- rownames(chisq.test.residuals)
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

