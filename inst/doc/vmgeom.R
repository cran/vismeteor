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
magnitudes <- merge(
  observations,
  as.data.frame(PER_2015_magn$magnitudes),
  by = 'magn.id'
)
magnitudes$magn <- as.integer(as.vector(magnitudes$magn))
head(magnitudes[magnitudes$Freq>0,], 5) # Example values

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(head(magnitudes[magnitudes$Freq>0,], 5))

## ----echo=TRUE, results='hide'------------------------------------------------
# maximum likelihood estimation (MLE) of r
result <- with(subset(magnitudes, (magnitudes$lim.magn - magnitudes$magn) > -0.5), {
    # log likelihood function
    ll <- function(r) -sum(Freq * dvmgeom(magn, lim.magn, r, log=TRUE))
    r.start <- 2.0 # starting value
    r.lower <- 1.2 # lowest expected value
    r.upper <- 4.0 # highest expected value
    # find minimum
    optim(r.start, ll, method='Brent', lower=r.lower, upper=r.upper, hessian=TRUE)
})

## ----echo=TRUE----------------------------------------------------------------
r.mean <- result$par # mean of r
print(r.mean)
r.var <- 1/result$hessian[1][1] # variance of r
print(r.var)

## ----echo=TRUE----------------------------------------------------------------
m.mean <- with(magnitudes, sum((lim.magn - magn) * Freq)/sum(Freq))
print(m.mean)

## ----echo=TRUE----------------------------------------------------------------
m.var <- with(magnitudes, {
    n <- sum(Freq)
    sum((lim.magn - magn - m.mean)^2 * Freq)/((n-1) * n)
})
print(m.var)

## ----echo=TRUE, results='hide'------------------------------------------------
r.mean.fun <- with(new.env(), {
    r <- seq(1.3, 3.5, 0.1)
    s <- log(r)
    m.mean <- -vmperception.l(s, deriv.degree = 1L)/vmperception.l(s)
    splinefun(m.mean, r)
})

## ----echo=TRUE----------------------------------------------------------------
r.mean <- r.mean.fun(m.mean)
print(r.mean)

## ----echo=TRUE----------------------------------------------------------------
r.var <- r.mean.fun(m.mean, deriv = 1L)^2 * m.var
print(r.var)

## ----echo=TRUE, results='hide'------------------------------------------------
magnitudes$p <- with(magnitudes, dvmgeom(m = magn, lm = lim.magn, r.mean))

## ----echo=TRUE, results='hide'------------------------------------------------
magn.min <- min(magnitudes$magn)

## ----echo=TRUE, results='asis'------------------------------------------------
idx <- magnitudes$magn == magn.min
magnitudes$p[idx] <- with(
    magnitudes[idx,],
    pvmgeom(m = magn + 1L, lm = lim.magn, r.mean, lower.tail = TRUE)
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

