## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(vismeteor)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE----------------------------------------------------------------
plot(
  vmperception,
  -0.5, 8,
  main = "Perception probability of visual meteor magnitudes",
  col = "blue",
  xlab = "limiting magnitude - meteor magnitude",
  ylab = "p"
)

## ----echo=TRUE, results='hide'------------------------------------------------
m <- seq(6, -4, -1)
limmag <- 6.5
r <- 2.0
p <- vismeteor::dvmgeom(m, limmag, r)
barplot(
    p,
    names.arg = m,
    main = paste0('Density (r = ', r, ', limmag = ', limmag, ')'),
    col = "blue",
    xlab = 'm',
    ylab = 'p',
    border = "blue",
    space = 0.5
)
axis(side = 2, at = pretty(p))

## ----echo=TRUE, results='hide'------------------------------------------------
m <- seq(6, -4, -1)
psi <- 5.0
limmag <- 6.5
p <- vismeteor::dvmideal(m, limmag, psi)
barplot(
    p,
    names.arg = m,
    main = paste0('Density (psi = ', psi, ', limmag = ', limmag, ')'),
    col = "blue",
    xlab = 'm',
    ylab = 'p',
    border = "blue",
    space = 0.5
)
axis(side = 2, at = pretty(p))

## ----echo=TRUE----------------------------------------------------------------
mt <- as.table(matrix(
    c(
        0.0, 0.0, 2.5, 0.5, 0.0, 1.0,
        0.0, 1.5, 2.0, 0.5, 0.0, 0.0,
        1.0, 0.0, 0.0, 3.0, 2.5, 0.5
    ), nrow = 3, ncol = 6, byrow = TRUE
))
colnames(mt) <- seq(6)
rownames(mt) <- c('A', 'B', 'C')
margin.table(mt, 1)
margin.table(mt, 2)

# contingency table with integer values
(mt.int <- vmtable(mt))
margin.table(mt.int, 1)
margin.table(mt.int, 2)

## ----echo=TRUE----------------------------------------------------------------
freq <- c(1,8,3,3,4,9,5,0,0,2,7,8,2,6,4)
f <- freq.quantile(freq, 10)
print(f)
print(tapply(freq, f, sum))

