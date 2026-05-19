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
set.seed(1)
x <- seq(0, 10, length.out = 200)
y <- rpois(length(x), lambda = 50 + 30 * sin(x))
dat <- data.frame(x = x, y = y)

fit <- \(d, knots) {
    f <- if (length(knots) == 0L) {
        y ~ splines::ns(x)
    } else {
        y ~ splines::ns(x, knots = knots)
    }
    glm(f, family = poisson(), data = d)
}
score_bic <- \(d, knots) fit(d, knots) |> BIC()

result <- select_knots(
    dat,
    knot_candidates = seq(1, 9, by = 0.5),
    score_fun = score_bic
)
final_knots <- sort(c(result$knots, result$fixed_knots))
final_fit <- fit(dat, final_knots)

op <- par(mar = c(7, 4, 4, 2) + 0.1)
plot(
    x, y,
    pch = 20, col = "gray",
    main = "Poisson spline fit with selected knots",
    xlab = "x", ylab = "count"
)
lines(x, 50 + 30 * sin(x), col = "darkgreen", lwd = 2, lty = 2)
lines(x, final_fit |> predict(type = "response"), col = "blue", lwd = 2)
abline(v = final_knots, col = "red", lty = 3)
legend(
    "bottom",
    inset = c(0, -0.35),
    legend = c("true mean", "spline fit", "selected knots"),
    col = c("darkgreen", "blue", "red"),
    lty = c(2, 1, 3),
    lwd = c(2, 2, 1),
    horiz = TRUE,
    bty = "n",
    xpd = TRUE
)
par(op)

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
    main = paste0("Density (r = ", r, ", limmag = ", limmag, ")"),
    col = "blue",
    xlab = "m",
    ylab = "p",
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
    main = paste0("Density (psi = ", psi, ", limmag = ", limmag, ")"),
    col = "blue",
    xlab = "m",
    ylab = "p",
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
    ),
    nrow = 3, ncol = 6, byrow = TRUE
))
colnames(mt) <- seq(6)
rownames(mt) <- c("A", "B", "C")
margin.table(mt, 1)
margin.table(mt, 2)

# contingency table with integer values
(mt_int <- vmtable(mt))
margin.table(mt_int, 1)
margin.table(mt_int, 2)

## ----echo=TRUE----------------------------------------------------------------
freq <- c(1, 8, 3, 3, 4, 9, 5, 0, 0, 2, 7, 8, 2, 6, 4)
f <- freq_quantile(freq, 10)
print(f)
print(tapply(freq, f, sum))

