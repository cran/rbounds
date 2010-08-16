## ------------------------------------------------------------------ ##
## Authors: Luke Keele, Ohio State University
##          Jason Morgan, Ohio State University 
## ------------------------------------------------------------------ ##

## ------------------------------------------------------------------ ##
## Main function.
## Rt/Rc:   Response for treatment and control.
## Dt/Dc:   Dose received.
## Zt/Zc:   encouragement received (i.e., the instrumental variable).
## b:       Null hypothesis for beta (default = 0).
## alpha:   confidence level of a one-sided hypothesis test.
## Beta:    defines the range over which to test sensitivity.
## BetaInc: value to increment over the defined range.
## ------------------------------------------------------------------ ##
iv_sens <- function(Rt, Rc, Dt, Dc, Zt = rep(1, length(Rt)),
                    Zc = rep(0, length(Rt)), b = 0, alpha = 0.025,
                    Beta = 2, BetaInc = 0.01) {

  ## Length of the response vector for the treated.
  n <- length(Rt)
  
  ## Check that the response and dose vectors have the same length.
  if (length(Rc) != n | length(Dt) != n | length(Dc) != n |
      length(Zt) != n | length(Zc) != n) {
    stop("Length of response, dose, and encouragement vectors not equal.")
  }

  ## Estimate the p-value for all null hypothesis in range as
  ## defined by Beta incremented by BetaInc.
  rng <- seq(-abs(Beta), abs(Beta), by = BetaInc)
  sens <- rep(NA, length(rng))
  for (i in 1:length(rng))
    sens[i] <- iv_pval(Rt, Rc, Dt, Dc, Zt, Zc, b = rng[i])
  Sens <- data.frame("Beta" = rng, "pval" = sens)
  
  ## Select some reasonable intervals from which to begin the
  ## search. Supress warnings, which will be dealt with cleanly below.
  suppressWarnings(x2 <- min(Sens$Beta[Sens$pval > alpha]))
  x3 <- x2 + 2*BetaInc; x1 <- x2 - 2*BetaInc
  suppressWarnings(y2 <- max(Sens$Beta[Sens$pval > alpha]))
  y3 <- y2 + 2*BetaInc; y1 <- y2 - 2*BetaInc

  ## Find the lower and upper bounds for beta only if boundaries were
  ## found from the initial sensitivity estimates.
  if (x2 == -Beta | any(is.infinite(c(x1, x2, x3)))) {
    cat("Lower boundary undefined. Try increasing Beta.\n")
    b1 <- -Inf
  }
  else {
    b1 <- iv_gssearch(x1, x2, x3, maim_dist, Rt, Rc, Dt, Dc, Zt, Zc)
  }

  if (y2 == Beta | any(is.infinite(c(y1, y2, y3)))) {
    cat("Upper boundary undefined. Try increasing Beta.\n")
    b2 <- Inf
  }
  else {
    b2 <- iv_gssearch(y1, y2, y3, maim_dist, Rt, Rc, Dt, Dc, Zt, Zc)
  }

  ## Calculate Hodges-Lehmann point estimate. This should be the Beta
  ## for which the p-value is maximized. Where more than one point is
  ## returned, the mean midpoint is chosen (see Rosenbaum 2010,
  ## p. 136). If one of the bounds is not defined, a warning is
  ## returned and the HL estimate is set to Inf.
  if (any(is.infinite(c(b1, b2)))) {
    cat("Confidence interval not defined. Hodges-Lehmann estimate not calculated.\n")
    HL <- Inf
  }
  else {
    T  <- Sens[Sens$pval == max(Sens$pval),]
    HL <- mean(min(T$Beta), max(T$Beta))
  }
  
  ## Construct and return result.
  ConfInt <- data.frame(lower = b1, upper = b2)
  rownames(ConfInt) <- c("Interval")

  Obj <- list("HL" = HL, "Alpha" = alpha, "Beta" = Beta,
              "BetaInc" = BetaInc, "Sens" = Sens,
              "Interval" = ConfInt)

  class(Obj) <- c("ivsens", class(Obj))
  Obj
}

## This is a cost function that will be passed to iv_gssearch() to
## minimize. It simply calculates the absolute distance between an IV
## p-value and the level of alpha chosen.
maim_dist <- function(Rt, Rc, Dt, Dc, Zt, Zc, b = 0, alpha = 0.025) {
  abs(iv_pval(Rt, Rc, Dt, Dc, Zt, Zc, b = b) - alpha)
}

## Calculate the p-value for the IV estimate given some hypothesized
## beta, b.
iv_pval <- function(Rt, Rc, Dt, Dc, Zt, Zc, b = 0) {
  n <- length(Rt)
  y <- (Zt - Zc)*((Rt - b*Dt) - (Rc - b*Dc))
  y <- y[Rt - Rc != 0]
  r <- rank(abs(y))
  T <- sum(r[y > b])
  lower <- psignrank(T, n, lower.tail = FALSE)
  upper <- psignrank(T, n, lower.tail = TRUE)
  min(lower, upper)
}

## Golden Section Search algorithm. This is a recursive function that
## minimizes some function FUN in the interval [x1, x3]. Including x2,
## is supposed to speed things up, but I can't see how it makes a
## large difference in this case.
##
## Source: Kiefer, J. (1953), "Sequential minimax search for a
## maximum", Proceedings of the American Mathematical Society 4:
## 502–506.
iv_gssearch <- function(x1, x2, x3, FUN, Rt, Rc, Dt, Dc, Zt, Zc,
                      tolerance = 1e-15) {
  tau <- 1.618033988749895  # golden ratio: (1 + sqrt(5))/2

  if (abs(x1 - x3) < tolerance) {
    return((x1+x3)/2)
  }
  else {
    x4 <- x2 + (2-tau)*(x3 - x2)
    f2 <- FUN(Rt, Rc, Dt, Dc, Zt, Zc, b = x2)
    f4 <- FUN(Rt, Rc, Dt, Dc, Zt, Zc, b = x4)
    
    if (f2 < f4) iv_gssearch(x1, x2, x4, FUN, Rt, Rc, Dt, Dc, Zt, Zc,
                           tolerance = tolerance)
    else iv_gssearch(x2, x4, x3, FUN, Rt, Rc, Dt, Dc, Zt, Zc,
                   tolerance = tolerance)
  }
}


## Print and summary functions.
print.ivsens <- function(x, ...) {
  HL <- round(x$HL, digits = 4)
  Interval <- round(x$Interval, digits = 4)
  Sens <- round(x$Sens, digits = 4)

  if (nrow(Sens) > 12) {
    s <- nrow(Sens)
    n <- c(1, floor((s / 12)) * 1:12, s)
    Sens <- Sens[n,]
  }

  cat("\nHodges-Lehmann point estimate for beta ...", HL, "\n")
  cat("\nConfidence interval for beta\n\n")
  print(Interval)
  cat("\nSensitivity matrix (subset)\n\n")
  print(Sens, digits = 4, row.names = FALSE)
  cat("\n")
}

summary.ivsens <- function(object, ...) {
  print.ivsens(object, ...)
}


## Plot function.
plot.ivsens <- function(x, ...) {
  Sens  <- x[["Sens"]]
  Int   <- x[["Interval"]]
  Alpha <- x[["Alpha"]]
  d <- (2*x[["Beta"]]) / 8

  Full <- xyplot(pval ~ Sens$Beta, data = Sens,
                 panel = function(x, y, ...) {
                   panel.abline(h = Alpha, lty = 3, col = "#CB181D")
                   panel.xyplot(x, y, pch = ".", cex = 2.0, ...)
                   panel.abline(v = Int[1,], lty = 2)
                 },
                 main = expression(paste("Full range of ", beta)),
                 xlab = expression(beta),
                 ylab = expression(paste(italic(p),"-value")))

  if (all(is.finite(as.numeric(Int)))) {
    Lower <- xyplot(pval ~ Sens$Beta, data = Sens,
                    subset = Sens$Beta >= (Int[1,1] - d) &
                    Sens$Beta <= (Int[1,1] + d),
                    panel = function(x, y, ...) {
                      panel.abline(h = Alpha, lty = 3, col = "#CB181D")
                      panel.xyplot(x, y, pch = ".", cex = 2.0, ...)
                      if (is.finite(Int[1,1]))
                        panel.abline(v = Int[1,1], lty = 2)
                    },
                    main = "Lower bound", xlab = expression(beta),
                    ylab = expression(paste(italic(p),"-value")))
    
    Upper <- xyplot(pval ~ Sens$Beta, data = Sens,
                    subset = Sens$Beta >= (Int[1,2] - d) &
                    Sens$Beta <= (Int[1,2] + d),
                    panel = function(x, y, ...) {
                      panel.abline(h = Alpha, lty = 3, col = "#CB181D")
                      panel.xyplot(x, y, pch = ".", cex = 2.0, ...)
                      panel.abline(v = Int[1,2], lty = 2)
                    },
                    main = "Upper bound", xlab = expression(beta),
                    ylab = expression(paste(italic(p),"-value")))

    print(Full, split = c(1,1,1,2), more = TRUE)
    print(Lower, position = c(0, 0, 0.5, 1), split = c(1,2,1,2), more = TRUE)
    print(Upper, position = c(0.5, 0, 1, 1), split = c(1,2,1,2), more = FALSE)  
  }
  else {
    print(Full, more = FALSE)
  }
}
