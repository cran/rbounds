binarysens <- function(x, y, Gamma=6, GammaInc=1) {
  D1 <- x 
  D0 <- y
  gamma <- seq(1, Gamma, by=GammaInc)
  mx <- D1 + D0
  if(D1 > D0) { D <- D1
  } else { D <- D0 }
  up <- c()
  lo <- c()
  series <- seq(D, mx, by=1)
  n.it <- length(gamma)
  for(i in 1:n.it) {
    p.plus <- gamma[i]/(1 + gamma[i])
    p.minus <- 1/(1 + gamma[i])
    up.tmp <- 1 - pbinom(D, mx, prob=p.plus)
    lo.tmp <- 1 - pbinom(D, mx, prob=p.minus)
    
    up <- c(up, up.tmp)
    lo <- c(lo, lo.tmp)
  }
  
  pval <- lo[1]
  bounds <- data.frame(gamma, round(lo, 5), round(up, 5))
  colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
  
  msg <- "Rosenbaum Sensitivity Test \n"
  note <- "Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors \n"
  
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  
  Obj
}

