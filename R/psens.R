
psens <- function(x,y=NULL, Gamma=6, GammaInc=1){
  
  if (is.numeric(x)){
  trt <- x
  ctrl <- y
  } else if(x$est > 0){
  ctrl <- Y[x$index.control]
  trt <- Y[x$index.treated]
	} else {
	ctrl <- Y[x$index.treated]
    trt <- Y[x$index.control]	
		}

  gamma <- seq(1, Gamma, by=GammaInc)
  m <- length(gamma)
  S <- length(ctrl)
  pvals <- matrix(NA, m, 2)
  diff <- trt - ctrl
  ranks <- rank(abs(diff), ties.method="average")
  psi <- as.numeric(diff > 0)
  T <- sum(psi * ranks)
  
  if(length(unique(ranks)) == length(ranks)){
    
    for(i in 1:m){
      p.plus <- gamma[i]/(1 + gamma[i])
      p.minus <- 1/(1+gamma[i])
      E.T.plus <- (p.plus*S*(S+1))/2
      V.T.plus <- p.plus*(1 - p.plus)*(S*(S+1)*(2*(S+1)))/6
      V.T.minus <- p.minus*(1 - p.minus)*(S*(S+1)*(2*(S+1)))/6
      E.T.minus <- (p.minus*S*(S+1))/2
      #Normal Approximation
      z.plus <- (T - E.T.plus)/sqrt(V.T.plus)
      z.minus <- (T - E.T.minus)/sqrt(V.T.minus)
      p.val.up <- 1 - pnorm(z.plus)
      p.val.low <- 1 - pnorm(z.minus)
      pvals[i,1] <- round(p.val.low, digits=4)
      pvals[i,2] <- round(p.val.up, digits=4)
    }
  } 
  else for(i in 1:m){
    p.plus <- gamma[i]/(1 + gamma[i])
    p.minus <- 1/(1+gamma[i])
    E.T.plus <- sum(ranks*p.plus)
    V.T <- sum(ranks^2 * p.plus*(1-p.plus))
    E.T.minus <- sum(ranks*p.minus)
	          #Normal Approximation
    z.plus <- (T - E.T.plus)/sqrt(V.T)
    z.minus <- (T - E.T.minus)/sqrt(V.T)
    p.val.up <- 1 - pnorm(z.plus)
    p.val.low <- 1 - pnorm(z.minus)
    pvals[i,1] <- round(p.val.low, digits=4)
    pvals[i,2] <- round(p.val.up, digits=4)
  }
  out <- cbind(gamma, pvals)
  colnames(out) <- c("Gamma", "L. Bound P-Value", "U. Bound P-Value")
  cat("Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value \n")
  print(out)
  cat("\n")
  cat("Note: Gamma is Log Odds of Differential Assignment To Treatment Due to Unobserved Factors \n")
}
	
		
