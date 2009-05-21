
hlsens <- function(x, y=NULL, pr=.1, Gamma=6, GammaInc=1)
{
  
  if (is.numeric(x)){
  trt <- x
  ctrl <- y
  }
  else {
        ctrl <- x$mdata$Y[x$mdata$Tr==0]
        trt <- x$mdata$Y[x$mdata$Tr==1]
      }
      
        gamma <- seq(1, Gamma, by=GammaInc)
        k <- length(gamma)
  
  ttau <- function(x) {
    tau <- x
    adj.trt <- trt - tau
    diff.2 <- adj.trt - ctrl
    ranks <- rank(abs(diff.2), ties.method="average")
    psi <- as.numeric(diff.2 > 0)
    sum(psi * ranks)
  }
  
  tau.up <- tau.l <- wilcox.test(trt, ctrl, paired=TRUE, conf.int=TRUE, exact=FALSE)$estimate
  #base <- c(1,tau.up,tau.l)
  eps <- 1.e-8 
  c.int <- matrix(0,k,2)
  s <- length(trt)
  
  for(i in 1:k){
    p.minus = 1/(1+gamma[i])
    p.plus = gamma[i]/(gamma[i]+1)
    t.min <- p.minus*(s*(s+1)/2)
    t.max <- p.plus*(s*(s+1)/2)
    lb <- t.min
    ub <- t.max

    while(abs(ub - lb) > eps){
      if (lb < ub){
        tau.old <- tau.up
        tau.up <-  tau.old  + pr
        ub <- ttau(tau.up)
      }
      else break
    }
    c.int[i,2] <- tau.up

    ub <- t.max
#    lb <- ttau(tau.l)
    lb <- t.min

    while (abs(ub - lb) > eps){
      if (lb <= ub){
        tau.old <- tau.l
        tau.l <-  tau.old  - pr
        lb <- ttau(tau.l)
      }
      else break
    }
    c.int[i,1] <- tau.l
  }
  
  out <- cbind(gamma, signif(c.int, digits=5))
  #out <- rbind(base, out)
  colnames(out) <- c("Gamma", "L. Bound HL Est.", "U. Bound HL Est.")
  cat("Rosenbaum Sensitivity Test for Hodges-Lehmann Point Estimate \n")
  print(out, scientific=FALSE)
  cat("\n")
  cat("Note: Gamma is Log Odds of Differential Assignment To Treatment Due To Unobserved Factors \n")
  
}

