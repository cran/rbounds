binarysens <- function(x, y=NULL, Gamma=6, GammaInc=1)
	{
	 if (length(x)==1){
  		ctrl <- x
  		trt <- y
  		}
  	else
  	{
		ctrl <- sum(Y[x$index.control])
        trt <- sum(Y[x$index.treated])
    }

  gamma <- seq(1, Gamma, by=GammaInc)
  mx <- ctrl + trt
  up <- c()
  lo <- c()
  series <- seq(trt,mx, by=1)
  n.it <- length(gamma)
  for(i in 1:n.it)
    {
      p.plus <- gamma[i]/(1 + gamma[i])
      p.minus <- 1/(1 + gamma[i])
      up.tmp <- sum(dbinom(series, mx, p=p.plus))
      lo.tmp <- sum(dbinom(series, mx, p=p.minus))
      up <- c(up, up.tmp)
      lo <- c(lo, lo.tmp)
    }
  
  out <- cbind(gamma, round(lo, 5), round(up, 5))
  colnames(out) <- c("Gamma", "Lower Bound P-Value", "Upper Bound P-Value")
  cat("Rosenbaum Sensitivity Test\n")
  print(out, scientific=FALSE)
  cat("\n")
  cat("Note: Gamma is Log Odds of Differential Assignment To Treatment Due To Unobserved Factors \n")
} #end of binary.sens function

