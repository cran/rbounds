\name{iv_sens}
\alias{iv_sens}

\title{
Function to calculate Rosenbaum bounds for IV Estimator based on Wilcoxon sign rank test. 
}

\description{\code{iv_sens} performs a non-parametric, instrumental
  variable sensitivity analysis on matched pairs following the logic of
  the Neyman-Rubin framework for causal inference. The function supports
  a variable-valued instrument.}

\usage{iv_sens(Rt, Rc, Dt, Dc, Gamma = 6, GammaInc = 1)}

\arguments{  
  \item{Rt,Rc}{Vectors of observed response outcomes for matched
    treatment and control observations, respectively.}

  \item{Dt,Dc}{Vectors of observed doses for matched observations, respectively.
  This is level of dose encouraged by the instrument.}

\item{Gamma}{Upper-bound on gamma parameter.} 

\item{GammaInc}{To set user specified increments for gamma parameter.} 
}

\details{Given matched pairs of observations on an instrument \code{Z}, 
 which encourages dose \code{D}, this function performs a Rosenbaum's bounds sensitivity analysis. Note  that matching is done on levels of the instrument.}

\value{Returns an object of class \code{rbounds}.}

\references{
  Angrist, Joshua D., Imbens, Guido W., and Rubin, Donald B. (1996).
  "Identification of Causal Effects Using Instrumental Variables."
  \emph{Journal of the American Statistical Association} 91/434, pp.
  444--455.

  Rosenbaum, Paul R. (1996). "Comment." \emph{Journal of the American 
  Statistical Association} 91/434, pp. 465--468.

  Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
  Springer-Verlag.

  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
}

\author{Luke Keele, University of Pennsylvania, \email{luke.keele@gmail.com}}
  
\seealso{See also \code{\link{binarysens}},  \code{\link{hlsens}}, \code{\link{mcontrol}}}

\examples{

data(angristlavy)

# Example from Ch 5 of Design of Observational Studies
	
#Extract Matched Outome Data
Rt <- angristlavy$avgmath[angristlavy$z==1]
Rc <- angristlavy$avgmath[angristlavy$z==0]

#Extract Matched Doses
#Doses Encouraged By Instrument - Here Class Size
Dt <- angristlavy$clasz[angristlavy$z==1]
Dc <- angristlavy$clasz[angristlavy$z==0]


#Run Sensitivity Analsyis - Rank Sum Test
iv_sens(Rc, Rt, Dc, Dt, Gamma=1.5, GammaInc=.01)


}
