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
 which encourages dose \code{D}, this function performs a Rosenbaum's bounds sensitivity analysis. Note   that matching is done on levels of the instrument. See example below.}

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

\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}
  
  Jason W. Morgan, Ohio State University, \email{morgan.746@osu.edu}
}

\seealso{See also \code{\link{data.prep}}, \code{\link{binarysens}},  \code{\link{hlsens}}, \code{\link[Matching]{Match}}, \code{\link{mcontrol}}}

\examples{
## Example from Rosenbaum (2010, ch. 5).

data(AngristLavy)

#Match on Economic Status Across Levels of the Instrument
rr <- Match(Y=AngristLavy$avgmath, Tr=AngristLavy$z, X=AngristLavy$pct_disadv, 
estimand ="ATC", M=2,  replace=FALSE)

#Extract Matched Outome Data
ctrl <- AngristLavy$avgmath[rr$index.control]
trt <- AngristLavy$avgmath[rr$index.treated]

#Extract Matched Doses
#Doses Encouraged By Instrument - Here Class Size
csize.trt <- AngristLavy$classize[rr$index.treated]
csize.ctrl <- AngristLavy$classize[rr$index.control]

#Run Sensitivity Analsyis
iv_sens(trt, ctrl, csize.trt, csize.ctrl, Gamma=2, GammaInc=.1)

}
