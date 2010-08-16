\name{iv_sens}
\alias{iv_sens}
\alias{iv_gssearch}
\alias{iv_pval}
\alias{maim_dist}
\alias{plot.ivsens}
\alias{print.ivsens}
\alias{summary.ivsens}

\title{Non-parametric Instrumental Variable Sensitivity Analysis for
  Causal Inference}

\description{\code{iv_sens} performs a non-parametric, instrumental
  variable sensitivity analysis on matched pairs following the logic of
  the Neyman-Rubin framework for causal inference. The function supports
  a variable-valued instrument.}

\usage{iv_sens(Rt, Rc, Dt, Dc, Zt = rep(1, length(Rt)), Zc = rep(0,
       length(Rt)), b = 0, alpha = 0.025, Beta = 2, BetaInc = 0.01)}

\arguments{  
  \item{Rt,Rc}{Vectors of observed response outcomes for matched
    treatment and control observations, respectively.}

  \item{Dt,Dc}{Vectors of observed doses for matched treatment and
    control observations, respectively.}
  
  \item{Zt,Zc}{Vectors of instrumental variables for matched treatment
    and control observations, respectively. Defaults to \code{Zt = 1}
    and \code{Zc = 0}, though variable values are supported.}

  \item{b}{Default value is 0.}

  \item{alpha}{Confidence level defining the boundary of the sensitivity
    analysis. Default value is \code{alpha = 0.025}.}
  
  \item{Beta}{An integer defining the range over which to calculate the
    sensitivity. The lower and upper bounds will be defined internally
    as \code{[-abs(Beta), abs(Beta)]}.}
  
  \item{BetaInc}{A floating point value indicating the values in
    \code{[-abs(Beta), abs(Beta)]} for which a \emph{p}-value will be
    calculated. The default is 0.01. Decreasing \code{BetaInc} will
    increase the accuracy of the resulting confidence interval but will
    also increase the calulation time.}
}

\details{Given matched pairs of observations, where the instrument for
  dose \code{D}, \code{Z}, defines the treated and control units, this
  function performs a non-parametric sensitivity analysis. \code{iv_pval},
  \code{iv_gssearch}, and \code{maim_dist} are functions used internally
  by \code{iv_sens}.}

\value{Returns an object of class \code{ivsens}.}

\references{
  Angrist, Joshua D., Imbens, Guido W., and Rubin, Donald B. (1996).
  "Identification of Causal Effects Using Instrumental Variables",
  \emph{Journal of the American Statistical Association} 91/434, pp.
  444--455.

  Rosenbaum, Paul R. (1996). "Comment", \emph{Journal of the American 
  Statistical Association} 91/434, pp. 465--468.

  Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
  Springer-Verlag.

  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
}

\author{
  Luke Keele, Ohio State University, \email{keele.4@osu.edu}
  
  Jason W. Morgan, Ohio State University, \email{morgan.746@osu.edu}
}

\seealso{
  Also see \code{\link{plot.ivsens}}.
}

\examples{
## Example from Rosenbaum (2010, ch. 5).

data(AngristLavy)
Trt  <- subset(AngristLavy, z == 0) 
Ctrl <- subset(AngristLavy, z == 1)

## Set treatment and control groups.
Rt1 <- Trt$avgmath
Rt2 <- Trt$avgverb
Dt  <- Trt$clasz
Rc1 <- Ctrl$avgmath
Rc2 <- Ctrl$avgverb
Dc  <- Ctrl$clasz

## Estimates.
(Math <- iv_sens(Rt1, Rc1, Dt, Dc, Beta = 1, BetaInc = .001))
(Verb <- iv_sens(Rt2, Rc2, Dt, Dc, Beta = 1, BetaInc = .001))

## Plot.
plot(Math)
plot(Verb)
}
