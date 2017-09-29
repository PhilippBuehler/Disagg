\name{optFun}
\alias{optFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to be optimized in the disaggregation via Weibull-function
}
\description{
calculates difference of volume and peak for given observation and approximated Weibull-function
}
\usage{
optFun(par, qStart, TM, qEnd, HQ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
parameters that are evaluated
}
  \item{qStart}{
disaggregated discharge of the time step before
}
  \item{TM}{
daily mean discharge of current day
}
  \item{qEnd}{
daily mean discharge of the next day
}
  \item{HQ}{
peak-discharge
}
}

\value{
difference of peak and volume of observations and fitted Weibull-function
}

\author{
Svenja Fischer
}

\examples{
following
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
