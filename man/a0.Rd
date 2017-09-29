\name{a0}
\alias{a0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parameter \eqn{a_0} of the polynomial
}
\description{
analytical solution for the parameter \eqn{a_0} of the polynomial
}
\usage{
a0(t, qs, q1, q2, q3)
}
\arguments{
  \item{qs}{
last disaggregated discharge of the day before
}
  \item{q1}{
daily mean discharge of the current day
}
  \item{q2}{
daily mean discharge of the next day
}
  \item{q3}{
daily mean discharge of the day after the next day
}
}

\value{
\eqn{Q_i(t)}
}

\author{
Svenja Fischer
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (t, qs, q1, q2, q3)
{
    a0t <- 1/72 * ((48 * t^3 + 144 * t^2 + 84 * t - 12) * qs +
        t * (8 * q3 - 46 * q2 - 46 * q1) + q3 + t^2 * (-12 *
        q3 + 96 * q2 - 228 * q1) + t^3 * (-16 * q3 + 56 * q2 -
        88 * q1) - 8 * q2 + 91 * q1)
    return(a0t)
  }
}
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
