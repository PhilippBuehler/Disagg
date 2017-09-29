\name{a1}
\alias{a1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parameter \eqn{a_1} of the polynomial
}
\description{
analytical solution for the parameter \eqn{a_1} of the polynomial
}
\usage{
a1(t, qs, q1, q2, q3)
}
%- maybe also 'usage' for other objects documented here.
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
    a1t <- -1/36 * ((72 * t^2 + 144 * t + 42) * qs + 4 * q3 +
        t * (-12 * q3 + 96 * q2 - 228 * q1) + t^2 * (-24 * q3 +
        84 * q2 - 132 * q1) - 23 * q2 - 23 * q1)
    return(a1t)
  }
}
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
