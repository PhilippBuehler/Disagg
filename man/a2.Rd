\name{a2}
\alias{a2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parameter \eqn{a_2} of the polynomial
}
\description{
analytical solution for the parameter \eqn{a_2} of the polynomial
}
\usage{
a2(t, qs, q1, q2, q3)
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
    a2t <- 1/6 * ((12 * t + 12) * qs - q3 + t * (-4 * q3 + 14 *
        q2 - 22 * q1) + 8 * q2 - 19 * q1)
    return(a2t)
  }
}
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
