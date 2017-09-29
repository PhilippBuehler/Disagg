\name{a3}
\alias{a3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parameter \eqn{a_3} of the polynomial
}
\description{
analytical solution for the parameter \eqn{a_3} of the polynomial
}
\usage{
a3(qs, q1, q2, q3)
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
function (qs, q1, q2, q3)
{
    a3t <- -1/9 * (6 * qs - 2 * q3 + 7 * q2 - 11 * q1)
    return(a3t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")

