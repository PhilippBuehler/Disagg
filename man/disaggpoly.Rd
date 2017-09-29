\name{disaggpoly}
\alias{disaggpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial for days without peak
}
\description{
Polynomial for days without peak
}
\usage{
disaggpoly(t, a0, a1, a2, a3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
time parameter
}
  \item{a0}{
parameter
}
  \item{a1}{
parameter
}
  \item{a2}{
parameter
}
  \item{a3}{
parameter
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
function (t, a0, a1, a2, a3)
{
    return(a3 * t^3 + a2 * t^2 + a1 * t + a0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")

