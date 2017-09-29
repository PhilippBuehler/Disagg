\name{lorentz}
\alias{lorentz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lorentz-function
}
\description{
The Lorentz-function given by
\eqn{Q_i(t)=\alpha_{1_i}+2\frac{\alpha_{2_i}}{\pi}\frac{\alpha_{3_i}}{4(t-\alpha_{4_i})^2+\alpha_{3_i}^2}}

with parameters \eqn{\alpha_{j_i}}
}
\usage{
lorentz(t, par = runif(4, 0, 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
function argument, here: time step
}
  \item{par}{
parameters \eqn{\alpha_{j_i}}, with starting values
}
}
\value{
returns \eqn{Q_i(t)}
}

\author{
Svenja Fischer
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (t, par = runif(4, 0, 10))
{
    par[1] + 2 * par[2]/pi * par[3]/(4 * (t - par[4])^2 + par[3]^2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
