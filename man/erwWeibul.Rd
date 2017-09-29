\name{erwWeibul}
\alias{erwWeibul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weibull-function
}
\description{
The Weibull-function given by
\eqn{{{Q}_{i}}(t)={{\beta }_{{{1}_{i}}}}+{{\beta }_{{{2}_{i}}}}{{\beta }_{{{3}_{i}}}}{{\beta }_{{{4}_{i}}}}{{((t-{{\beta }_{{{5}_{i}}}})\beta_{3_i})}^{{{\beta }_{{{4}_{i}}}}-1}}\exp (-{{\beta }_{{{3}_{i}}}}{{(t+{{\beta }_{{{5}_{i}}}})}^{{{\beta }_{{{4}_{i}}}}}})}
with parameters \eqn{\beta_{j_i}}
}
\usage{
erwWeibul(t, par = c(0, 1, -1, 5, -0.2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
function argument, here: time step
}
  \item{par}{
parameters \eqn{\beta_{j_i}}, with starting values
}
}

\value{
returns \eqn{Q_i(t)}
}

\author{
Svenja Fischer
}

\examples{

## The function is currently defined as
function (t, par = c(0, 1, -1, 5, -0.2))
{
    par[1] + par[2] * par[3] * par[4] * ((t - par[5]) * par[3])^(par[4] -
        1) * exp(-par[3] * (t - par[5])^par[4])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
