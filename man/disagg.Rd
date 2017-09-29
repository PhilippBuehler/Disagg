\name{disagg}
\alias{disagg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disaggregation of daily mean discharges into hourly discharges with respect to the peak
}
\description{
In case of no peak data a third degree polynomial is fitted to the mean discharges. If the peak is known, a Lorentz- or Weibull-function is fitted. In both cases, the volume balance and the continuity have to be respected. In the case of a peak also the maximum of the diaggregated discharges has to be equal to the peak value. The procedure is based on Wagner (2012).
}
\usage{
disagg(q, scheitelTag = c(), decision = "wagner", diagnose = FALSE, param = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
A dataframe with the first column being the daily date (as character) in the format dd.mm.YY and the second column being the daily discharges (numeric).
}
  \item{scheitelTag}{
A dataframe with the first column being the peak date (as character) in the format dd.mm.YY and the second column being the peak discharges (numeric). Can be of any time-resolution (daily,monthly, yearly,...). If not specified, only daily discharges will be disaggregated without peak-knowledge.
}
  \item{decision}{
specifies decision for Lorentz- or Weibull-function, either "wagner" or "bestfit". Default value: "wagner". If "wagner" is chosen, the decision is made according to Wagner (2012), unless any of the fits aborted. In this case, the decision is changed to "bestfit", where the best optimisation (smaller deviation in peak and volume) is chosen, and a warning is given.
}
  \item{diagnose}{
Logical: shall the results of the optimisation steps be given. Default: FALSE
}
  \item{param}{
Logical: shall the parameters of the polynomial for every time intervall be given. Default: FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a dataframe with the following entries:
\item{thour}{the new time steps (in hours)}
\item{dissAgg}{the disaggregated hourly discharge values}
\item{diagnose}{the optimisation steps, if diagnose=TRUE}
\item{param}{the polynomial parameters for every time step, if param=TRUE}
}
\references{
Wagner, M. (2012): Regionalisierung von Hochwasserscheiteln auf Basis einer gekoppelten Niederschlag-Abfluss-Statistik mit besonderer Beachtung von Extremereignissen. Dissertation. Technische Universität Dresden.
}
\author{
Svenja Fischer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
following...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{manip}
