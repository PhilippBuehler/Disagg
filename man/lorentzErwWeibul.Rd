\name{lorentzErwWeibul}
\alias{lorentzErwWeibul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disaggregation-Decision
}
\description{
decides, wether the Weibull or the Lorentz-function should be used as approximation for the discharges.
}
\usage{
lorentzErwWeibul(qStart, TM, qEnd, HQ, qPrec, decision)
}
\arguments{
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
  \item{qPrec}{
daily mean discharge of the day before
}
  \item{decision}{
which aspect should be used to make  the decision: either "wagner" or "bestfit"
}
}
\details{
If the decision is "wagner", the function is chosen according to the decisions given in Wagner (2012). If the decision is "bestfit" the function with the smallest peak and volume difference compared to the observations is used.
}
\value{
disaggregated hourly discharges according to the decision made
}
\references{
Wagner, M. (2012): Regionalisierung von Hochwasserscheiteln auf Basis einer gekoppel-ten Niederschlag-Abfluss-Statistik mit besonderer Beachtung von Extremereignissen. Dissertation. Technische Universität Dresden.
}
\author{
Svenja Fischer
}

\examples{
following

## The function is currently defined as
function (qStart, TM, qEnd, HQ, qPrec, decision)
{
    stopifnot(decision \%in\% c("wagner", "bestfit"))
    if (any(is.na(c(qStart, TM, qEnd, HQ, qPrec))) | TM > HQ) {
        diagRow <- matrix(c(qStart, TM, qEnd, HQ, rep(NA, 9)),
            nrow = 1)
        colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ",
            "wei1", "wei2", "wei3", "wei4", "wei5", "lor1", "lor2",
            "lor3", "lor4")
        res <- rep(NA, 24)
        attr(res, "diagnose") <- diagRow
        return(res)
    }
    parFitLor <- cma_es(par = c(0, 1, 0.5), fn = optFunLor, qStart = qStart,
        TM = TM, qEnd = qEnd, HQ = HQ, control = list(maxit = 1000))
    count <- 0
    while (parFitLor$value > (HQ - TM)/2 & count < 50) {
        parFitLor <- cma_es(par = c(0, 1, runif(1)), fn = optFunLor,
            qStart = qStart, TM = TM, qEnd = qEnd, HQ = HQ, control = list(maxit = 1000))
        count <- count + 1
    }
    if (count > 0 & count < 50)
        warning("Retried Lorentz optimisation ", count, " times.")
    if (count == 50) {
        decision <- "bestfit"
        warning("Lorentz optimisation aborted for", paste(c("qStart",
            "TM", "qEnd", "HQ"), round(c(qStart, TM, qEnd, HQ),
            2), sep = ": "), "with error:", parFitLor$value)
    }
    par <- parFitLor$par
    p1 <- qStart - 2 * par[1]/pi * par[2]/(4 * (0 - par[3])^2 +
        par[2]^2)
    dayDisLor <- lorentz(0:23/24, c(p1, par))
    p1Lor <- p1
    parFit <- optim(par = c(1, 1, 1.2, 0), fn = optFun, qStart = qStart,
        TM = TM, qEnd = qEnd, HQ = HQ, method = "SANN", control = list(maxit = 10000))
    count <- 0
    while (parFit$value > (HQ - TM)/2 & count < 50) {
        parFit <- optim(par = c(1, 1, 1.2, runif(1)), fn = optFun,
            qStart = qStart, TM = TM, qEnd = qEnd, HQ = HQ, method = "SANN",
            control = list(maxit = 10000))
        count <- count + 1
    }
    if (count > 0 & count < 50)
        warning("Retried Weibull optimisation ", count, " times.")
    if (count == 50) {
        decision <- "bestfit"
        warning("Weibull optimisation aborted for", paste(c("qStart",
            "TM", "qEnd", "HQ"), round(c(qStart, TM, qEnd, HQ),
            2), sep = ": "), "with error:", parFit$value)
    }
    par <- parFit$par
    p1 <- qStart - par[1] * par[2] * par[3] * ((0 - par[4]) *
        par[2])^(par[3] - 1) * exp(-par[2] * (0 - par[4])^par[3])
    dayDisWei <- erwWeibul(0:23/24, c(p1, par))
    diagRow <- matrix(c(qStart, TM, qEnd, HQ, p1, parFit$par,
        p1Lor, parFitLor$par), nrow = 1)
    colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ", "wei1",
        "wei2", "wei3", "wei4", "wei5", "lor1", "lor2", "lor3",
        "lor4")
    attr(dayDisWei, "diagnose") <- diagRow
    attr(dayDisLor, "diagnose") <- diagRow
    if (decision == "wagner") {
        A <- TM > qPrec | TM > qEnd
        B <- (TM > dayDisLor[24] & dayDisLor[24] > qEnd) | (TM <
            dayDisLor[24] & dayDisLor[24] < qEnd)
        if (A) {
            if (B)
                return(dayDisLor)
            else return(dayDisWei)
        }
        C <- qPrec <= TM & TM <= qEnd
        if (C)
            return(dayDisWei)
        D <- qPrec >= TM & TM >= qEnd
        if (D) {
            if (B)
                return(dayDisLor)
            else return(dayDisWei)
        }
        E <- qPrec > TM & TM < qEnd
        if (E)
            return(dayDisLor)
    }
    if (decision == "bestfit") {
        if (parFitLor$value < parFit$value)
            return(dayDisLor)
        else return(dayDisWei)
    }
  }
}
\keyword{manip}
