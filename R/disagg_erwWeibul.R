library(cmaes)

erwWeibul <- function(t, par=c(0, 1, -1, 5, -0.2)) {
  par[1] + par[2] * par[3] * par[4] * ((t - par[5])*par[3])^(par[4]-1) * exp(-par[3]*(t-par[5])^par[4])
}

optFun <- function(par, qStart, TM, qEnd, HQ) {
  eVal <- 1e6+runif(1,-1e2,1e2)

  if(abs(par[2]) < 1e-3)
    return(eVal)

  if(prod(par[1:3]) < 0)
    return(eVal)

  p1 <- qStart - par[1] * par[2] * par[3] * ((0 - par[4])*par[2])^(par[3]-1) * exp(-par[2]*(0-par[4])^par[3])

  if(is.nan(p1) | is.na(p1))
    return(eVal)

  if(par[4] > 0)
    return(eVal)

  dayDis <- erwWeibul(0:23/24, c(p1,par))

  if(any(is.nan(dayDis)) | any(is.infinite(dayDis)) | any(is.na(dayDis)))
    return(eVal)
  if (any(dayDis < 0))
    return(eVal)
  if (dayDis[24] <= qEnd)
    return(eVal)
  if (dayDis[24] >= HQ)
    return(eVal)

  if (((dayDis[2]-p1+1e-4)/(max(dayDis)-p1)) < 0.01)
    return(eVal)

  absDiffHQ <- abs(max(dayDis)-HQ)
  absDiffTM <- abs(sum(erwWeibul(0:99/100, c(p1,par)))/100 - TM)

  return(absDiffHQ+absDiffTM)
}

lorentz <- function(t, par=runif(4,0,10)) {
  par[1] + 2 * par[2]/pi*par[3]/(4*(t-par[4])^2+par[3]^2)
}

optFunLor <- function(par, qStart, TM, qEnd, HQ) {
  eVal <- 1e6+runif(1,-1e2,1e2)

  if(abs(par[2]) < 1e-3 & abs(par[3]) < 1e-3)
    return(eVal)

  if(par[1] * par[2] < 0)
    return(eVal)

  p1 <- qStart - 2 * par[1]/pi*par[2]/(4*(0-par[3])^2+par[2]^2)

  if(is.nan(p1))
    return(eVal)

  # ensure minimum width of Lorentz function
  if(abs(par[2]) < 0.2)
    return(eVal)

  dayDis <- lorentz(0:23/24, c(p1,par))

  if(any(is.nan(dayDis)) | any(is.infinite(dayDis)))
    return(eVal)
  if (any(dayDis < 0))
    return(eVal)
  if (dayDis[24] < qEnd)
    return(eVal)
  if (dayDis[24] >= HQ)
    return(eVal)

  absDiffHQ <- abs(max(dayDis)-HQ)
  absDiffTM <- abs(sum(lorentz(0:99/100, c(p1,par)))/100 - TM)

  return(absDiffHQ+absDiffTM)
}

#############
lorentzErwWeibul <- function(qStart, TM, qEnd, HQ, qPrec, decision) {
  library(cmaes)
  stopifnot(decision %in% c("wagner", "bestfit"))

  if(any(is.na(c(qStart, TM, qEnd, HQ, qPrec))) | TM > HQ ) {
    diagRow <- matrix(c(qStart, TM, qEnd, HQ, rep(NA, 9)), nrow = 1)
    colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ",
                           "wei1", "wei2", "wei3", "wei4", "wei5",
                           "lor1", "lor2", "lor3", "lor4")
    res <- rep(NA,24)
    attr(res, "diagnose") <- diagRow
    return(res)
  }

  # Lorentzfunktion
  parFitLor <- cma_es(par=c(0,1,0.5), fn=optFunLor,
                     qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                    control = list(maxit=1e3))
  count <- 0
  while (parFitLor$value > (HQ-TM)/2 & count < 50) {
    parFitLor <- cma_es(par=c(0,1,runif(1)), fn=optFunLor,
                       qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                       control = list(maxit=1e3))
    count <- count + 1
  }

  if (count > 0 & count < 50 )
    warning("Retried Lorentz optimisation ", count, " times.")

  if(count == 50) {
    decision <- "bestfit"
    warning("Lorentz optimisation aborted for",
            paste(c("qStart", "TM", "qEnd", "HQ"), round(c(qStart, TM, qEnd, HQ),2), sep=": "),
            "with error:", parFitLor$value)
  }

  par <- parFitLor$par
  p1 <- qStart - 2 * par[1]/pi*par[2]/(4*(0-par[3])^2+par[2]^2)
  dayDisLor <- lorentz(0:23/24, c(p1,par))

  p1Lor <- p1

  ## erweiterte Weibull
  parFit <- optim(par=c(1,1,1.2,0), fn=optFun,
                  qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                  method = "SANN", control = list(maxit=1e4))
  count <- 0
  while (parFit$value > (HQ-TM)/2 & count < 50) {
    parFit <- optim(par=c(1,1,1.2,runif(1)), fn=optFun,
                    qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                    method = "SANN", control = list(maxit=1e4))
    count <- count + 1
  }

  if (count > 0 & count < 50)
    warning("Retried Weibull optimisation ", count, " times.")

  if(count == 50) {
    decision <- "bestfit"
    warning("Weibull optimisation aborted for",
            paste(c("qStart", "TM", "qEnd", "HQ"), round(c(qStart, TM, qEnd, HQ),2), sep=": "),
            "with error:", parFit$value)
  }

  par <- parFit$par
  p1 <- qStart - par[1] * par[2] * par[3] * ((0 - par[4])*par[2])^(par[3]-1) * exp(-par[2]*(0-par[4])^par[3])
  dayDisWei <- erwWeibul(0:23/24, c(p1,par))

  diagRow <- matrix(c(qStart, TM, qEnd, HQ, p1, parFit$par, p1Lor, parFitLor$par), nrow = 1)
  colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ",
                         "wei1", "wei2", "wei3", "wei4", "wei5",
                         "lor1", "lor2", "lor3", "lor4")

  attr(dayDisWei, "diagnose") <- diagRow
  attr(dayDisLor, "diagnose") <- diagRow

  if (decision == "wagner") {
    A <- TM > qPrec | TM > qEnd # wachsend - fallend
    B <- (TM > dayDisLor[24] & dayDisLor[24] > qEnd) | (TM < dayDisLor[24] & dayDisLor[24] < qEnd) # zwischen aktuell und Folgetag

    if (A) {
      if (B)
        return(dayDisLor)
      else
        return(dayDisWei)
    }

    C <- qPrec <= TM & TM <= qEnd # monoton steigend

    if (C)
      return(dayDisWei)

    D <- qPrec >= TM & TM >= qEnd # monoton fallend

    if (D) {
      if (B)
        return(dayDisLor)
      else
        return(dayDisWei)
    }

    E <- qPrec > TM & TM < qEnd # fallend - wachsend

    if (E)
      return(dayDisLor)
  }
  if (decision == "bestfit") {
    if (parFitLor$value < parFit$value)
      return(dayDisLor)
    else
      return(dayDisWei)
  }
}


###

a3<-function(qs,q1,q2,q3){
  a3t<- -1/9*(6*qs-2*q3+7*q2-11*q1)
  return(a3t)
}

a2<-function(t,qs,q1,q2,q3){
  a2t<-1/6*((12*t+12)*qs-q3+t*(-4*q3+14*q2-22*q1)+8*q2-19*q1)
  return(a2t)
}

a1<-function(t,qs,q1,q2,q3){
  a1t<- -1/36*((72*t^2+144*t+42)*qs+4*q3+t*(-12*q3+96*q2-228*q1)+t^2*(-24*q3+84*q2-132*q1)-23*q2-23*q1)
  return(a1t)
}

a0<-function(t,qs,q1,q2,q3){
  a0t<- 1/72*((48*t^3+144*t^2+84*t-12)*qs+t*(8*q3-46*q2-46*q1)+q3+t^2*(-12*q3+96*q2-228*q1)+t^3*(-16*q3+56*q2-88*q1)-8*q2+91*q1)
  return(a0t)
}

disaggpoly<-function(t,a0,a1,a2,a3){
  return(a3*t^3+a2*t^2+a1*t+a0)
}

disagg<-function(q, scheitelTag=c(), decision="wagner", diagnose=FALSE, param=FALSE){
  library(cmaes)

  q[,1]<-as.Date(q[,1], "%d.%m.%Y")
  if(length(scheitelTag) > 0)
    scheitelTag[,1]<-as.Date(scheitelTag[,1], "%d.%m.%Y")

  n <- length(q[,1])
  outv<-numeric()

  ###Startwerte berechnen mit t=1###
  a30<-a3(q[1,2],q[1,2],q[2,2],q[3,2])
  a20<-a2(0.5,q[1,2],q[1,2],q[2,2],q[3,2])
  a10<-a1(0.5,q[1,2],q[1,2],q[2,2],q[3,2])
  a00<-a0(0.5,q[1,2],q[1,2],q[2,2],q[3,2])

  polyParam <- NULL

  if(param)
    polyParam <- rbind(polyParam, c(a00, a10, a20, a30))

  outv <- disaggpoly(1:24/24,a00,a10,a20,a30)

  #für jeden Ursprungszeitschritt neue Koeffizienten

  pb <- txtProgressBar(0, n-3, style = 3)

  diagList <- NULL

  for (j in 1:(n-3)) {
    setTxtProgressBar(pb, j)

    qS <- outv[24*j]
    if(is.na(qS))
      qS <- q[j+1,2]

    if (q[j+1,1] %in% scheitelTag[,1]) {
      sInd <- which(scheitelTag[,1] == q[j+1,1])
      resDisAg <- lorentzErwWeibul(qStart=qS, TM=q[j+1,2],
                                   qEnd=q[j+2,2], HQ=scheitelTag[sInd,2],
                                   qPrec=q[j,2], decision)
      if(param)
        polyParam <- rbind(polyParam, c(NA, NA, NA, NA))

      outv <- c(outv, resDisAg)
      if(diagnose)
        diagList <- rbind(diagList, attr(resDisAg, "diagnose"))
    } else {
      a3j<-a3(qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a2j<-a2(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a1j<-a1(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a0j<-a0(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])

      if(param)
        polyParam <- rbind(polyParam, c(a0j, a1j, a2j, a3j))

      ##nutze Startwerte als Startwertbedingung##
      outv <- c(outv, disaggpoly(1:24/24,a0j,a1j,a2j,a3j))
    }
  }

  close(pb)

  thour <- seq(from=ISOdate(format(q[1,1], "%Y"),format(q[1,1], "%m"),format(q[1,1], "%d"), hour=0), length.out=length(outv), by="hour")

  res <- data.frame(thour=thour,
                    dissAgg = outv)

  attr(res, "diagnose") <- diagList
  attr(res, "polyParam") <- polyParam

  return(res)
}
