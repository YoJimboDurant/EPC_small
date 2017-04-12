getDetLim <- function(input, groupno="censQ1"){

  if(input$distribution =="Log-Normal"){
    dist="qlnorm"
    qargs = list(meanlog=input$meanlog, sdlog=input$sdlog)
  }else{
    if(input$distribution=="Gamma"){
      dist="qgamma"
      qargs = list(shape=input$shape, scale=input$scale)
    }else{
      if(input$distribution=="Normal"){
        dist="qnorm"
        qargs = list(input$mu, sd =input$sd)

      }
    }
  }
  input1 <- c(input$censQ1, input$censQ2, input$censQ3)
  qargs$p <- input1[c("censQ1", "censQ2", "censQ3") %in% groupno]
  x <- do.call(dist, qargs)
  return(x)
}


createData <- function(input,...) {
  if(input$distribution =="Log-Normal"){
    dist="lnorm"
    qargs = list(meanlog=input$meanlog, sdlog=input$sdlog)
    pop.mean <- exp(input$meanlog + 0.5 *input$sdlog ^2)
    pop.var <- (exp(input$sdlog^2)-1) * exp(2 *input$meanlog +input$sdlog ^2)

  }else{
    if(input$distribution=="Gamma"){
      dist="gamma"
      qargs = list(shape=input$shape, scale=input$scale)
      pop.mean <- input$shape * input$scale
      pop.var <- input$shape * input$scale^2

    }else{
      if(input$distribution=="Normal"){
        dist="norm"
        qargs = list(input$mu, sd =input$sd)
        pop.mean <- input$mu
        pop.var <- input$sd^2


      }
    }
  }

  totalObs <- input$obs1 + input$obs2 +input$obs3
  rdist <- paste("r", dist, sep="")
  qdist <- paste("q", dist, sep="")

  qargs$p <- c(input$censQ1, input$censQ2, input$censQ3)
  cenValues <- do.call(qdist, qargs)
  cenValues <- rep(cenValues, times = c(input$obs1, input$obs2, input$obs3))

  cenValue <- sample(cenValues, totalObs, replace=FALSE)
  cenValue <- rep(cenValue, times = input$simNum)

  # create sample from distribution
  qargs$p <- NULL
  qargs$n <- totalObs * input$simNum
  dist <- do.call(rdist, qargs)
  cen = cenValue > dist


  dfx <- data.frame(value=dist, cen=cenValue>dist, cenValue=cenValue)
  dfx$obs <- dfx$value
  dfx$obs[dfx$cen] <- cenValue[dfx$cen]
  dfx$obs_zero <- dfx$value
  dfx$obs_zero[dfx$cen] <- 0
  dfx$obs_half <-dfx$value
  dfx$obs_half[dfx$cen] <- 0.5 * cenValue[dfx$cen]
  dfx$obs_full <- dfx$value
  dfx$obs_full[dfx$cen] <- cenValue[dfx$cen]
  dfx$obs_sqrt2 <- dfx$value
  dfx$obs_sqrt2[dfx$cen] <- cenValue[dfx$cen]/sqrt(2)
  dfx$run <- input$myValue
  dfx$sampleNum <- rep(1:input$simNum, each=totalObs)
  return(list(dfx = dfx, pop.mean = pop.mean, pop.var = pop.var))
}

createInput <- function(distribution, censQ1, censQ2, censQ3, obs1,
                        obs2, obs3, simNum=100,...){

  input <-list(distribution=distribution, censQ1=censQ1, censQ2=censQ2,
               censQ3=censQ3, obs1=obs1, obs2=obs2, obs3=obs3, simNum=simNum,...)

  dataObj <- createData(input)
  dfx <- dataObj$dfx
  dfx$pop.mean <- dataObj$pop.mean
  dfx$pop.var <- dataObj$pop.var

  return(dfx)

}

#create test dataframe to generate data



mean.boot <- function(x,ind)
{
  return(c(mean(x[ind]),var(x[ind])/length(ind)))
}

bootRosMean <- function(dfx,i){
  dfx <- dfx[i,]
  c(mean(dfx$modeled), (var(dfx$modeled)/length(dfx$modeled)))

}

bootMean <- function(x, i){
  x <- x[i]
  c(mean(x), (var(x)/length(x)))

}

rosSimpleBoot <- function(myros, conf = 0.95, R=1999, M=50){
  require(NADA)
  require(boot)
  mydfx <- data.frame(myros)
  myboot <- boot(mydfx, bootRosMean, R=R)
  bootci <- boot.ci(myboot)
  bootCI <- data.frame(matrix(c(bootci$norm[2:3], bootci$basic[4:5],
                                bootci$perc[4:5], bootci$bca[4:5],
                                bootci$student[4:5]), nrow=1))
  names(bootCI) <- c("ros.normal.LCL", "ros.normal.UCL", "ros.basic.LCL", "ros.basic.UCL",
                     "ros.percentile.LCL", "ros.percentile.UCL", "ros.BCa.LCL", "ros.BCa.UCL",
                     "ros.studentized.LCL", "ros.studentized.UCL")
  return(bootCI)
}



lmbootRos_MI <- function(rosObj, R.model=10, R.boot=4999, p=c(0.025,0.975)){

  require(NADA)
  require(simpleboot)
  stopifnot(c("ros", "lm") %in% class(rosObj))
  stopifnot(is.numeric(R.model))
  stopifnot(is.numeric(R.boot))
  stopifnot(is.numeric(p))
  stopifnot(length(p)==2)





  # identify censored values
  cenValues <- rosObj$modeled[rosObj$censored]
  ppValues <- rosObj$pp[rosObj$censored]
  sampleSize <- length(rosObj$modeled)

  lboot <- lm.boot(rosObj, R=R.model, new.xpts=data.frame(pp.nq=qnorm(ppValues)))

  # generate imputed values
  impDF <- ldply(lboot$boot.list, function(lx) data.frame(lx$fitted))

  # downweight the observations

  impDF$weight <- length(ppValues)/length(impDF$lx.fitted)
  impDF$lx.fitted <- exp(impDF$lx.fitted)

  # extract the observed data
  obs <- rosObj$modeled[!rosObj$censored]

  #upweight the observed data
  weight <- rep(1, length(obs))

  sampleData <- c(obs, impDF$lx.fitted)
  sampleWeight <- c(weight, impDF$weight) # weight should sum to original sample size

  myboot <- boot(sampleData, bootMean, R= R.boot, w=sampleWeight)
  bootci <- boot.ci(myboot, t0=mean(rosObj),
                    t=myboot$t[,1], var.t0=sd(rosObj)^2, var.t=myboot$t[,2])
  bootCI <- data.frame(matrix(c(bootci$normal[2:3], bootci$basic[4:5],
                                bootci$perc[4:5], bootci$bca[4:5],
                                bootci$student[4:5]), nrow=1))
  names(bootCI) <- c("ros.normal.LCL", "ros.normal.UCL", "ros.basic.LCL", "ros.basic.UCL",
                     "ros.percentile.LCL", "ros.percentile.UCL", "ros.BCa.LCL", "ros.BCa.UCL",
                     "ros.studentized.LCL", "ros.studentized.UCL")

  return(bootCI)

}


# use Land's exact interval with ROS
landRos <- function(rosObj){
  stopifnot(c("ros", "lm") %in% class(rosObj))
  rosDf <- data.frame(rosObj)
  mean(cenmle(rosDf$modeled, rep(FALSE, length(rosDf$modeled))))
}

## biasRMSE calculates Root Mean Square Error and bias for a mean and standard deviation
biasRMSE <- function(dfx, pop.mean=tesit$pop.mean, pop.var=tesit$pop.var, variable="sample"){

    if(any(grepl("[.]sd", names(dfx)))){
        variable.mean <- paste(variable,"mean", sep=".")
        variable.sd <- paste(variable,"sd", sep=".")
        variable.mean.rmse =rmse(pop.mean, dfx[variable.mean])
        variable.sd.rmse =rmse(sqrt(pop.var), dfx[variable.sd])
        variable.mean.bias = sum((dfx[variable.mean]-pop.mean)/pop.mean)/dim(dfx)[1]
        variable.sd.bias = sum((dfx[variable.sd]-sqrt(pop.var))/sqrt(pop.var))/dim(dfx)[1]
        rmseBias <- c(variable.mean.rmse, variable.sd.rmse, variable.mean.bias, variable.sd.bias)
        names(rmseBias) <- paste0(variable, c(".mean.rmse", ".sd.rmse",".mean.bias",".sd.bias"))
    } else {
        variable.mean <- paste(variable,"mean", sep=".")
        variable.mean.rmse =rmse(pop.mean, dfx[variable.mean])
        variable.mean.bias = sum((dfx[variable.mean]-pop.mean)/pop.mean)/dim(dfx)[1]
        rmseBias <- c(variable.mean.rmse, variable.mean.bias)
        names(rmseBias) <- paste0(variable, c(".mean.rmse", ".mean.bias"))

    }
    return(rmseBias)
}



smallKM = function(obs,cen){
# this function takes NADA format data and returns KM mean
#   method is based off of D. Helsel 2012 Book. I wrote it for bootstrapping
#   Function does not require the overhead of the cenfit and runs approximately 10
#   times faster.

#   Args:
#     obs: observed value or censoring limit
#     cen: logical value if obs values is censored or not
#
#   Returns:
#     mean value using KM statistics
#
  num.risk <- cumsum(table(obs))
  num.det <- table(obs * !cen)
  num.det <-num.det[rownames(num.det)!=0]
  det.values <-as.numeric(rownames(num.det))
  min.stime=min(0,min(det.values))
  num.risk <- num.risk[rownames(num.det)]
  S <- cumprod(((num.risk -num.det)/num.risk)[length(num.risk):1])

  sum(c(diff(c(min.stime, det.values[length(det.values):1])),0) *
    c(1, S))
}

# bootstrapping functions ----

kmmean_2 = function(x,d) {
  x <- x[d,]
  smallKM(x$obs, x$censored)
}



mymeans2 <-function(mydata, R=2000){
  boot(mydata, kmmean_2, R, weights=mydata$weight)
}



bci <- function(myboot) {
  unname(boot.ci(myboot))
}


kmmean <-  function(dfx, d) {
    dfx <-dfx[d,]
    smallKM(dfx$obs, dfx$cen)
}

smallKMSimpleBoot <- function(mydfx, R=1999, ...){

  myboot <- boot(mydfx, kmmean, R=R)
  bootci <- boot.ci(myboot, type=c("norm","basic", "perc", "bca"))
  bootCI <- data.frame(matrix(c(bootci$norm[2:3], bootci$basic[4:5],
                                bootci$perc[4:5], bootci$bca[4:5]
                                ), nrow=1))
  names(bootCI) <- c("smallkm.normal.LCL", "smallkm.normal.UCL", "smallkm.basic.LCL", "smallkm.basic.UCL",
                     "smallkm.percentile.LCL", "smallkm.percentile.UCL", "smallkm.BCa.LCL", "smallkm.BCa.UCL"
                     )
  return(bootCI)
}

#from  Hadley's github
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name)=="guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


gatherResults <- function(path.glob, header=TRUE) {
    return(
        rbind(
            read.table(
                normalizePath(Sys.glob(path.glob)),
                header)))
}
