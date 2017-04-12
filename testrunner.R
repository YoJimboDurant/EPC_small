# load needed libraries:
library(shiny)
library(lattice)
library(plyr)
library(EnvStats)
library(reshape2)

source("./utilityFUnctions.R")



set.seed(1)

input <- list(distribution="Log-Normal", sdlog=1, meanlog=0, censQ1=0.25, 
              censQ2=0.0, censQ3=0, obs1=10,obs2=0,obs3=0, simNum=100)

theData <- createData(input)

dfx <- theData$dfx
pop.mean <- theData$pop.mean
pop.var <- theData$pop.var


dfx <- theData$dfx

#ID too highly censored data:
dfx <- ddply(dfx, .(sampleNum), transform, nUncen = sum(!cen))
valid.samples <- unique(dfx$sampleNum)

pop.mean <- theData$pop.mean
pop.var <- theData$pop.var


dfx <- ddply(dfx, .(sampleNum), transform, nUncen = sum(!cen))
dfx <- ddply(dfx, .(sampleNum), function(x){
  x$nCenLevels <- length(unique(x$obs[x$cen]))
  x
})

fix1Cen <- function(dfx){
  if(dfx$nCenLevels[1]==1){
    minCen <- min(dfx$obs[dfx$cen])
    dfx$cen <- dfx$obs < minCen| dfx$cen
    dfx$obs[dfx$cen] <- minCen
  }
  
  return(dfx)
}

dfx <- ddply(dfx, .(sampleNum), fix1Cen)

lambdas <- dlply(dfx, .(sampleNum), function(dfx){
  if(dfx$nUncen[1] <3 ) { 
    results <- NA
  }else{
    if (any(dfx$cen)){
      results <- boxcoxCensored(dfx$obs, dfx$cen, lambda=seq(0,1,0.1))
    } else {
      results <- boxcox(dfx$obs, lambda=seq(0,1,0.1))
    }
  }
  
  return(results)
}
)


rec.dist <- function(lx) {
  lambdas = lx
  if(!is.na(lx[1])){
    
    best.lambda <- lambdas$lambda[lambdas$objective==max(lambdas$objective)]
    if(best.lambda <0.15) rec.dist <- "Log Normal"
    if(0.15 < best.lambda & best.lambda < 0.6) rec.dist <- "Gamma"
    if(best.lambda > 0.6) rec.dist <- "Normal"
  }else{rec.dist = "<3 uncensored"}
  return(rec.dist)
}



rec.dist.fit <- ldply(lambdas, function(lx) rec.dist(lx))



# normal distribution
norm.envStat <- dlply(dfx, .(sampleNum), function(dfx) {
  if(dfx$nUncen[1] <3){
    return(NA)
    
    
  }else{
    if(any(dfx$cen)){
      enormCensored(dfx$obs, dfx$cen, ci=TRUE, ci.type="upper", 
                    ci.method="normal.approx")
    }else{
      enorm(dfx$obs, ci=TRUE, ci.type="upper", 
            ci.method="exact") 
    }
  }
}
)

norm.envStat.UCL <- ldply(norm.envStat, function(lx){
  if(!is.na(lx[[1]])){
    unname(lx$interval$limits[2])
  }else{
    return(NA)
  }
})

ldply(norm.envStat, function(lx){
  if(!is.na(lx[[1]])){
    qnorm(pnorm(-3), mean = lx$parameters[[1]], sd   = lx$parameters[[2]])
  }else{
    return(NA)
  }
})

}

norm.envStat.UCL.coverage <- mean(norm.envStat.UCL$V1 >=pop.mean) 


# gamma distribution
gamma.envStat <- dlply(dfx, .(sampleNum), function(dfx) {
  if(dfx$nUncen[1] <3){
    return(NA)
    
    
  }else{
    if(any(dfx$cen)){
      egammaAltCensored(dfx$obs, dfx$cen, ci=TRUE, ci.type="upper", 
                        ci.method="normal.approx")
    }else{
      egammaAlt(dfx$obs, ci=TRUE, ci.type="upper", 
                ci.method="normal.approx") 
    }
  }
}

)

gamma.envStat.UCL <- ldply(gamma.envStat, function(lx){
  if(!is.na(lx[[1]])){
    unname(lx$interval$limits[2])
  }else{
    return(NA)
  }
})

# lognormal distribution
ln.envStat <- dlply(dfx, .(sampleNum), function(dfx) {
  if(dfx$nUncen[1] <3){
    return(NA)
    
    
  }else{
    if(any(dfx$cen)){
      elnormAltCensored(dfx$obs, dfx$cen, ci=TRUE, ci.type="upper", 
                        ci.method="cox")
    }else{
      elnormAlt(dfx$obs, ci=TRUE, ci.type="upper", 
                ci.method="cox") 
    }
  }
}
)

ln.envStat.UCL <- ldply(ln.envStat, function(lx){
  if(!is.na(lx[[1]])){
    unname(lx$interval$limits[2])
  }else{
    return(NA)
  }
})



#data.frame for UCL's
rec.df <- data.frame(sampNum = rec.dist.fit$sampleNum,
                     rec.dist=rec.dist.fit$V1, 
                     Gamma = gamma.envStat.UCL$V1,
                     LogNormal=ln.envStat.UCL$V1,
                     Normal=norm.envStat.UCL$V1
)


myChoice <- function(dfx){
  if(dfx$rec.dist=="Gamma") dfx$rec.UCL <- dfx$Gamma
  if(dfx$rec.dist=="Log Normal") dfx$rec.UCL <- dfx$LogNormal
  if(dfx$rec.dist=="Normal") dfx$rec.UCL <- dfx$Normal
  dfx
}

rec.df <- ddply(rec.df, .(sampNum), function(dfx) myChoice(dfx))


norm.envStat.UCL.coverage <- mean(norm.envStat.UCL$V1 >=pop.mean, na.rm=TRUE) 
gamma.envStat.UCL.coverage <- mean(gamma.envStat.UCL$V1 >=pop.mean, na.rm=TRUE) 
ln.envStat.UCL.coverage <- mean(ln.envStat.UCL$V1 >=pop.mean, na.rm=TRUE)

rec.UCL.coverage <- mean(rec.df$rec.UCL >=pop.mean, na.rm=TRUE)

envStat.UCL.coverage = c(gamma.envStat.UCL.coverage,
                         ln.envStat.UCL.coverage,
                         norm.envStat.UCL.coverage,
                         rec.UCL.coverage)
names(envStat.UCL.coverage) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL")
print(envStat.UCL.coverage)
results = list(
  data = dfx,
  dist = lambdas,
  rec.dist.fit = rec.dist.fit,
  envStat.UCL.coverage = envStat.UCL.coverage)

xyplot(value ~ sampNum |variable,
  data = melt(rec.df, id.vars = c("sampNum", "rec.dist")),
  layout=c(1,4),
  as.table=TRUE)



# Bayes -------------------------------------------------------------------
Bayes_lognormal <- ddply(dfx, .(sampleNum), function(dfx){
  if(dfx$nUncen[1] <3){
    return(NA)
    
    
  }else{
  EH_Cen_Bayes(dfx$obs, dfx$cen, conf=c(0.05,0.95))[[2]]
}
}

)


