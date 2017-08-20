
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {
  
  datasetInput <- reactive({
    isolate({thedata <- createData(input)
             dfx <- thedata$dfx
             
             #ID too highly censored data:
             dfx <- ddply(dfx, .(sampleNum), transform, nUncen = sum(!cen))
       
             valid.samples <- unique(dfx$sampleNum)
             
             pop.mean <- thedata$pop.mean
             pop.var <- thedata$pop.var
   
            
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
             
             rec.dist.alt <- function(lx) {
               lambdas = lx
               if(!is.na(lx[1])){
                 lambdas$objective <- lambdas$objective[lambdas$lambda<0.7]
                 lambdas$lambda <- lambdas$lambda[lambdas$lambda<0.7] 
                 best.lambda <- lambdas$lambda[lambdas$objective==max(lambdas$objective)]
                 if(best.lambda <0.15) rec.dist <- "Log Normal"
                 if(0.15 < best.lambda & best.lambda <= 0.7) rec.dist <- "Gamma"
                 
               }else{rec.dist = "<3 uncensored"}
               return(rec.dist)
             }
             
             
             
             rec.dist.fit <- ldply(lambdas, function(lx) rec.dist(lx))
             rec.dist.alt.fit <- ldply(lambdas, function(lx) rec.dist.alt(lx))
         
             
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
             
            envStat.qnorm3 <- ldply(norm.envStat, function(lx){
              if(!is.na(lx[[1]])){
                qnorm(pnorm(-3), mean = lx$parameters[[1]], sd   = lx$parameters[[2]])
              }else{
                return(NA)
              }
            })
            
            testZero <- envStat.qnorm3$V1 < 0
            testZero[is.na(testZero)] <- FALSE
           
            rec.dist.fit$V1[testZero] <- rec.dist.alt.fit$V1[testZero]
            
            
            
             norm.envStat.UCL <- ldply(norm.envStat, function(lx){
               if(!is.na(lx[[1]])){
                unname(lx$interval$limits[2])
                  }else{
                    return(NA)
               }
              })
             
            
            
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
             
             chby95 <- ddply(dfx, .(sampleNum), function(dfx) {
               if(dfx$nUncen[1] <3){
                 return(NA)
                 
                 
               }else{
                 chebyUCL(dfx$obs, dfx$cen)
               }
             })
             
             # Durant's method. Uses ROS + a limiting factor for GSD:
            
            dur.ros <- dlply(dfx, .(sampleNum), function(dfx){
              ros(dfx$obs,dfx$cen)
            })
            
            dur.mean <- ldply(dur.ros, mean)
            
            CoxLand.UCL <- ldply(dur.ros, function(lx){
              ros.df <- as.data.frame(lx)
              X <-  ros.df$modeled
              y <- log(X)
              n <- length(y)
              sd_y = sd(y)
              mean_y <- mean(y)
              
              xbar = exp(mean(y) + 0.5 * (sd_y)^2)
              ucl <- exp(mean_y + (0.5 * (sd_y)^2) + qnorm(0.95) * sqrt(
                ((sd_y^2)/2) +
                  ((sd_y^4)/(2*(n-1)))
              )
              )
              
              # if UCL > 5x mean, set to 5x mean
              
              if(ucl > 5 * mean(X)) ucl <- 5 * mean(X)
              return(ucl)
                      
            })
            
           
             
             #data.frame for UCL's
             rec.df <- data.frame(sampNum = rec.dist.fit$sampleNum,
                                  rec.dist=rec.dist.fit$V1, 
                                  Gamma = gamma.envStat.UCL$V1,
                                  LogNormal=ln.envStat.UCL$V1,
                                  Normal=norm.envStat.UCL$V1,
                                  CoxLand = CoxLand.UCL$V1,
                                  cheby95 = chby95$V1
             )
             
             
             myChoice <- function(dfx){
               if(dfx$rec.dist=="Gamma") dfx$rec.UCL <- dfx$Gamma
               if(dfx$rec.dist=="Log Normal") dfx$rec.UCL <- dfx$LogNormal
               if(dfx$rec.dist=="Normal") dfx$rec.UCL <- dfx$Normal
               dfx
             }
             
             rec.df <- ddply(rec.df, .(sampNum), function(dfx) myChoice(dfx))
          
            norm.envStat.UCL.coverage <- mean(norm.envStat.UCL$V1 >=pop.mean, na.rm=TRUE)
            norm.envStat.UCL.rope.coverage <- mean(norm.envStat.UCL$V1 >=pop.mean * input$ROPE, na.rm=TRUE)
            norm.envStat.UCL.miss <- 1 - mean(norm.envStat.UCL$V1[norm.envStat.UCL$V1 < pop.mean]/pop.mean, na.rm=TRUE)
            norm.envStat.UCL.bias <- sum((norm.envStat.UCL$V1 - pop.mean) / 
                                           pop.mean, na.rm=TRUE) / 
              length(norm.envStat.UCL$V1[!is.na(norm.envStat.UCL$V1)])
            norm.envStat.UCL.rmse <-  rmse(pop.mean, 
                                           norm.envStat.UCL$V1[!is.na(norm.envStat.UCL$V1)])
            norm.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                               norm.envStat.UCL$V1[!is.na(norm.envStat.UCL$V1)])))
            
            
            gamma.envStat.UCL.coverage <- mean(gamma.envStat.UCL$V1 >=pop.mean, na.rm=TRUE)
            gamma.envStat.UCL.rope.coverage <- mean(gamma.envStat.UCL$V1 >=pop.mean * input$ROPE, na.rm=TRUE)
            gamma.envStat.UCL.miss <- 1 - mean(gamma.envStat.UCL$V1[gamma.envStat.UCL$V1 < pop.mean]/pop.mean, na.rm=TRUE)
            gamma.envStat.UCL.bias <- sum((gamma.envStat.UCL$V1 - pop.mean) / 
                                           pop.mean, na.rm=TRUE) / 
              length(gamma.envStat.UCL$V1[!is.na(gamma.envStat.UCL$V1)])
            gamma.envStat.UCL.rmse <-  rmse(pop.mean, 
                                           gamma.envStat.UCL$V1[!is.na(gamma.envStat.UCL$V1)])
            gamma.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                            gamma.envStat.UCL$V1[!is.na(gamma.envStat.UCL$V1)])))
            
            
            
            ln.envStat.UCL.coverage <- mean(ln.envStat.UCL$V1 >=pop.mean, na.rm=TRUE)
            ln.envStat.UCL.rope.coverage <- mean(ln.envStat.UCL$V1 >=pop.mean * input$ROPE, na.rm=TRUE)
            ln.envStat.UCL.miss <- 1 - mean(ln.envStat.UCL$V1[ln.envStat.UCL$V1 < pop.mean]/pop.mean, na.rm=TRUE)
            ln.envStat.UCL.bias <- sum((ln.envStat.UCL$V1 - pop.mean) / 
                                            pop.mean, na.rm=TRUE) / 
              length(ln.envStat.UCL$V1[!is.na(ln.envStat.UCL$V1)])
            ln.envStat.UCL.rmse <-  rmse(pop.mean, 
                                            ln.envStat.UCL$V1[!is.na(ln.envStat.UCL$V1)])
            ln.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                            ln.envStat.UCL$V1[!is.na(ln.envStat.UCL$V1)])))
            
            cheby95.UCL.coverage <- mean(chby95$V1 >=pop.mean, na.rm=TRUE)
            cheby95.UCL.rope.coverage <- mean(chby95$V1 >=pop.mean * input$ROPE, na.rm=TRUE)
            cheby95.UCL.miss <- 1 - mean(chby95$V1[chby95$V1 < pop.mean]/pop.mean, na.rm=TRUE)
            cheby95.UCL.bias <- sum((chby95$V1 - pop.mean) / 
                                         pop.mean, na.rm=TRUE) / 
              length(chby95$V1[!is.na(chby95$V1)])
            cheby95.UCL.rmse <-  rmse(pop.mean, 
                                      chby95$V1[!is.na(chby95$V1)])
            cheby95.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                            chby95$V1[!is.na(chby95$V1)])))
            
            
            rec.UCL.coverage <- mean(rec.df$rec.UCL >=pop.mean, na.rm=TRUE)
            rec.UCL.rope.coverage <- mean(rec.df$rec.UCL >=pop.mean * input$ROPE, na.rm=TRUE)
            rec.UCL.miss <- 1 - mean(rec.df$rec.UCL[rec.df$rec.UCL < pop.mean]/pop.mean, na.rm=TRUE)
            rec.UCL.bias <- sum((rec.df$rec.UCL - pop.mean) / 
                                            pop.mean, na.rm=TRUE) / 
              length(rec.df$rec.UCL[!is.na(rec.df$rec.UCL)])
            rec.UCL.rmse <-  rmse(pop.mean, 
                                  rec.df$rec.UCL[!is.na(rec.df$rec.UCL)])
            rec.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                             rec.df$rec.UCL[!is.na(rec.df$rec.UCL)])))
            
            
            
            
            CoxLand.UCL.coverage <- mean(CoxLand.UCL$V1 >=pop.mean, na.rm=TRUE)
            CoxLand.UCL.rope.coverage <- mean(CoxLand.UCL$V1 >=pop.mean * input$ROPE, na.rm=TRUE)
            CoxLand.UCL.miss <- 1 - mean(CoxLand.UCL$V1[CoxLand.UCL$V1 < pop.mean]/pop.mean, na.rm=TRUE)
            CoxLand.UCL.bias <- sum((CoxLand.UCL$V1 - pop.mean) / 
                                         pop.mean, na.rm=TRUE) / 
              length(CoxLand.UCL$V1[!is.na(CoxLand.UCL$V1)])
            CoxLand.UCL.rmse <-  rmse(pop.mean, 
                                         CoxLand.UCL$V1[!is.na(CoxLand.UCL$V1)])
            
            CoxLand.UCL.rMedSe <- sqrt(median(se(pop.mean, 
                                                 CoxLand.UCL$V1[!is.na(CoxLand.UCL$V1)])))
            CoxLand.UCL.bias <- sum((CoxLand.UCL$V1 - pop.mean) / 
                                      pop.mean, na.rm=TRUE) / 
              length(CoxLand.UCL$V1[!is.na(CoxLand.UCL$V1)])
            
             envStat.UCL.coverage = c(gamma.envStat.UCL.coverage,
                                      ln.envStat.UCL.coverage,
                                      norm.envStat.UCL.coverage,
                                      rec.UCL.coverage,
                                      CoxLand.UCL.coverage,
                                      cheby95.UCL.coverage)
             names(envStat.UCL.coverage) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", "Cox.Lands", 
                                              "Cheby95")
             print(envStat.UCL.coverage)
            
            
            envStat.UCL.rope.coverage = c(gamma.envStat.UCL.rope.coverage,
                                     ln.envStat.UCL.rope.coverage,
                                     norm.envStat.UCL.rope.coverage,
                                     rec.UCL.rope.coverage,
                                     CoxLand.UCL.rope.coverage,
                                     cheby95.UCL.rope.coverage)
            names(envStat.UCL.rope.coverage) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", 
                                                  "Cox.Lands", "Cheby95")
            print(envStat.UCL.rope.coverage)
            
            
            envStat.UCL.bias <- c(gamma.envStat.UCL.bias,
                                  ln.envStat.UCL.bias,
                                  norm.envStat.UCL.bias,
                                  rec.UCL.bias,
                                  CoxLand.UCL.bias,
                                  cheby95.UCL.bias)
            
            names(envStat.UCL.bias) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", "Cox.Lands",
                                         "Cheby95")
            print(envStat.UCL.bias)
            
            envStat.UCL.rmse <- c(gamma.envStat.UCL.rmse/pop.mean,
                                  ln.envStat.UCL.rmse/pop.mean,
                                  norm.envStat.UCL.rmse/pop.mean,
                                  rec.UCL.rmse/pop.mean,
                                  CoxLand.UCL.rmse/pop.mean,
                                  cheby95.UCL.rmse/pop.mean)
            names(envStat.UCL.rmse) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", "CoxLands",
                                         "Cheby95")
            print(envStat.UCL.rmse)
            
            envStat.UCL.miss <- c(gamma.envStat.UCL.miss,
                                  ln.envStat.UCL.miss,
                                  norm.envStat.UCL.miss,
                                  rec.UCL.miss,
                                  CoxLand.UCL.miss,
                                  cheby95.UCL.miss)
            names(envStat.UCL.miss) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", "COxLands",
                                         "Cheby95")
            print(envStat.UCL.miss)
            
            RMedSE <- c(gamma.UCL.rMedSe/pop.mean,
                        ln.UCL.rMedSe/pop.mean,
                        norm.UCL.rMedSe/pop.mean,
                        rec.UCL.rMedSe/pop.mean,
                        CoxLand.UCL.rMedSe/pop.mean,
                        cheby95.UCL.rMedSe/pop.mean)
            names(RMedSE) <- c("Gamma", "Log-Normal", "Normal", "Rec. UCL", "COxLands",
                               "Cheby95")
            print(RMedSE)

# Bayesian ----------------------------------------------------------------


            
            
# Bayes_lognormal <- ddply(dfx, .(sampleNum), function(dfx){
#   if(dfx$nUncen[1] <3){
#     return(NA)
#     
#     
#   }else{
#     EH_Cen_Bayes(dfx$obs, dfx$cen, restriction = 5, conf=c(0.05,0.95))[[2]]
#   }
# }
# 
# )
              
# Bayes_coverage <- mean(Bayes_lognormal >= pop.mean, na.rm=TRUE)

             results = list(
             data = dfx,
             dist = lambdas,
             rec.dist.fit = rec.dist.fit,
             envStat.UCL.coverage = envStat.UCL.coverage,
             envStat.UCL.rope.coverage = envStat.UCL.rope.coverage,
             envStat.UCL.bias = envStat.UCL.bias,
             envStat.UCL.rmse = envStat.UCL.rmse,
             pop.mean = pop.mean,
            
            
             rec.df = rec.df
                                      
             )

})

myValue <- input$myValue1
             
             
             return(results)
   
    
  })
    
             
             
             
             
             
             

             

  
  output$distGraph <-renderPlot({
    if(input$distribution=='Log-Normal'){
      .x <- seq(qlnorm(0.01, input$meanlog, input$sdlog), qlnorm(0.99, input$meanlog, input$sdlog), length.out=100)
      plot(.x, dlnorm(.x, meanlog=input$meanlog, sdlog=input$sdlog), xlab="x", ylab="Density",
           type="l")
      abline(h=0, col="gray")
    }
    
    if(input$distribution=='Gamma'){
      .x <- seq(qgamma(0.01, shape=input$shape, scale=input$scale),qgamma(0.99, shape=input$shape, scale=input$scale) , length.out=100)
      plot(.x, dgamma(.x, shape=input$shape, scale=input$scale), xlab="x", ylab="Density",
           type="l")
      abline(h=0, col="gray")
    }
    
    if(input$distribution=='Normal'){
      .x <- seq(qnorm(0.01, mean=input$mu, sd=input$sd),qnorm(0.99, mean=input$mu, sd=input$sd), length.out=100)
      plot(.x, dnorm(.x, mean=input$mu, sd=input$sd), xlab="x", ylab="Density",
           type="l")
      abline(h=0, col="gray")
    }
    
  })
  
  output$Lambdaview <- renderPlot({
    results <- datasetInput()
    dist <- results$dist
    dist <- dist[!is.na(dist)]
    lambdas <- ldply(dist, function(lx) data.frame(objective = lx$objective, lambda=lx$lambda) )
    boxplot(lambdas$objective ~ lambdas$lambda, col="light green", xlab=expression(lambda), 
            ylab="PPCC", main="boxcoxCensored Output", ylim=c(0.6,1))
  })
  


 

output$judgementPlot <- renderPlot({
  results <- datasetInput()
  myDf <- data.frame(
    dist = names(results$envStat.UCL.coverage),
    rope = results$envStat.UCL.rope.coverage,
    bias = results$envStat.UCL.bias,
    rmse = results$envStat.UCL.rmse)
    xyplot(rmse ~ bias, data=myDf, groups=dist, auto.key=TRUE, type=c("g","p"), 
           par.settings=list(superpose.symbol = list(pch =c(1,2,3,4,5))),
           main = "Bias and RMSE of UCL", xlab="Bias", ylab="RMSE/Population Mean")
})

  output$UCLgraph <- renderPlot({
    results <- datasetInput()
    rec.df <- results$rec.df
    pop.mean <- results$pop.mean
    x1 <- xyplot(value ~ variable, 
           data = melt(rec.df, id.vars = c("sampNum", "rec.dist")),
           panel = function(...){
             panel.xyplot(..., jitter.x=TRUE,
                          factor=1)
             panel.abline(h=log10(pop.mean), lty=3)
             print(pop.mean)
           },
           #layout=c(1,4),
           #as.table=TRUE,
           scales=list(y=list(log=10)),
           xlab="Distribution",
           ylab="UCL",
           main="Plot of UCL's",
           
           pch=19, col="black", cex=0.3)
    print(x1)
  })


output$UCLcoverage <- renderDataTable({
  results <- datasetInput()
  data.frame(Distribution = names(results$envStat.UCL.coverage), Overall.Coverage = signif(results$envStat.UCL.coverage,2), 
             ROPE.coverage = signif(results$envStat.UCL.rope.coverage,2),
             Bias = signif(results$envStat.UCL.bias,2),
             RMSE = signif(results$envStat.UCL.rmse,2))
  
  
  
}, options=list(paging=FALSE, ordering=FALSE, info=FALSE, searching=FALSE))

output$UCLrecommend <- renderDataTable({
  results <- datasetInput()
  
  rec.dist.fit <- table(factor(results$rec.dist.fit$V1))
  dfx <- data.frame(rec.dist.fit)
  names(dfx) <- c("Distribution", "Frequency")
    
  return(dfx)
}, options=list(paging=FALSE, ordering=FALSE, info=FALSE, searching=FALSE))

  
# output$Bayesgraph <- renderPlot({
#   results <- datasetInput()
#   Bayes_lognormal <- results$Bayes_lognormal
#   x1 <- xyplot(V1 ~ sampleNum, type=c("g","h"),
#                data = Bayes_lognormal,
#                as.table=TRUE,
#                xlab="Simulation",
#                ylab="UCL",
#                main="Plot of UCL's")
#   print(x1)
# })
# 
# 
# output$bayesView <- renderPlot({
#   results <- datasetInput()
#   barplot( results$Bayes_coverage, col="light blue",
#            main="Bayes Coverage", ylim=c(0,1))
#   abline(h=0.95, lty=3)
#   
# })
#   

})
