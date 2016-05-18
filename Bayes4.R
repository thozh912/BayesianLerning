library(stats)
library(mvtnorm)
library(coda)
ebaydata <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
head(ebaydata)
ebaydata <- ebaydata[,-2]
glmresult <- glm(nBids ~.-1 ,  family = "poisson", ebaydata)
summary(glmresult)

#It looks like VerifyID, Sealed, MajBlem, LogBook and MinBidShare
#are significant predictors in this model.

str(glmresult)
paste("The maximum liklihood estimator of beta coefficients:")
glmresult$coefficients
LogPostPoisson <- function(betaVect,y,X){
  X <- as.matrix(X)
  nPara <- length(betaVect)
  linPred <- X%*% betaVect
  Loglik <- sum( y * linPred  - exp(linPred))
  Logprior <- dmvnorm(t(betaVect), rep(0,8), 
                      sigma = 100 * solve( t(X)%*%X ), log=TRUE)
  return(Loglik + Logprior)
}

initVal <- rep(0,8)
OptimResults<-optim(initVal,LogPostPoisson,gr=NULL,
                    y = ebaydata$nBids,X = ebaydata[,-1],method=c("BFGS"),
                    control=list(fnscale=-1),hessian=TRUE)
Jay <- -OptimResults$hessian
betatilde <- OptimResults$par
posteriordraw <- rmvnorm(1, mean = betatilde, sigma = solve(Jay))

#Metropolis-Hastings Algorithm for any logpostfunc
outerfunc <- function(propdensity,c,innerfunc,...){
  nIter <-10000
  startmeanvector <- as.vector(propdensity[[1]])
  sigma <- as.matrix(propdensity[[2]])
  draws <- matrix(0,nrow = nIter,ncol = length(startmeanvector))
  draws[1,] <- startmeanvector
  for(i in 2:nIter){
    proposaldraw <- rmvnorm(1, mean = draws[i-1,],sigma = sigma)
    unifdraw <- runif(1)
    alpha <- min(1,exp(innerfunc(t(proposaldraw),...) - innerfunc(draws[i-1,],...)))
    if(unifdraw < alpha){
      draws[i,] <- proposaldraw
    }else{
      draws[i,] <- draws[i-1,]
    }
  }
  
  return(draws)
}

prior <- list(betatilde, solve(Jay)) 
c <- 2.4/sqrt(8)

mhresult <- outerfunc(prior,c,LogPostPoisson,
                      y = ebaydata$nBids , X = ebaydata[,-1])

par(mfrow = c(2,2))
for(i in 1:dim(mhresult)[2]){
  plot(mhresult[,i],type="l",col=i,
       ylab = c("beta coef. no: ",i),xlab = "iteration")
}
effsizes <- effectiveSize(as.mcmc(mhresult))

burnins <- mhresult[1:1000,]
withoutburnins <- mhresult[1001:dim(mhresult)[1],]
phis <- exp(withoutburnins)
for(i in 1:dim(mhresult)[2]){
  hist(phis[,i], breaks = 100,col=i, main = c("Posterior dist. of phi coef. no: ",i),
       xlab = "phi value",border = "white")
}

newX <- c(1, 1, 1, 0, 0, 0, 1, 0.5)

predictresult <- withoutburnins %*% newX
predictpoissonpar <- exp(predictresult)
predictprob <- ppois(0.5,lambda = predictpoissonpar)
finalres <- mean(predictprob)
paste("Mean probability for no bids on the new auction:",finalres)

