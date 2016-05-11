###################################################################################
# Author: Mattias Villani, Sveriges Riksbank and Stockholm University. 
#         E-mail: mattias.villani@gmail.se
# Script to illustrate numerical maximization of the Logistic or Probit regression
###################################################################################

###########   BEGIN USER INPUTS   ################
Probit <- 1 # If Probit <-0, then logistic model is used.
chooseCov <- c(1:16) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################


# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

library("msm")

# Loading data from file
Data<-read.table("SpamReduced.dat",header=TRUE)  # Spam data from Hastie et al.
y <- as.vector(Data[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- as.matrix(Data[,2:17]);
covNames <- names(Data)[2:length(names(Data))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

# Defining the functions that returns the log posterior (Logistic and Probit models). Note that the first input argument of

# this function must be the one that we optimize on, i.e. the regression coefficients.

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
                                      
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

LogPostProbit <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
                                      
  # MQ change:
  # instead of logLik <- sum(y*log(pnorm(linPred)) + (1-y)*log(1-pnorm(linPred)) ) type in the equivalent and
  # much more numerically stable; 
  logLik <- sum(y*pnorm(linPred, log.p = TRUE) + (1-y)*pnorm(linPred, log.p = TRUE, lower.tail = FALSE))
  # The old expression: logLik2 <- sum(y*log(pnorm(linPred)) + (1-y)*log(1-pnorm(linPred)) )
  abs(logLik) == Inf
  #print('-----------------')
  #print(logLik)
  #print(logLik2)
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

# Calling the optimization routine Optim. Note the auxilliary arguments that are passed to the function logPost
# Note how I pass all other arguments of the function logPost (i.e. all arguments except betaVect which is the one that we are trying to optimize over) to the R optimizer.
# The argument control is a list of options to the optimizer. Here I am telling the optimizer to multiply the objective function (i.e. logPost) by -1. This is because
# Optim finds a minimum, and I want to find a maximum. By reversing the sign of logPost I can use Optim for my maximization problem.

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,dim(X)[2])); 
# Or a random starting vector: as.vector(rnorm(dim(X)[2]))
# Set as OLS estimate: as.vector(solve(crossprod(X,X))%*%t(X)%*%y); # Initial values by OLS

if (Probit==1){
  logPost = LogPostProbit;
} else{
  logPost = LogPostLogistic;
}
  
OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# # Printing the results to the screen
# names(OptimResults$par) <- covNames # Naming the coefficient by covariates
# approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
# names(approxPostStd) <- covNames # Naming the coefficient by covariates
# print('The posterior mode is:')
# print(OptimResults$par)
# print('The approximate posterior standard deviation is:')
# approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
# print(approxPostStd)
drawnbetaposterior <- rmvnorm(1, mean = OptimResults$par, sigma = solve(-OptimResults$hessian))


betanow <- rmvnorm(1, sigma = 100 * diag(16))
betanow <- as.vector(betanow)
betanows <- matrix(0,nrow = 100, ncol = 16)
betanows[1,] <- betanow
j <- 2
while(j <  101){
  ueyes <- matrix(0, nrow = 1,ncol = length(y))
  for( i in 1:length(y)){
    if(y[i] == 0){
      ueyes[i] <- rtnorm(1, mean = X[i,] %*% betanows[j-1,], upper = 0)    
    } else{
      ueyes[i] <- rtnorm(1, mean = X[i,] %*% betanows[j-1,], lower = 0)
    }
  }
  ueyes <- as.vector(ueyes)
  #themvmodel <- as.data.frame(cbind(resp = ueyes, X))
  betanow <- lm(ueyes ~ X - 1)$coef
  betanows[j,] <- betanow
  j <-j + 1
}
plot(betanows[,1],type="l",ylim = c(-50,50),main = "beta coefficents vs # iterations",
     xlab = "# iterations", ylab = "beta coefficient values")
for(i in 2:16){
  lines(betanows[,i],col=i)
}

spamindicators <- matrix(0,nrow = length(y),ncol = 2)
for(i in 1:length(y)){
  if(pnorm(X[i,] %*% betanows[dim(betanows)[1],]) >= 0.5){
    spamindicators[i,2] <- 1    
  }
  if(pnorm(X[i,] %*% as.vector(drawnbetaposterior)) >= 0.5){
    spamindicators[i,1] <- 1    
  } 
}

paste("Classification rate for posterior mode coefficients: ",
      signif(100 * sum(y == spamindicators[,1]) /length(y),3), "%")
paste("Classification rate for gibbs sampler coefficients: ",
      signif(100 * sum(y == spamindicators[,2]) /length(y),3), "%")

