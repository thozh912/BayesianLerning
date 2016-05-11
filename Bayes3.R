raindata <- read.table(file = "https://raw.githubusercontent.com/STIMALiU/BayesLearnCourse/master/Labs/rainfall.dat")
head(raindata)
n <- dim(raindata)[1]

#This is a draw from Inv-chisq(n-1,ssquared)
invchisqaccordingtoslide <- function(n,ssquared){
  eks <- rchisq(1, df = n - 1)
  sigmasquared <- (n - 1) * ssquared / eks
  return(sigmasquared)
}


sigmasquared <- rep(0,1000)
mu <- rep(30,1000)
mu[1] <- 30
sigmasquared[1] <- 1000

ybar <- mean(raindata[,1])
ssquared <- var(raindata[,1])


tausquared_0 <- 300

mu_0 <- 20

nu_0 <- 4
sigmasquared_0 <- 40

j <- 2
while(j <  1001){
    w <- (n / sigmasquared[j-1]) / ((n / sigmasquared[j-1]) + 1 / tausquared_0) 
    tausquared_n <- solve(n / sigmasquared[j-1] + 1 / tausquared_0)
    mu_n <- w * ybar + ( 1 - w) * mu_0
    mu[j] <- rnorm(1, mean = mu_n, sd = sqrt( tausquared_n))
    secondtermbelow <- (nu_0 * sigmasquared_0 + sum((raindata - mu[j])^2)) / (n + nu_0)
    sigmasquared[j] <- invchisqaccordingtoslide(n,secondtermbelow)
    j <- j + 1
}

plot(mu, type = "l")
plot(sigmasquared, type = "l")

muac <- acf(mu)
(sum(muac$acf[2:31])*2)+1 #the IF

sigcf <- acf(sigmasquared)
(sum(sigcf$acf[2:31])*2)+1 #the IF

x <- as.matrix(raindata['V1'])
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
ndistr <- dnorm(xGrid , mean = mu[1000], sd = sqrt(sigmasquared[1000]))

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax))
lines(xGrid, ndistr, col = "red") 
legend("topright", box.lty = 1, legend = c("Data histogram",'Normal Model'), 
                                          col = c("black", 'red'), lwd = 2)


#1.b
x <- as.matrix(raindata['V1'])

j <- NULL
# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(30,nComp) # Prior mean of theta
tau2Prior <- rep(300,nComp) # Prior std theta
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 120 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  thetaDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    thetaDraws[j] <- rgamma(1,param[j],1)
  }
  thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
  return(thetaDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
theta <- vector("list", 1000)
theta[[1]] <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- vector("list", 1000)
sigma2[[1]] <- rep(var(x),nComp)
probObsInComp <- vector("list", 1000)
probObsInComp[[1]] <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))


for (k in 1:nIter){
  #message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  #print(nAlloc)
  # Update components probabilities
  w <- rDirichlet(alpha + nAlloc)
  
  # Update theta's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[[k]][[j]]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    theta[[k+1]][[j]] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[[k+1]][[j]] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - theta[[k+1]][[j]])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[[k]][j] <- w[j]*dnorm(x[i], mean = theta[[k]][[j]], sd = sqrt(sigma2[[k]][[j]]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp[[k]]/sum(probObsInComp[[k]])))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k == 120)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,theta[[k]][[j]],sd = sqrt(sigma2[[k]][[j]]))
      mixDens <- mixDens + w[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}



theta1 <- c()
for(i in 1:120){
 theta1[i] <- theta[[i]][[1]]
}
plot(theta1, type ="l", ylim = c(0, 90))
t1ac <- acf(theta1[30:120])
(sum(t1ac$acf[2:19])*2)+1 #the IF


theta2 <- c()
for(i in 1:120){
  theta2[i] <- theta[[i]][[2]]
}
plot(theta2, type = "l")
ac <- acf(theta2)
(sum(ac$acf[2:20])*2)+1 #the IF

sigma1 <- c()
for(i in 1:120){
  sigma1[i] <- sigma2[[i]][[1]]
}
plot(sigma1, type ="l")
sig1 <- acf(sigma1)
(sum(sig1$acf[2:20])*2)+1 #the IF

sigmaSecond <- c()
for(i in 1:120){
  sigmaSecond[i] <- sigma2[[i]][[2]]
}
plot(sigmaSecond, type = "l")
sig2 <- acf(sigmaSecond)
(sum(sig2$acf[2:20])*2)+1 #the IF

probObsInComp1 <- c()
for(i in 1:120){
  probObsInComp1[i] <- probObsInComp[[i]][[1]]
}
plot(probObsInComp1, type ="l")
p1 <- acf(probObsInComp1)
(sum(p1$acf[2:20])*2)+1 #the IF

probObsInComp2 <- c()
for(i in 1:120){
  probObsInComp2[i] <- probObsInComp[[i]][[2]]
}
plot(probObsInComp2, type="l")
p2 <- acf(probObsInComp2)
(sum(p2$acf[2:20])*2)+1 #the IF

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density"), col=c("black","red"), lwd = 2)


#1.c
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, ndistr, type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)
#the mixture seems to be the best model in terms of fitting





#2
###########   BEGIN USER INPUTS   ################
Probit <- 1 # If Probit <-0, then logistic model is used.
chooseCov <- c(1:16) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################


# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

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
  
  logLik <- sum(y*log(pnorm(linPred)) + (1-y)*log(1-pnorm(linPred)) );
  abs(logLik) == Inf
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}



# Calling the optimization routine Optim. Note the auxilliary arguments that are passed to the function logPost
# Note how I pass all other arguments of the function logPost (i.e. all arguments except betaVect which is the one that we are trying to optimize over) to the R optimizer.
# The argument control is a list of options to the optimizer. Here I am telling the optimizer to multiply the objective function (i.e. logPost) by -1. This is because
# Optim finds a minimum, and I want to find a maximum. By reversing the sign of logPost I can use Optim for my maximization problem.
initVal <- as.vector(solve(crossprod(X,X))%*%t(X)%*%y); # Initial values by OLS

if (Probit==1){
  logPost = LogPostProbit;
} else{
  logPost = LogPostLogistic;
}

OptimResults <-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
names(OptimResults$par) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- covNames # Naming the coefficient by covariates
print('The posterior mode is:')
print(OptimResults$par)
print('The approximate posterior standard deviation is:')
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
print(approxPostStd)

