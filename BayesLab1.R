library(geoR)
library(emdbook)

#1
s = 14
n  = 20
f = n - s
alph = 2
bet = 2

postbetanums <- rbeta(10000, alph + s, bet + f)
paste("Theoretical posterior beta distr. mean:",
      round((alph + s) / (alph + s + bet + f),3))
paste("Theoretical posterior beta distr. sd:",
      round(sqrt((alph + s) * (bet + f) / 
                   ((alph + s + bet + f)^2 * (alph + s + bet + f + 1))),3))
hist(postbetanums,breaks = 50)
paste("Posterior beta sample mean:",round(mean(postbetanums),3))
paste("Posterior beta sample sd:",round(sd(postbetanums),3))
#It checks out to 3 decimal places

postbetalessthan <- postbetanums < 0.4
ratio <- sum(postbetalessthan) / 10000
paste("Theoretical prob. posterior theta less than 0.4:",
      round(pbeta(0.4,alph + s,bet + f),5))
paste("Posterior beta sample prob. theta less than 0.4:",round(ratio,5))
# similar numbers

postlogit <- density(log(postbetanums / (1 - postbetanums)))
postlogit
plot(postlogit, main="Posterior distr. of logit(theta)")
hist(log(postbetanums / (1 - postbetanums)), breaks = 50)

#2
lognormmu <- 3.5

y <- c(14,
       25, 45, 25, 30, 33, 19, 50, 34 , 67)

gini <- function(sigma){
  return( 2*pnorm(sigma / sqrt(2)) - 1)
}
#sigma must be a number
lognormPosterior <- function(sigma,y){
  n <- length(y)
  
  liklihood <- prod( dlnorm(y,meanlog = lognormmu,sdlog = sigma))
  posterior <- liklihood * 1/sigma^2
  return(posterior)
}

sigmavals <- seq(0.001,1,0.01)

lognormpostvals <-c()
for( i in 1:length(sigmavals)){
  lognormpostvals[i] <- lognormPosterior(sigmavals[i],y)
}
plot(sigmavals,lognormpostvals,
     main = "Posterior dist for sigmas given non-informative prior",
     ylab = "prop. to posterior prob dist. of sigma")

invChidensity <- function(x,y,konstant = 1){
  n <- length(y)
  tausquared <- sum((log(y) - 3.5)^2) / n
  res <- konstant * x ^ (-2 * ( n / 2 + 1)) * exp( -1/(2 * x^2) * n * tausquared)
  return(res)
}

x <- seq(0,1,0.000025)
x <- x[-1]
maxval <- max(invChidensity(x,y,konstant = 1), na.rm = TRUE)
plot(x,invChidensity(x,y,konstant = 1/maxval),type ="l",
     main = "inverse chi-squared distribution of sigma squared")
randvals <- runif(40000)
pickvals <-c()
for(i in 1:40000){
  if(randvals[i] < invChidensity(x[i],y,konstant = 1/maxval)){
   pickvals <- c(pickvals,x[i]) 
  }
}
hist(pickvals, breaks = 50,main = "Draws from the posterior inverse chi-squared distr.")
hist(gini(sqrt(pickvals)),breaks = 50, main = "Posterior Gini coef. distr.")

#3
radianobs <- c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
#kappa must be a constant
vonMisesPosterior <- function(kappa,y){
  n <- length(y)
  
  liklihood <- prod( exp( kappa * cos(y - 2.39)) / 
                       (2 * pi * besselI(x = kappa, nu = 0)))
  posterior <- liklihood * dexp(x = kappa,rate = 1)
  return(posterior)
}
kappavals <- seq(0.001,5.001,0.1)
#kappavals <- lseq(0.001,5,100)
vonMisespostvals <-c()
for( i in 1:length(kappavals)){
  vonMisespostvals[i] <- vonMisesPosterior(kappavals[i],radianobs)
}

plot(kappavals,vonMisespostvals,main="Posterior distr. of Kappa",
     xlab = "Prior Kappa",
     ylab = "posterior Kappa prob. proportion")
paste("Posterior mode:",signif(kappavals[which.max(vonMisespostvals)]))

#4
diseasedata <- data.frame(pop = c(120342, 235967,
           243745,
           197452, 
           276935, 
           157222))
diseasedata <- cbind(diseasedata,cases = c(2,5,3,5,3,1))
#searchforbeta
beta <- seq(1,5,0.1) + 0.001
#gammanums <- rgamma(10000, 4*beta, beta)

plot(beta,pgamma(5,4*beta, beta) - pgamma(3,4*beta, beta),
     main = "Search for beta",ylab="P(3<lambda<5)")
abline(h = 0.5,col ="red")
paste("optimal beta value:",beta[9] - 0.001)

x <- seq(0.001,10,0.01)
plot(x, dgamma(x,1.8 * 4,1.8), ,type ="l",main="Prior and Posterior distr. for Lambda",
     ylim = c(0,1.3),ylab="Prob. distr.")

newshape <- 1.8 * 4 + sum(diseasedata$cases)
newrate <- 1.8 + sum(diseasedata$pop) / 100000
lines(x, dgamma(x,newshape,newrate), col = "red")
legend("topright",legend = c("Prior","Posterior"),lty = c(1,1), col =c("black","red"))
paste("posterior probability lambda between 3 and 5:",
      signif(pgamma(5,newshape,newrate) - pgamma(3,newshape,newrate),3))
