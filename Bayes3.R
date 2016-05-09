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
mu[1] <- 300
sigmasquared[1] <- 40

ybar <- mean(raindata[,1])
ssquared <- var(raindata[,1])


tausquared_0 <- 300
mu_0 <- 300
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

