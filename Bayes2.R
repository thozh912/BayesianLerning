library(plotrix)
library(emdbook)
votredat <- matrix( c(208, 45, 46, 35, 110, 189, 34, 53, 88, 808,
                       403, 58, 74, 42, 146, 413, 127, 93, 57, 1413,
                       370, 51, 60, 47, 67, 401, 59, 61, 15, 1131,
                      383, 89, 86, 65, 45, 567, 74, 79, 17, 1405,
                       1364, 243, 266, 189, 368, 1570, 294, 286, 177, 4757),
                      nrow = 5, byrow = TRUE)
row.names(votredat) <- c("18-29","30-49","50-64","65+","Total")
colnames(votredat) <- c("M", "C", "FP", "KD", "MP", "S", "V", "SD", "Others", "Total")
votedata <- as.data.frame(votredat)
prioralphas <- c( 30, 6, 7, 6, 7, 30, 6, 6, 2)
j <- 1
thetavalues <- array(0, dim= c(4,9,1000))
while(j < 1001){
  for(i in 1:4){
    thetavalues[i,,j] <- rdirichlet(1,alpha = prioralphas + as.numeric(votedata[i,1:9]))
  }
  j <- j + 1
}
themean <- matrix(0,nrow = 4, ncol = 9)
thesd <- matrix(0,nrow = 4, ncol = 9)
for(i in 1:4){
  for(k in 1:9){
    themean[i,k] <- mean(thetavalues[i,k,])
    thesd[i,k] <- sd(thetavalues[i,k,])
  }
}
par(mfrow = c(2,2))
for( i in 1:4){
  plotCI(1:9,themean[i,],ui = themean[i,] + 2 *thesd[i,], li = themean[i,] - 2 *thesd[i,],
         xlab = "Political party #", ylab = "vote probability",
         main = c("posterior vote probabilities",paste("age group",i)))
}

counter <- rep(0,4)
for(i in 1:4){
  for(j in 1:1000){
    if( sum(thetavalues[i,1:4,j]) < sum(thetavalues[i,5:7,j]) ){
      counter[i] <- counter[i] + 1
    }
    
  }
}
counter <- counter / 1000
paste("Posterior probability of Red-Greens winning in age group",1,"is",counter[1],"percent")
paste("Posterior probability of Red-Greens winning in age group",2,"is",counter[2],"percent")
paste("Posterior probability of Red-Greens winning in age group",3,"is",counter[3],"percent")
paste("Posterior probability of Red-Greens winning in age group",4,"is",counter[4],"percent")

counter <- 0
for(k in 1:1000){
  Yais <- rmultinom(1,6300000, prob = c(0.2, 0.3, 0.3, 0.2))
  thetaais <- matrix(0,4,9)
  actualvotes <- matrix(0,4,9)
  for( i in 1:4){
    thetaais[i,] <- rdirichlet(1,alpha = prioralphas + as.numeric(votedata[i,1:9]))
  }
  
  for( i in 1:4){
    actualvotes[i,] <- rmultinom(1,Yais[i],prob = thetaais[i,])
  }
  
  allagevotes <- colSums(actualvotes)
  if(sum(allagevotes[1:4]) < sum(allagevotes[5:7])){
    counter <- counter + 1
  }
}
paste("posterior probability of red-green election win:",counter / 1000)

par(mfrow = c(1,1))
Japandata <- read.delim("JapanTemp.dat", header = TRUE, sep = "",stringsAsFactors = FALSE)

plot(Japandata$time, Japandata$temp)
lines(Japandata$time,predict(lm(temp ~ time + I(time^2) , data = Japandata)))
sd(Japandata$temp - predict(lm(temp ~ time + I(time^2) , data = Japandata)))
#its about 2.17
# we gonna pick s_0^2 = 25 and nu_0 = 50  andbeta_0 = c(11.58, 57.83, -50.82) and
#omega_0 is diag(3)

xes <- seq(0,10,0.001)
nus <- 1:100
ssquareds <- lseq(1,1000,100)

invchisq <- function(x,nu,ssquared){
  return( x^(-2 * (nu / 2 + 1)) * exp( - nu * ssquared / (2 * x^2)))
}

degrees <- matrix(0,100,5001)
for( i in 1:100){
  degrees[i,] <- invchisq(xes,nus[i],1) 
}

sdependence <- matrix(0,100,5001)
for(i in 1:100){
  sdependence[i,] <- invchisq(xes,1,ssquareds[i])
}
plot(xes,sdependence[1,],type = "l")
for( i in 1:100){
  lines(xes,sdependence[i,])
}

library(MASS)

thet <- matrix(ncol = 3, nrow = 100)
for(i in 1:100){
  X <- rchisq(1, 50)
  sigsquare <- 50 * 25 / X
  thet[i, ] <- mvrnorm(1, c(11.58, 57.83, -50.82), sigsquare*(1/2)*diag(3)) 
}


x <- cbind(1, Japandata$time)
x <- cbind(x, Japandata$time^2)
y <- list()

  for(j in 1:nrow(thet)){
    y[[j]] <- x%*%thet[j,]
  }
 
plot(Japandata$time,y[[1]],type = "l", ylim = c(10, 40))
for( i in 2:100){
  lines(Japandata$time,y[[i]])
}
points(Japandata$time, Japandata$temp, col ="orange")


bzero <- c(11.58, 57.83, -50.82)
betamle <- solve( t(x) %*% x) %*% t(x) %*% as.matrix(Japandata$temp)
bn <- solve(t(x)%*%x + 2*diag(3))%*% (t(x)%*%x%*% betamle + 2*diag(3) %*% bzero)


vn <- 365 + 50
snsquare <- (50*25 + t(Japandata$temp)%*%Japandata$temp + bzero%*%(2*diag(3))%*% bzero -
  t(bn)%*%(2*diag(3))%*%bn) / vn

omegan <- t(x)%*%x + 2*diag(3)

thet1 <- matrix(ncol = 3, nrow = 100)
for(i in 1:100){
  X <- rchisq(1, vn)
  sigsquare <- (vn *snsquare) / X
  thet1[i, ] <- mvrnorm(1, t(bn), as.numeric(sigsquare)*solve(omegan) )
}


y1 <- list()

for(j in 1:nrow(thet1)){
  y1[[j]] <- x%*%thet1[j,]
}

plot(Japandata$time,y1[[1]],type = "l", ylim = c(10, 40))
for( i in 2:100){
  lines(Japandata$time,y1[[i]])
}
points(Japandata$time, Japandata$temp, col ="orange")

class(thet1)
xtildes <- -thet1[,2] / (2 * thet1[,3])

pickthisday <-c()

vect <- 1:365 /365
for(i in 1:100){
  pickthisday[i] <- which.min(abs(vect-xtildes[i]))
}
hist(pickthisday,breaks = 365, xlim =c(0,365),
     main = "Distribution of hottest expected day of year")
