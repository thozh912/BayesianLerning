s = 14
n  = 20
f = n - s
alph = 2
bet = 2

postbetanums <- rbeta(10000, alph + s, bet + f)
paste("Theoretical posterior beta distr. mean:",round((alph + s) / (alph + s + bet + f),3))
paste("Theoretical posterior beta distr. sd:",round(sqrt((alph + s) * (bet + f) / 
                                                           ((alph + s + bet + f)^2 * (alph + s + bet + f + 1))),3))
hist(postbetanums,breaks = 50)
paste("Posterior beta sample mean:",round(mean(postbetanums),3))
paste("Posterior beta sample sd:",round(sd(postbetanums),3))
#It checks out to 3 decimal places

postbetalessthan <- postbetanums < 0.4
ratio <- sum(postbetalessthan) / 10000
paste("Theoretical prob. posterior theta less than 0.4:",round(pbeta(0.4,alph + s,bet + f),5))
paste("Posterior beta sample prob. theta less than 0.4:",round(ratio,5))

#similar