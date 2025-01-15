# Reproduce results from the paper from the paper "A Bayesian error model for synthesis and sequencing of oligonucleotides", Marrs, FW, Gratz, D, and Erkkila, TH.
#
#
# © 2024. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

setwd("INSERT_WD")
source("functions.R")

################################################################################
#### Simulation study 2
#
set.seed(100)

# Generate data
n <- 1e5   # number of reads, used 1e6 in paper
dat <- simulate_errors(simnum=2, nsim=n)  # second simulation in paper

names(dat)
dim(dat$ne)
xx <- 0:75
p2 <- sapply(xx, function(x) mean(dat$ne[,2]==x))
p4 <- sapply(xx, function(x) mean(dat$ne[,4]==x))
p5 <- sapply(xx, function(x) mean(dat$ne[,5]==x))

ylim = log10(c(p2, p4, p5))
ylim[is.infinite(ylim)] <- NA
plot(NA, NA, xlim=range(xx), ylim=range(ylim, na.rm=TRUE), 
     xlab="Number of errors", ylab="log10 distribution")
colors <- RColorBrewer::brewer.pal(8, "Set2")
lines(xx, log10(p2), col=colors[1], lwd=2)
lines(xx, log10(p4), col=colors[2], lwd=1, lty=2)
lines(xx, log10(p5), col=colors[3], lwd=1, lty=3)
legend("topright", 
       c("Synthesis errors", "Forward read errors", "Reverse read errors"), 
       lwd=c(2,1,1), lty=c(1:3), col=colors, bty='n'
       )



#################################################################
# Estimate model parameters
#


# Fit, Bayesian
writedir <- "./results/results_bayes_repro_1M"
dir.create(writedir)
mcmc.error(dat$ne[,1:5], writedir=writedir, nmcmc=1e3)   # need at least 1e4 for good mixing
params.noq <- read.table(file.path(writedir, "params.txt"), header=TRUE)
e.noq <- process_bayes(dat$ne, params.noq[-1,], 75, withq=FALSE)

# Fit, Bayesian with q
writedir <- "./results_bayes_with_q"
dir.create(writedir)
mcmc.error.qscore(dat$ne, writedir=writedir, nmcmc=1e3)   # need at least 1e4 for good mixing
params.q <- read.table(file.path(writedir, "params.txt"), header=TRUE)
e.q <- process_bayes(dat$ne, params.q[-1,], 75, withq=TRUE)

# Table of LLs for Bayesian fit (a la Table 1)
ii <- tail(1:nrow(params.noq),round(nrow(params.noq)/2))
ll.noq <- mean(params.noq[ii,8])
ll.q <- mean(params.q[ii,10])
dic.noq <- dic2(params.noq[ii,8])
dic.q <- dic2(params.q[ii,10])
cat("LLs (no Q and Q): \n")
cat(c(round(ll.noq/1e6,2), round(ll.q/1e6,2)), "\n")
cat("DIC (no Q and Q): \n")
cat(c(round(dic.noq/1e6,2), round(dic.q/1e6,2)), "\n")


# Frequentist fits
ests <- c("naive", "one-step", "converged")
runmat <- expand.grid(c(FALSE, TRUE), ests)
resfreq <- NULL
for(i in 1:nrow(runmat)){
  mod <- frequ_est(dat$ne, npl=75, est=runmat[i,2], boot=as.logical(runmat[i,1]), 
                   nboot=1000)   # used 1000 in paper
  resfreq <- cbind(resfreq, mod$vsim)
}




# Table of error estimates
phat.freq <- apply(resfreq, 2, mean)
se.freq <- apply(resfreq, 2, sd)



################################################################################
# make plots from paper
#   Note: may not exactly match due to sample sizes and shorter Bayesian MCMC chain
#

# Synthesis error rate
colors <- RColorBrewer::brewer.pal(8, "Set2")
par(mar=c(3,3,1,1), cex.axis=1.4, cex.lab=1.4, mgp=c(1.7,.55,0)) #, mfrow=c(4,1))
boxes <- cbind(resfreq[,seq(1, 18, by=3)], e.noq[,1],  e.q[,1])*100
boxplot(boxes, 
        ylab="Synthesis error rate (%)", col="white", xaxt="n",
        border=c(colors[c(1,2,1,2,1,2,3,4)]))
abline(h=dat$s1*100, 
       lwd=2, lty=2, 
       col="gray50")
axis(1, at=seq(1.5, 9, by=2), labels = c("norm.", "1step", "conv.", "bayes"))
legend("topleft", 
       c("Analytical", "Bootstrapped", "Bayesian (no q)",  "Bayesian (w/ q)"), col=colors, bty='n',
       lty=1, pch=NA
       )

# Forward read error rate
boxes <- cbind(resfreq[,seq(1, 18, by=3)+1], e.noq[,2],  e.q[,2])*100
boxplot(boxes, 
        ylab="Fwd. read error rate (%)", col="white", xaxt="n",
        border=c(colors[c(1,2,1,2,1,2,3,4)]))
abline(h=dat$r1*100, 
       lwd=2, lty=2, 
       col="gray50")
axis(1, at=seq(1.5, 9, by=2), labels = c("norm.", "1step", "conv.", "bayes"))

# Reverse read error rate
boxes <- cbind(resfreq[,seq(1, 18, by=3)+2], e.noq[,3],  e.q[,3])*100
boxplot(boxes, 
        ylab="Rev. read error rate (%)", col="white", xaxt="n",
        border=c(colors[c(1,2,1,2,1,2,3,4)]))
abline(h=dat$r2*100, 
       lwd=2, lty=2, 
       col="gray50")
axis(1, at=seq(1.5, 9, by=2), labels = c("norm.", "1step", "conv.", "bayes"))




# Bayesian estimated distribution and truth NO Q
# simulate distributions
nsim <- 1e3
npl <- 75
n.avg <- nrow(params)
x=0:npl
d0 <- d1 <- d2 <- matrix(0, length(x), nsim)
params <- params.noq
for(jj in 1:nsim){
  k <- sample(tail(1:nrow(params), n.avg), nsim, replace=TRUE)[jj]
  p0.hat <- params[k,1]
  alpha.hat <- params[k,2]
  q.hat <- params[k,3]
  f2 <- sapply(x, int.pi, npl=npl, q=q.hat, kmax=100)
  
  d0[,jj] <- log( dbinom(x, npl, p0.hat)*alpha.hat + (1-alpha.hat)*f2)
  
  mu1.hat <- params[k,4]
  eta1.hat <- params[k,5]
  d1[,jj] <- my.dbb(x, npl, mu1.hat, eta1.hat, log=TRUE)
  
  mu2.hat <- params[k,6]
  eta2.hat <- params[k,7]
  d2[,jj] <- my.dbb(x, npl, mu2.hat, eta2.hat, log=TRUE)
}

y.range <- c(-15,0)
par(mar=c(3,3,2,1), cex.axis=1.2, cex.lab=1.2, mgp=c(1.7,.55,0)) #, mfrow=c(4,1))
plot(NA, NA, 
     xlab="Number of errors",
     ylab="Logarithmic probability distribution", 
     pch=1, 
     col=colors[1], 
     xlim=range(x),
     ylim=y.range, main="No Q"
)

ub <- apply(d0, 1, quantile, .975)
lb <- apply(d0, 1, quantile, .025)
polygon(c(x, rev(x)), c(ub, rev(lb)), border = NA, col=paste0(colors[1], 60))
lines(x, apply(d0, 1, median), lty=1, lwd=1, col=colors[1])


ub <- apply(d1, 1, quantile, .975)
lb <- apply(d1, 1, quantile, .025)
polygon(c(x, rev(x)), c(ub, rev(lb)), border = NA, col=paste0(colors[2], 60))
lines(x, apply(d1, 1, median), lty=1, lwd=1, col=colors[2])


ub <- apply(d2, 1, quantile, .975)
lb <- apply(d2, 1, quantile, .025)
polygon(c(x, rev(x)), c(ub, rev(lb)), border = NA, col=paste0(colors[3], 60))
lines(x, apply(d2, 1, median), lty=1, lwd=1, col=colors[3])

xx <- 0:75
pts1 <- sapply(xx, function(x) mean(dat$ne[,2] == x))
points(xx, log(pts1), pch=1,
       col=colors[1])

pts1 <- sapply(xx, function(x) mean(dat$ne[,4] == x))
points(xx, log(pts1), pch=2,
       col=colors[2])

pts1 <- sapply(xx, function(x) mean(dat$ne[,5] == x))
points(xx, log(pts1), pch=3,
       col=colors[3])


legend("topright",
       c("Synthesis", "Read 1", "Read 2"),
       pch=1:3, col=colors, bty="n", lty=1)




