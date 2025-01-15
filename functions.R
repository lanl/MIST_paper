# Functions to estimate hierarchical Bayesian error model
#  from the paper "A Bayesian error model for synthesis and sequencing of oligonucleotides", Marrs, FW, Gratz, D, and Erkkila, TH.
#
# © 2024. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


################################################################################
#
#  Bayesian error model
#

# MCMC for Bayesian model estimation,
mcmc.error <- function(ne,   # error counts of five types
                       writedir=getwd(),  # where to write results to
                       nmcmc=1e4,   # number of mcmc samples
                       nburn=round(nmcmc*3/4),    # number of samples to omit
                       thin=10,    # save interval
                       save.latent.params=FALSE,   # should latent parameters be saved?
                       save.burn=TRUE,   # should the burned iterations be saved?
                       verbose=FALSE,   # should progress be printed out?
                       print.interval=round(nmcmc/10)   # how often should progress be printed out?
)
{
  joint.prop=FALSE
  
  nsd <- round(nmcmc/20)
  nn <- nrow(ne)
  
  
  ne <- as.matrix(ne)
  ne[is.na(ne)] <- 0  
  if(ncol(ne) != 5){
    stop("wrong input size of ne")
  }
  ngen <- nrow(ne)
  n.read <- rowSums(ne[, 1:5])
  npl <- my.mode(n.read)
  if(npl != min(n.read)){
    warning(paste0("Minimal rowsum ", min(n.read), " is not equal current payload of ", npl))
  }
  
  e0 <- ne[,2]  #rowSums(ne.sub[,2:3]) 
  e0 <- e0[e0<npl]
  
  p0 <- mean(e0[e0 < npl/10] / npl)
  p0 <- phat0 <- mean(p0, 1 - mean(e0 == 0)^(1/npl))
  
  x = round(npl/10):round(npl*2/3)
  mini.distr <- sapply(x, function(z) sum(e0 == z))
  mini.distr <- mini.distr/sum(mini.distr)
  y = log(mini.distr)
  q0 <- 1 - exp(coef(lm(y[!is.infinite(y)] ~ x[!is.infinite(y)]))[2])
  
  xx <- c(0:1e5)/1e5
  dg <- (1-q0)^c(xx*npl) # / pgeom(npl, q.hat)
  y.temp <- sample(xx, nn, replace=TRUE, prob=dg)
  qbar <- mean(y.temp)
  
  
  alpha.temp <- (mean(e0/npl) - qbar) / (p0 - qbar)  # match expectation
  alpha.temp <- min(alpha.temp, .9999)
  
  alphas <- 10^seq(log10(alpha.temp) - .1, min(log10(alpha.temp) + .1, log10(.99999)), length.out=100)
  ll <- 0*alphas
  for(i in 1:length(alphas)){
    dg <- dgeom(e0, q0) / pgeom(npl, q0)
    ll[i] <- sum(log(alphas[i]*dbinom(e0, npl, p0) + (1-alphas[i])*dg))
  }
  # plot(alphas, ll)
  alpha0 <- alphas[which.max(ll)]
  
  ## MLE step to really hone in on p0, q0, alpha0
  params0 <- c(p0, alpha0, q0)
  res <- optim(qnorm(params0), neg.ll.mix, x=e0, npl=npl)
  if(res$convergence == 0){
    params0 <- pnorm(res$par)
  }
  p0 <- params0[1]
  alpha0 <- params0[2]
  q0 <- params0[3]
  
  
  mu1 <- mean(ne[,4]/npl)
  eta1 <- mu1*(1-mu1) / var(ne[,4]/npl) - 1
  res = optim(c(qnorm(mu1), log(eta1)),
              beta.opt,
              x=ne[,4]/npl)
  if(res$convergence == 0){
    mu1 <- pnorm(res$par[1])
    eta1 <- exp(res$par[2])
  }
  
  mu2 <- mean(ne[,5]/npl)
  eta2 <- mu2*(1-mu2) / var(ne[,5]/npl) - 1
  res = optim(c(qnorm(mu2), log(eta2)),
              beta.opt,
              x=ne[,5]/npl)
  if(res$convergence == 0){
    mu2 <- pnorm(res$par[1])
    eta2 <- exp(res$par[2])
  }
  
  a.mu1 <- mu1
  b.mu1 <- (1-mu1)
  a.eta1 <- eta1
  b.eta1 <- 1
  
  a.mu2 <- mu2
  b.mu2 <- (1-mu2)
  a.eta2 <- eta2
  b.eta2 <- 1
  
  
  ## Initialize p1i and p2i
  pi <- matrix(0, ngen, 3)  # latent probability of high dispersion errors for synthesis, read 1, and read 2
  pi[,2] <- rbeta(ngen, ne[,4] + eta1*mu1, n.read - ne[,4] + eta1*(1-mu1))
  pi[,3] <- rbeta(ngen, ne[,5] + eta2*mu2, n.read - ne[,5] + eta2*(1-mu2))
  
  
  # initialize p0i and zi
  y <- ne[,2] #rowSums(n.sim[,2:3])
  temp <- appx.q.beta(q0, npl)
  pi[,1] <- rbeta(length(y), y + 10*temp[1], n.read-y+10*temp[2])
  i=1
  temp1 <- log(alpha0) + dbinom(y, n.read, p0, log=TRUE)
  temp2 <- log(1-alpha0) + dbinom(y, n.read, pi, log=TRUE)  # + dg
  pz <- 1 / (1 + exp(temp2 - temp1))
  zi <- rbinom(ngen, 1, pz)
  
  xx <- c(0:1e5)/1e5
  dg <- (1-q0)^(xx*npl)
  # pi=sample(xx, size=ngen, replace=TRUE, prob=dg)
  pi[zi == 1,1] <- sample(xx, size=sum(zi==1), replace=TRUE, prob=dg)
  temp <- appx.q.beta(q0, npl)
  a0 <- y + 1*temp[1]
  b0 <- n.read-y+1*temp[2]
  pi[zi==0,1] <- rbeta(sum(zi==0), a0[zi==0], b0[zi==0])
  
  
  # MCMC proposal parameters
  tau.pi <- c(1, 1, 1)
  tau.p0 <- 1e-4
  tau.q <- .005
  tau.mu.eta <- c(.001, .001, .001)
  
  
  # Prior parameters
  a.p0 <- npl*p0
  b.p0 <- npl*(1-p0)
  a.alpha <- npl*alpha0
  b.alpha <- npl*(1-alpha0)
  a.q <- npl*q0
  b.q <- npl*(1-q0)
  
  
  ll.data <- ll.compute(p0, pi, zi, ne)
  ll.latent <- compute.ll.latent(zi, pi, 
                                 p0, alpha0, q0, 
                                 mu1, eta1, mu2, eta2, 
                                 n.read)
  
  # Setup to sample
  params <- matrix(c(p0, alpha0, q0, mu1, eta1, mu2, eta2, sum(ll.data), sum(ll.data) + ll.latent), nrow=1)
  colnames(params) <- c("p0", "alpha", "q", "mu1", "eta1", "mu2", "eta2", "ll", "ll.w.latent")
  
  accept.p0 <- accept.q <- 0
  accept.mu <- accept.eta <- rep(0,3)
  accept.pi <- matrix(0, nrow(pi), 3)
  accept.pi1.z0 <- accept.pi1.z1 <- 0
  
  # Initialize adaptive proposals
  params.save <- matrix(0, nsd, length(params))
  params.save[1,] <- params
  params.save[2,] <- params
  
  
  # verbose <- TRUE
  
  # write.interval <- round(nmcmc/10)
  
  
  write.table(params, file=file.path(writedir, "params.txt"), col.names=TRUE, 
              row.names=FALSE, append=FALSE)
  temp <- matrix(c(accept.p0, accept.q, accept.mu[2:3]), nrow=1)
  colnames(temp) <- c("accept.p0", "accept.q", "accept.mu1", "accept.mu2")
  write.table(temp, file=file.path(writedir, "accept_rates.txt"),col.names=TRUE, 
              row.names=FALSE, append=FALSE)
  
  if(verbose){
    cat("Chain initialized \n")
  }
  
  if(save.latent.params){
    write.table(accept.pi, file=file.path(writedir, "accept_pi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
    write.table(zi, file=file.path(writedir, "zi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
    write.table(pi, file=file.path(writedir, "pi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
  }
  
  
  
  for(i in 2:nmcmc){
    #####
    ##  Update p0 directly
    if(sum(zi ==1)>0){
      p0.old <- params[1]
      
      # adaptive proposal
      if(i <= nsd + 1){
        tau.prop <- tau.p0
      } else {
        tau.prop <- sd(qnorm(params.save[,1]))
        tau.prop <- max(tau.prop, .Machine$double.eps)
      }
      # tau.prop <- tau.p0
      p0.prop <- pnorm(rnorm(1, qnorm(p0.old), tau.prop))
      d.prop <- 1/2*(qnorm(p0.old)^2 - qnorm(p0.prop)^2)
      
      # p0.prop <- rbeta(1, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2]) + 1)
      # d.prop <- dbeta(p0.old, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2])  + 1, log=TRUE) -  dbeta(p0.prop, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2])  + 1, log=TRUE)
      
      
      # Compute likelihoods, starting with priors (FLAT PRIOR)
      dnew <- dbeta(p0.prop, a.p0, b.p0, log=TRUE)
      dold <- dbeta(p0.old, a.p0, b.p0, log=TRUE)
      
      ll.new <- ll.compute(p0.prop, pi[zi==1,], rep(1, sum(zi==1)), ne[zi==1,])
      
      alpha <- sum(ll.new - ll.data[zi==1]) + dnew - dold + d.prop
      if(log(runif(1,0,1)) < alpha){
        accept <- TRUE
        p.out <- p0.prop
        ll.data[zi==1] <- ll.new   # save data IF accepted sample
      } else {
        accept <- FALSE
        p.out <- p0.old
      }
      params[1] <- p.out
      accept.p0 <- accept.p0*(i-1)/i + accept/i
      
    } else {
      accept.p0 <- accept.p0*(i-1)/i 
    }
    
    
    #####
    ## Udpate pi[,1]
    y <- ne[,2] #rowSums(n.sim[,2:3])
    p <- params[1]
    alpha <- params[2]
    mu0 <- q <- params[3]
    
    temp <- appx.q.beta(mu0, npl)
    a0 <- y + tau.pi[1]*temp[1]
    # b0 <- (npl - y) + tau.pi[1]*temp[2]
    b0 <- (n.read - y) + tau.pi[1]*temp[2]
    pi.prop <- pi
    pi.prop[zi==0,1] <- rbeta(sum(zi==0), a0[zi==0], b0[zi==0])
    pi.prop[zi==1,1] <- rbeta(sum(zi==1), temp[1], temp[2])
    
    d.prop.new <- d.prop.old <- rep(0, nrow(pi))
    d.prop.new[zi==0] <- dbeta(pi.prop[zi==0,1], a0[zi==0], b0[zi==0], log=TRUE)
    d.prop.old[zi==0] <- dbeta(pi[zi==0,1], a0[zi==0], b0[zi==0], log=TRUE)
    d.prop.new[zi==1] <- dbeta(pi.prop[zi==1,1], temp[1], temp[2], log=TRUE)
    d.prop.old[zi==1] <- dbeta(pi[zi==1,1], temp[1], temp[2], log=TRUE)
    
    # prior.new <- npl*pi.prop[,1]*log(1-mu0) # geom prior
    # prior.old <- npl*pi[,1]*log(1-mu0) # geom prior
    prior.new <- n.read*pi.prop[,1]*log(1-mu0) # geom prior
    prior.old <- n.read*pi[,1]*log(1-mu0) # geom prior
    
    # log.phi.new <- t(log(apply(pi.prop[zi==0,], 1, compute.phi.direct)))
    # j.new <- t(apply(pi.prop[zi==0,], 1, jacobian.direct, log=TRUE))
    # j.old <- t(apply(pi, 1, jacobian.direct, log=TRUE))
    # j.0 <- t(apply(pi0, 1, jacobian.direct, log=TRUE))
    
    ll.new <- 0*ll.data  # apply(n.sim*log.phi.new, 1, sum) + j.new
    ll.old <- ll.data #apply(n.sim*log.phi.old, 1, sum) + j.old
    ll.new[zi==0] <- ll.compute(params[1], pi.prop[zi==0,], rep(0, sum(zi==0)), ne[zi==0,])
    # rowSums(n.sim[zi==0,]*log.phi.new) + j.new
    ll.old[zi==1] <- 0
    
    # Accept/reject
    a <- ll.new - ll.old + prior.new - prior.old - d.prop.new + d.prop.old
    accept.vec <- log(runif(length(a),0,1)) < a
    pi[accept.vec,1] <- pi.prop[accept.vec,1]   # update p0i
    i.new <- as.logical(accept.vec*c(zi==0))   
    ll.data[which(i.new)] <- ll.new[which(i.new)]   # update ll.data, only for high dispersion process
    
    # Acceptance rates
    accept.pi[,1] <- accept.pi[,1]*(i-1)/i + accept.vec/i
    accept.pi1.z0 <- accept.pi1.z0*(i-1)/i + mean(accept.vec[zi==0])/i
    accept.pi1.z1 <- accept.pi1.z1*(i-1)/i + mean(accept.vec[zi==1])/i
    
    
    
    #####
    ## Update q
    x <- pi[,1]
    q <- params[3]
    # adaptive proposal for q
    if(i <= nsd + 1){
      tau.prop <- tau.q
    } else {
      tau.prop <- sd(qnorm(params.save[,3]))
      tau.prop <- max(tau.prop, .Machine$double.eps)
    }
    q.prop <- pnorm(rnorm(1, qnorm(q), tau.prop))
    d.prop <- 1/2*(qnorm(q)^2 - qnorm(q.prop)^2)
    # ll.new <- npl*sum(x)*log(1-q.prop) +  length(x)*log(npl*log(1-q.prop) / ((1-q.prop)^npl - 1))
    # ll.old <- npl*sum(x)*log(1-q) +  length(x)*log(npl*log(1-q) / ((1-q)^npl - 1))
    ll.new <- sum(n.read*x*log(1-q.prop) +  log(n.read*log(1-q.prop) / ((1-q.prop)^n.read - 1)))
    ll.old <- sum(n.read*x*log(1-q) +  log(n.read*log(1-q) / ((1-q)^n.read - 1)))
    a <- ll.new - ll.old + dbeta(q.prop, a.q, b.q, log=TRUE) - dbeta(q, a.q, b.q, log=TRUE) + d.prop
    accept <- log(runif(1,0,1)) < a
    if(accept){
      params[3] <- q.prop
    } else {
      params[3] <- q
    }
    accept.q <- accept.q*(i-1)/i + accept/i
    
    
    #####
    ## Sample zi
    ll.z0 <- ll.z1 <- ll.data
    ll.z0[zi==1] <- ll.compute(params[1], pi[zi==1,], rep(0, sum(zi==1)), ne[zi==1,])
    ll.z1[zi==0] <- ll.compute(params[1], pi[zi==0,], rep(1, sum(zi==0)), ne[zi==0,])
    
    temp1 <- log(params[2]) + ll.z1
    temp2 <- log(1-params[2]) + ll.z0
    
    # pz <- exp(temp1) / (exp(temp1) + exp(temp2))
    pz <- 1 / (1 + exp(temp2-temp1))
    zi <- rbinom(length(zi), 1, pz)
    
    # update ll
    ll.data[zi==1] <- ll.z1[zi==1]
    ll.data[zi==0] <- ll.z0[zi==0]
    
    
    #####
    ## sample alpha
    params[2] <- update.alpha(zi, a.alpha, b.alpha)
    
    
    #####
    ## Udpate pi[,2] and pi[,3] by column
    for(j in sample(2:3)){
      mu0 <- params[4 + 2*(j-2)]
      nu0 <- params[5 + 2*(j-2)]
      
      x <- ne[,2 + j]
      
      # pi conjugacy for proposal
      a0 = mu0*nu0 + x*tau.pi[j]
      b0 <- nu0*(1-mu0) + (n.read - x)*tau.pi[j]
      pi.prop <- rbeta(length(x), a0, b0)
      pi.old <- pi[,j]
      
      d.prop.new <- dbeta(pi.prop, a0, b0, log=TRUE)
      d.prop.old <- dbeta(pi.old, a0, b0, log=TRUE)
      
      prior.new <- dbeta(pi.prop, mu0*nu0, nu0*(1-mu0), log=TRUE)
      prior.old <- dbeta(pi.old, mu0*nu0, nu0*(1-mu0), log=TRUE)
      
      p.prop <- pi
      p.prop[,j] <- pi.prop
      ll.new <- ll.compute(params[1], p.prop, zi, ne)
      
      alpha <- ll.new - ll.data + prior.new - prior.old + d.prop.old - d.prop.new
      accept.vec <- log(runif(length(alpha))) < alpha
      
      pi[accept.vec,j] <- pi.prop[accept.vec]
      ll.data[accept.vec] <- ll.new[accept.vec]    # update ll for accepted entries
      
      accept.pi[,j] <- accept.pi[,j]*(i-1)/i + accept.vec/i
    }
    
    
    
    # Update mu and eta jointly
    if(i <= nsd + 1){
      tau.prop <- tau.mu.eta
    } else {
      tau.prop <- c(sd(qnorm(params.save[,4])), sd(log(params.save[,5])))
      tau.prop[1] <- max(tau.prop[1], .Machine$double.eps)
      tau.prop[2] <- max(tau.prop[2], .Machine$double.eps)
    }
    res <- update.mu.eta.joint.jacob(pi[,2],  # data
                                     params[4],
                                     params[5], # old parameters
                                     a.mu1, b.mu1,  # beta prior for mu
                                     a.eta1, b.eta1,  # gamma prior for eta
                                     tau.prop)
    accept.mu[2] <- accept.mu[2]*(i-1)/i + res$accept/i
    params[4] <- res$mu
    params[5] <- res$eta
    
    
    # Update mu and eta jointly
    if(i <= nsd + 1){
      tau.prop <- tau.mu.eta
    } else {
      tau.prop <- c(sd(qnorm(params.save[,6])), sd(log(params.save[,7])))
      tau.prop[1] <- max(tau.prop[1], .Machine$double.eps)
      tau.prop[2] <- max(tau.prop[2], .Machine$double.eps)
    }
    res <- update.mu.eta.joint.jacob(pi[,3],  # data
                                     params[6],
                                     params[7], # old parameters
                                     a.mu2, b.mu2,  # beta prior for mu
                                     a.eta2, b.eta2,  # gamma prior for eta
                                     tau.prop)
    accept.mu[3] <- accept.mu[3]*(i-1)/i + res$accept/i
    params[6] <- res$mu
    params[7] <- res$eta
    
    params[8] <- sum(ll.data)
    ll.latent <- compute.ll.latent(zi, pi, 
                                   params[1], params[2], params[3], 
                                   params[4], params[5], 
                                   params[6], params[7], 
                                   n.read)
    params[9] <- params[8] + (ll.latent)
    
    # Save params for adaptive proposals
    params.save[i %% nsd + 1, ] <- params
    
    
    
    if(verbose & (i %% print.interval == 0)){
      cat("###################################################\n")
      cat("iteration: ", i, "of", nmcmc, "\n")
      cat("accept p0: ", round(accept.p0,4), "\n")
      cat("accept q: ", round(accept.q,4), "\n")
      cat("accept mus: ", round(accept.mu,4), "\n")
      cat("min accept pi: ", round(apply(accept.pi, 2, min),4), "\n")
      cat("Avg. accept pi: ", round(apply(accept.pi, 2, mean),4), "\n")
      cat("Fraction zero accept pi: ", round(apply(accept.pi==0, 2, mean),4), "\n")
      cat("Average accept pi[,1] by zi==0:", round(accept.pi1.z0, 3), "and zi==1: ", round(accept.pi1.z1, 3), "\n")
      cat("###################################################\n")
    }
    
    writeout <- i %% thin == 0 
    if(!save.burn){
      writeout <- i %% thin == 0 & i > nburn
    }
    
    if(writeout){
      write.table(params, file=file.path(writedir, "params.txt"), col.names = FALSE, 
                  row.names=FALSE, append=TRUE)
      temp <- matrix(c(accept.p0, accept.q, accept.mu[2:3]), nrow=1)
      colnames(temp) <- c("accept.p0", "accept.q", "accept.mu1", "accept.mu2")
      write.table(temp, file=file.path(writedir, "accept_rates.txt"), col.names=FALSE, 
                  row.names=FALSE, append=TRUE)
      
      if(save.latent.params){
        write.table(accept.pi, file=file.path(writedir, "accept_pi.txt"), col.names=FALSE,  
                    row.names=FALSE, append=FALSE)
        write.table(zi, file=file.path(writedir, "zi.txt"), col.names=FALSE, 
                    row.names=FALSE, append=FALSE)
        write.table(pi, file=file.path(writedir, "pi.txt"), col.names=FALSE,  
                    row.names=FALSE, append=FALSE)
      }
    }
  }
  
  return()
}



################################################################################
#
#  Bayesian error model
#

# MCMC for Bayesian model estimation, including Q-scores
mcmc.error.qscore <- function(ne,  # error counts
                             writedir=getwd(),  # where to write results to
                             nmcmc=1e4,   # number of mcmc samples
                             nburn=round(nmcmc*3/4),    # number of samples to omit
                             thin=10,    # save interval
                             save.latent.params=FALSE,
                             save.burn=TRUE,   # should the burned iterations be saved?
                             verbose=FALSE,   # should progress be printed out?
                             print.interval=round(nmcmc/10)  # how often should progress be printed out?

)
{
  joint.prop=FALSE
  init=NULL
  
  ngen <- nrow(ne)
  n.read <- rowSums(ne[, 1:5])
  npl <- my.mode(n.read)
  if(npl != min(n.read)){
    warning(paste0("Minimal rowsum ", min(n.read), " is not equal current payload of ", npl))
  }
  
  
  # joint.prop <- TRUE
  q1.i <-  scale(qnorm(ne[,6]/npl))
  q2.i <-  scale(qnorm(ne[,7]/npl))
  
  
  nsd <- round(nmcmc/20)
  nn <- nrow(ne)
  
  
  ne <- as.matrix(ne)
  ne[is.na(ne)] <- 0  
  if(ncol(ne) != 7){
    stop("wrong input size of ne")
  }
  ne <- ne[,1:5]
  
  
  if(is.null(init)){
    e0 <- ne[,2]  #rowSums(ne.sub[,2:3]) 
    e0 <- e0[e0<npl]
    
    p0 <- mean(e0[e0 < npl/10] / npl)
    p0 <- phat0 <- mean(p0, 1 - mean(e0 == 0)^(1/npl))
    
    x = round(npl/10):round(npl*2/3)
    mini.distr <- sapply(x, function(z) sum(e0 == z))
    mini.distr <- mini.distr/sum(mini.distr)
    y = log(mini.distr)
    q0 <- 1 - exp(coef(lm(y[!is.infinite(y)] ~ x[!is.infinite(y)]))[2])
    
    xx <- c(0:1e5)/1e5
    dg <- (1-q0)^c(xx*npl) # / pgeom(npl, q.hat)
    y.temp <- sample(xx, nn, replace=TRUE, prob=dg)
    qbar <- mean(y.temp)
    
    
    alpha.temp <- (mean(e0/npl) - qbar) / (p0 - qbar)  # match expectation
    alpha.temp <- min(alpha.temp, .9999)
    
    # cat(c(qbar, q0, p0, alpha.temp))
    
    alphas <- 10^seq(log10(alpha.temp) - .1, min(log10(alpha.temp) + .1, log10(.99999)), length.out=100)
    ll <- 0*alphas
    for(i in 1:length(alphas)){
      dg <- dgeom(e0, q0) / pgeom(npl, q0)
      ll[i] <- sum(log(alphas[i]*dbinom(e0, npl, p0) + (1-alphas[i])*dg))
    }
    # plot(alphas, ll)
    alpha0 <- alphas[which.max(ll)]
    
    ## MLE step to really hone in on p0, q0, alpha0
    params0 <- c(p0, alpha0, q0)
    res <- optim(qnorm(params0), neg.ll.mix, x=e0, npl=npl)
    if(res$convergence == 0){
      params0 <- pnorm(res$par)
    }
    p0 <- params0[1]
    alpha0 <- params0[2]
    q0 <- params0[3]
    
    
    
    
    # Method of moments -> MLE estimates
    #    start with zero slope
    x=ne[,4]
    mu1 <-mean(x/npl)
    b1.init <- qnorm(mu1)
    a1.init <- 0
    eta1 <- mu1*(1-mu1) / var(x/npl) - 1
    res = optim(c(a1.init, b1.init, log(eta1)),
                beta.opt.slope.v2,
                x=x, q=q1.i, npl=n.read)
    if(res$convergence == 0){
      a1.init <- res$par[1]
      b1.init <- res$par[2]
      eta1 <- exp(res$par[3])
    }
    
    # Method of moments -> MLE estimates
    #    start with zero slope
    x=ne[,5]
    mu2 <-mean(x/npl)
    b2.init <- qnorm(mu1)
    a2.init <- 0
    eta2 <- mu2*(1-mu2) / var(x/npl) - 1
    res = optim(c(a2.init, b2.init, log(eta2)),
                beta.opt.slope.v2,
                x=x, q=q2.i, npl=n.read)
    if(res$convergence == 0){
      a2.init <- res$par[1]
      b2.init <- res$par[2]
      eta2 <- exp(res$par[3])
    }
    
    mu.a1 <- a1 <- a1.init
    mu.b1 <- b1 <- b1.init
    sd.a1 <- sd.b1 <- 1e3
    a.eta1 <- eta1
    b.eta1 <- 1
    
    mu.a2 <- a2 <- a2.init
    mu.b2 <- b2 <- b2.init
    sd.a2 <- sd.b2 <- 1e3
    a.eta2 <- eta2
    b.eta2 <- 1
    
    
    ## Initialize p1i and p2i
    mu1 <- pnorm(a1.init*q1.i + b1.init)
    mu2 <- pnorm(a2.init*q2.i + b2.init)
    pi <- matrix(0, ngen, 3)  # latent probability of high dispersion errors for synthesis, read 1, and read 2
    pi[,2] <- rbeta(ngen, ne[,4] + eta1*mu1, n.read - ne[,4] + eta1*(1-mu1))
    pi[,3] <- rbeta(ngen, ne[,5] + eta2*mu2, n.read - ne[,5] + eta2*(1-mu2))
    
    
    # initialize p0i and zi
    y <- ne[,2] #rowSums(n.sim[,2:3])
    temp <- appx.q.beta(q0, npl)
    pi[,1] <- rbeta(length(y), y + 10*temp[1], n.read-y+10*temp[2])
    i=1
    temp1 <- log(alpha0) + dbinom(y, n.read, p0, log=TRUE)
    temp2 <- log(1-alpha0) + dbinom(y, n.read, pi, log=TRUE)  # + dg
    pz <- 1 / (1 + exp(temp2 - temp1))
    zi <- rbinom(ngen, 1, pz)
    
    xx <- c(0:1e5)/1e5
    dg <- (1-q0)^(xx*npl)
    # pi=sample(xx, size=ngen, replace=TRUE, prob=dg)
    pi[zi == 1,1] <- sample(xx, size=sum(zi==1), replace=TRUE, prob=dg)
    temp <- appx.q.beta(q0, npl)
    a0 <- y + 1*temp[1]
    b0 <- n.read-y+1*temp[2]
    pi[zi==0,1] <- rbeta(sum(zi==0), a0[zi==0], b0[zi==0])
  
    
    ll.data <- ll.compute(p0, pi, zi, ne)
    mu1 <- pnorm(a1*q1.i + b1)
    mu2 <- pnorm(a2*q2.i + b2)
    ll.latent <- compute.ll.latent(zi, pi, 
                                   p0, alpha0, q0, 
                                   mu1, eta1, mu2, eta2, 
                                   n.read)
    # ll.latent <- ll.latent + apply(pi, 2, )
    params0 <- c(p0, alpha0, q0, a1, b1, eta1, a2, b2, eta2, sum(ll.data), 
                 sum(ll.data)+(ll.latent))
  } else {
    params.1 <- init[1,]
    params0 <- init[nrow(init),]
    p0 <- c(unlist(params.1[1]))
    alpha0 <- c(unlist(params.1[2]))
    q0 <- c(unlist(params.1[3]))
    
    mu.a1 <- a1 <- c(unlist(params.1[4]))
    mu.b1 <- b1 <- c(unlist(params.1[5]))
    sd.a1 <- sd.b1 <- 1e3
    a.eta1 <- eta1 <- c(unlist(params.1[6]))
    b.eta1 <- 1
    
    mu.a2 <- a2 <- c(unlist(params.1[7]))
    mu.b2 <- b2 <- c(unlist(params.1[8]))
    sd.a2 <- sd.b2 <- 1e3
    a.eta2 <- eta2 <- c(unlist(params.1[9]))
    b.eta2 <- 1
    
    mu1 <- pnorm(a1*q1.i + b1)
    mu2 <- pnorm(a2*q2.i + b2)
    pi <- matrix(0, ngen, 3)  # latent probability of high dispersion errors for synthesis, read 1, and read 2
    pi[,2] <- rbeta(ngen, ne[,4] + eta1*mu1, n.read - ne[,4] + eta1*(1-mu1))
    pi[,3] <- rbeta(ngen, ne[,5] + eta2*mu2, n.read - ne[,5] + eta2*(1-mu2))
    
    
    # initialize p0i and zi
    y <- ne[,2] #rowSums(n.sim[,2:3])
    temp <- appx.q.beta(q0, npl)
    
    temp1 <- log(alpha0) + dbinom(y, n.read, p0, log=TRUE)
    temp2 <- log(1-alpha0) + dbinom(y, n.read, pi, log=TRUE)  # + dg
    pz <- 1 / (1 + exp(temp2 - temp1))
    zi <- rbinom(ngen, 1, pz)
    
    xx <- c(0:1e5)/1e5
    dg <- (1-q0)^(xx*npl)
    # pi=sample(xx, size=ngen, replace=TRUE, prob=dg)
    pi[zi == 1,1] <- sample(xx, size=sum(zi==1), replace=TRUE, prob=dg)
    temp <- appx.q.beta(q0, npl)
    a0 <- y + 1*temp[1]
    b0 <- n.read-y+1*temp[2]
    pi[zi==0,1] <- rbeta(sum(zi==0), a0[zi==0], b0[zi==0])
    
    
    ll.data <- ll.compute(p0, pi, zi, ne)
    mu1 <- pnorm(a1*q1.i + b1)
    mu2 <- pnorm(a2*q2.i + b2)
    ll.latent <- compute.ll.latent(zi, pi, 
                                   p0, alpha0, q0, 
                                   mu1, eta1, mu2, eta2, 
                                   n.read)
    
    params0[10:11] <- c(sum(ll.data), sum(ll.data)+(ll.latent))
  }

  params0 <- c(unlist(params0))
  
  # Setup to sample
  # MCMC proposal parameters
  tau.pi <- c(1, 1, 1)
  tau.p0 <- 1e-4
  tau.q <- .005
  tau.mu.eta <- c(.001, .001)
  
  
  # Prior parameters
  a.p0 <- npl*p0
  b.p0 <- npl*(1-p0)
  a.alpha <- npl*alpha0
  b.alpha <- npl*(1-alpha0)
  a.q <- npl*q0
  b.q <- npl*(1-q0)
  
  params <- matrix(params0, nrow=1)
  colnames(params) <- c("p0", "alpha", "q", 
                     "a1", "b1", "eta1", 
                     "a2", "b2", "eta2", "ll", "ll.w.latent")
  
  accept.p0 <- accept.q <- 0
  accept.mu <- accept.eta <- rep(0,3)
  accept.pi <- matrix(0, nrow(pi), 3)
  accept.pi1.z0 <- accept.pi1.z1 <- 0
  
  
  
  # Initialize adaptive proposals
  params.save <- matrix(0, nsd, length(params))
  params.save[1,] <- params0
  params.save[2,] <- params0
  colnames(params.save) <- c("p0", "alpha", "q", 
                        "a1", "b1", "eta1", 
                        "a2", "b2", "eta2", "ll", "ll.w.latent")
  
  # verbose <- TRUE
  
  # write.interval <- round(nmcmc/10)
  
  
  write.table(params, file=file.path(writedir, "params.txt"), col.names=TRUE, 
              row.names=FALSE, append=FALSE)
  temp <- matrix(c(accept.p0, accept.q, accept.mu[2:3]), nrow=1)
  colnames(temp) <- c("accept.p0", "accept.q", "accept.mu1", "accept.mu2")
  write.table(temp, file=file.path(writedir, "accept_rates.txt"),col.names=TRUE, 
              row.names=FALSE, append=FALSE)
  
  if(verbose){
    cat("Chain initialized \n")
  }
  
  if(save.latent.params){
    write.table(accept.pi, file=file.path(writedir, "accept_pi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
    write.table(zi, file=file.path(writedir, "zi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
    write.table(pi, file=file.path(writedir, "pi.txt"), col.names=FALSE,  
                row.names=FALSE, append=FALSE)
  }
  
  
  
  for(i in 2:nmcmc){
    #####
    ##  Update p0 directly
    if(sum(zi ==1)>0){
      p0.old <- params[1]
      
      # adaptive proposal
      if(i <= nsd + 1){
        tau.prop <- tau.p0
      } else {
        tau.prop <- sd(qnorm(params.save[,1]))
        tau.prop <- max(tau.prop, .Machine$double.eps)
      }
      # tau.prop <- tau.p0
      p0.prop <- pnorm(rnorm(1, qnorm(p0.old), tau.prop))
      d.prop <- 1/2*(qnorm(p0.old)^2 - qnorm(p0.prop)^2)
      
      # p0.prop <- rbeta(1, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2]) + 1)
      # d.prop <- dbeta(p0.old, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2])  + 1, log=TRUE) -  dbeta(p0.prop, sum(n.sim[zi==1,2]) + p0, sum(zi==1)*npl - sum(n.sim[zi==1,2])  + 1, log=TRUE)
      
      
      # Compute likelihoods, starting with priors (FLAT PRIOR)
      dnew <- dbeta(p0.prop, a.p0, b.p0, log=TRUE)
      dold <- dbeta(p0.old, a.p0, b.p0, log=TRUE)
      
      ll.new <- ll.compute(p0.prop, pi[zi==1,], rep(1, sum(zi==1)), ne[zi==1,])
      
      alpha <- sum(ll.new - ll.data[zi==1]) + dnew - dold + d.prop
      if(log(runif(1,0,1)) < alpha){
        accept <- TRUE
        p.out <- p0.prop
        ll.data[zi==1] <- ll.new   # save data IF accepted sample
      } else {
        accept <- FALSE
        p.out <- p0.old
      }
      params[1] <- p.out
      accept.p0 <- accept.p0*(i-1)/i + accept/i
      
    } else {
      accept.p0 <- accept.p0*(i-1)/i 
    }
    
    
    #####
    ## Udpate pi[,1]
    y <- ne[,2] #rowSums(n.sim[,2:3])
    p <- params[1]
    alpha <- params[2]
    mu0 <- q <- params[3]
    
    temp <- appx.q.beta(mu0, npl)
    a0 <- y + tau.pi[1]*temp[1]
    # b0 <- (npl - y) + tau.pi[1]*temp[2]
    b0 <- (n.read - y) + tau.pi[1]*temp[2]
    pi.prop <- pi
    pi.prop[zi==0,1] <- rbeta(sum(zi==0), a0[zi==0], b0[zi==0])
    pi.prop[zi==1,1] <- rbeta(sum(zi==1), temp[1], temp[2])
    
    d.prop.new <- d.prop.old <- rep(0, nrow(pi))
    d.prop.new[zi==0] <- dbeta(pi.prop[zi==0,1], a0[zi==0], b0[zi==0], log=TRUE)
    d.prop.old[zi==0] <- dbeta(pi[zi==0,1], a0[zi==0], b0[zi==0], log=TRUE)
    d.prop.new[zi==1] <- dbeta(pi.prop[zi==1,1], temp[1], temp[2], log=TRUE)
    d.prop.old[zi==1] <- dbeta(pi[zi==1,1], temp[1], temp[2], log=TRUE)
    
    # prior.new <- npl*pi.prop[,1]*log(1-mu0) # geom prior
    # prior.old <- npl*pi[,1]*log(1-mu0) # geom prior
    prior.new <- n.read*pi.prop[,1]*log(1-mu0) # geom prior
    prior.old <- n.read*pi[,1]*log(1-mu0) # geom prior
    
    # log.phi.new <- t(log(apply(pi.prop[zi==0,], 1, compute.phi.direct)))
    # j.new <- t(apply(pi.prop[zi==0,], 1, jacobian.direct, log=TRUE))
    # j.old <- t(apply(pi, 1, jacobian.direct, log=TRUE))
    # j.0 <- t(apply(pi0, 1, jacobian.direct, log=TRUE))
    
    ll.new <- 0*ll.data  # apply(n.sim*log.phi.new, 1, sum) + j.new
    ll.old <- ll.data #apply(n.sim*log.phi.old, 1, sum) + j.old
    ll.new[zi==0] <- ll.compute(params[1], pi.prop[zi==0,], rep(0, sum(zi==0)), ne[zi==0,])
    # rowSums(n.sim[zi==0,]*log.phi.new) + j.new
    ll.old[zi==1] <- 0
    
    # Accept/reject
    a <- ll.new - ll.old + prior.new - prior.old - d.prop.new + d.prop.old
    accept.vec <- log(runif(length(a),0,1)) < a
    pi[accept.vec,1] <- pi.prop[accept.vec,1]   # update p0i
    i.new <- as.logical(accept.vec*c(zi==0))   
    ll.data[which(i.new)] <- ll.new[which(i.new)]   # update ll.data, only for high dispersion process
    
    # Acceptance rates
    accept.pi[,1] <- accept.pi[,1]*(i-1)/i + accept.vec/i
    accept.pi1.z0 <- accept.pi1.z0*(i-1)/i + mean(accept.vec[zi==0])/i
    accept.pi1.z1 <- accept.pi1.z1*(i-1)/i + mean(accept.vec[zi==1])/i
    
    
    
    #####
    ## Update q
    x <- pi[,1]
    q <- params[3]
    # adaptive proposal for q
    if(i <= nsd + 1){
      tau.prop <- tau.q
    } else {
      tau.prop <- sd(qnorm(params.save[,3]))
      tau.prop <- max(tau.prop, .Machine$double.eps)
    }
    q.prop <- pnorm(rnorm(1, qnorm(q), tau.prop))
    d.prop <- 1/2*(qnorm(q)^2 - qnorm(q.prop)^2)
    # ll.new <- npl*sum(x)*log(1-q.prop) +  length(x)*log(npl*log(1-q.prop) / ((1-q.prop)^npl - 1))
    # ll.old <- npl*sum(x)*log(1-q) +  length(x)*log(npl*log(1-q) / ((1-q)^npl - 1))
    ll.new <- sum(n.read*x*log(1-q.prop) +  log(n.read*log(1-q.prop) / ((1-q.prop)^n.read - 1)))
    ll.old <- sum(n.read*x*log(1-q) +  log(n.read*log(1-q) / ((1-q)^n.read - 1)))
    a <- ll.new - ll.old + dbeta(q.prop, a.q, b.q, log=TRUE) - dbeta(q, a.q, b.q, log=TRUE) + d.prop
    accept <- log(runif(1,0,1)) < a
    if(accept){
      params[3] <- q.prop
    } else {
      params[3] <- q
    }
    accept.q <- accept.q*(i-1)/i + accept/i
    
    
    #####
    ## Sample zi
    ll.z0 <- ll.z1 <- ll.data
    ll.z0[zi==1] <- ll.compute(params[1], pi[zi==1,], rep(0, sum(zi==1)), ne[zi==1,])
    ll.z1[zi==0] <- ll.compute(params[1], pi[zi==0,], rep(1, sum(zi==0)), ne[zi==0,])
    
    temp1 <- log(params[2]) + ll.z1
    temp2 <- log(1-params[2]) + ll.z0
    
    # pz <- exp(temp1) / (exp(temp1) + exp(temp2))
    pz <- 1 / (1 + exp(temp2-temp1))
    zi <- rbinom(length(zi), 1, pz)
    
    # update ll
    ll.data[zi==1] <- ll.z1[zi==1]
    ll.data[zi==0] <- ll.z0[zi==0]
    
    
    #####
    ## sample alpha
    params[2] <- update.alpha(zi, a.alpha, b.alpha)
    
    
    #####
    ## Udpate pi[,2] and pi[,3] by column
    for(j in sample(2:3)){
      # now this is slope-determined mu0
      a0 <- params[4 + 3*(j-2)]
      b0 <- params[5 + 3*(j-2)]
      if(j == 2){
        q.temp <- q1.i
      } else {
        q.temp <- q2.i
      }
      mu0 <- pnorm(a0*q.temp + b0)
      nu0 <- params[6 + 3*(j-2)]
      
      x <- ne[,2 + j]
      
      # pi conjugacy for proposal
      a0 = mu0*nu0 + x*tau.pi[j]
      b0 <- nu0*(1-mu0) + (n.read - x)*tau.pi[j]
      pi.prop <- rbeta(length(x), a0, b0)
      pi.old <- pi[,j]
      
      d.prop.new <- dbeta(pi.prop, a0, b0, log=TRUE)
      d.prop.old <- dbeta(pi.old, a0, b0, log=TRUE)
      
      prior.new <- dbeta(pi.prop, mu0*nu0, nu0*(1-mu0), log=TRUE)
      prior.old <- dbeta(pi.old, mu0*nu0, nu0*(1-mu0), log=TRUE)
      
      p.prop <- pi
      p.prop[,j] <- pi.prop
      ll.new <- ll.compute(params[1], p.prop, zi, ne)
      
      alpha <- ll.new - ll.data + prior.new - prior.old + d.prop.old - d.prop.new
      accept.vec <- log(runif(length(alpha))) < alpha
      
      pi[accept.vec,j] <- pi.prop[accept.vec]
      ll.data[accept.vec] <- ll.new[accept.vec]    # update ll for accepted entries
      
      accept.pi[,j] <- accept.pi[,j]*(i-1)/i + accept.vec/i
    }
    
    
    
    # Update mu and eta jointly
    if(i <= nsd + 1){
      tau.prop <- tau.mu.eta
      joint.prop.temp <- FALSE
      sigma.temp <- diag(3)
    } else {
      if(joint.prop){
        sigma.temp <- cov(cbind(params.save[,4:5], params.save[,6]))
        tau.prop <- rep(1, 3)
      } else {
        tau.prop <- c(sd(params.save[,4]), sd(params.save[,5]), sd(log(params.save[,6])))/2
        tau.prop[1] <- max(tau.prop[1], .Machine$double.eps)
        tau.prop[2] <- max(tau.prop[2], .Machine$double.eps)
        tau.prop[3] <- max(tau.prop[3], .Machine$double.eps)
        sigma.temp <- diag(3)
      }
      
      joint.prop.temp <- joint.prop
    }

    # Update a, b, and eta jointly for beta distribution
    res <- update.mu.eta.slope(pi[,2], q1.i, # data, p_i for appropriate dimension
                                    params[4], params[5], params[6], # old parameters
                                    mu.a1, sd.a1,  # normal prior for a
                                    mu.b1, sd.b1,  # normal prior for a
                                    a.eta1, b.eta1,  # gamma prior for eta
                                    tau.prop, # 3-vector of independent proposal sds
                                    sigma=sigma.temp,   # optional multivariate
                                    joint.prop=joint.prop.temp   # joint proposal??
    ) 
    accept.mu[2] <- accept.mu[2]*(i-1)/i + res$accept/i
    params[4] <- res$a
    params[5] <- res$b
    params[6] <- res$eta
    
    
    # Update mu and eta jointly
    if(i <= nsd + 1){
      tau.prop <- tau.mu.eta
      joint.prop.temp <- FALSE
      sigma.temp <- diag(3)
    } else {
      if(joint.prop){
        sigma.temp <- cov(cbind(params.save[,4:5], params.save[,6]))
        tau.prop <- rep(1, 3)
      } else {
        tau.prop <- c(sd(params.save[,4]), sd(params.save[,5]), sd(log(params.save[,6])))/2
        tau.prop[1] <- max(tau.prop[1], .Machine$double.eps)
        tau.prop[2] <- max(tau.prop[2], .Machine$double.eps)
        tau.prop[3] <- max(tau.prop[3], .Machine$double.eps)
        sigma.temp <- diag(3)
      }
      
      joint.prop.temp <- joint.prop
    }
    
    # Update a, b, and eta jointly for beta distribution
    res <- update.mu.eta.slope(pi[,3], q2.i, # data, p_i for appropriate dimension
                               params[7], params[8], params[9], # old parameters
                               mu.a2, sd.a2,  # normal prior for a
                               mu.b2, sd.b2,  # normal prior for b
                               a.eta2, b.eta2,  # gamma prior for eta
                               tau.prop, # 3-vector of independent proposal sds
                               sigma=sigma.temp,   # optional multivariate
                               joint.prop=joint.prop.temp   # joint proposal??
    ) 
    accept.mu[3] <- accept.mu[3]*(i-1)/i + res$accept/i
    params[7] <- res$a
    params[8] <- res$b
    params[9] <- res$eta
    
    params[10] <- sum(ll.data)
    
    # ll.latent <- dbeta(zi)
    ll.latent <- compute.ll.latent(zi, pi, 
                                   params[1], params[2], params[3], 
                                   pnorm(params[4]*q1.i + params[5]), params[6], 
                                   pnorm(params[7]*q2.i + params[8]), params[9], 
                                   n.read)
    params[11] <- params[10] + (ll.latent)
    
    # Save params for adaptive proposals
    params.save[i %% nsd + 1, ] <- params
    
    
    
    if(verbose & (i %% print.interval == 0)){
      cat("###################################################\n")
      cat("iteration: ", i, "of", nmcmc, "\n")
      cat("accept p0: ", round(accept.p0,4), "\n")
      cat("accept q: ", round(accept.q,4), "\n")
      cat("accept mus: ", round(accept.mu,4), "\n")
      cat("min accept pi: ", round(apply(accept.pi, 2, min),4), "\n")
      cat("Avg. accept pi: ", round(apply(accept.pi, 2, mean),4), "\n")
      cat("Fraction zero accept pi: ", round(apply(accept.pi==0, 2, mean),4), "\n")
      cat("Average accept pi[,1] by zi==0:", round(accept.pi1.z0, 3), "and zi==1: ", round(accept.pi1.z1, 3), "\n")
      cat("###################################################\n")
    }
    
    writeout <- i %% thin == 0 
    if(!save.burn){
      writeout <- i %% thin == 0 & i > nburn
    }
    
    if(writeout){
      write.table(params, file=file.path(writedir, "params.txt"), col.names = FALSE, 
                  row.names=FALSE, append=TRUE)
      temp <- matrix(c(accept.p0, accept.q, accept.mu[2:3]), nrow=1)
      colnames(temp) <- c("accept.p0", "accept.q", "accept.mu1", "accept.mu2")
      write.table(temp, file=file.path(writedir, "accept_rates.txt"), col.names=FALSE, 
                  row.names=FALSE, append=TRUE)
      
      if(save.latent.params){
        write.table(accept.pi, file=file.path(writedir, "accept_pi.txt"), col.names=FALSE,  
                    row.names=FALSE, append=FALSE)
        write.table(zi, file=file.path(writedir, "zi.txt"), col.names=FALSE, 
                    row.names=FALSE, append=FALSE)
        write.table(pi, file=file.path(writedir, "pi.txt"), col.names=FALSE,  
                    row.names=FALSE, append=FALSE)
      }
    }
  }
  
  return()
}



################################################################################
#
# Simulate data
#

simulate_errors = function(simnum=1, # which simulation settings, 1 or 2?
                           nsim=1e6, # number of oligos
                           npl=75,  # number of payload bases
                           include_q=TRUE,   # include q-score information in output
                           seed=100)   # random seed
{
  set.seed(seed)
  
  # Set parameters
  if(simnum==1){
    
    # npl <- 75
    p0 <- .01
    alpha0 <- .9
    q0 <- .1
    a1 <- .4
    b1 <- -3.2
    eta1 <- 200
    a2 <- .2
    b2 <- -2.8
    eta2 <- 100
    
    p1.i.true <- rbeta(nsim, .0015*500, 500*(1-.0015))
    q.temp = round(-10*log10(p1.i.true))
    q.temp[q.temp <= 1] <- 1
    q.temp[q.temp >= 40] <- 40
    p1.i <- 10^(-q.temp)/10*npl
    q1.i <- scale(qnorm(p1.i/npl))
    q1.true <- q.temp
    
    p2.i.true <- rbeta(nsim, .002*400, 400*(1-.002))
    q.temp = round(-10*log10(p2.i.true))
    q.temp[q.temp <= 1] <- 1
    q.temp[q.temp >= 40] <- 40
    p2.i <- 10^(-q.temp)/10*npl
    q2.i <- scale(qnorm(p2.i/npl))
    q2.true <- q.temp
    
    
  } else if (simnum == 2){
    # npl <- 75
    p0 <- .05
    alpha0 <- .91
    q0 <- .05
    a1 <- .6
    b1 <- -3
    eta1 <- 100
    a2 <- .3
    b2 <- -2.
    eta2 <- 60
    
    p1.i.true <- rbeta(nsim, .002*400, 400*(1-.002))
    q.temp = round(-10*log10(p1.i.true))
    q.temp[q.temp <= 1] <- 1
    q.temp[q.temp >= 40] <- 40
    p1.i <- 10^(-q.temp)/10*npl
    q1.i <- scale(qnorm(p1.i/npl))
    q1.true <- q.temp
    
    
    p2.i.true <- rbeta(nsim, .003*350, 350*(1-.003))
    q.temp = round(-10*log10(p2.i.true))
    q.temp[q.temp <= 1] <- 1
    q.temp[q.temp >= 40] <- 40
    p2.i <- 10^(-q.temp)/10*npl
    q2.i <- scale(qnorm(p2.i/npl))
    q2.true <- q.temp
    
  }
  
  # True error rates
  s1 <- mean(p0*alpha0 + (1-alpha0)*mu.q(q0, npl))
  r1 <- mean(pnorm(a1*q1.i + b1))
  r2 <- mean(pnorm(a2*q2.i + b2))
  
  # simulate data
  
  
  
  params.true <- c(p0, alpha0, q0, a1, b1, eta1, a2, b2, eta2)
  
  
  zi <- zi.true <- rbinom(nsim, 1, alpha0)
  # y <- rep(0, ngen)
  
  
  xx <- c(0:1e5)/1e5
  dg <- (1-q0)^(xx*npl)
  pi <- sample(xx, nsim, replace=TRUE, prob=dg)
  pi.true <- pi
  
  pi.sim <- matrix(NA, nsim, 3)
  pi.sim[,1] <- pi
  pi.sim[zi==1,1] <- p0
  
  mu1 <- pnorm(a1*q1.i + b1)
  mu2 <- pnorm(a2*q2.i + b2)
  pi.sim[,2] <- rbeta(nsim, eta1*mu1, eta1*(1-mu1))
  pi.sim[,3] <- rbeta(nsim, eta2*mu2, eta2*(1-mu2))
  
  phi.sim <- t(apply(pi.sim, 1, compute.phi.direct))
  ne <- t(apply(phi.sim, 1, function(z) rmultinom(1, npl, z)))
  
  if (include_q){
    ne <- cbind(ne, p1.i, p2.i)
  }
  
  return(list(ne=ne, params=params.true, 
              s1=s1, r1=r1, r2=r2))
}




################################################################################
#
# Frequentist estimators
#

frequ_est = function(ne,    # 5-column matrix of error observations, each row is a read
                     npl=NA,  # if number of payload bases is NA, autodetects
                        est="naive", # estimator, one of naive, one step, or converged estimator
                        boot=TRUE,   # should uncertainty be bootstrapped (TRUE) or analytical (FALSE)?
                        nboot=1e3,  # number of bootstrap samples
                     maxit=100,   # maximum iterations for converged estimator
                     tol=1e-6,   # convergence tolerance
                     seed=1,   # random seed
                     simulate=TRUE  # simulate analytical values for boxplot (nboot times)
                     )   
{
  set.seed(as.numeric(seed))
  est <- substr(tolower(est),1,2)
  
  # Input processing
  n <- nrow(ne)
  n.read <- rowSums(ne[, 1:5])
  if(is.na(npl)){
    npl <- my.mode(n.read)
    if(npl != min(n.read)){
      warning(paste0("Minimal rowsum ", min(n.read), " is not equal current payload of ", npl))
    }
  }
  
  # Bootstrap if desired
  if(boot == TRUE){
    phat <- matrix(NA, nboot, 3)
    for(i in 1:nboot){
      ii <- sample(1:n, n, replace = TRUE)
      if(est == "na"){
        phat[i,] <- colMeans(ne[ii,c(2,4,5)]/npl)
      } else if (est == "on"){
        phat[i,] <- one.step.est(ne[ii,], npl)$v
      } else if (est == "co"){
        phat[i,] <- conv.est(ne[ii,], npl, 
                        maxit=maxit, tol=tol)$v
      } else {
        stop("Estimator must be one of naive, one-step, or converged")
      }
    }
    
    v <- apply(phat, 2, mean)
    vcov <- cov(phat)
    
  } else {  #analytical non-bootstrap uncertainties
    if(est == "na"){
      v <- colMeans(ne[,c(2,4,5)]/npl)
      cov.temp <- cov(ne[,c(2,4,5)]/npl)
      vcov <- 1/(n)*cov.temp
    } else if (est == "on"){
      res <- one.step.est(ne, npl)
      v <- res$v
      vcov <- 1/(n)*solve(res$H)
    } else if (est == "co"){
      res <- conv.est(ne, npl, 
                      maxit=maxit, tol=tol)
      v <- res$v
      vcov <- 1/(n)*solve(res$H)
    } else {
      stop("Estimator must be one of naive, one-step, or converged")
    }
    
    phat <- NA
    if(as.logical(simulate)){
      require(mvtnorm)
      phat <- rmvnorm(nboot, v, vcov)
    }
  }
  
  return(list(v=v, vcov=vcov, 
              vsim=phat))
}




################################################################################
#
#  Helpers and subfunctions
#


# Transform generative parameters to global mean error rate estimates
process_bayes <- function(ne, # nobs x 5 matrix of observations
                          params,   # nmcmc x 7 matrix of parameters corresponding to mcmc chain
                          npl,  # specify payload size 
                          withq=FALSE, 
                          nmax=500)   # max number of samples to consider
{
  n <- nrow(ne)
  nmcmc <- nrow(params)
  if(nmcmc > nmax){
    params <- params[tail(1:nrow(params),nmax),]
    nmcmc <- nrow(params)
  }
  means <- matrix(NA, nmcmc, 3)
  colnames(means) <- c("s1", "r1", "r2")
  
  # Global synthesis error rate
  mq <- sapply(params[,3], mu.q, npl)
  means[,1] <- params[,1]*params[,2] + (1-params[,2])*mq
  
  if(!withq){
    if(ncol(params) < 6){stop("Parmeter matrix must have at least 6 columns")}
    if(ncol(ne) < 5){stop("ne matrix must have at least 5 columns")}
    means[,2] <- params[,4]
    means[,3] <- params[,6]
    
  } else {
    if(ncol(params) < 8){stop("Parmeter matrix must have at least 8 columns")}
    if(ncol(ne) < 7){stop("ne matrix must have at least 7 columns")}
    #
    q1.i <-scale(qnorm(ne[,6]/npl))
    means[,2] <- sapply(1:nmcmc, function(x) mean(pnorm(params[x,4]*q1.i + params[x,5])))
    
    q2.i <- scale(qnorm(ne[,7]/npl))
    means[,3] <- sapply(1:nmcmc, function(x) mean(pnorm(params[x,7]*q2.i + params[x,8])))
  }
  
  return(means)
}




# Compute expectation for rv x \in [0,1], f(x) \propto (1-q)^(n*x)
mu.q <- function(q, n)
{
  # # old sim version, same as below
  # xx <- c(0:nsim)/nsim
  # dg <- (1-q)^c(xx*npl) # / pgeom(npl, q.hat)
  # # y.temp <- sample(xx, nn, replace=TRUE, prob=dg)
  # # return(mean(y.temp))
  # dg <- dg/sum(dg)
  # return(sum(dg*xx))
  
  # Analytical version
  L <- log(1-q)
  c <- (1-q)^n
  b <- n*L/(c-1)
  e <- b*(n*c*L -c + 1) / (n^2*L^2)
  e
}


# Function to compute phi directly as derived
compute.phi.direct <- function(v)
{
  p0=v[1]
  p1=v[2]
  p2=v[3]
  
  phi <- rep(0, 5)
  phi[1] <- (1-p0)*(1-p1)*(1-p2) + 1/9*p0*p1*p2
  phi[2] <- 1/3*(1-p0)*p1*p2 + p0*(1-p1)*(1-p2) + 2/9*p0*p1*p2
  phi[3] <- 2/3*(1-p0)*p1*p2 + 2/9*p0*p1*p2 + 2/3*p0*(p1*(1-p2) + p2*(1-p1))
  phi[4] <- (1-p0)*p1*(1-p2) + 2/9*p0*p1*p2 + 1/3*p0*p2*(1-p1)
  phi[5] <- (1-p0)*p2*(1-p1) + 2/9*p0*p1*p2 + 1/3*p0*p1*(1-p2)
  
  phi
}

# My coding of density of beta-binomial
#  N is number of sample in binomial (p integrated out)
#  a and b are canonical beta parameters
my.dbb <- function(x, N, a, b, log=TRUE)
{
  temp <- lchoose(N, x) + lbeta(x + a, N - x + b) - lbeta(a, b)
  if(!log){
    temp <- exp(temp)
  } 
  return(temp)
}

# DIC
dic2 <- function(ll)
{
  -2*mean(ll) + 2*var(ll)
}

# Integrate out pi, for x \sim Binom(npl, pi)
#  f(pi) \propto (1-q)^(npl*x) on [0,1]
int.pi <- function(y, npl, q, kmax=100, A=NULL)
{
  y <- as.integer(y)
  if(is.null(A)){
    A <- lchoose(npl,y) + lbeta(y+1, npl+1-y) + log(npl*log(1-q)/((1-q)^npl - 1))
  }
  ratio <- (y + 1 + 1:kmax - 1)/(npl + 2 + 1:kmax - 1)
  coef <- rep(0, kmax)
  for(i in 1:kmax){
    coef[i] <- prod(ratio[1:i])
  }
  terms <- rep(0, kmax)
  for(i in 1:kmax){
    terms[i] <- coef[i]*(npl*log(1-q))^i/factorial(i)
  }
  (1+sum(terms))*exp(A)
}

# Analytical density for rv x \in [0,1]
#   f(x) \propto (1-q)^(n*x)
fq <- function(x, q, npl)
{
  (1-q)^(x*npl)*npl*log(1-q) / ((1-q)^npl - 1)
}


# Approximate dg(p0_i;q) distribution with a beta distribution
appx.q.beta <- function(q, npl, nx=1e5)
{
  xx <- c(0:nx)/nx
  dg <- (1-q)^c(xx*npl)
  dg <- dg/sum(dg)  # skip integral
  
  mean <- sum(dg*xx)
  ex2 <- sum(dg*xx^2)
  var <- ex2 - mean^2
  
  eta <- mean*(1-mean) / var - 1
  # eta <- 10 #sqrt(.Machine$double.eps)
  
  a <- eta*mean
  b <- eta*(1-mean)
  
  return(c(a,b))
}


# MLE for beta-binomial with slope parameter
beta.opt.slope.v2 <- function(v, x, q, npl)
{
  mu <- pnorm(v[1]*q + v[2])
  eta <- exp(v[3])
  a <- eta*mu
  b <- eta*(1-mu)
  
  # eps <- 1e-10
  # x[x<eps] <- eps
  # x[x>1-eps] <- 1-eps
  
  # -sum(sapply(1:length(x), function(j) my.dbb(x[j], npl, a[j], b[j])))
  -sum(my.dbb(x, npl, a, b))
}



# Update mixture probability
update.alpha <- function(x,
                         e1, f1)
{
  m <- length(x)
  rbeta(1, e1 + sum(x), f1 + m - sum(x))
}


# Update mu and eta jointly for beta distribution
update.mu.eta.joint.jacob <- function(x,  # data
                                      mu, eta, # old parameters
                                      a1, a2,  # beta prior for mu
                                      b1, b2,  # gamma prior for eta
                                      tau.vec) # 2-vector of 
{
  # propose, props cancel
  theta <- c(qnorm(mu), log(eta))
  prop <- rnorm(2, theta, tau.vec)
  mu.new <- pnorm(prop[1])
  eta.new <- exp(prop[2])
  
  d.prop.mu <- 1/2*(qnorm(mu)^2 - qnorm(mu.new)^2)
  d.prop.eta <- log(eta.new) - log(eta)
  
  # LL evaluation
  ll.new <- sum(dbeta(x, mu.new*eta.new, eta.new*(1-mu.new), log=TRUE))
  ll.old <- sum(dbeta(x, mu*eta, eta*(1-mu), log=TRUE))
  
  prior.new <- dbeta(mu.new, a1, a2, log=TRUE) + dgamma(eta.new, b1, b2, log=TRUE)
  prior.old <- dbeta(mu, a1, a2, log=TRUE) + dgamma(eta, b1, b2, log=TRUE)
  
  alpha <- ll.new + prior.new - ll.old - prior.old + d.prop.mu + d.prop.eta
  accept <- log(runif(1,0,1)) < alpha
  mu.out <- mu
  eta.out <- eta
  if(accept){
    mu.out <- mu.new
    eta.out <- eta.new
  }
  
  return(list(mu=mu.out, eta=eta.out, accept=accept))
}



# Update a, b, and eta jointly for beta distribution
update.mu.eta.slope <- function(x, q.i, # data, p_i for appropriate dimension
                                a, b, eta, # old parameters
                                mu.a, sd.a,  # normal prior for a
                                mu.b, sd.b,  # normal prior for a
                                b1, b2,  # gamma prior for eta
                                tau.vec, # 3-vector of independent proposal sds
                                sigma=diag(3),   # optional multivariate
                                joint.prop=TRUE   # joint proposal??
                                ) 
{
  # propose, props cancel
  theta <- c(a, b, log(eta))
  if(joint.prop){
    require(mvtnorm)
    prop <- rmvnorm(2, theta, sigma)
  } else {
    prop <- rnorm(3, theta, tau.vec)
  }
  a.prop <- prop[1]
  b.prop <- prop[2]
  eta.new <- exp(prop[3])
  
  # d.prop.mu <- 1/2*(qnorm(mu)^2 - qnorm(mu.new)^2)
  d.prop.eta <- log(eta.new) - log(eta)
  
  # LL evaluation
  mu.new <- pnorm(a.prop*q.i + b.prop)
  mu <- pnorm(a*q.i + b)
  ll.new <- sum(dbeta(x, mu.new*eta.new, eta.new*(1-mu.new), log=TRUE))
  ll.old <- sum(dbeta(x, mu*eta, eta*(1-mu), log=TRUE))
  
  prior.new <- dnorm(a.prop, mu.a, sd.a, log=TRUE) + dnorm(b.prop, mu.b, sd.b, log=TRUE) + dgamma(eta.new, b1, b2, log=TRUE)
  prior.old <- dnorm(a, mu.a, sd.a, log=TRUE) + dnorm(b, mu.b, sd.b, log=TRUE) + dgamma(eta, b1, b2, log=TRUE)
  
  alpha <- ll.new + prior.new - ll.old - prior.old + d.prop.eta
  accept <- log(runif(1,0,1)) < alpha
  a.out <- a
  b.out <- b
  eta.out <- eta
  if(accept){
    a.out <- a.prop
    b.out <- b.prop
    eta.out <- eta.new
  }
  
  return(list(a=a.out, b=b.out, eta=eta.out, accept=accept))
}




# Reverse complement of a string; pretty slow
#  (just to check what is implemented in available packages)
my.rc <- function(v)
{
  vv <- unlist(strsplit(as.character(v), ""))
  w <- vv <- vv[length(vv):1]  # reverse
  # i.a <- 
  w[which(vv == "A")] <- "T"
  w[which(vv == "T")] <- "A"
  w[which(vv == "C")] <- "G"
  w[which(vv == "G")] <- "C"
  u <- paste0(w, collapse="")
  
  return(u)
}



# MLE step to really hone in on p0, q0, alpha0
neg.ll.mix <- function(w, x, npl)
{
  v <- pnorm(w)
  p0 <- v[1]
  alpha <- v[2]
  q <- v[3]
  
  -sum(log(alpha*dbinom(x, npl, p0) + (1-alpha)*dgeom(x, q)/pgeom(npl, q)))
}



# Function to optimize beta parameters (MLE)
beta.opt <- function(w,x)
{
  v=c(pnorm(w[1]), exp(w[2]))
  a=v[2]*v[1]
  b=v[2]*(1-v[1])
  
  eps <- 1e-10
  x[x<eps] <- eps
  x[x>1-eps] <- 1-eps
  
  -sum(dbeta(x,a,b,log=TRUE))
}
#  EXAMPLE USE:
# res = optim(c(qnorm(mean(x)), log(eta.true)), 
#             beta.opt, 
#             x=x)
# res
# pnorm(res$par[1])
# exp(res$par[2])


# Fast 3x3 determinant
my.det <- function(A, log=FALSE)
{
  temp <- A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2]) - A[1,2]*(A[2,1]*A[3,3] - A[2,3]*A[3,1]) + A[1,3]*(A[2,1]*A[3,2] - A[2,2]*A[3,1])
  
  if(log){
    return(log(abs(temp)))
  } else {
    return(temp)
  }
}


# Compute jacobian to work in phi-space (instead of probability space)
jacobian.direct <- function(v, log=FALSE)
{
  # p0=v[1]
  # p1=v[2]
  # p2=v[3]
  # 
  # J <- matrix(0, 3, 3)
  # J[,1] <- c(1 - p1 - p2 + 8/9*p1*p2,
  #            -p1 +1/3*p2 + 8/9*p1*p2,
  #            -p2 +1/3*p1 + 8/9*p1*p2
  # )
  # J[,2] <- c(1/3*p2 - p0 + 8/9*p0*p2,
  #            1 - p2 - p0 + 8/9*p0*p2,
  #            -p2 + 1/3*p0 + 8/9*p0*p2
  # )
  # J[,3] <- c(1/3*p1 - p0 + 8/9*p0*p1,
  #            -p1 + 1/3*p0 + 8/9*p0*p1,
  #            1 - p1 - p0 + 8/9*p0*p1
  # )
  # 
  
  # determinant(J, logarithm = log)$mod
  # return(my.det(J, log=log))
  J <- 0
  if(!log){
    J <- 1
  }
  return(J)
}




# Synthesis error rate parameter sampling
update.p0.vec <- function(nmat,  # matrix for zi==1
                          ps,    # read error parameters, subset to zi ==1
                          p.old, # previous global probability for low dispersion process
                          tau,   # number of observations for beta proposal
                          mu0, nu0  # beta priors for p0
)
{
  # beta priors
  a0 <- mu0*nu0
  b0 <- (1-mu0)*nu0
  tau.p0=tau
  eps <- 1e-300
  
  # # proposal
  # p.prop <- rbeta(1, p.old*tau.p0, (1-p.old)*tau.p0)
  # prior.prop <- dbeta(p.prop, p.old*tau.p0, (1-p.old)*tau.p0, log=TRUE)
  # prior.reverse <- dbeta(p.old, p.prop*tau.p0, (1-p.prop)*tau.p0, log=TRUE) 
  
  # Proposal
  # x <- sum(nmat[,2:3])
  # a1 <- a0 + tau*x
  # b1 <- b0 + tau*(npl*nrow(nmat) - x)
  # a1 <- tau*x
  # b1 <- tau*(npl*nrow(nmat) - x)
  # p.prop <- rbeta(1, a1, b1)
  # prior.prop <- dbeta(p.prop, a1, b1, log=TRUE)
  # prior.reverse <- dbeta(p.old, a1, b1, log=TRUE)
  
  p.prop <- pnorm(rnorm(1, qnorm(p.old), tau))
  prior.prop <- prior.reverse <- 0
  
  # Compute likelihoods, starting with priors
  dnew <- dbeta(p.prop, a0, b0, log=TRUE)
  dold <- dbeta(p.old, a0, b0, log=TRUE)
  
  ps.new <- cbind(p.prop, ps[,-1,drop=FALSE])
  ps.old <- cbind(p.old, ps[,-1,drop=FALSE])
  
  phi.new <- t(apply(ps.new, 1, compute.phi.direct))
  phi.old <- t(apply(ps.old, 1, compute.phi.direct))
  
  # j.new <- sum(apply(ps.new, 1, jacobian.direct,log=TRUE))
  # j.old <- sum(apply(ps.old, 1, jacobian.direct, log=TRUE))
  # 
  
  
  # Handle zeros
  if(sum(phi.new == 0) > 0){
    if(sum(phi.old == 0) > 0){
      i1 <- intersect(which(phi.new == 0), which(phi.old == 0))
      phi.new[i1] <- phi.old[i1] <- 1   # these will cancel
      
      i2 <- setdiff(which(phi.new == 0), which(phi.old == 0))
      phi.new[i2] <- eps
      
      i3 <- setdiff(which(phi.old == 0), which(phi.new == 0))
      phi.old[i3] <- eps
    } else {
      phi.new[phi.new == 0] <- eps
    }
  } else {
    if(sum(phi.old == 0) > 0){
      phi.old[phi.old == 0] <- eps
    }
  }
  
  dnew <- dnew + sum(log(phi.new)*nmat)  #+ j.new
  dold <- dold + sum(log(phi.old)*nmat)  #+ j.old
  # 
  # for(i in 1:nrow(ps)){
  #   
  #   phi.old <- compute.phi(c(p.old, ps[i,-1]))
  #   
  #   dnew <- dnew + sum(log(phi.new)*nmat[i,])
  #   dold <- dold + sum(log(phi.old)*nmat[i,])
  # }
  
  # accept/reject
  alpha <- dnew - dold + prior.reverse - prior.prop
  if(log(runif(1,0,1)) < alpha){
    accept <- TRUE
    p.out <- p.prop
  } else {
    accept <- FALSE
    p.out <- p.old
  }
  
  return(list(p=p.out, accept=accept))
  
}


# compute data log likelihood from params and ns
ll.compute <- function(p0, pi, zi, ns)
{
  pi[zi==1, 1] <- p0
  log.phi <- t(log(apply(pi, 1, compute.phi.direct)))
  # log.j <- t(apply(pi, 1, jacobian.direct, log=TRUE))
  rowSums(ns[,1:5]*log.phi)  #+ log.j
}



# Update latent probability (p.read, p. write, etc.) for oligo i
#  NOT IN USE TO MY KNOWLEDGE
sample.pi.vec <- function(ns,   # 5-column matrix, number of errors 1:5, and number of payload = sum(ns)
                          p.old,  # 3-column matrix, p.synth error, p.read1.error, p.read2.error
                          j,   # column to update
                          npl,   # payload length
                          tau0,  # tuning parameter, beta proposal number of equiv obs
                          mu0, nu0,   # beta prior parameters
                          zin, p0.in)   # other params
  
{
  if(sum(zin == 1) > 0 & j > 1){
    p.old[zin==1, 1] <- p0.in
  }
  
  nn <- nrow(p.old)
  pi.old <- p.old[,j]
  # eps <- 1e-300
  # eps <- (.Machine$double.eps)^15
  # pi.old[pi.old < eps] <- eps
  # pi.old[pi.old > 1-eps] <- 1-eps
  # 
  if(j > 1){
    x <- ns[,2 + j]
    
    # pi conjugacy for proposal
    a0 = mu0*nu0 + x*tau0
    b0 <- nu0*(1-mu0) + (npl - x)*tau0
    pi.prop <- rbeta(length(x), a0, b0)
    
    d.prop.new <- dbeta(pi.prop, a0, b0, log=TRUE)
    d.prop.old <- dbeta(pi.old, a0, b0, log=TRUE)
    
    prior.new <- dbeta(pi.prop, mu0*nu0, nu0*(1-mu0), log=TRUE)
    prior.old <- dbeta(pi.old, mu0*nu0, nu0*(1-mu0), log=TRUE)
    
  } else {
    x <- rowSums(ns[,2:3])
    # appx beta distr for appx proposa
    temp <- appx.q.beta(mu0, npl)
    
    
    # pi conjugacy for proposal, uniform prior (for proposal)
    a0 =x + tau0*temp[1]
    b0 <- (npl - x) + tau0*temp[2]
    pi.prop <- rbeta(length(x), a0, b0)
    
    d.prop.new <- dbeta(pi.prop, a0, b0, log=TRUE)
    d.prop.old <- dbeta(pi.old, a0, b0, log=TRUE)
    
    # pi.old <- p.old[,j]
    # pi.prop <- pi.old
    
    # Only update when zi == 0 (high dispersion process)
    # pi.prop[zi==0] <- pnorm(qnorm(pi.old[zi==0]) + rnorm(length(pi.old[zi==0]), 0, tau0))
    # pi.prop <- pnorm(qnorm(pi.old) + rnorm(length(pi.old), 0, tau0))
    # d.prop.new <- d.prop.old <- 0
    
    prior.old <- npl*pi.old*log(1-mu0) # geom prior
    prior.new <- npl*pi.prop*log(1-mu0) # geom prior
  }
  
  # compute likelihoods
  p.prop <- p.old
  p.prop[,j] <- pi.prop   # put in the right position
  
  phi.old <- t(apply(p.old, 1, compute.phi.direct)) # compute.phi(p.old)
  phi.new <-t(apply(p.prop, 1, compute.phi.direct)) # compute.phi(p.prop)
  
  # j.old <- t(apply(p.old, 1, jacobian.direct, log=TRUE)) # compute.phi(p.old)
  # j.new <- t(apply(p.prop, 1, jacobian.direct, log=TRUE)) # compute.phi(p.old)
  
  
  # # Dealw ith Infs
  # d.prop.new[d.prop.new == Inf] <- .Machine$double.xmax
  # d.prop.new[d.prop.new == -Inf] <- -.Machine$double.xmax
  # 
  # d.prop.old[d.prop.old == Inf] <- .Machine$double.xmax
  # d.prop.old[d.prop.old == -Inf] <- -.Machine$double.xmax
  # 
  # j.new[j.new == Inf] <- .Machine$double.xmax
  # j.new[j.new == -Inf] <- -.Machine$double.xmax
  # 
  # j.old[j.old == Inf] <- .Machine$double.xmax
  # j.old[j.old == -Inf] <- -.Machine$double.xmax
  # 
  # # Handle zeros
  # if(sum(phi.new == 0) > 0){
  #   if(sum(phi.old == 0) > 0){
  #     i1 <- intersect(which(phi.new == 0), which(phi.old == 0))
  #     phi.new[i1] <- phi.old[i1] <- 1   # these will cancel
  #     
  #     i2 <- setdiff(which(phi.new == 0), which(phi.old == 0))
  #     phi.new[i2] <- eps
  #     
  #     i3 <- setdiff(which(phi.old == 0), which(phi.new == 0))
  #     phi.old[i3] <- eps
  #   } else {
  #     phi.new[phi.new == 0] <- eps
  #   }
  # } else {
  #   if(sum(phi.old == 0) > 0){
  #     phi.old[phi.old == 0] <- eps
  #   }
  # }
  
  
  
  ll.old <- apply(log(phi.old)*ns, 1, sum) #+ j.old   # data LL 
  ll.new <- apply(log(phi.new)*ns, 1, sum) #+ j.new    #+ dbeta(p0.prop,a0,b0,log=TRUE)
  
  alpha <- ll.new - ll.old + prior.new - prior.old + d.prop.old - d.prop.new
  accept.vec <- log(runif(length(alpha))) < alpha
  p.out <- pi.old
  p.out[accept.vec] <- pi.prop[accept.vec]
  
  
  # nn <- nrow(p.old)
  # eps <- sqrt(.Machine$double.eps)
  # p.old[p.old < eps] <- eps
  # p.old[p.old > 1-eps] <- 1-eps
  # 
  # # prior params
  # a0 <- mu0*nu0
  # b0 <- (1-mu0)*nu0
  # 
  # 
  # # if(sum(zin == 1) > 0 & j > 1){
  # #   p.old[zin==1, 1] <- p0.in
  # # }
  # # 
  # 
  # # beta proposal
  # p0.prop <- rbeta(nn, p.old[,j]*tau0, (1-p.old[,j])*tau0)
  # p0.prop[p0.prop < eps] <- eps
  # p0.prop[p0.prop > 1 - eps] <- 1 - eps
  # 
  # 
  # prop.new <- dbeta(p0.prop, p.old[,j]*tau0, (1-p.old[,j])*tau0, log=TRUE)
  # prop.reverse <- dbeta(p.old[,j], p0.prop*tau0, (1-p0.prop)*tau0, log=TRUE)
  # 
  # 
  # if(j != 1){
  #   prior.old <- dbeta(p.old[,j],a0,b0,log=TRUE)  ## prior
  #   prior.prop <- dbeta(p0.prop,a0,b0,log=TRUE)  ## prior
  # } else {
  #   prior.old <- npl*p.old[,j]*log(1-mu0) # dbeta(,a0,b0,log=TRUE)  ## prior
  #   prior.prop <- npl*p0.prop*log(1-mu0) # dbeta(,a0,b0,log=TRUE)  ## LL + prior
  # }
  # 
  # 
  # # compute likelihoods
  # p.prop <- p.old
  # p.prop[,j] <- p0.prop   # put in the right position
  # 
  # phi.old <- t(apply(p.old, 1, compute.phi)) # compute.phi(p.old)
  # phi.new <-t(apply(p.prop, 1, compute.phi)) # compute.phi(p.prop)
  # 
  # dold <- apply(log(phi.old)*ns, 1, sum)  # data LL 
  # dnew <- apply(log(phi.new)*ns, 1, sum)    #+ dbeta(p0.prop,a0,b0,log=TRUE)
  # 
  # # accept/reject
  # alpha <- dnew - dold + prior.prop - prior.old - prop.new + prop.reverse
  # accept.vec <- log(runif(nn,0,1)) < alpha
  # 
  # p.out <- p.old
  # p.out[accept.vec,j] <- p0.prop[accept.vec]
  
  return(list(pi=p.out, accept=accept.vec))
}


# Compute numeric mode
my.mode <- function(x)
{
  tab <- table(x)
  return(as.numeric(names(tab)[which.max(tab)]))
}



# Compute contribution to complete LL of latent variables zi and pi
compute.ll.latent <- function(zi, pi, p0, alpha, q, mu1, eta1, mu2, eta2, n.read)
{
  ll.latent <- dbinom(zi, 1, alpha, log=TRUE)
  
  # ll.latent[zi==1] <- ll.latent[zi==1]
  ll.latent[zi==0] <- ll.latent[zi==0] + n.read[zi==0]*pi[zi==0,1]*log(1-q) + log(n.read[zi==0]*log(1-q) / ((1-q)^n.read[zi==0] - 1))
  
  ll.latent <- ll.latent + dbeta(pi[,2], mu1*eta1, eta1*(1-mu1), log=TRUE)
  ll.latent <- ll.latent + dbeta(pi[,3], mu2*eta2, eta2*(1-mu2), log=TRUE)
  
  sum(ll.latent)
}




# Given p (vector v), compute gradient of phi wrt p
compute.grad <- function(v)
{
  mat <- matrix(NA,3,5)
  mat[1,1] <- -(1-v[2])*(1-v[3]) + v[2]*v[3]/9
  mat[2,1] <- -(1-v[1])*(1-v[3]) + v[1]*v[3]/9
  mat[3,1] <- -(1-v[1])*(1-v[2]) + v[1]*v[2]/9
  
  mat[1,2] <- (1-v[2])*(1-v[3])-v[2]*v[3]/3
  mat[2,2] <- 1/3*(1-v[1])*v[3] - v[1]*(1-v[3])+2/9*v[1]*v[3]
  mat[3,2] <- 1/3*(1-v[1])*v[2] - v[1]*(1-v[2])+2/9*v[1]*v[2]
  
  mat[1,3] <- -4/9*v[2]*v[3] + 2/3*v[2]*(1-v[3]) + 2/3*v[3]*(1-v[2])
  mat[2,3] <- 2/3*(1-v[1])*v[3] + 2/9*v[1]*v[3] + 2/3*v[1]*(1 - 2*v[3])
  mat[3,3] <- 2/3*(1-v[1])*v[2] + 2/9*v[1]*v[2] + 2/3*v[1]*(1 - 2*v[2])
  
  mat[1,4] <- -v[2]*(1-v[3]) + 2/9*v[2]*v[3] + 1/3*v[3]*(1-v[2])
  mat[2,4] <- (1-v[1])*(1-v[3]) - 1/9*v[1]*v[3]
  mat[3,4] <- -(1-v[1])*v[2] + 2/9*v[1]*v[2] +1/3*v[1]*(1-v[2])
  
  mat[1,5] <- -v[3]*(1-v[2]) + 2/9*v[2]*v[3] + 1/3*v[2]*(1-v[3])
  mat[2,5] <-  -(1-v[1])*v[3] + 2/9*v[1]*v[3] +1/3*v[1]*(1-v[3])
  mat[3,5] <- (1-v[1])*(1-v[2]) - 1/9*v[1]*v[2]
  
  mat
}


compute.info <- function(v, ne)
{
  phi <- compute.phi.direct(v)
  A <- compute.grad(v)   # matrix of dphi/dp, 3x5
  dl.dp <- rep(0, 3)   #(a.vec-1)/v - (b.vec-1)/(1-v)
  M.mat <- matrix(0, nrow(ne), 3)
  # array(nrow(ne), 3, 5)
  for(j in 1:3){
    M.mat[,j] <- as.matrix(ne[, 1:5]) %*% c(c(1/phi)*A[j,])
    dl.dp[j] <- mean(M.mat[,j])
    # for(k in 1:3){
    #   Mk <- ne %*% (c(1/phi)*A[k,])
    # }
  }
  
  H <- matrix(0, 3, 3)
  for(j in 1:3){
    for(k in (j):3){
      H[j,k] <- H[k,j] <- mean(M.mat[,j]*M.mat[,k])
    }
  }
  return(list(H=H, S=dl.dp))
}


one.step.est <- function(ne, npl)
{
  v <- colMeans(ne[,c(2,4,5)]/npl)
  res <- compute.info(v, ne[,1:5])
  vnew <- v + solve(res$H, res$S)
  
  return(list(v=vnew, H=res$H, S=res$S))
}


conv.est <- function(ne, npl, 
                     maxit=100, tol=1e-6)
{
  v <- colMeans(ne[,c(2,4,5)]/npl)
  res <- compute.info(v, ne[,1:5])
  eps <- solve(res$H, res$S)
  delta <- sqrt(mean((eps^2)))
  vnew <- v + eps
  v <- vnew
  # vnew
  
  count <- 1
  while(count < maxit & delta > tol){
    count <- count + 1
    
    res <- compute.info(v, ne[,1:5])
    eps <- solve(res$H, res$S)
    vnew <- v + eps
    
    delta <- sqrt(mean((eps^2)))
    v <- vnew
  }
  
  return(list(v=v, 
              nit=count,
              delta=delta, tol=tol, 
              S=res$S, H=res$H))
}
