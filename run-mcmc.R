#helper functions
lu <- function(x) {
  if (is.factor(x)){
    length(levels(x))
  } else {
    length(unique(x))
  }
}

metropolis <- function(alpha, gamma, sigma, jump_sd, xb, reps, y, n, t, sample.alpha=TRUE) {
  #rep alpha & gamma values to match the length of the election observations (for likelihood calculations)
  alpha_long <- rep(alpha, reps)
  gamma_long <- rep(gamma, reps)
  #calculate denominator of acceptance ratio with the "old" values of alpha & gamma
  #browser()
  old <- alpha_long + gamma_long*t
  olde <- 1/(1+exp(old))
  cum.reps <- cumsum(reps)
  denom <- diff(c(0, cumsum(y*log(1 - olde) + (n-y)*log(olde))[cum.reps]))
  #proposed alpha & gamma (\alpha* & \gamma*)
  n.players <- length(alpha)
  if (sample.alpha) {
    alpha_star <- rnorm(n.players, alpha, jump_sd)
    #overwrite alpha_long with proposals
    alpha_long <- rep(alpha_star, reps)
    prior.star <- dnorm(alpha_star, xb, sigma, log=TRUE)
    prior.old <- dnorm(alpha, xb, sigma, log=TRUE)
  } else { #sample gamma
    gamma_star <- rnorm(n.players, gamma, jump_sd)
    #overwrite alpha_long with proposals
    gamma_long <- rep(gamma_star, reps)
    prior.star <- dnorm(gamma_star, xb, sigma, log=TRUE)
    prior.old <- dnorm(gamma, xb, sigma, log=TRUE)
  }
  #repeated calculations in the acceptance ratio
  star <- alpha_long + gamma_long*t
  stary <- 1/(1+exp(star))
  num <- diff(c(0, cumsum(y*log(1 - stary) + (n-y)*log(stary))[cum.reps]))
  ratio <- exp(num + prior.star - denom - prior.old)
  # impose acceptance probability
  if (sample.alpha) {
    alpha <- ifelse(runif(n.players) < ratio, alpha_star, alpha)
    return(list(alpha=alpha, ratio=ratio))
  } else{
    gamma <- ifelse(runif(n.players) < ratio, gamma_star, gamma)
    return(list(gamma=gamma, ratio=ratio))
  }
}

#MCMC workhorse
#Data on the year level (ie, number of ballots) should be input as dat1
#Data on the player level (ie, career statistics) should be input as dat2
#Both dat1 and dat2 should include player names as "Name" in the first column
#dat2 should include one or more 'predictors'
run_mcmc <- function(dat1, dat2, reps, alpha0 = 0, sigma_alpha0 = 1, gamma0 = 0, sigma_gamma0 = 1, 
                     beta_int0 = 0, mu_int = 0, sigma_int = 5, beta_slope0 = 0, mu_slope = 0, 
                     sigma_slope = 5, a=10, b=10, n.reps = 5000, adapt = 1000, tune = TRUE){
  require(MASS)
  X <- cbind(1, dat2[,-1])
  n.predictors <- dim(X)[2]
  n.players <- dim(X)[1]
  X <- as.matrix(X, ncol=n.predictors)
  XtX <- t(X) %*% X
  stopifnot(dat2$Name == unique(dat2$Name))
  #rename relevant player data
  #browser()
  y <- dat1$Votes
  t <- dat1$YoB
  n <- dat1$NumBallots
  alpha_keep <- matrix(numeric(n.reps*n.players), 
                       nrow=n.reps, ncol=n.players)
  alpha_keep[1,] <- alpha <- if (length(alpha0)==1) rep(alpha0, n.players) else alpha0
  gamma_keep <- matrix(numeric(n.reps*n.players), 
                       nrow=n.reps, ncol=n.players)
  gamma_keep[1,] <- gamma <- if (length(gamma0)==1) rep(gamma0, n.players) else gamma0
  beta_int_keep <- matrix(numeric(n.reps*n.predictors), 
                          nrow=n.reps, ncol=n.predictors)
  beta_int_keep[1,] <- beta_int <- beta_int0
  beta_slope_keep <- matrix(numeric(n.reps*n.predictors), 
                            nrow=n.reps, ncol=n.predictors)
  beta_slope_keep[1,] <- beta_slope <- beta_slope0
  sigma_alpha_keep <- numeric(n.reps)
  sigma_alpha_keep[1] <- sigma_alpha <- sigma_alpha0
  sigma_gamma_keep <- numeric(n.reps)
  sigma_gamma_keep[1] <- sigma_gamma <- sigma_gamma0
  alpha_long <- rep(alpha, reps)
  gamma_long <- rep(gamma, reps)
  star <- alpha_long + gamma_long*t
  stary <- 1/(1+exp(star))
  loglik <- sum(y*log(1-stary) + (n-y)*log(stary))
  loglik_keep <- numeric(n.reps)
  loglik_keep[1] <- loglik  
  #jumping std devs for metropolis (include these as an option?)
  jump_alpha <- rep(0.10, n.players)
  jump_gamma <- rep(0.05, n.players)
  jump1 <- 1
  jump2 <- 1
  A <- 1.1 
  B <- 1.1^(-44/56)
  for (k in 2:n.reps) {
    if (k%%1000==0) cat(k, "\n")
    
    #sample \beta^{int}
    sigma.int.inv <- (1/sigma_int^2)*diag(nrow=n.predictors) + XtX/sigma_alpha^2
    sigma.int <- solve(sigma.int.inv)
    mu.int <- sigma.int %*% (sigma.int.inv %*% matrix(mu_int, nrow=n.predictors, ncol=1) + (t(X) %*% alpha)/sigma_alpha^2)
    beta_int <- mvrnorm(n=1, mu=mu.int, Sigma=sigma.int)
    
    #sample \beta^{slope}
    sigma.slope.inv <- (1/sigma_slope^2)*diag(nrow=n.predictors) + XtX/sigma_gamma^2
    sigma.slope <- solve(sigma.slope.inv)
    mu.slope <- sigma.slope %*% (sigma.slope.inv %*% matrix(mu_slope, nrow=n.predictors, ncol=1) + (t(X) %*% gamma)/sigma_gamma^2)
    beta_slope <- mvrnorm(n=1, mu=mu.slope, Sigma=sigma.slope)
    
    beta_int_keep[k,] <- beta_int
    beta_slope_keep[k,] <- beta_slope
    
    xb_int <- X %*% beta_int
    #metropolis step for \alpha
    stuff <- metropolis(alpha=alpha, gamma=gamma, sigma=sigma_alpha, jump_sd=jump_alpha, xb=xb_int, reps=reps, y=y, n=n, t=t)
    alpha_keep[k,] <- alpha <- stuff$alpha
    if (tune & k < adapt) jump_alpha <- ifelse(stuff$ratio > 0.44, jump_alpha*A, jump_alpha*B)
    
    xb_slope <- X %*% beta_slope
    #metropolis step for \gamma
    stuff <- metropolis(alpha=alpha, gamma=gamma, sigma=sigma_gamma, jump_sd=jump_gamma, xb=xb_slope, reps=reps, y=y, n=n, t=t, sample.alpha=FALSE)
    gamma_keep[k,] <- gamma <- stuff$gamma
    if (tune & k < adapt) jump_gamma <- ifelse(stuff$ratio > 0.44, jump_gamma*A, jump_gamma*B)
    
    #metropolis step for \sigma_{\alpha}
    sa_star <- rnorm(1, sigma_alpha, jump1)
    if (sa_star > a | sa_star < 0) {
      ratio <- 0
    } else {
      num <- sum(dnorm(alpha, mean=xb_int, sd=sa_star, log=TRUE))
      denom <- sum(dnorm(alpha, mean=xb_int, sd=sigma_alpha, log=TRUE))
      ratio <- exp(num - denom)
    }
    # impose acceptance probability
    sigma_alpha <- ifelse(runif(1) < ratio, sa_star, sigma_alpha)
    sigma_alpha_keep[k] <- sigma_alpha
    if (tune & k < adapt) jump1 <- ifelse(ratio > 0.44, jump1*A, jump1*B)
    
    #metropolis step for \sigma_{\gamma}
    sg_star <- rnorm(1, sigma_gamma, jump2)
    if (sg_star > b | sg_star < 0) {
      ratio <- 0
    } else {
      num <- sum(dnorm(gamma, mean=xb_slope, sd=sg_star, log=TRUE))
      denom <- sum(dnorm(gamma, mean=xb_slope, sd=sigma_gamma, log=TRUE))
      ratio <- exp(num - denom)
    }
    # impose acceptance probability
    sigma_gamma <- ifelse(runif(1) < ratio, sg_star, sigma_gamma)
    sigma_gamma_keep[k] <- sigma_gamma
    if (tune & k < adapt) jump2 <- ifelse(ratio > 0.44, jump2*A, jump2*B)
    
    #calculate log-likelihood (note that we need the 'accepted' values of \alpha and \gamma)
    alpha_long <- rep(alpha, reps)
    gamma_long <- rep(gamma, reps)
    star <- alpha_long + gamma_long*t
    stary <- 1/(1+exp(star))
    loglik <- sum(y*log(1-stary) + (n-y)*log(stary))
    loglik_keep[k] <- loglik
  }
  colnames(alpha_keep) <- dat2$Name
  colnames(gamma_keep) <- dat2$Name
  colnames(beta_int_keep) <- paste0("beta_int", 1:n.predictors)
  colnames(beta_slope_keep) <- paste0("beta_slope", 1:n.predictors)
  return(list(alphas=alpha_keep, 
              gammas=gamma_keep, 
              beta_ints=beta_int_keep, 
              beta_slopes=beta_slope_keep, 
              sigma_alphas=sigma_alpha_keep, 
              sigma_gammas=sigma_gamma_keep, 
              logliks=loglik_keep))
}