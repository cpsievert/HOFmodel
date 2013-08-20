#setwd("~/Stats/HOF/")
setwd("~/Desktop/HOFmodel")

#players <- read.csv("player_data.csv")
#elections <- read.csv("election_data.csv")
HOFdat <- read.csv("HOFvotingdata.csv")
#we only consider elections after 1967
election_dat <- subset(HOFdat, Year >= 1967)
#disregard players that get no votes?
player_dat <- subset(election_dat, Votes != 0)
#format years on ballot
player_dat$YoB <- as.integer(gsub("[a-z]+", "", player_dat$YoB))
#grab columns of interest and sort them by player name and years on ballot
election_vars <- c("Year", "Votes", "YoB", "NumBallots")
stats <- c("WAR", "WAR7", "JAWS")
datAll <- player_dat[c("Name", election_vars, stats)]
datAll <- plyr::arrange(datAll, Name, YoB)
#break up variables into their inherit levels of info
election_dat <- datAll[c("Name", election_vars)]
election_dat$prop <- with(election_dat, Votes/NumBallots)
player_dat <- unique(datAll[c("Name", stats)])
#standardize the covariates
#player_dat <- reshape::rescaler(player_dat)
for (i in 2:4) player_dat[,i] <- scale(player_dat[,i])



#helper functions
lu <- function(x) {
  if (is.factor(x)){
    length(levels(x))
  } else {
    length(unique(x))
  }
}


  require(MASS)
  require(coda)


  X <- cbind(1, player_dat[,-1])
  n.predictors <- dim(X)[2]
  n.players <- dim(X)[1]
  X <- as.matrix(X, ncol=n.predictors)
  for (i in 2:4) X[,i] <- runif(n.players, -1, 1)
  XtX <- t(X) %*% X
  stopifnot(player_dat$Name == unique(player_dat$Name))
  #this ensures nothing funky happens because of factors (with table, for instance) if a player is removed
  if (is.factor(election_dat$Name)) election_dat$Name <- as.character(election_dat$Name)  
  if (is.factor(player_dat$Name)) player_dat$Name <- as.character(player_dat$Name)
  repeats <- table(election_dat$Name)


# simulate fake values for effects of covariates:
beta_int.sim <- c(0.2, -0.3, 0.5, 0.2)
beta_slope.sim <- c(0, 0.03, 0.09, -0.02)

sigma_alpha.true <- 0.3
sigma_gamma.true <- 0.1

# simulate slopes and intercepts for each player:
alpha.sim <- rnorm(n.players, X %*% beta_int.sim, sigma_alpha)
gamma.sim <- rnorm(n.players, X %*% beta_slope.sim, sigma_gamma)

alpha_long <- rep(alpha.sim, repeats)
gamma_long <- rep(gamma.sim, repeats)

t <- election_dat$YoB
n <- election_dat$NumBallots
temp <- exp(alpha_long + gamma_long*t)
p.sim <- temp/(1 + temp)

y.sim <- rbinom(sum(repeats), size=n, prob=p.sim)


#options to be included in run_mcmc
alpha0 <- 0
gamma0 <- 0
beta_int0 <- 0
beta_slope0 <- c(-0.3, 0.3, -0.2, 0.4)
mu_int <- numeric(n.predictors)
sigma_int <- 5
mu_slope <- numeric(n.predictors)
sigma_slope <- 5
sigma_alpha0 <- 1
sigma_gamma0 <- 1
n.reps <- 5000
adapt <- 1000
burnin <- 500
tune <- TRUE
A <- 1.1 
B <- 1.1^(-44/56)


  #rename relevant player data
  #y <- election_dat$Votes
  y <- y.sim
  alpha_keep <- matrix(numeric(n.reps*n.players), 
                       nrow=n.reps, ncol=n.players)
  alpha_keep[1,] <- alpha <- rep(alpha0, n.players)     #set first alpha vector to alpha0
  gamma_keep <- matrix(numeric(n.reps*n.players), 
                       nrow=n.reps, ncol=n.players)
  gamma_keep[1,] <- gamma <- rep(gamma0, n.players)       #set first gamma vector to gamma0
  beta_int_keep <- matrix(numeric(n.reps*n.predictors), 
                          nrow=n.reps, ncol=n.predictors)
  beta_int_keep[1,] <- beta_int <- rep(beta_int0, n.predictors)
  beta_slope_keep <- matrix(numeric(n.reps*n.predictors), 
                            nrow=n.reps, ncol=n.predictors)
  beta_slope_keep[1,] <- beta_slope <- beta_slope0
  sigma_alpha_keep <- numeric(n.reps)
  sigma_alpha_keep[1] <- sigma_alpha <- sigma_alpha0
  sigma_gamma_keep <- numeric(n.reps)
  sigma_gamma_keep[1] <- sigma_gamma <- sigma_gamma0
  alpha_long <- rep(alpha, repeats)
  gamma_long <- rep(gamma, repeats)
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
  
  for (k in 2:n.reps) {
    if (k%%1000==0) print(k)
    
    #sample \beta^{int}
    sigma.int.inv <- (1/sigma_int^2)*diag(nrow=n.predictors) + XtX/sigma_alpha^2
    sigma.int <- solve(sigma.int.inv)
    mu.int <- sigma.int %*% (sigma.int.inv %*% matrix(mu_int, ncol=1) + (t(X) %*% alpha)/sigma_alpha^2)
    beta_int <- mvrnorm(n=1, mu=mu.int, Sigma=sigma.int)
    
    #sample \beta^{slope}
    sigma.slope.inv <- (1/sigma_slope^2)*diag(nrow=n.predictors) + XtX/sigma_gamma^2
    sigma.slope <- solve(sigma.slope.inv)
    mu.slope <- sigma.slope %*% (sigma.slope.inv %*% matrix(mu_slope, ncol=1) + (t(X) %*% gamma)/sigma_gamma^2)
    beta_slope <- mvrnorm(n=1, mu=mu.slope, Sigma=sigma.slope)

    beta_int_keep[k,] <- beta_int
    beta_slope_keep[k,] <- beta_slope
    
    #metropolis step for \alpha
    alpha_old <- alpha
    #rep alpha & gamma values to match the length of the election observations (for likelihood calculations)
    alpha_long <- rep(alpha, repeats)
    gamma_long <- rep(gamma, repeats)
    #proposed alpha & gamma (\alpha* & \gamma*)
    alpha_star <- rnorm(n.players, alpha, jump_alpha)
    alpha_star_long <- rep(alpha_star, repeats)
    #repeated calculations in the likelihood ratio (for alpha acceptance)
    xb_int <- X %*% beta_int
    a.star <- alpha_star_long + gamma_long*t
    a.stary <- 1/(1+exp(a.star))
    a.num <- diff(c(0, cumsum(y*log(1 - a.stary) + (n-y)*log(a.stary))[cumsum(repeats)]))
    a.old <- alpha_long + gamma_long*t
    a.olde <- 1/(1+exp(a.old))
    a.denom <- diff(c(0, cumsum(y*log(1 - a.olde) + (n-y)*log(a.olde))[cumsum(repeats)]))
    a.prior.star <- dnorm(alpha_star, xb_int, sigma_alpha, log=TRUE)
    a.prior.old <- dnorm(alpha, xb_int, sigma_alpha, log=TRUE)
    #lrho <- ((alpha - xb_int)^2 - (alpha_star - xb_int)^2)/2 + a.num - a.denom
    ratio <- exp(a.num + a.prior.star - a.denom - a.prior.old)
    # impose acceptance probability
    alpha <- ifelse(runif(n.players) < ratio, alpha_star, alpha)
    alpha_keep[k,] <- alpha
    if (tune & k < adapt) jump_alpha <- ifelse(ratio > 0.44, jump_alpha*A, jump_alpha*B)
    #metropolis step for \gamma
    gamma_old <- gamma
    alpha_long <- rep(alpha, repeats)
    gamma_long <- rep(gamma, repeats)
    gamma_star <- rnorm(n.players, gamma, jump_gamma)
    gamma_star_long <- rep(gamma_star, repeats)
    #repeated calculations in the likelihood ratio (for alpha acceptance)
    xb_slope <- X %*% beta_slope
    b.star <- alpha_long + gamma_star_long*t
    b.stary <- 1/(1+exp(b.star))
    b.num <- diff(c(0, cumsum(y*log(1 - b.stary) + (n-y)*log(b.stary))[cumsum(repeats)]))
    b.old <- alpha_long + gamma_long*t
    b.olde <- 1/(1+exp(b.old))
    b.denom <- diff(c(0, cumsum(y*log(1 - b.olde) + (n-y)*log(b.olde))[cumsum(repeats)]))
    b.prior.star <- dnorm(gamma_star, xb_slope, sigma_gamma, log=TRUE)
    b.prior.old <- dnorm(gamma, xb_slope, sigma_gamma, log=TRUE)
    ratio <- exp(b.num + b.prior.star - b.denom - b.prior.old)
    # impose acceptance probability
    gamma <- ifelse(runif(n.players) < ratio, gamma_star, gamma)
    gamma_keep[k,] <- gamma
    if (tune & k < adapt) jump_gamma <- ifelse(ratio > 0.44, jump_gamma*A, jump_gamma*B)
    #metropolis step for \sigma_{\alpha}
    sa_star <- abs(rnorm(1, sigma_alpha, jump1))
    num <- sum(dnorm(alpha, mean=xb_int, sd=sa_star, log=TRUE))
    denom <- prod(dnorm(alpha, mean=xb_int, sd=sigma_alpha))
    ratio <- exp(num/denom)
    # impose acceptance probability
    sigma_alpha <- ifelse(runif(1) < ratio, sa_star, sigma_alpha)
    sigma_alpha_keep[k] <- sigma_alpha
    if (tune & k < adapt) jump1 <- ifelse(ratio > 0.44, jump1*A, jump1*B)
    
    #metropolis step for \sigma_{\gamma}
    sg_star <- abs(rnorm(1, sigma_gamma, jump2))
    num <- prod(dnorm(gamma, mean=xb_slope, sd=sg_star))
    denom <- prod(dnorm(gamma, mean=xb_slope, sd=sigma_gamma))
    ratio <- num/denom
    # impose acceptance probability
    sigma_gamma <- ifelse(runif(1) < ratio, sg_star, sigma_gamma)
    sigma_gamma_keep[k] <- sigma_gamma
    if (tune & k < adapt) jump2 <- ifelse(ratio > 0.44, jump2*A, jump2*B)
    
    #calculate log-likelihood (note that we need the 'accepted' values of \alpha and \gamma)
    gamma_long <- rep(gamma, repeats)
    star <- alpha_long + gamma_long*t
    stary <- 1/(1+exp(star))
    loglik <- sum(y*log(1-stary) + (n-y)*log(stary))
    loglik_keep[k] <- loglik
  }



par(mfrow=c(1,1))
plot(1:n.reps, loglik_keep, type="l")

plot((adapt + 1):n.reps, loglik_keep[(adapt + 1):n.reps])

par(mfrow=c(1,2))
plot(1:n.reps, sigma_alpha_keep, col=1, type="l")
abline(h = sigma_alpha.true, col=2, lwd=2, lty=2)
plot(1:n.reps, sigma_gamma_keep, col=1, type="l")
abline(h = sigma_gamma.true, col=2, lwd=2, lty=2)
#acceptance probabilites (on the low side?)
lu(sigma_alpha_keep)/n.reps
lu(sigma_gamma_keep)/n.reps

par(mfrow=c(2,2))
for (i in 1:4){
  plot(1:n.reps, beta_int_keep[,i], col=1, type="l")
  abline(h = beta_int.sim[i], col=2, lwd=2, lty=2)
}

par(mfrow=c(2,2))
for (i in 1:4){
  plot(1:n.reps, beta_slope_keep[,i], col=1, type="l", ylim=c(-0.45, 0.45))
  abline(h = beta_slope.sim[i], col=2, lwd=2, lty=2)
}

plot(1:n.reps, alpha_keep[,1])
lu(alpha_keep[,1])
acc.prop.alpha <- apply(alpha_keep[(adapt + 1):n.reps,], 2, lu)/n.reps
hist(acc.prop.alpha)


samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  plot(1:n.reps, alpha_keep[,samp[i]], type="l", ylim=c(-2, 2), main=paste("Player", samp[i], repeats[samp[i]]))
  abline(h = alpha.sim[samp[i]], col=2, lwd=2, lty=2)
}

pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  pct[i,] <- quantile(alpha_keep[(adapt + 1):n.reps, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}

cover.50 <- sum(pct[,2] < alpha.sim & pct[,4] > alpha.sim)/n.players
cover.95 <- sum(pct[,1] < alpha.sim & pct[,5] > alpha.sim)/n.players
cover.50
cover.95


plot(1:n.reps, gamma_keep[,1])
lu(gamma_keep[,1])
acc.prop.gamma <- apply(gamma_keep[(adapt + 1):n.reps,], 2, lu)/n.reps
hist(acc.prop.gamma)


samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  plot(1:n.reps, gamma_keep[,samp[i]], type="l", ylim=c(-0.5, 0.5), main=paste("Player", samp[i], repeats[samp[i]]))
  abline(h = gamma.sim[samp[i]], col=2, lwd=2, lty=2)
}

pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  pct[i,] <- quantile(gamma_keep[(adapt + 1):n.reps, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}

cover.50 <- sum(pct[,2] < gamma.sim & pct[,4] > gamma.sim)/n.players
cover.95 <- sum(pct[,1] < gamma.sim & pct[,5] > gamma.sim)/n.players  
cover.50
cover.95  
  
