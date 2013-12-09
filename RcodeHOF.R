# Load some libraries:
#library(reshape)
#library(ggplot2)
library(coda)

head(player_dat)
idx <- NULL
for (i in seq_along(player_dat)){
  l <- is.numeric(player_dat[,i])
  idx <- c(idx, l)
}

#integer columns that we don't want to consider in PC analysis
vars1 <- c("Year", "Votes", "NumBallots")

#should we consider a seperate model for pitchers and batters????
vars2 <- c("")
names(player_dat) %in% summary(player_dat)

source("getDat.R")
#source("simDat.R")
@

<<realDat, cache=TRUE, results='hide'>>=
#real data
n.iter <- 5000
adapt <- 1000
burnin <- 1000
#estimate sigmas, different starting values
source("run-mcmc.R") #load run_mcmc
mcmc_chain <- function(df1, df2, reps, n.players, n.predictors, n.chains=3, ...) {
  res <- list(NULL)
  for (i in 1:n.chains) {
    sigma_alpha0 <- runif(1, 0.2, 0.5)
    sigma_gamma0 <- runif(1, 0.2, 0.5)
    beta_int0 <- rnorm(n.predictors, 0, 1)
    beta_slope0 <- rnorm(n.predictors, 0, 1)
    X <- cbind(1, as.matrix(df2[,-1], ncol=n.predictors-1))
    alpha0 <- rnorm(n.players, X %*% beta_int0, sigma_alpha0)
    gamma0 <- rnorm(n.players, X %*% beta_slope0, sigma_gamma0)
    res[[i]] <- run_mcmc(dat1=df1, dat2=df2, reps, alpha0=alpha0, gamma0=gamma0, sigma_alpha0=sigma_alpha0, sigma_gamma0=sigma_gamma0, beta_int0=beta_int0, beta_slope0=beta_slope0, ...)
  }
  return(res)
}
repeats <- table(election_dat$Name)
res <- mcmc_chain(df1=election_dat, df2=player_dat, reps=repeats, n.players, n.predictors, n.chains=3, n.reps=n.iter)
@

<<likelihood-real>>=
par(mfrow=c(1,2))
#all iterations
result <- res[[1]]
plot(1:n.iter, result$logliks[1:n.iter], col=1, type="l", ylab="Likelihood", xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$logliks[1:n.iter], col=j, type="l")
  }
}
#second-half iterations
result <- res[[1]]
plot(500:n.iter, result$logliks[500:n.iter], col=1, type="l")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(500:n.iter, result$logliks[500:n.iter], col=j, type="l")
  }
}
@

<<slopes-real>>=
par(mfrow=c(2,2))
for (i in 1:n.predictors) {
  result <- res[[1]]
  plot(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=1, type="l", ylab=paste0("beta_slope", eval(i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=j, type="l")
    }
  }
}
@

<<ints-real>>=
par(mfrow=c(2,2))
for (i in 1:n.predictors) {
  result <- res[[1]]
  plot(1:n.iter, result$beta_ints[,paste0("beta_int", i)], col=1, type="l", ylab=paste0("beta_int", eval(i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$beta_ints[,paste0("beta_int", i)], col=j, type="l")
    }
  }
}
@

<<alphas-real>>=
samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  result <- res[[1]]
  alpha.i <- samp[i]
  plot(1:n.iter, result$alphas[,paste0("alpha", alpha.i)], col=1, type="l", ylim=c(-10, 10), , main=paste0(player_dat$Name[alpha.i], ", n = ", repeats[alpha.i]), ylab=paste0("alpha", eval(alpha.i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$alphas[,paste0("alpha", alpha.i)], col=j, type="l")
    } 
  }
}
@

<<gammas-real>>=
samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  result <- res[[1]]
  gamma.i <- samp[i]
  plot(1:n.iter, result$gammas[,paste0("gamma", samp[i])], col=1, type="l", ylim=c(-15, 15), main=paste0(player_dat$Name[gamma.i], ", n = ", repeats[gamma.i]),  ylab=paste0("gamma", eval(gamma.i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$gammas[,paste0("gamma", samp[i])], col=j, type="l")
    }
  }
}
@


<<sigmas-real>>=
par(mfrow=c(1,2))
result <- res[[1]]
plot(1:n.iter, result$sigma_alphas, col=1, type="l", ylab=expression(sigma[alpha]), xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$sigma_alphas, col=j, type="l")
  }
}
result <- res[[1]]
plot(1:n.iter, result$sigma_gammas, col=1, type="l", ylab=expression(sigma[gamma]), xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$sigma_gammas, col=j, type="l")
  }
}
@


<<R-real, cache=TRUE>>=
loglik.list <- lapply(res, function(x) mcmc(x$logliks))
#gelman.diag(loglik.list)
gelman.plot(loglik.list, main="Gelman plot of log-likelihood")

par(mfrow=c(1,2))
alpha.list <- lapply(res, function(x) mcmc(x$alphas))
alphaR <- gelman.diag(alpha.list)
hist(alphaR[[1]][,2], main=expression(paste("Histogram of ", hat(R), " values for ", alpha)), xlab=expression(hat(R)))
gamma.list <- lapply(res, function(x) mcmc(x$gammas))
gammaR <- gelman.diag(gamma.list)
hist(gammaR[[1]][,2], main=expression(paste("Histogram of ", hat(R), " values for ", gamma)), xlab=expression(hat(R)))
@

<<R-betas-real, fig.height=9>>=
beta.list <- lapply(res, function(x) mcmc(cbind(x$beta_ints, x$beta_slopes)))
gelman.diag(beta.list)
gelman.plot(beta.list)
@

<<R-sigmas-real>>=
sigma.list <- lapply(res, function(x) mcmc(cbind(x$sigma_alphas, x$sigma_gammas)))
gelman.diag(sigma.list)
gelman.plot(sigma.list)
plot(mcmc.list(sigma.list))
@

<<segamma-real, cache=TRUE>>=
gamma.sum <- summary(mcmc.list(gamma.list))
sds <- gamma.sum[[1]][,2]
plot(as.numeric(repeats), sds)
@

% maybe
% <<coda-plots-real>>=
% samp <- paste0("alpha", sample(n.players, 3))
% alpha.sublist <- lapply(res, function(x) mcmc(x$alphas[,samp]))
% plot(mcmc.list(alpha.sublist))
% 
% samp <- paste0("gamma", sample(n.players, 3))
% gamma.sublist <- lapply(res, function(x) mcmc(x$gammas[,samp]))
% plot(mcmc.list(gamma.sublist))
% @


<<find-ps>>=
#Compute $p_{it}$ values for simulating $y_{it}$.
alphas <- res[[1]]$alphas[(burnin+adapt+1):n.iter,]
gammas <- res[[1]]$gammas[(burnin+adapt+1):n.iter,]
stuff <- sapply((1:15-1)/14, function(x) exp(alphas+gammas*x))
p <- stuff/(1+stuff)
ps <- array(p, dim=c(dim(alphas), 15))
players <- as.character(player_dat$Name)
dimnames(ps)[[2]] <- players
dimnames(ps)[[3]] <- paste0("YoB", 1:15)
@

<<residual-plots>>=
election_dat$YoB_orig <- election_dat$YoB*14 + 1
y.hats <- NULL
for (i in 1:dim(election_dat)[1]){
  obs <- election_dat[i,]
  n <- obs$NumBallots
  t <- obs$YoB_orig
  p.sim <- ps[,obs$Name,t]
  n.sim <- length(p.sim)
  y.sim <- rbinom(n=n.sim, size=rep(n, n.sim), prob=p.sim)
  y.hat <- median(y.sim)
  y.hats <- c(y.hats, y.hat)
}
election_dat$y_hat <- y.hats
election_dat$resid <- with(election_dat, Votes - y_hat)
qplot(x=YoB_orig, y=resid, data=election_dat)
qplot(x=NumBallots, y=resid, data=election_dat)+facet_wrap(~YoB_orig)
qplot(x=Year, y=resid, data=election_dat)+facet_wrap(~YoB_orig)
qplot(x=prop, y=resid/NumBallots, data=election_dat)+facet_wrap(~YoB_orig)
@


<<rmse, eval=FALSE, results='hide'>>=
yrs <- 2010:max(election_dat$Year)
rmses <- NULL
#compute the rmse for each year
for (i in yrs) {
  election_train <- subset(election_dat, Year < i)
  election_test <- subset(election_dat, Year == i)
  n <- unique(election_test$NumBallots)
  stopifnot(all(election_test$NumBallots %in% n))
  player_train <- subset(player_dat, Name %in% unique(election_train$Name))
  player_test <- subset(player_dat, Name %in% unique(election_test$Name))
  repeats <- table(election_train$Name)
  fit <- run_mcmc(dat1=election_train, dat2=player_train, reps=repeats, n.reps=n.iter)
  #get parameters for players that have been on ballot before
  idx <- names(repeats) %in% election_test$Name 
  alphas <- fit$alphas[(adapt+burnin+1):n.iter, idx]
  gammas <- fit$gammas[(adapt+burnin+1):n.iter, idx]
  old <- names(repeats)[idx]
  new.players <- election_test$Name[!election_test$Name %in% old]
  #get data for just the new players to simulate their values of alpha/beta/p/y
  #election_new <- subset(election_test, Name %in% new.players)
  player_new <- subset(player_test, Name %in% new.players)
  X <- as.matrix(cbind(1, player_new[,-1]))
  #get point estimates for betas
  beta_ints <- colMeans(as.matrix(fit$beta_ints))
  beta_slope <- colMeans(as.matrix(fit$beta_slopes))
  Xb_int <- X %*% beta_ints
  Xb_slope <- X %*% beta_slope
  sigma_alpha <- mean(fit$sigma_alphas)
  sigma_gamma <- mean(fit$sigma_gammas)
  n.sim <- dim(alphas)[1]
  new_alphas <- matrix(0, nrow=n.sim, ncol=length(Xb_int))
  colnames(new_alphas) <- player_new$Name
  for (i in Xb_int){
    new_alphas <- rnorm(n.sim, mean=i, sd=sigma_alpha)
  }
  new_gammas <- matrix(0, nrow=n.sim, ncol=length(Xb_slope))
  colnames(new_gammas) <- player_new$Name
  for (i in Xb_slope) {
    new_gammas <- rnorm(n.sim, mean=Xb_slope, sd=sigma_gamma)
  }
  all_alphas <- cbind(alphas, new_alphas)
  all_gammas <- cbind(gammas, new_gammas)
  election <- election_test[match(colnames(all_alphas), election_test$Name),]
  t <- election$YoB
  num <- exp(all_alphas + t(t(all_gammas) * t))
  p <- num/(1+num)
  #summary(c(p)) #shouldn't have any probabilities over 1....
  ys <- rbinom(n=length(p), size=n, prob=p)
  preds <- colMeans(matrix(ys, nrow=n.sim, ncol=dim(all_alphas)[2]))
  actual <- election$Votes
  rmse <- sqrt(sum((actual - preds)^2)/length(preds))
  rmses <- c(rmse, rmses)
}

qplot(rmses, geom="line")#+facet_wrap(~model)

@



<<compare, cache=TRUE>>=
#ballots.yr <- plyr::arrange(unique(HOFdat[c("Year", "NumBallots")]), Year)
#function to compare posterior distribution of p to the observed values
#this function assumes that ps & election_dat exist in the global environment :\
compare <- function(player="Jack Morris", geom="line", alpha=.01) { 
  edat <- subset(election_dat, Name %in% player)
  if (dim(edat)[2] == 0) stop("No election data was found. Are you sure you spelled the name right?")
  if (dim(edat)[2] == 1) stop("Only only election year exists for this player. Why bother?")
  edat$YoB <- round(14*edat$YoB+1)
  YoB <- paste0("YoB", edat$YoB)
  n <- edat$NumBallots
  p <- ps[,player,]
  N <- dim(p)[1]
  bins <- matrix(0, nrow=N, ncol=length(YoB))
  ctr <- 1
  for (i in YoB) {
    pi <- p[,i]
    bins[,ctr] <- rbinom(N, size=n[ctr], prob=pi)
    ctr <- ctr+1
  }
  bins2 <- t(t(bins)/n)
  propDat <- melt(data.frame(bins2, sim=1:dim(bins2)[1]), id.vars="sim")
  propDat$variable <- as.numeric(gsub("X", "", propDat$variable))
  if (geom=="boxplot") simlayer <- geom_boxplot(data=propDat, aes(y=value, x=factor(variable)))
  if (geom=="line") simlayer <- geom_line(data=propDat, aes(y=value, x=variable, group=sim), alpha=alpha)
  ggplot()+simlayer+ylab("Proportion of Votes")+xlab("Years on Ballot")+
    geom_line(data=edat, aes(y=prop, x=YoB), color="red")+geom_point(data=edat, aes(y=prop, x=YoB), color="red", size=10, shape="*")+
    ggtitle(paste("Observed versus the posterior distribution of proportion of votes for\n", player))+
    geom_hline(yint=.75, color="blue")
}
compare() #alpha probably should be a function of the number of iterations
compare("Mark McGwire")
compare("Gary Carter")
compare("Ryne Sandberg")
compare("Gil Hodges")
compare()
@



<<prop-inside>>=
#calculate the percentage of ALL observed values that fall within their respective 50% credible interval
qs <- apply(ps, c(2, 3), function(x) quantile(x, c(0.25, 0.75)))
df <- data.frame(qs[1,,])
df$Name <- rownames(df)
m.df <- melt(df, variable_name="YoB")
m.df$YoB <- (as.numeric(gsub("YoB", "", m.df$YoB))-1)/14
names(m.df) <- gsub("value", "lower", names(m.df))
df2 <- data.frame(qs[2,,])
df2$Name <- rownames(df2)
m.df$upper <- melt(df2)$value
eldat <- plyr::join(election_dat, m.df, by=c("Name", "YoB"))
#https://stat.ethz.ch/pipermail/r-help/2008-August/170749.html
eldat$inside <- as.numeric(with(eldat, (prop-lower) * (upper-prop) > 0))
with(eldat, sum(inside)/length(inside))
@


<<interest, eval=FALSE, fig.cap="Posterior distribution of the average percentage of ballots as a function of the number of years on the ballot. The percentage of ballots received increases with the number of years on ballot.">>=
#pf <- function(paramz, t) with(paramz, alphas+gammas*t)
#l <- lapply(1:10, function(x) pf(params, x))
means <- lapply(l, function(x) rowMeans(x))
names(means) <- as.character(1:10)
df <- data.frame(means)
m.df <- melt(df)
ggplot(data=m.df, aes(x=factor(variable), y=value))+geom_boxplot()+ylab("Avg percentage of Ballots")+xlab("Years on Ballot")
@


\section{Verification of Sampler using Simulated Data}

<<verify1, eval=FALSE>>=
#######  This chunk is mainly for debugging purposes
#known sigmas WORKS
source("run-mcmc2.R")
res <- list(NULL)
res[[1]] <- run_mcmc(dat1=election_dat.sim, dat2=player_dat.sim, sigma_alpha0=sigma_alpha.true, sigma_gamma0=sigma_gamma.true)
res[[2]] <- run_mcmc(dat1=election_dat.sim, dat2=player_dat.sim, sigma_alpha0=sigma_alpha.true, sigma_gamma0=sigma_gamma.true)

#known sigmas, but different starting values - WORKS (but low acceptance rates??)
mcmc_chain <- function(n.players, n.predictors, n.chains=3, ...) {
  res <- list(NULL)
  for (i in 1:n.chains) {
    sigma_alpha0 <- sigma_alpha.true
    sigma_gamma0 <- sigma_gamma.true
    beta_int0 <- rnorm(4, 0, 1)
    beta_slope0 <- rnorm(4, 0, 1)
    alpha0 <- rnorm(n.players, X %*% beta_int0, sigma_alpha0)
    gamma0 <- rnorm(n.players, X %*% beta_slope0, sigma_gamma0)
    res[[i]] <- run_mcmc(dat1=election_dat.sim, dat2=player_dat.sim, alpha0=alpha0, gamma0=gamma0, sigma_alpha0=sigma_alpha0, sigma_gamma0=sigma_gamma0, beta_int0=beta_int0, beta_slope0=beta_slope0, ...)
  }
  return(res)
}
res <- mcmc_chain(n.players, n.predictors, n.chains=2, n.reps=n.iter)
@

<<mcmc_chain, cache=TRUE, results='hide'>>=
repeats <- table(election_dat.sim$Name)
res <- mcmc_chain(df1=election_dat.sim, df2=player_dat.sim, reps=repeats, n.players=n.players, n.predictors=n.predictors, n.chains=2, n.reps=n.iter)
# #without dispersed starting values
# res <- list(NULL)
# res[[1]] <- run_mcmc(dat1=election_dat.sim, dat2=player_dat.sim)
# res[[2]] <- run_mcmc(dat1=election_dat.sim, dat2=player_dat.sim)
@

First we compute the Gelman-Rubin convergence diagnostic to monitor convergence.

<<likelihood>>=
par(mfrow=c(1,2))
#all iterations
result <- res[[1]]
plot(1:n.iter, result$logliks[1:n.iter], col=1, type="l", ylab="Likelihood", xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$logliks[1:n.iter], col=j+1, type="l")
  }
}
abline(h=loglik.true, col=2, lwd=2, lty=2)
#second-half iterations
result <- res[[1]]
plot(500:n.iter, result$logliks[500:n.iter], col=1, type="l", ylab="Likelihood", xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(500:n.iter, result$logliks[500:n.iter], col=j+1, type="l")
  }
}
abline(h=loglik.true, col=2, lwd=2, lty=2)
@

<<slopes>>=
par(mfrow=c(2,2))
for (i in 1:n.predictors) {
  result <- res[[1]]
  plot(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=1, type="l", ylab=paste0("beta_slope", eval(i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=j+1, type="l")
    }
  }
  abline(h=beta_slope.true[i], col=2, lwd=2, lty=2)
}
@

<<ints>>=
par(mfrow=c(2,2))
for (i in 1:n.predictors) {
  result <- res[[1]]
  plot(1:n.iter, result$beta_ints[,paste0("beta_int", i)], col=1, type="l", ylab=paste0("beta_int", eval(i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$beta_ints[,paste0("beta_int", i)], col=j+1, type="l")
    }
  }
  abline(h=beta_int.true[i], col=2, lwd=2, lty=2)
}
@

50 and 95 percent coverage rates

<<alpha_cover>>=
pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  pct[i,] <- quantile(res[[1]]$alphas[(adapt + 1):n.iter, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}
# 50% coverage rate
sum(pct[,2] < alpha.true & pct[,4] > alpha.true)/n.players
# 95% coverage rate
sum(pct[,1] < alpha.true & pct[,5] > alpha.true)/n.players
@

<<alphas_rates>>=
par(mfrow=c(length(res),1))
alpha.rates <- list(NULL)
for (i in 1:length(res)) {
  alpha.rates[[i]] <- apply(res[[i]]$alphas[(adapt + 1):n.iter,], 2, lu)/n.iter
  hist(alpha.rates[[i]])
}
@

<<alphas>>=
samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  result <- res[[1]]
  alpha.i <- samp[i]
  plot(1:n.iter, result$alphas[,paste0("alpha", alpha.i)], col=1, type="l", ylim=c(-10, 10), ylab=paste0("alpha", eval(alpha.i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$alphas[,paste0("alpha", alpha.i)], col=j+1, type="l")
    } 
  }
  abline(h = alpha.true[alpha.i], col=2, lwd=2, lty=2)
}
@

50 and 95 percent coverage rates

<<gamma_cover>>=
g.pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  g.pct[i,] <- quantile(res[[1]]$gammas[(adapt + 1):n.iter, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}
# 50% coverage rate (for first chain)
sum(g.pct[,2] < gamma.true & g.pct[,4] > gamma.true)/n.players
# 95% coverage rate (for first chain)
sum(g.pct[,1] < gamma.true & g.pct[,5] > gamma.true)/n.players
@

<<gamma_rates>>=
par(mfrow=c(length(res),1))
gamma.rates <- list(NULL)
for (i in 1:length(res)) {
  gamma.rates[[i]] <- apply(res[[i]]$gammas[(adapt + 1):n.iter,], 2, lu)/n.iter
  hist(gamma.rates[[i]])
}
@

<<gammas>>=
samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  result <- res[[1]]
  gamma.i <- samp[i]
  plot(1:n.iter, result$gammas[,paste0("gamma", gamma.i)], col=1, type="l", ylim=c(-15, 15), ylab=paste0("gamma", eval(gamma.i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$gammas[,paste0("gamma", gamma.i)], col=j+1, type="l")
    }
  }
  abline(h = gamma.true[gamma.i], col=2, lwd=2, lty=2)
}
@


<<sigmas>>=
par(mfrow=c(1,2))
result <- res[[1]]
plot(1:n.iter, result$sigma_alphas, col=1, type="l")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$sigma_alphas, col=j+1, type="l")
  }
}
abline(h=sigma_alpha.true, col=2, lwd=2, lty=2)
result <- res[[1]]
plot(1:n.iter, result$sigma_gammas, col=1, type="l")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$sigma_gammas, col=j+1, type="l")
  }
}
abline(h=sigma_gamma.true, col=2, lwd=2, lty=2)
@


<<R, cache=TRUE>>=
loglik.list <- lapply(res, function(x) mcmc(x$logliks))
#gelman.diag(loglik.list)
gelman.plot(loglik.list, main="Gelman plot of log-likelihood")

par(mfrow=c(1,2))
alpha.list <- lapply(res, function(x) mcmc(x$alphas))
alphaR <- gelman.diag(alpha.list)
hist(alphaR[[1]][,2], main=expression(paste("Histogram of ", hat(R), " values for ", alpha)), xlab=expression(hat(R)))
gamma.list <- lapply(res, function(x) mcmc(x$gammas))
gammaR <- gelman.diag(gamma.list)
hist(gammaR[[1]][,2], main=expression(paste("Histogram of ", hat(R), " values for ", gamma)), xlab=expression(hat(R)))
@

<<R-betas, fig.height=9>>=
beta.list <- lapply(res, function(x) mcmc(cbind(x$beta_ints, x$beta_slopes)))
gelman.diag(beta.list)
gelman.plot(beta.list)
@

<<R-sigmas>>=
sigma.list <- lapply(res, function(x) mcmc(cbind(x$sigma_alphas, x$sigma_gammas)))
gelman.diag(sigma.list)
gelman.plot(sigma.list)
plot(mcmc.list(sigma.list))
@


% maybe
% <<coda-plots>>=
% samp <- paste0("alpha", sample(n.players, 3))
% alpha.sublist <- lapply(res, function(x) mcmc(x$alphas[,samp]))
% plot(mcmc.list(alpha.sublist))
% 
% samp <- paste0("gamma", sample(n.players, 3))
% gamma.sublist <- lapply(res, function(x) mcmc(x$gammas[,samp]))
% plot(mcmc.list(gamma.sublist))
% @
