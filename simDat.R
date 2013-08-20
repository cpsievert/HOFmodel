player_dat.sim <- player_dat
election_dat.sim <- election_dat

# simulate fake values for covariates:
for (i in names(player_dat.sim)[-1]) player_dat.sim[,i] <- runif(n.players, -1, 1)

# simulate fake values for effects of covariates:
# note that too high of a variance will break run_mcmc
beta_int.true <- rnorm(n.predictors, 0, 2)
beta_slope.true <- rnorm(n.predictors, 0, 2)
#beta_int.true <- c(0.2, -0.3, 0.5, 0.2)
#beta_slope.true <- c(0, 0.03, 0.09, -0.02)

sigma_alpha.true <- 1
sigma_gamma.true <- 1

# simulate slopes and intercepts for each player:
X <- cbind(1, as.matrix(player_dat.sim[,-1], ncol=n.predictors-1))
alpha.true <- rnorm(n.players, X %*% beta_int.true, sigma_alpha.true)
gamma.true <- rnorm(n.players, X %*% beta_slope.true, sigma_gamma.true)

repeats <- table(election_dat.sim$Name)
alpha_long.true <- rep(alpha.true, repeats)
gamma_long.true <- rep(gamma.true, repeats)

t <- election_dat.sim$YoB
n <- election_dat.sim$NumBallots
temp <- exp(alpha_long.true + gamma_long.true*t)
p.sim <- temp/(1 + temp)

y.sim <- rbinom(sum(repeats), size=n, prob=p.sim)
election_dat.sim$Votes <- y.sim
election_dat.sim$prop <- election_dat.sim$Votes/election_dat.sim$NumBallots

loglik.true <- sum(y.sim*log(p.sim) + (n-y.sim)*log(1-p.sim))