# Load a few libraries and set working directory:
rm(list=ls())
setwd("~/Stats/HOFmodel/")
library(coda)
lu <- function(x) length(unique(x))

# Read in the data:
HOFdat <- read.csv("HOFvotingdata.csv", stringsAsFactors=FALSE)

# Discard Warren Spahn's weird 1958 vote (cast while he was currently playing)
HOFdat <- HOFdat[-2213,]

# Append a variable that is an indicator of whether the player was a pitcher:
HOFdat$P <- as.numeric(HOFdat$position == "P")

# Select players for the analysis:
included.players <- HOFdat[HOFdat$YoB == "1st" & HOFdat$Year > 1966, "Name"]

# Select only those whose first ballot appearance was in 1967 or after:
#election_dat <- HOFdat[HOFdat[, "Year"] - yob >= 1966, ]
election_dat <- HOFdat[HOFdat$Name %in% included.players, ]

# Look at year of election and number of years on ballot:
election_dat$YoB <- as.integer(gsub("[a-z]+", "", election_dat$YoB))

# Analysis below of the number of players retained on ballot from year to year
# Should be non-increasing, but a few special cases.
# Looks like in 1985 there were a handful of players added back to the ballot
# who had previously been left off ballot.
# We'll treat their ballot years as contiguous, even though there were gaps.

# Disregard players that get no votes:
#player_dat <- subset(election_dat, Votes != 0)
player_dat <- election_dat

# Format years on ballot
#player_dat$YoB <- as.integer(gsub("[a-z]+", "", player_dat$YoB))

# Grab columns of interest and sort them by player name and years on ballot:
election_vars <- c("Year", "Votes", "YoB", "NumBallots")
stats <- c("WAR", "P")




dim(unique(election_dat[, c(3, 7:42, 44:46)]))


stats <- colnames(election_dat)[c(9:11, 14:39, 43)]


dim(unique(election_dat[, -c(1, 2, 4, 5, 6, 43)]))

st <- c(3, 7:40)


dim(unique(election_dat[, c(3, 7:40)]))

standard <- c("WAR")  # Variables that need to be standardized

# data:
datAll <- player_dat[c("Name", election_vars, stats)]
datAll <- plyr::arrange(datAll, Name, YoB)

# Break up variables into their inherit levels of info
election_dat <- datAll[c("Name", election_vars)]
election_dat$prop <- with(election_dat, Votes/NumBallots)
player_dat <- unique(datAll[c("Name", stats)])

# Standardize the covariates:
for (i in standard) player_dat[,i] <- scale(player_dat[,i])

# Rescale time to [0, 1]
election_dat$YoB <- (election_dat$YoB - 1)/14

n.players <- dim(player_dat)[1]
n.predictors <- dim(player_dat)[2]



### After the data is set up:
n.iter <- 5000
adapt <- 1000
burnin <- 1000


source("run-mcmc.R") # load run_mcmc

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
    res[[i]] <- run_mcmc(dat1=df1, dat2=df2, reps, alpha0=alpha0, gamma0=gamma0, sigma_alpha0=sigma_alpha0,
                         sigma_gamma0=sigma_gamma0, beta_int0=beta_int0, beta_slope0=beta_slope0, ...)
  }
  return(res)
}

repeats <- table(election_dat$Name)

t1 <- Sys.time()
res <- mcmc_chain(df1=election_dat, df2=player_dat, reps=repeats, n.players, n.predictors, n.chains=3, n.reps=n.iter)
t2 <- Sys.time()
t2 - t1
# 38 seconds for 3 chains, 5000 iterations each




# Plot the log-likelihood:
par(mfrow=c(1,2))

# All iterations
result <- res[[1]]
plot(1:n.iter, result$logliks[1:n.iter], col=1, type="l", ylab="Likelihood", xlab="Iteration Number")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(1:n.iter, result$logliks[1:n.iter], col=j, type="l")
  }
}

# Just zoom in on second-half iterations
result <- res[[1]]
plot(500:n.iter, result$logliks[500:n.iter], col=1, type="l")
if (length(res) > 1) {
  for (j in 2:length(res)) {
    result <- res[[j]]
    lines(500:n.iter, result$logliks[500:n.iter], col=j, type="l")
  }
}

par(mfrow=c(2,2))
for (i in 1:n.predictors) {
  result <- res[[1]]
  plot(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=1, type="l",  
       ylab=paste0("beta_slope", eval(i)), xlab="Iteration Number")
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- res[[j]]
      lines(1:n.iter, result$beta_slopes[,paste0("beta_slope", i)], col=j, type="l")
    }
  }
}




















### Additional analysis of years on ballot:

t <- table(election_dat[, "Year"], election_dat[, "YoB"])
nt <- dim(t)[1]

pdf(file="fig_ballot_dropoff.pdf", width=10, height=10)
par(mfrow=c(3, 3))
for (i in 1:nt) {
  plot(i:min(i + 14, nt) + 1966, t[cbind(i:min(i + 14, nt), 1:min(15, nt - i + 1))], xlab="Year", main=i + 1966)
  lines(i:min(i + 14, nt) + 1966, t[cbind(i:min(i + 14, nt), 1:min(15, nt - i + 1))])
}
dev.off()

for (i in 1:nt) {
  vec <- t[cbind(i:min(i + 14, nt), 1:min(15, nt - i + 1))]
  if (sum(diff(vec) > 0) > 0) {
    print("")
    print("")
    print("")
  	for (j in 1:1:min(15, nt - i + 1)) {
  	  print(election_dat[election_dat$Year == 1966 + i + j - 1 & election_dat$YoB == j, 1:5])
  	}
  }
}


