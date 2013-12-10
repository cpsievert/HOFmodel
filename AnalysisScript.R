# Load a few libraries and set working directory:
rm(list=ls())
setwd("~/Stats/HOFmodel/")
library(coda)
lu <- function(x) length(unique(x))

# Read in the data:
HOFdat <- read.csv("HOFvotingdata.csv", stringsAsFactors=FALSE, as.is=TRUE)

# Discard Warren Spahn's weird 1958 vote (cast while he was currently playing)
HOFdat <- HOFdat[-2213,]

# Append a variable that is an indicator of whether the player was a pitcher:
#HOFdat$P <- as.numeric(HOFdat$position == "P")

# Select players for the analysis:
included.players <- HOFdat[HOFdat$YoB == "1st" & HOFdat$Year > 1966, "Name"]

# Select only those whose first ballot appearance was in 1967 or after:
#election_dat <- HOFdat[HOFdat[, "Year"] - yob >= 1966, ]
election_dat <- HOFdat[HOFdat$Name %in% included.players, ]

# Look at year of election and number of years on ballot:
election_dat$YoB <- as.integer(gsub("[a-z]+", "", election_dat$YoB))

# For pitchers, I'll enter 0 for NA batting stats:
election_dat[is.na(election_dat[, "BA"]), c("BA", "OBP", "SLG", "OPS", "OPS.Plus")] <- 0

# Also, note that Jpos has only 10 unique values --- one for each position.
# If we use this, we should divide JAWS by Jpos, at the least
# For the first iteration, we'll probably not include these stats.


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
#election_vars <- c("Year", "Votes", "YoB", "NumBallots")
#stats <- c("WAR", "P")

# 636 unique players, and at least 40 variables that are player level stats (not player-year level stats):
lu(election_dat[, 3])
dim(unique(election_dat[, c(3, 7:42, 44:46)]))

# start a new data.frame called data with one row per player:
data <- unique(election_dat[, c(3, 7:42, 44:46)])
rownames(data) <- 1:dim(data)[1]

# pick off the 39 statistics that are common for every player:
stats <- colnames(data)[9:39]

# create a position matrix, with "P" as the baseline position:
pos.names <- sort(unique(data[,"position"]))[order(table(data[, "position"]), decreasing=TRUE)]
lp <- length(pos.names)
pos.mat <- matrix(0, nrow(data), lp - 1)
for (i in 1:10) pos.mat[, i] <- as.numeric(data[, "position"] == pos.names[i + 1])
colnames(pos.mat) <- pos.names[2:lp]

# logical values for pitchers and batters:
pitchers <- data[, "position"] == "P"
batters <- !pitchers

# Look at pitching stats of batters:
data[batters, ]

# Look at batting stats of pitchers:
data[pitchers, ]

# Look at NA values:
count.na <- function(x) sum(is.na(x))
apply(data, 2, count.na)
# OK, for batters, we can ignore all pitching stats.




### Try the baseline model for batters, starting with voting in 1997 (30 years training):

# Compute the voting percentages:
V <- election_dat[, "Votes"]/election_dat[, "NumBallots"]

# set up batting and pitching statistics:
bs <- colnames(data)[c(4, 9:21)]
ps <- colnames(data)[c(4, 22:34)]



# year to predict:
t <- 1997

# select type:
type <- 1  # 1 = batters, 2 = pitchers

# select either batters or pitchers:
sel <- election_dat[, "Year"] < t & election_dat[, "YoB"] == 1 & (election_dat[, "position"]=="P") == c(FALSE, TRUE)[type]

# Get the matrix for batters up to this year:
sel.players <- batters == c(TRUE, FALSE)[type]
m <- match(election_dat[sel, "Name"], data[sel.players, "Name"])
X.mat <- as.matrix(data[sel.players, ][m, bs])
X.names <- data[sel.players, "Name"][m]

pc <- prcomp(X.mat, scale=TRUE)



y.vec <- V[sel]

x.mat <- pc$x[, 1:6]

fit <- lm(y.vec ~ x.mat)





















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





### Principal Components Analysis on the whole data set:

### Batting:
bs <- colnames(data)[c(4, 9:21)]
Xb <- as.matrix(data[batters, bs])

b <- prcomp(Xb, scale=TRUE)
e <- eigen(cor(Xb))
e$values

plot(sqrt(e$values))
points(b$sdev, col=2)
# same. good.

plot(cumsum(e$values/sum(e$values)), las=1)
points(cumsum(b$sdev^2/sum(b$sdev^2)), col=2)
# looks like the first 6 PCs are useful




# Plot Batting PCs:
pdf(file="fig_pc_batting.pdf", width=6, height=6)
par(mfrow=c(1, 1))
for (i in 1:6) {
  plot(0, 0, type="n", xlim=c(-1.5, 1.5), ylim=c(1, 15), bty="n", xlab="PC Loading", yaxt="n", ylab="")
  for (j in 1:14) {
  	lines(c(0, b$rotation[j, i]), c(j, j), lwd=2)
  	text(b$rotation[j, i], j, colnames(Xb)[j], pos=c(2, 4)[as.numeric(b$rotation[j, i] > 0) + 1])
  }
  abline(v=0, lty=2, col=gray(0.7))
  text(-1.6, 15, "Top Ten Scores", pos=4, cex=0.6)
  text(-1.6, 14:5, paste(data[batters, "Name"][order(b$x[, i], decreasing=TRUE)[1:10]], round(sort(b$x[,i], decreasing=TRUE)[1:10], 1), sep=", "), pos=4, cex=0.6)
  text(1.6, 15, "Bottom Ten Scores", pos=2, cex=0.6)
  text(1.6, 14:5, paste(data[batters, "Name"][order(b$x[, i], decreasing=FALSE)[1:10]], round(sort(b$x[,i], decreasing=FALSE)[1:10], 1), sep=", "), pos=2, cex=0.6)
  title(main=paste0("Principal Component #", i, ": ", round(b$sdev[i]^2/sum(b$sdev^2)*100, 1), "% of Variation", sep=""))
}
dev.off()


### Pitching:
ps <- colnames(data)[c(4, 22:34)]
Xp <- as.matrix(data[pitchers, ps])

p <- prcomp(Xp, scale=TRUE)
e <- eigen(cor(Xp))
e$values

plot(sqrt(e$values))
points(p$sdev, col=2)
# same. good.

plot(cumsum(e$values/sum(e$values)), las=1)
points(cumsum(p$sdev^2/sum(p$sdev^2)), col=2)
# looks like the first 7 PCs are useful for pitchers:


# Plot them:
pdf(file="fig_pc_pitching.pdf", width=6, height=6)
par(mfrow=c(1, 1))
for (i in 1:7) {
  plot(0, 0, type="n", xlim=c(-1.5, 1.5), ylim=c(1, 15), bty="n", xlab="PC Loading", yaxt="n", ylab="")
  for (j in 1:14) {
  	lines(c(0, p$rotation[j, i]), c(j, j), lwd=2)
  	text(p$rotation[j, i], j, colnames(Xp)[j], pos=c(2, 4)[as.numeric(p$rotation[j, i] > 0) + 1])
  }
  abline(v=0, lty=2, col=gray(0.7))
  text(-1.6, 15, "Top Ten Scores", pos=4, cex=0.6)
  text(-1.6, 14:5, paste(data[pitchers, "Name"][order(p$x[, i], decreasing=TRUE)[1:10]], round(sort(p$x[,i], decreasing=TRUE)[1:10], 1), sep=", "), pos=4, cex=0.6)
  text(1.6, 15, "Bottom Ten Scores", pos=2, cex=0.6)
  text(1.6, 14:5, paste(data[pitchers, "Name"][order(p$x[, i], decreasing=FALSE)[1:10]], round(sort(p$x[,i], decreasing=FALSE)[1:10], 1), sep=", "), pos=2, cex=0.6)
  title(main=paste0("Principal Component #", i, ": ", round(p$sdev[i]^2/sum(p$sdev^2)*100, 1), "% of Variation", sep=""))
}
dev.off()



