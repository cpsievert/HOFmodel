# Load a few libraries and set working directory:
rm(list=ls())
setwd("~/Stats/HOFmodel/")
lu <- function(x) length(unique(x))
su <- function(x) sort(unique(x))
count.na <- function(x) sum(is.na(x))
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))
library(arm)

# Read in the data:
data <- read.csv(file="HOFregression.csv", as.is=TRUE)
n <- dim(data)[1]

# get number of players:
np <- lu(data[, "Name"])

# Create variables for previous year's vote:
prev <- matrix(0, n, 14)
for (i in 1:n) {
  if (data[i, "YoB"] > 1) {
    sel <- data[, "Name"] == data[i, "Name"] & data[, "YoB"] < data[i, "YoB"]
    prev[i, 1:sum(sel)] <- data[sel, "p"][order(data[sel, "Year"], decreasing=TRUE)]
  }
}
colnames(prev) <- paste0("prev", 1:14)
data <- cbind(data, prev)

# Create some per-year statistics:
#batter.py <- c("R", "H", "HR", "RBI", "SB", "BB")
#batter.avg <- c("BA", "OBP", "SLG")
#pitcher.py <- c("W", "L", "GS", "SV", "IP", "SO")
#pitcher.avg <- c("ERA", "WHIP")
#lpy <- length(c(batter.py, pitcher.py))

#per.year <- matrix(NA, n, lpy)
#colnames(per.year) <- paste0(c(batter.py, pitcher.py), "py")
#for (i in 1:lpy) per.year[, i] <- data[, c(batter.py, pitcher.py)[i]]/data[, "Yrs"]

# Add the per-year stats to the data.frame:
#data <- cbind(data, per.year)

# Add in awards data:
AllStarpy <- data[, "all.star"]/data[, "Yrs"]
data <- cbind(data, AllStarpy)

# Add indicators for 8 batting positions (DH is the baseline)
for (i in c("C", "1B", "2B", "3B", "SS", "LF", "CF", "RF")) assign(paste0("pos", i), as.numeric(data[, "position"] == i))
data <- cbind(data, posC, pos1B, pos2B, pos3B, posSS, posLF, posCF, posRF)

# Add in previous year's vote squared
prev1.squared <- data[, "prev1"]^2

# Add the mean vote percentage of the top-k (k = 1, 2, 3, 4, 5) first-year ballot players in each year:
ny <- lu(data[, "Year"])
first.ballot.crowd <- matrix(NA, ny, 5)
for (i in 1:ny) {
  for (k in 1:5) {
  	sel <- data[, "Year"] == su(data[, "Year"])[i] & data[, "YoB"] == 1
  	first.ballot.crowd[i, k] <- mean(data[sel, "p"][1:k], na.rm=TRUE)
  }
}
rownames(first.ballot.crowd) <- 1967:2014

fb <- matrix(NA, n, 5)
for (k in 1:5) fb[, k] <- first.ballot.crowd[(1:ny)[data[, "Year"] - 1966], k]
colnames(fb) <- paste0("top", 1:5)

# Append these to the data:
data <- cbind(data, prev1.squared, fb)

# one-team players:
one <- read.csv("oneteam_raw.csv", as.is=TRUE)
oneteam <- as.numeric(data[, "Name"] %in% one[1:35, 1])
data <- cbind(data, oneteam)

# add milestones:
hr500 <- as.numeric(data[, "HR"] >= 500)
h3000 <- as.numeric(data[, "H"] >= 3000)
w300 <- as.numeric(data[, "W"] >= 300)
k3000 <- as.numeric(data[, "SO"] >= 3000)
data <- cbind(data, hr500, h3000, w300, k3000)


# Try removing Bonds and Clemens first year:
data <- data[-which(data[, "Name"] %in% c("Barry Bonds", "Roger Clemens") & data[, "Year"] == 2013), ]
n <- dim(data)[1]

# Get number of ballots each year:
nb <- aggregate(data[, "NumBallots"], list(year=data[, "Year"]), median)[, 2]
nb[48] <- 569  # insert 569 for 2014 (just use number from last year)

# Set up type of model:
type <- numeric(n)
#for (i in 1:10) {
#  sel <- data[, "YoB"] == 1 & data[, "position"] == c("C", "1B", "2B", "3B", "SS", "LF", "CF", "RF", "DH", "P")[i]
#  type[sel] <- i
#}
#type[data[, "YoB"] > 1] <- 11
type[data[, "position"] != "P" & data[, "YoB"] == 1] <- 1     # batters
type[data[, "position"] == "P" & data[, "YoB"] == 1] <- 2     # pitchers
type[data[, "YoB"] > 1] <- 3                                  # returning
#type[data[, "YoB"] > 1] <- data[data[, "YoB"] > 1, "YoB"] + 1 # returning
nt <- lu(type)
var.names <- as.list(rep(NA, nt))


# R0: get base statistics for batters and pitchers:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs")
var.names[[3]] <- c("prev1")

# R1 (with per-year stats):
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs")
var.names[[3]] <- c("prev1")
# Not better at all!

# R2 with All-Stars and Awards:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young")
var.names[[3]] <- c("prev1")

# R3 with indicators for positions:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp", "posC", "pos1B", "pos2B", "pos3B", "posSS", "posLF", "posCF", "posRF")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young")
var.names[[3]] <- c("prev1")

# R4 with different models for each position:
for (i in 1:9) {
  var.names[[i]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                      "rookie", "gold.gloves", "mvp")
}
var.names[[10]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young")
var.names[[11]] <- c("prev1")


# R5: different type for each year on ballot:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young")
for (j in 3:16) var.names[[j]] <- paste0("prev", 1:min(k, j-2))



# R6: Add in crowded ballot, and quadratic term for previous year:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp", "oneteam", "hr500", "h3000")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young", "w300", "k3000")
var.names[[3]] <- c("prev1", "prev1.squared", "top3")





### Set up a vector to store all predictions:
pred <- rep(NA, n)
pred.mat <- matrix(NA, n, 5)
qbounds <- function(x) quantile(x, c(0.025, 0.250, 0.500, 0.750, 0.975))

coef <- as.list(rep(NA, nt))
lt <- length(1997:2014)
for (j in 1:nt){
  coef[[j]] <- matrix(NA, lt, length(var.names[[j]]) + 1)
  colnames(coef[[j]]) <- c("Intercept", var.names[[j]])
  rownames(coef[[j]]) <- 1997:2014
}

in.samp <- matrix(NA, lt, nt)

# Loop through years and positions (batter vs. pitcher):
for (t in 1997:2014) {
  print(t)
  for (j in 1:nt) {

    # Set up the design matrix for this type of prediction:
    if (j %in% 1:2) sel <- type == j & data[, "Year"] < t
    if (j > 2) sel <- type == j & data[, "Year"] < t & data[, "prev1"] >= 0.05
    X.mat <- as.matrix(data[sel, var.names[[j]]])

    # Scale the inputs, keeping the means and sds:
    x.mean <- apply(X.mat, 2, mean)
    x.sd <- apply(X.mat, 2, sd)
    X.scale <- X.mat
    for (i in 1:dim(X.mat)[2]) {
      if (x.sd[i] != 0) X.scale[, i] <- (X.mat[, i] - x.mean[i])/x.sd[i]
    }
    
    # Fit the model using weak priors:
    fit <- bayesglm(data[sel, "p"] ~ X.scale, weights=data[sel, "NumBallots"], family=binomial(link = "logit"), 
                    prior.mean=0, prior.scale=2.5)
    in.samp[t - 1996, j] <- sd(fit$fitted.values - data[sel, "p"])
    
    # Store the coefficients:
    coef[[j]][t - 1996, ] <- coef(fit)

    # predict this type for the year of interest:
    sel.test <- type == j & data[, "Year"] == t

    if (sum(sel.test) > 0) {
      X.mat <- as.matrix(data[sel.test, var.names[[j]]])
      if (t == 2014 & j == 3) {
      	X.mat[, "top3"] <- mean(pred[data[, "Year"] == 2014 & type %in% 1:2][1:3], na.rm=TRUE)
      }
      X.scale <- X.mat
      for (i in 1:dim(X.mat)[2]) {
        if (x.sd[i] != 0) X.scale[, i] <- (X.mat[, i] - x.mean[i])/x.sd[i]
      }
      beta <- mvrnorm(1000, mu=coef(fit), Sigma=summary(fit)$cov.scaled)
      pred[sel.test] <- expit(coef(fit)[1] + X.scale %*% matrix(coef(fit)[-1], ncol=1))
      pred.sim <- expit(beta[,1] + X.scale %*% t(beta[, -1]))
      votes.sim <- matrix(rbinom(sum(sel.test)*1000, size=nb[t - 1966], prob=pred.sim), sum(sel.test), 1000)/nb[t - 1966]
      pred.mat[sel.test, ] <- t(apply(votes.sim, 1, qbounds))
    }
  }
}


# Now, boxplots:
sel.pred <- data[, "Year"] > 1996 & data[, "Year"] < 2014
rmse <- sqrt(mean((pred[sel.pred] - data[sel.pred, "p"])^2))
resids <- data[sel.pred, "p"] - pred[sel.pred]


rmse.vec <- numeric(nt)
sel.vec <- as.list(rep(NA, nt))
for (j in 1:3) sel.vec[[j]] <- sel.pred & type == j
for (i in 1:nt) rmse.vec[i] <- sqrt(mean((pred[sel.vec[[i]]] - data[sel.vec[[i]], "p"])^2))


out.samp <- matrix(NA, lt, 3)
n.out <- matrix(NA, lt, 3)
for (i in 1:lt) {
  for (j in 1:3) {
  	sel <- data[, "Year"] == i + 1996 & type == j
  	out.samp[i, j] <- sd(pred[sel] - data[sel, "p"])
  	n.out[i, j] <- sum(sel)
  }
}

par(mfrow=c(1, 3))
for (j in 1:3) {
  plot(1997:2014, in.samp[, j], type="l", ylim=range(c(in.samp[, j], out.samp[-18, j])))
  lines(1997:2013, out.samp[-18, j], lty=2)
}


### R0: Baseline Model:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs")
var.names[[3]] <- c("prev1")
# 11.47 %
# 0.18021465 0.09683615 0.05688873

# R2: With All-Stars and Awards:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young")
var.names[[3]] <- c("prev1")
# 9.81%
# 0.14494421 0.09511777 0.05688873


# R6: Add in crowded ballot, and quadratic term for previous year:
var.names[[1]] <- c("Yrs", "G", "AB", "R", "H", "HR", "RBI", "SB", "BB", "BA", "OBP", "SLG", "drugs", "AllStarpy", 
                    "rookie", "gold.gloves", "mvp", "oneteam", "hr500", "h3000")
var.names[[2]] <- c("Yrs", "W", "L", "ERA", "WHIP", "G.1", "GS", "SV", "IP", "H.1", "HR.1", "BB.1", "SO", "drugs", 
                    "AllStarpy", "rookie", "gold.gloves", "mvp", "cy.young", "w300", "k3000")
var.names[[3]] <- c("prev1", "prev1.squared", "top3")
# 8.82%
# 0.12923726 0.09577417 0.04677941


# Excluding Bonds and Clemens 2013:

# 7.90%
# 0.11904610 0.06842750 0.04677941



cover50 <- numeric(nt)
cover95 <- numeric(nt)
for (j in 1:nt) {
  s <- sel.pred & type == j
  cover50[j] <- sum(data[s, "p"] > pred.mat[s, 2] & data[s, "p"] < pred.mat[s, 4])/sum(s)
  cover95[j] <- sum(data[s, "p"] > pred.mat[s, 1] & data[s, "p"] < pred.mat[s, 5])/sum(s)
}


z <- round(pred.mat, 3)*100
colnames(z) <- c("lower95", "lower50", "prediction", "upper50", "upper95")
cbind(z, data)[data[, "Year"] == 2014, 1:15]




# Look at residuals with nametags:
sel.big <- abs(resids) > 0.1  # select big residuals
xl <- "Predicted Percentage"
pdf(file="fig_residuals_R2.pdf", width=11, height=6.5)
plot(pred[sel.pred], resids, type="n", las=1, ylab="Actual Vote % - Predicted Vote %", yaxt="n", xlim=c(0, 1.1), xlab=xl, xaxt="n")
axis(2, at=seq(-1, 1, 0.2), labels=paste(seq(-100, 100, 20), "%", sep=""), las=1)
axis(1, at=seq(0, 1, 0.2), labels=paste0(seq(0, 100, 20), "%"))
abline(h=seq(-1, 1, 0.1), col=gray(0.8))
text(pred[sel.pred][sel.big], resids[sel.big], paste(data[sel.pred, "Name"], 
     data[sel.pred, "Year"], sep="-")[sel.big], cex=0.6, col=as.numeric(data[sel.pred, "YoB"][sel.big] == 1) + 1)
points(pred[sel.pred][!sel.big], resids[!sel.big], col=as.numeric(data[sel.pred, "YoB"][!sel.big] == 1) + 1)
legend("topright", inset=0.01, col=c(1, 2), pch=19, legend=c("Returning Player", "First Ballot"))
#title(main="M3: Including model for non-first ballot (Residuals with absolute value > 10%)")
dev.off()


z <- data.frame(data, pred)[data[, "Year"] == 2014, ]
zz <- data.frame(Name=z[, "Name"], Predicted=round(z[, "pred"], 3)*100)



top <- data.frame(election_dat[sel.pred, c("Year", "Name", "X.vote")], Predicted=pred[sel.pred], Residual=resids)[order(resids, decreasing=TRUE), ][1:5,]

bottom <- data.frame(election_dat[sel.pred, c("Year", "Name", "X.vote")], Predicted=pred[sel.pred], Residual=resids)[order(resids, decreasing=FALSE), ][1:5,]

rownames(top) <- 1:5
rownames(bottom) <- 1:5
top[, 4:5] <- round(top[, 4:5], 3)*100
bottom[, 4:5] <- round(bottom[, 4:5], 3)*100
top[, 3] <- as.numeric(substr(top[, 3], 1, nchar(top[, 3]) - 1))
bottom[, 3] <- as.numeric(substr(bottom[, 3], 1, nchar(bottom[, 3]) - 1))
colnames(top)[3] <- "Actual"
colnames(bottom)[3] <- "Actual"




sd(resids[election_dat[sel.pred, "YoB"] == 1 & election_dat[sel.pred, "position"] != "P"]) # batters
sd(resids[election_dat[sel.pred, "YoB"] == 1 & election_dat[sel.pred, "position"] == "P"]) # pitchers
sd(resids[election_dat[sel.pred, "YoB"] > 1]) # second ballot or later

sum(election_dat[sel.pred, "YoB"] == 1 & election_dat[sel.pred, "position"] != "P")
sum(election_dat[sel.pred, "YoB"] == 1 & election_dat[sel.pred, "position"] == "P")
sum(election_dat[sel.pred, "YoB"] > 1)



# Look at residuals with nametags, only for the second or later time on ballot:
sel.return <- election_dat[sel.pred, "YoB"] > 1
xl <- "Predicted Percentage"
pdf(file="fig_residuals_M3_repeats.pdf", width=11, height=6.5)
plot(pred[sel.pred][sel.return], resids[sel.return], type="n", las=1, ylab="Actual Vote % - Predicted Vote %", 
     yaxt="n", xlim=c(0, 1.1), xlab=xl, xaxt="n")
axis(2, at=seq(-1, 1, 0.2), labels=paste(seq(-100, 100, 20), "%", sep=""), las=1)
axis(1, at=seq(0, 1, 0.2), labels=paste0(seq(0, 100, 20), "%"))
abline(h=seq(-1, 1, 0.1), col=gray(0.8))
text(pred[sel.pred][sel.return], resids[sel.return], paste(election_dat[sel.pred, "Name"][sel.return], 
     election_dat[sel.pred, "Year"][sel.return], sep="-"), cex=0.6)
title(main="M3: Everything (Returning Players to Ballot)")
dev.off()


this.year <- data[, "Year"] == 2014
type <- ifelse(data[this.year, "position"] == "P" & data[this.year, "YoB"] == 1, 2, 
               ifelse(data[this.year, "position"] != "P" & data[this.year, "YoB"] == 1, 1, 3))
rmse.vec <- c(13.5, 11.3, 4.8)
Lower <- round(pmax(pred[this.year]*100 - 2*rmse.vec[type], 0), 1)
Upper <- round(pmin(pred[this.year]*100 + 2*rmse.vec[type], 100), 1)

zzz <- data.frame(zz, Lower, Upper)




n.bins <- 20
rmse.bin <- numeric(n.bins)
qq <- quantile(prev.fit$fitted, seq(0, 1, length=n.bins + 1))
for (i in 1:n.bins) {
  sel.bin <- prev.fit$fitted > qq[i] & prev.fit$fitted <= qq[i + 1]
  rmse.bin[i] <- sd(prev.fit$fitted[sel.bin] - V[sel][sel.bin])
}

rmse.match <- numeric(17)
for (i in 1:17) rmse.match[i] <- max(which(x2014[i, 1] > qq))
x2014[1:17, 2] <- x2014[1:17, 1] - 1.96*rmse.bin[rmse.match]
x2014[1:17, 3] <- x2014[1:17, 1] + 1.96*rmse.bin[rmse.match]

zzz <- data.frame(Name=zz[, 1], Prediction=round(x2014[, 1]*100, 1), Lower=round(x2014[, 2]*100, 1), Upper=round(x2014[, 3]*100, 1))
zzz <- zzz[order(zzz[, 2], decreasing=TRUE), ]
rownames(zzz) <- 1:36

predict(prev.fit, newdata=data.frame(v.prev, trend, bal.year, top3firstballot), type="response", se.fit=T)

# Fill in Jose Rijo's prediction:
pred[is.na(pred) & election_dat[, "Name"] == "Jose Rijo"] <- V[election_dat[, "Year"] == 2001 & 
     election_dat[, "Name"] == "Jose Rijo"]

# Now, boxplots:
sel.pred <- election_dat[, "Year"] > 1996 & election_dat[, "Year"] < 2014
rmse <- sqrt(mean((pred[sel.pred] - V[sel.pred])^2))
resids <- V[sel.pred] - pred[sel.pred]


sqrt(mean((pred[sel.pred] - V[sel.pred])^2))

# M0: baseline batting, baseline pitching, previous vote: 0.1101
# M1: batting and pitching + AS.percent, Drugs, GG, and baseline previous vote: 0.993
# M2: batting and pitching everything, baseline previous vote: 0.984
# M3: batting and pitching everything, regression for previous vote: 0.0949








# Make a plot of mean difference from one year to the next based on year on ballot:
vdiff <- matrix(NA, 14, 3)
for (i in 1:14) {
  sel <- data[, "YoB"] == i + 1 & data[, "Year"] != 2014
  vdiff[i, 1] <- sum(sel)
  vdiff[i, 2] <- mean(data[sel, "p"] - data[sel, "prev1"])
  vdiff[i, 3] <- sd(data[sel, "p"] - data[sel, "prev1"])
}

lower <- vdiff[, 2] - 2*vdiff[, 3]/sqrt(vdiff[, 1])
upper <- vdiff[, 2] + 2*vdiff[, 3]/sqrt(vdiff[, 1])

par(mfrow=c(1, 1))
plot(2:15, vdiff[, 2], ylim=range(c(lower, upper)), las=1)
for (i in 1:14) lines(c(i+1, i+1), c(lower[i], upper[i]))
abline(h = 0, lty=2)









