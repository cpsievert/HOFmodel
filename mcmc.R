#explore results
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
  abline(h = beta_int.true[i], col=2, lwd=2, lty=2)
}

par(mfrow=c(2,2))
for (i in 1:4){
  plot(1:n.reps, beta_slope_keep[,i], col=1, type="l", ylim=c(-0.45, 0.45))
  abline(h = beta_slope.true[i], col=2, lwd=2, lty=2)
}

plot(1:n.reps, alpha_keep[,1])
lu(alpha_keep[,1])
acc.prop.alpha <- apply(alpha_keep[(adapt + 1):n.reps,], 2, lu)/n.reps
hist(acc.prop.alpha)


samp <- sample(n.players, 6)
par(mfrow=c(2,3))
for (i in 1:6){
  plot(1:n.reps, alpha_keep[,samp[i]], type="l", ylim=c(-2, 2), main=paste("Player", samp[i], repeats[samp[i]]))
  abline(h = alpha.true[samp[i]], col=2, lwd=2, lty=2)
}

pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  pct[i,] <- quantile(alpha_keep[(adapt + 1):n.reps, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}

cover.50 <- sum(pct[,2] < alpha.true & pct[,4] > alpha.true)/n.players
cover.95 <- sum(pct[,1] < alpha.true & pct[,5] > alpha.true)/n.players
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
  abline(h = gamma.true[samp[i]], col=2, lwd=2, lty=2)
}

pct <- matrix(NA, n.players, 5)
for (i in 1:n.players){
  pct[i,] <- quantile(gamma_keep[(adapt + 1):n.reps, i], p=c(0.025, 0.250, 0.500, 0.750, 0.975))
}

cover.50 <- sum(pct[,2] < gamma.true & pct[,4] > gamma.true)/n.players
cover.95 <- sum(pct[,1] < gamma.true & pct[,5] > gamma.true)/n.players  
cover.50
cover.95  
