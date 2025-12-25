############################################################
# Example: Minimal demonstration of GameTalker
############################################################
source("main.R")

# Time points (same as original experiment)
times <- 0:10

# Toy abundance trajectories 
a1 <- read.csv("FT-LN18.csv")[1:11,-1]
a2 <- read.csv("ST-LN18.csv")[1:11,-1]

b1 <- read.csv("FT-3T3.csv")[1:11,-1]
b2 <- read.csv("ST-3T3.csv")[1:11,-1]


a <- (a1 + a2)/2
b <- (b1 + b2)/2

X <-  matrix(as.numeric(unlist(c(a))),nrow=11)
Y <-  matrix(rep(as.numeric(a[1,]),11),nrow=11,byrow=T)
Z <- X - Y + 1

X2 <-  matrix(as.numeric(unlist(c(b))),nrow=11)
Y2 <-  matrix(rep(as.numeric(b[1,]),11),nrow=11,byrow=T)
Z2 <- X2 - Y2 + 1

# Pack data 
data_obs <- list(
  X = X,
  Y = Y,
  Z = Z,
  times = times
)

# Initial guess (same structure as original)
par.init <- c(6, 0.32)

# Fit intrinsic growth for species X
fit_X <- optim_BFGS.ind(
  par.curve = par.init,
  z = data_obs$X,
  times = data_obs$times
)

# Estimated intrinsic curve
mu_X <- LC_get_mu1(fit_X$par, data_obs$times)

# Time-varying drivers 
inter_X <- splinefun(times, data_obs$X[,6])
inter_Y <- splinefun(times, data_obs$Y[,6])
inter_Z <- splinefun(times, data_obs$Z[,6])


# Pairwise interaction parameters
para1 <- c(1, 10, 1,  0.15, 1)
para2 <- c(1, 10, 1, -0.10, 1)
para3 <- c(1, 10, 1,  0.05, 1)

# Run GameTalker pairwise model
predict.function <- function(object, newdata, ...) {
  list(y = object(newdata))
}

res <- ode.sovle12(
  para1 = para1,
  para2 = para2,
  para3 = para3,
  times = times,
  inter1 = inter_X,
  inter2 = inter_Y,
  inter3 = inter_Z
)

# Visualize simulated dynamics
matplot(
  res$nt,
  res$NG,
  type = "l",
  lwd = 2,
  lty = 1,
  col = c("darkred", "darkgreen", "darkblue"),
  xlab = "Time",
  ylab = "Abundance",
  main = "GameTalker: toy inference example"
)

legend(
  "topleft",
  legend = c("Species X", "Species Y", "Species Z"),
  col = c("darkred", "darkgreen", "darkblue"),
  lwd = 2,
  bty = "n"
)





