fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
optim(c(-1.2,1), fr)
(res <- optim(c(-1.2,1), fr, grr, method = "BFGS"))
optimHess(res$par, fr, grr)
optim(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)

#1 a)
#############
library(MASS)
data(leuk)
leuk.pos <- leuk[leuk$ag == "present", ]

# Log-likelihood
loglik <- function(theta, y) {
  n <- nrow(y)
  t <- y[ , 3]
  x <- log(y[ , 1] / 10000)
  - (- n * log(theta[1]) + theta[2] * sum(x) - sum(t * exp( theta[2] * x)) / theta[1])
}

# OLS and Initial values ?????
theta0 <- c(50, .5)
loglik(theta0, leuk.pos)
loglik(c(56.85, 0.482), leuk.pos)

# MLE
myMLE <- optim(theta0, loglik, y = leuk.pos)
MLE.e <- myMLE$par
MLE.e

# Plot loglik function
theta1 <- seq(30, 80, by = 0.01)
theta2 <- seq(0, 1, by = 0.01)
res1 <- rep(NA, length(theta1))
res2 <- rep(NA, length(theta2))

# minimization w.r.t theta1
for( i in 1 : length(theta1)) {
  res1[i] <- loglik(c(theta1[i], MLE.e[2]), leuk.pos)
}

plot(x = theta1, res1, type = 'l', ylab = "log likelihood", 
     main = "log likelihood with respect to theta1")
abline(v = theta1[which.min(res1)], col = 'red')

# minimization w.r.t. theta2
for( i in 1 : length(theta2)) {
  res2[i] <- loglik(c(MLE.e[1], theta2[i]), leuk.pos)
}

plot(x = theta2, res2, type = 'l', ylab = "log likelihood", 
     main = "log likelihood with respect to theta2")
abline(v = theta2[which.min(res2)], col = 'red')

### Problem 1 b)
post.den <- function(theta2, y) {
  n <- nrow(y)
  t <- y[, 3]
  x <- log(y[, 1] / 10000)
  
  factorial(n + 1) * ((0.1) ^ 6) * exp(theta2 * sum(x)) *
    (0.001 + sum(t * exp(theta2 * x))) ^ (- n - 1) 
}

myres <- c()
theta2 <- seq(- 1.5, 1.5, by = 0.01)
for(i in 1 : length(theta2)){
  myres[i] <- post.den(theta2[i], leuk.pos)
}

plot(theta2, myres, type = 'l', 
     main = "the unnormalized marginal posterior density for theta2")

### Problem 1 c)
myv <- function(theta2, y) {
  myres2 <- c()
  for( i in 1 : length(theta2)) {
    myres2[i] <- sqrt(post.den(theta2[i], y)) * theta2[i]
  }
  v1 <- myres2[which.min(myres2)]
  v2 <- myres2[which.max(myres2)]
  c(v1, v2)
}
theta2 <- seq(-1.5, 1.5, by=0.01)
v <- myv(theta2, leuk.pos)

ustar <- sqrt(myres[which.max(myres)])
vstar1 <- v[1]
vstar2 <- v[2]

myru <- function(n, u, v1, v2) {
  count <- 0
  x <- c()
  
  while(count <= n) {
    n1 <- 1.5 * n
    x0 <- rep(NA, n1)
    u <- runif(n1, 0, u)
    v <- runif(n1, v1, v2)
    x <- v / u
    for(i in 1 : n1) {
      x0[i] <- post.den(x[i], leuk.pos)
    }
    ind <- (u ^ 2 <= x0)
    ind <- ind[- which(is.na(ind))]
    x <- c(x, x[ind])
    count <- count + length(ind)
  }
  
  x[1 : n]
}

mysample <- myru(10000, ustar, vstar1, vstar2)
plot(density(mysample))
mysample[mysample < 1.5 & mysample > -1.5]
plot(density(mysample[mysample < 1.5 & mysample > -1.5]))

