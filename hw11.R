num.sample <- 10000
n <- 9

# Stardard Normal
################
## Generate sample
set.seed(2)
sample <- matrix(rnorm(num.sample * n), nrow = num.sample)

trimean <- function(data) {
  data <- sort(data)
  (data[3] + 2 * data[5] + data[7]) / 4
}

trimean.est <- apply(sample, MARGIN = 1, trimean)
var.norm <- mean(trimean.est ^ 2)
se.norm <- sd(trimean.est ^ 2) / sqrt(length(trimean.est))

## Variance reduction
new.trimean <- function(data) {
  x.hat <- sum(data) / n
  s.hat <- sqrt( sum((data - x.hat) ^ 2) / (n - 1) )
  c <- ( data - x.hat ) / s.hat
  trimean(c)
}

new.trimean.est <- apply(sample, MARGIN = 1, new.trimean)
var.new.norm <- 1 / n + mean(new.trimean.est ^ 2)  
se.new.norm <- sd(new.trimean.est ^ 2) / sqrt( length(new.trimean.est) )

# t(4)
##################
## Generate sample
df <- 4
set.seed(87)
sample <- matrix(rt(num.sample * n, df = df), nrow = num.sample)
trimean <- function(data) {
  data <- sort(data)
  (data[3] + 2 * data[5] + data[7]) / 4
}

trimean.est <- apply(sample, MARGIN = 1, trimean)
var.t4 <- mean(trimean.est ^ 2)
se.t4 <- sd(trimean.est ^ 2) / sqrt(length(trimean.est))

## Variance reduction
set.seed(87)
V <- matrix(sqrt(rchisq(num.sample * n, df = df) / df), nrow = num.sample)
data <- cbind(sample, V)

new.trimean <- function(data) {
  x <- data[1 : n]
  v <- data[(n + 1): length(data)]
  x.hat <- sum(x * (v ^ 2)) / sum(v ^ 2)
  s.hat <- sqrt( sum(((x - x.hat) ^ 2) * (v ^ 2)) / (n - 1) )
  c <- (x - x.hat) / s.hat
  1 / sum(v ^ 2) + trimean(c) ^ 2
}

new.trimean.est <- apply(data, MARGIN = 1, new.trimean)
var.new.t4 <- mean(new.trimean.est)
se.new.t4 <- sd(new.trimean.est) / sqrt(length(new.trimean.est))

# t(10)
##################
## Generate sample
df <- 10
set.seed(56)
sample <- matrix(rt(num.sample * n, df = df), nrow = num.sample)
trimean <- function(data) {
  data <- sort(data)
  (data[3] + 2 * data[5] + data[7]) / 4
}

trimean.est <- apply(sample, MARGIN = 1, trimean)
var.t10 <- mean(trimean.est ^ 2)
se.t10 <- sd(trimean.est ^ 2) / sqrt(length(trimean.est))

## Variance reduction
set.seed(56)
V <- matrix(sqrt(rchisq(num.sample * n, df = df) / df), nrow = num.sample)
data <- cbind(sample, V)

new.trimean <- function(data) {
  x <- data[1 : n]
  v <- data[(n + 1): length(data)]
  x.hat <- sum(x * (v ^ 2)) / sum(v ^ 2)
  s.hat <- sqrt( sum(((x - x.hat) ^ 2) * (v ^ 2)) / (n - 1) )
  c <- (x - x.hat) / s.hat
  1 / sum(v ^ 2) + trimean(c) ^ 2
}

new.trimean.est <- apply(data, MARGIN = 1, new.trimean)
var.new.t10 <- mean(new.trimean.est)
se.new.t10 <- sd(new.trimean.est) / sqrt(length(new.trimean.est))

# xtable
###############
row1 <- c(var.norm, se.norm, var.new.norm, se.new.norm)
row2 <- c(var.t4, se.t4, var.new.t4, se.new.t4)
row3 <- c(var.t10, se.t10, var.new.t10, se.new.t10)
table <- data.frame(rbind(row1, row2, row3))
colnames(table) <- c("Crude Variance", "Crude SE", "Est. Var", "Est. SE")
rownames(table) <- c("Standard normal", "t(4)", "t(10)")
xtable(table, caption = "Crude Monte Carlo and estimated variances and SEs", label = "tab", digits = 4)
