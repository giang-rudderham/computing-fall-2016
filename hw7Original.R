longley
install.packages("gmp")
library("gmp")
#1a
coeff <- function(x, y) {
    x <- as.bigq(x)
    y <- as.bigq(y)
    coeff <- as.double(solve(t(x) %*% x) %*% t(x) %*% y)
    coeff
}

x <- as.matrix(cbind(1, longley[ , -7]))
y <- longley[ , 7]

exact <- coeff(x, y)
exact

#1b
cholesky <- function(x, y, center = FALSE) {
    if (center) {
        colmean <- colMeans(x)
        x <- sweep(x, 2, colmean)
        x <- as.matrix(cbind(1, x))
    }
    if (!center) {
        x <- as.matrix(cbind(1, x))
    }
    lT <- chol(t(x) %*% x)
    z <- forwardsolve(t(lT), t(x) %*% y)
    cholesky <- backsolve(lT, z)
    cholesky
}

x <- as.matrix(longley[ , -7])
y <- longley[ , 7]

cholesky1 <- cholesky(x, y)
cholesky2 <- cholesky(x, y, center = TRUE)


#1c
x <- as.matrix(cbind(1, longley[ , -7]))
y <- longley[ , 7]
qr <- lm.fit(x, y)$coefficients

exact - cholesky1
exact - cholesky2
exact - qr

#Trial code
###########################
fit <- lsfit(longley[ , -7], longley[ , 7])$coef #use built-in lsfit function.

y <- longley[ , 7]

x <- as.matrix(cbind(1, longley[ , -7]))
xT <- t(x)
xTy <- xT %*% y
xTx <- xT %*% x
coeff1 <- as.double(solve(xTx) %*% xTy) #compute beta-hat using floating point data.

u <- as.bigq(x)
v <- as.bigq(y)
coeff2 <- as.double(solve(t(u) %*% u) %*% t(u) %*% v) #compute beta-hat using
#arbitrary precision rational data.

coeff1 - coeff2 #There's a difference between using floating point data
#and arbitrary precision rational data

compute1 <- as.double(solve(xTx) %*% xTy) #convert result to double precision numbers.
compute2 <- solve(xTx) %*% xTy #not convert
compute1 - compute2 # All 0

compute1 - fit #There's a difference between using lsfit function and computing by hand.

amatrix <- array(1:9, c(3,3))
mean <- colMeans(amatrix)
sweep(amatrix, 2, mean) #subtract each column by column means

#Check function cholesky when center = T
colmean <- colMeans(x)
newx <- sweep(x, 2, colmean)
newx <- as.matrix(cbind(1, newx))
newx
lT <- chol(t(newx) %*% newx)
z <- forwardsolve(t(lT), t(newx) %*% y)
cholesky3 <- backsolve(lT, z)
cholesky2 <- cholesky(x, y, center = T)
cholesky3 - cholesky2 #All 0.

#2a
####################
set.seed(1)
y <- rnorm(1000)
a <- 0.5
lldense <- function(y, a) {
  dim <- length(y)
  c <- diag(dim)
  c[cbind(2 : dim, 1: (dim - 1))] <- a
  c[cbind(1 : (dim - 1), 2 : dim)] <- a
  u <- chol(c)
  diagonals <- diag(u)
  logdet <- sum(log(diagonals))
  
  l <- t(u)
  I <- diag(dim)
  inverse <- matrix(NA, nrow = dim, ncol = dim)
  
  for (i in 1: dim) {
    tmp <- forwardsolve(l, I[ , i])
    inverse[ , i] <- backsolve(u, tmp)
  }
  
  -logdet - (1/2) * t(y) %*% inverse %*% y
}

lldense(y, a)  

#Trial code
###################
##Find logdetC
dim <- length(y)
c <- diag(dim)
c[cbind(2 : dim, 1: (dim - 1))] <- a
c[cbind(1 : (dim - 1), 2 : dim)] <- a
c
u <- chol(c)
diagonals <- diag(u)
diagonals
sum(log(diagonals))

##Find inverse of C
l <- t(u)
l
I <- diag(dim)
I
inverse <- matrix(NA, nrow = dim, ncol = dim)

for (i in 1: dim) {
    tmp <- forwardsolve(l, I[ , i])
    inverse[ , i] <- backsolve(u, tmp)
}
inverse 

-sum(log(diagonals)) - (1/2) * t(y) %*% inverse %*% y
c %*% inverse - I #close
solve(c) - inverse #close. Inverse has similar print representation as solve(c)

#b
########################
library(Matrix)
llsparse <- function(y, a) {
  dim <- length(y)
  diags <- list(rep(1, dim), rep(a, dim))
  c <- bandSparse(dim, k = c(0 : 1), diag = diags, symm = TRUE)
  u <- chol(c)
  logdet <- sum(log(diag(u)))
  
  -logdet - (1/2) * t(y) %*% solve(c) %*% y
}

llsparse(y, a)  
#c
##########################
system.time(lldense(y, a))
system.time(llsparse(y, a))

#Generate correlated normal data
####################################
R = matrix(cbind(1, .50, 0,  .50, 1, .5,  0, .5, 1), nrow = 3)
R
U = t(chol(R))
U
nvars = dim(U)[1]
numobs = 100000
set.seed(1)
random.normal = matrix(rnorm(nvars*numobs, 0, 1), nrow = nvars, ncol = numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
orig.raw = as.data.frame(t(random.normal))
names(raw) = c("response","predictor1","predictor2")
cor(raw)
plot(head(raw, 100))
plot(head(orig.raw,100))
