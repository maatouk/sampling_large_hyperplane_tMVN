library(rSPDE) # matern cov fct
library(mvtnorm)
library(FastGP)
library(Rfast)
library(MASS)


#####################################################
## Matern family cov fct with \nu smooth para & \ell length-scale para
#####################################################
k <- function(h, nu, l) {
  matern.covariance(h, sqrt(2 * nu)/l, nu = nu, sigma = 1)
}

## Matern kernel covariance matrix
kMat <- function(knot, nu, l) {
  k(outer(knot, knot, '-'), nu, l)
}


#####################################################
######## Maatouk & Bay2017 Basis functions ##########
#####################################################
h <- function(x) {
  ifelse(x >= -1 & x <= 1, 1 - abs(x), 0)
}
hi <- function(x, u, i) {
  delta <- (max(u) - min(u)) / (length(u) - 1)
  h((x - u[i]) / delta)
}
#####################################################

#####################################################
########### function of design matrix ###############
#####################################################
fctv <- function(x, u, M, N1) {
  n <- length(x)
  v <- matrix(NA, nrow = n, ncol = (M * N1))
  for (j in 1 : (M * N1)) {
    v[, j] = hi(x, u, i = j)
  }
  return(v)
}
#####################################################z



#####################################################
########### Naive KLE one sample path ###############
#####################################################
KLE <- function(u, p, nu, l, tol) {
  N <- length(u)
  if (missing(tol)) {
    tol <- 1e-8
  }
  Gamma <- kMat(u, nu, l) + tol * diag(N)
  eig <- eigen(Gamma)
  value <- eig$values[1 : p]
  vector <- eig$vectors[, 1 : p]
  eta <- rnorm(p)
  return(as.vector(vector %*% (sqrt(value) * eta)))
}
#####################################################



#####################################################
####### Naive KLE more than one sample path #########
#####################################################
KLE_v <- function(nbsim, u, p, nu, l, tol) {
  N <- length(u)
  if (missing(tol)) {
    tol <- 1e-8
  }
  Gamma <- kMat(u, nu, l) + tol * diag(N)
  eig <- eigen(Gamma)
  value <- eig$values[1 : p]
  vector <- eig$vectors[, 1 : p]
  eta <- matrnorm(p, nbsim)
  return(vector %*% (sqrt(value) * eta))
}
#####################################################





#####################################################
######### naive Matheron's update rule ##############
#####################################################
MUR <- function(nbsim, u, A, y, nu, l, sigN, tol) {
  n <- length(y)
  N <- length(u)
  Gamma <- kMat(u, nu, l)
  f <- KLE_v(nbsim, u, p, nu, l, tol)
  # f <- t(mvtnorm::rmvnorm(nbsim, mean = rep(0, N), sigma = Gamma, method = 'eigen'))
  GA <- Gamma %*% t(A)
  return(f + GA %*% chol2inv(chol(A %*% GA + sigN^2 * diag(n))) %*% (y - A %*% f))
}
#####################################################


#####################################################
####### function for length-scale estimating ########
#####################################################
# function for uniroot:
fl <- function(l, para) { 
  # para[1] = x, para[2] = y and para[3] = nu of MK : Matern kernel function
  # para[4] = pre-specified value of the correlation
  a <- k(abs(para[1] - para[2]), para[3], l)
  return(a - para[4])
}
#####################################################




l_est <- function(nu, range, val) {
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para <- c(range[1], range[2], nu, val)
  rl <- uniroot(f = fl, interval = c(0.000001, 100000), para)
  return(rl$root)
}
#####################################################


#####################################################
################ LS KLE one sample ##################
#####################################################
LS.KLE <- function(u, N1, p, M, nu, l, tau, tol, sseedLS) {
  if (missing(sseedLS))
      set.seed(sseedLS)
  if (missing(tol)) 
    tol <- 1e-8
  if (missing(M))
    M <- 1
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M * N1)
    stop("Error: the length of the vector \'u\' must be \'M\' times \'N1\'")
      
  u1 <- u[1 : N1]
  u2 <- u[(N1 + 1) : (2 * N1)]
  Gamma11 <- tau * kMat(u1, nu, l) + tol * diag(N1)
      
  if (M == 1) 
    return(as.vector(mvtnorm::rmvnorm(n = 1, mean = rep(0, N1), sigma = Gamma11, method = 'chol')))
  else 
    Gamma12 <- tau * k(outer(u1, u2, '-'), nu, l)
      
  eig11 <- eigen(Gamma11)
  value11 <- eig11$values[1 : p]
  vector11 <- eig11$vectors[, 1 : p]
  K12 <- (((t(vector11) %*% (Gamma12)) %*%
             (vector11)) / 
            sqrt(tcrossprod(value11)))
  L12 <- t(chol(diag(p) - crossprod(K12)))
  eta <- matrnorm(M, p)
  etaT <- matrix(NA, nrow = M, ncol = p)
  f <- matrix(NA, nrow = M, ncol = N1)
  f[1,] <- vector11 %*% (sqrt(value11) * eta[1, ])  
  etaT[1, ] <- eta[1, ]
  for (i in 2 : M) {
    etaT[i, ] <- t(K12) %*% etaT[(i - 1), ] + L12 %*% eta[i, ]
    f[i, ] <- vector11 %*% (sqrt(value11) * etaT[i, ])  
  }
  return(as.vector(t(f)))
}
#####################################################





#####################################################
########## LS KLE more than one sample ##############
#####################################################

LS.KLE_v <- function(nbsim, u, N1, p, M, nu, l, tau, tol) {
  if (missing(nbsim)) 
    nbsim <- 10
  if (missing(tol)) 
    tol <- 1e-8
  if (missing(M)) 
    M <- 1
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M * N1)
    stop("The length of the vector u must be M times N1")
  
  u1 <- u[1 : N1]
  u2 <- u[(N1 + 1) : (2 * N1)]
  # delta <- 1 / (M * (N1 - 1))
  Gamma11 <- tau * kMat(u1, nu, l) + tol * diag(N1)
  
  if (M == 1) 
    return(as.vector(mvtnorm::rmvnorm(n = 1, mean = rep(0, N1), sigma = Gamma11, method = 'chol')))
  else 
    Gamma12 <- tau * k(outer(u1, u2, '-'), nu, l)
  
  eig11 <- eigen(Gamma11)
  value11 <- eig11$values[1 : p]
  vector11 <- eig11$vectors[, 1 : p]
  K12 <- (((t(vector11) %*% (Gamma12)) %*% 
             (vector11)) / sqrt(tcrossprod(value11)))
  L12 <- t(chol(diag(p) - crossprod(K12)))
  eta <- matrnorm(p, nbsim)
  etaT <- list()
  etaT[[1]] <- eta
  f <- list()
  f[[1]] <- vector11 %*% (sqrt(value11) * eta)
  for (i in 2 : M) {
    eta <- matrnorm(p, nbsim)
    etaT[[i]] <- t(K12) %*% etaT[[i - 1]] + L12 %*% eta
    f[[i]] <- vector11 %*% (sqrt(value11) * etaT[[i]])
  }
  return(do.call(rbind, f))
}
#####################################################



#####################################################
################### LS KLE MUR ######################
#####################################################
## one posterior sample path
LS.KLE_MUR <- function(u, f, A, y, N1, p, M, nu, l, tausq, sigN) {
  n <- length(y)
  # f <- LS.KLE(u, N1, p, M, nu, l, tol)
  Gamma <- tausq * kMat(u, nu, l)
  GA <- Gamma %*% t(A)
  invAGA <- chol2inv(cholesky(A %*% GA + sigN^2 * diag(n)))
  return(f + GA %*% invAGA %*% (y - A %*% f))
}
#####################################################



#####################################################
###### LS KLE MUR more than one sample path #########
#####################################################
LS.KLE_MUR_v <- function(nbsim, u, A, y, N1, p, M, nu, l, sigN, tausq, tol) {
  n <- length(y)
  f <- LS.KLE_v(nbsim, u, N1, p, M, nu, l, tau = tausq, tol)
  Gamma <- tausq * kMat(u, nu, l)
  GA <- Gamma %*% t(A)
  return(f + GA %*% chol2inv(cholesky(A %*% GA + sigN^2 * diag(n))) %*% (y - A %*% f))
}
#####################################################



#####################################################
################# samp WC function ##################
#####################################################
### Functions related to Wood and Chan algorithm of drawing samples
## Order of the circulant matrix:
## minimum value of g and m so that G can be embedded into C
min_g <- function(knot) {
  N <- length(knot)
  g <- ceiling(log(2 * N, 2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n-1), the condition is modified!
  return("g" = g)
}

## forming the circulant matrix:
circulant <- function(x) {
  n <- length(x)
  mat <- matrix(0, nrow = n, ncol = n)
  for (j in 1 : n) {
    mat[j, ] <- c(x[-(1 : (n + 1 - j))], x[1 : (n + 1 - j)])
  }
  return(mat)
}
#####################################################



## Function for forming the vector of circulant matrix:
circ_vec <- function(knot, g, nu, l, tausq) {
  delta_N <- 1/(length(knot)-1)
  m <- 2**g
  cj <- integer()
  for (j in 1 : m) {
    if (j <= (m/2))
      cj[j] <- (j - 1) * delta_N
    else
      cj[j] <- (m - (j - 1)) * delta_N
  }
  x <- (tausq * MK(cj, 0, l, nu))
  return(x)
}
#####################################################


#####################################################
## Function for finding a g such that C is nnd:
## without forming the circulant matrix and without computing eigen values:
#####################################################
C.eval <- function(knot, g, nu, l, tausq) {
  vec <- circ_vec(knot, g, nu, l, tausq)
  val <- fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
  # vector is by construction is symmetric!
  ev <- min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}
#####################################################


nnd_C <- function(knot, g, nu, l, tausq) {
  C.vec <- C.eval(knot, g, nu, l, tausq)$vec
  eval <- C.eval(knot, g, nu, l, tausq)$min.eig.val
  if (eval > 0)
    return(list("cj" = C.vec, "g" = g))
  else {
    g <- g + 1
    nnd_C(knot, g, nu, l, tausq)
  }
}
#####################################################



## computing the eigen values of C using FFT:
eigval <- function(knot, nu, l, tausq) {
  g <- min_g(knot)
  c.j <- nnd_C(knot, g, nu, l, tausq)$cj
  lambda <- Re(fft(c.j))
  if (min(lambda) > 0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}
#####################################################





#####################################################
### Samples drawn using Wood and Chan Algorithm #####
#####################################################
samp.WC <- function(knot, nu, l, tausq, sseedWC = 1) {
  N <- length(knot)
  lambda <- eigval(knot, nu, l, tausq)
  m <- length(lambda)
  samp.vec <- rep(0, N)
  a <- rep(0, m)
  set.seed(sseedWC)
  a[1] <- sqrt(lambda[1]) * rnorm(1) / sqrt(m)
  a[(m / 2) + 1] <- sqrt(lambda[(m / 2) + 1]) * rnorm(1) / sqrt(m)
  i <- sqrt(as.complex(-1))
  for (j in 2 : (m/2)) {
    uj <- rnorm(1)
    vj <- rnorm(1)
    a[j] <- (sqrt(lambda[j]) * (uj + i * vj)) / (sqrt(2 * m))
    a[m + 2 - j] <- (sqrt(lambda[j]) * (uj - i * vj)) / (sqrt(2 * m))
  }
  samp <- fft(a)
  samp.vec <- Re(samp[1 : N])
  return(samp.vec)
}
#####################################################




## end
