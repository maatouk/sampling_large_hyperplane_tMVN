## Numerical examples illustrated in COSTA's article
## Comparison between naive MUR and LS.KLE (computational running time)
## Comparison between LS.KLE and Wood and Chan FFT (computational running time)
## Numerical performance of LS.KLE with size of 10M
## Nonparametric fct estimation with synthetic data with fixed hyperparameters
## Two real data studies (age-income and fossil) with fixed hyperparameters


setwd("~/Documents/Recherche/Article_high-dimensional-hyperplaneMVN_2021/High dimension code/Official_codes/GitHub/sampling_large_hyperplane_tMVN")
source('all_models.R')
source('all-fcts.R')


## Choice of a numerical example
## put 'yes' between quotations

KLEvsLS.KLE = ''
WCvsLS.KLEN = ''
WCvsLS.KLEnu = ''
LS.KLE10M = ''
synthetic = 'yes'
ageincome = ''
fossil = ''




if (KLEvsLS.KLE == 'yes') {
  ####################################################
  ####### naive KLE versus Large scale KLE ###########
  ####################################################
  par(mar = c(3.1, 3.1, 1.1, 1.1)) # adapt margins
  nu <- 2.5 # smoothness parameter Matern Kernel (MK)
  l <- l_est(nu, c(0, 1), 0.05) # length-scale parameter MK
  N <- seq(from = 100, to = 1000, length = 10) # size of the MVN
  N1 <- 50
  M <- N / N1
  p <- 30
  trial <- 25
  timeKLE <- matrix(data = NA, nrow = trial, ncol = length(N))
  timeLS.KLE <- matrix(data = NA, nrow = trial, ncol = length(N))
  for (j in 1 : trial) {
    print(j)
    for (i in 1 : length(N)) {
      u <- seq(from = 0, to = 1, length = N[i])
      timeLS.KLE[j, i] <- system.time(LS.KLE(u = u, N1 = N1, p = p, M = M[i],
                                             nu = nu, l = l, tau = 1, tol = 1e-12, sseedLS = i))[3]
      timeKLE[j, i] <- system.time(KLE(u, p, nu, l))[3]
    }
  } 
  averageLS.KLE <- colMeans(timeLS.KLE)
  averageKLE <- colMeans(timeKLE)
  plot(N, averageLS.KLE, type = 'l', lwd = 2, ylim = range(averageLS.KLE, averageKLE), xlab = '', ylab = '')
  title(xlab = 'Dimension', ylab = 'Average Time (s)', line = 2)
  lines(N, averageKLE, type = 'l', lty = 2, lwd = 2)
  legend(200, 0.9, c('KLE', 'LS.KLE'),
         lty = c(2, 1), cex = 1,
         text.font = 2, box.lty = 0, lwd = 2)
  ##############################################################
}


if (WCvsLS.KLEN == 'yes') {
  ####################################################
  ####### Comparison between samp.WC & LS.KLE ########
  ####################################################
  ## The Matérn smoothness parameter nu is fixed
  par(mar = c(3.1, 3.1, 1.4, 1.1)) # adapt margins
  nu <- 0.5 # smoothness parameter Matern Kernel (MK)
  l <- 0.4 # length-scale parameter MK
  N <- seq(from = 1000, to = 50000, length = 11) # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N / N1
  p <- 30 # truncation expansion parameter
  trial <- 25
  timeWC <- matrix(data = NA, nrow = trial, ncol = length(N))
  timeLS <- matrix(data = NA, nrow = trial, ncol = length(N))
  for (j in 1 : trial) {
    print(j)
    for (i in 1 : length(N)) {
      u <- seq(from = 0, to = 1, length = N[i])
      timeLS[j, i] <- system.time(LS.KLE(u, N1 = N1, p, M = M[i], nu = nu, l = l, tau = 1, tol = 1e-12, sseedLS = i))[3]
      timeWC[j, i] <- system.time(samp.WC(knot = u, nu = nu, l = l, tausq = 1, sseedWC = i))[3]
    }
  } 
  averageLS <- colMeans(timeLS)
  averageWC <- colMeans(timeWC)
  # smooth <- paste("smoothness~parameter~", expression(nu), "~", nu)
  plot(N, averageLS, type = 'l', lwd = 2,
       ylim = range(averageLS, averageWC), xlab = '', ylab = '')
  # mtext(parse(text = smooth))
  title(xlab = 'Dimension', ylab = 'Average Time (s)', line = 2)
  mtext(text = c('smoothness parameter'), side = 3, line = 0.6, cex = 1)
  mtext(bquote(nu == .(nu)), side = 3, line = -0.1, cex = 1)
  lines(N, averageWC, type = 'l', lty = 2, lwd = 2)
  legend(1000, 0.33, c('FFT.WC', 'LS.KLE'),
         lty = c(2, 1), cex = 1,
         text.font = 2, box.lty = 0, lwd = 2, bg = 'transparent')
  ####################################################
}

if (WCvsLS.KLEnu == 'yes') {
  ### running time as a fct of the smoothness parameter nu 
  par(mar = c(3.1, 3.1, 1.1, 1.1)) # adapt margins
  nu <- seq(from = 0.5, to = 1.5, length = 10) # smoothness parameter
  N <- 10000 # sizes of MVN
  N1 <- 100 # size of 1st subdomain
  M <- N / N1 # nb of subdomains
  p <- 30 # expansion retains parameter
  trial <- 25 # nb of replicates
  l <- 0.4 # length-scale parameter MK
  timeWC <- matrix(data = NA, nrow = length(nu), ncol = trial)
  timeLS <- matrix(data = NA, nrow = length(nu), ncol = trial)
  u <- seq(from = 0, to = 1, length = N)
  for (i in 1 : length(nu)) {
    print(i)
    for (j in 1 : trial) {
      timeLS[i, j] <- system.time(LS.KLE(u, N1 = N1, p, M = M, nu = nu[i], l = l, tau = 1, tol = 1e-12, sseedLS = j))[3]
      timeWC[i, j] <- system.time(samp.WC(u, nu = nu[i], l = l, tausq = 1,sseedWC = j))[3]
    }
  } 
  averageLS <- rowMeans(timeLS)
  averageWC <- rowMeans(timeWC)
  plot(nu, averageLS, type = 'l', lwd = 2,
       ylim = range(averageLS, averageWC),
       xlab = '', ylab = '')
  title(ylab = 'Average Time (s)', line = 2)
  mtext(text = expression("smoothness parameter" ~ nu),
        side = 1, line = 2)
  mtext(text =  paste("Dimension N = ", N), side = 3, line = 0.2, cex = 1)
  lines(nu, averageWC, type = 'l', lty = 2, lwd = 2)
  legend(0.5, 0.4, c('FFT.WC', 'LS.KLE'),
         lty = c(2, 1), cex = 1, text.font = 2, 
         box.lty = 0, lwd = 2, bg = 'transparent')
  #####################################################
}


if (LS.KLE10M == 'yes') {
  ####################################################
  ##### Performance Fast.LS with N=10,000,000 #########
  ####################################################
  ## second test N1=100
  par(mar = c(3.1, 3.1, 1.1, 1.1)) # adapt margins
  nu <- 1.5 # smoothness parameter Matern Kernel (MK)
  l <- 0.2 # length-scale parameter MK
  N1 <- 100 # size of 1st subdomain
  p <- 30
  M <- seq(from = 1000, to = 100000, length = 10) # nbs of blocks
  N <- M * N1 # sizes of the MVN
  trial <- 5
  timeLS <- matrix(data = NA, nrow = trial, ncol = length(M))
  for (i in 1 : trial) {
    print(i)
    for (j in 1 : length(M)) {
      u <- seq(from = 0, to = 1, length = N[j])
      timeLS[i, j] <- system.time(LS.KLE(u, N1, p, M = M[j], nu = nu, l, tau = 1, tol = 1e-12, sseedLS = j))[3]
    }
  }
  averageLS <- colMeans(timeLS)
  plot(N, averageLS, type = 'l', lwd = 2, ylim = range(averageLS),
       xlab = '', ylab = '')
  title(xlab = 'Dimension', ylab = 'Average Time (s)', line = 2)
  legend(1000000, 1.2, c('LS.KLE'),
         lty = 1, cex = 1, text.font = 2, 
         box.lty = 0, lwd = 2, bg = 'transparent')
  ####################################################
}



if (synthetic == 'yes') {
  ####################################################
  ############ runing time nonparametric #############
  ################# fct estimation ###################
  ####################################################
  #### Large scale KLE and MUR
  N1 <- 300
  p <- 30 # KLE truncation parameter
  M <- 5 # nb of subdomains
  nbsim <- 5000 # nb of simulation
  trial <- 5
  nu <- 2.5 # smoothness parameter Matern Kernel MK 
  l <- l_est(nu, range = c(0, 1), 0.2) # length-scale
  u <- seq(from = 0, to = 1, length = (M * N1))
  tol <- 1e-12
  # delta <- 1/(M * N1 - 1)
  ### data
  ntot <- 50 # nb of total data
  ntr <- 30 # nb training data
  nte <- ntot-ntr # nb test data
  sigN <- 0.05 # sd noise
  f <- function(x) {
    x * cos(2 * x)
  }
  ## split data
  set.seed(123)
  xtot <- runif(n = ntot, min = 0, max = 1)
  ytot <- f(xtot) + rnorm(n = ntot, mean = 0, sd = sigN)
  timeLS <- rep(NA, trial) # run time Large scale
  timeMUR <- rep(NA, trial) # run time naive MUR
  print("MCMC replicates:")
  for (Q in 1 : trial) {
    print(paste0('replicate = ', Q))
    set.seed(2 * Q)
    ind <- sample.int(ntot, ntr)
    xtr <- xtot[ind]
    ytr <- ytot[ind]
    xte <- xtot[-ind]
    yte <- ytot[-ind]
    ytrue.te <- f(xte)
    A <- fctv(xtr, u, M, N1)
    timeLS[Q] <- system.time(LS.KLE_MUR_v(nbsim, u, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol))[3]
    timeMUR[Q] <- system.time(MUR(nbsim, u, A, y = ytr, nu, l, sigN, tol))[3]
  }
  ## Illustration
  par(mar = c(3.1, 3.1, 1.9, 1.1)) # adapt margins
  post_samp_LS <- LS.KLE_MUR_v(nbsim, u, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol)
  post_samp_MUR <- MUR(nbsim, u, A, y = ytr, nu, l, sigN, tol)
  t <- seq(from = 0, to = 1, length = 500)
  Y_LS <- fctv(t, u, M, N1) %*% post_samp_LS
  Y_MUR <- fctv(t, u, M, N1) %*% post_samp_MUR
  tmp_LS <- apply(Y_LS, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_LS <- tmp_LS[1, ]
  fmean_LS <- tmp_LS[2, ]
  f_upp_LS <- tmp_LS[3, ]
  tmp_MUR <- apply(Y_MUR, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_MUR <- tmp_MUR[1, ]
  fmean_MUR <- tmp_MUR[2, ]
  f_upp_MUR <- tmp_MUR[3, ]
  ## Illustration: Large-scale KLE_MUR
  plot(t, f(t), type = 'l', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_LS, f_upp_LS, f(u), ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeLS), 1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'x', ylab = 'y(x)', line = 2)
  polygon(c(t, rev(t)), y = c(f_low_LS, rev(f_upp_LS)), border = F, col = 'gray')
  lines(t, f(t), lwd = 2, lty = 1, col = 'black')
  lines(t, fmean_LS, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  legend(0.1, 0, c('true function', 'posterior mean'), lwd = 2,
         lty = c(1, 2), bg = 'transparent', box.lty = 0)
  ## Illustration: Naive Matheron's update rule
  plot(t, f(t), type = 'l', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_MUR, f_upp_MUR, f(u), ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeMUR), 1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'x', ylab = 'y(x)', line = 2)
  polygon(c(t, rev(t)), y = c(f_low_MUR, rev(f_upp_MUR)),border = F, col = 'gray')
  lines(t, f(t), lwd = 2, lty = 1, col = 'black')
  lines(t, fmean_MUR, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  legend(0.1, 0, c('true function', 'posterior mean'), lwd = 2,
         lty = c(1, 2), bg = 'transparent', box.lty = 0)
  ####################################################
}



if (ageincome == 'yes') {
  ####################################################
  ######## Age income real application ###############
  ####################################################
  library(SemiPar) #income real data semi parametric regression 
  data(age.income)
  attach(age.income)
  xtot <- age # input
  ytot <- log.income # total data
  ntot <- length(xtot) # nb of total data
  ntr <- floor(ntot * 0.8) # nb training data
  nte <- ntot - ntr # nb test data
  N1 <- 150 # size of the 1st subdomain
  p <- 30 # KLE truncation parameter
  M <- 10 # nb of subdomains
  nbsim <- 5000 # nb of simulation
  trial <- 25 # nb of replicates
  nu <- 2.5 # smoothness kernel parameter
  l <- 30 # length-scale parameter
  tol <- 1e-12
  sigN <- 1 # noise sd
  my_knots <- seq(from = min(xtot), to = max(xtot), length = (M*N1))
  ## split data
  timeLS <- rep(NA, trial) # run time Large scale
  timeMUR <- rep(NA, trial) # run time naive MUR
  for (Q in 1 : trial) {
    print(paste0('replicate = ', Q))
    set.seed(2 * Q)
    ind <- sample.int(ntot, ntr)
    xtr <- xtot[ind]
    ytr <- ytot[ind]
    xte <- xtot[-ind]
    yte <- ytot[-ind]
    A <- fctv(xtr, my_knots, M, N1)
    timeLS[Q] <- system.time(LS.KLE_MUR_v(nbsim, my_knots, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol))[3]
    timeMUR[Q] <- system.time(MUR(nbsim, my_knots, A, y = ytr, nu, l, sigN, tol))[3]
  }
  ## Illustration
  par(mar = c(3.1, 3.1, 1.9, 1.1)) # adapt margins
  t <- seq(from = min(xtot), to = max(xtot), length = 100)
  post_samp_LS <- LS.KLE_MUR_v(nbsim, u = my_knots, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol)
  post_samp_MUR <- MUR(nbsim, u = my_knots, A, y = ytr, nu, l, sigN, tol)
  Y_LS <- fctv(t, u = my_knots, M, N1) %*% post_samp_LS
  Y_MUR <- fctv(t, u = my_knots, M, N1) %*% post_samp_MUR
  tmp_LS <- apply(Y_LS, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_LS <- tmp_LS[1, ]
  fmean_LS <- tmp_LS[2, ]
  f_upp_LS <- tmp_LS[3, ]
  tmp_MUR <- apply(Y_MUR, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_MUR <- tmp_MUR[1, ]
  fmean_MUR <- tmp_MUR[2, ] 
  f_upp_MUR <- tmp_MUR[3, ]
  ## large-scale KLE and MUR
  plot(xtr, ytr, pch = '*', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_LS, f_upp_LS,ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeLS),1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'age', ylab = 'log income', line = 2)
  polygon(c(t, rev(t)), y = c(f_low_LS, rev(f_upp_LS)), border = F, col = 'gray')
  lines(t, fmean_LS, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  ## Naive Matheron's update rule
  plot(xtr, ytr, pch = '*', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_MUR, f_upp_MUR, ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeMUR),1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'age', ylab = 'log income', line = 2)
  polygon(c(t, rev(t)), y = c(f_low_MUR, rev(f_upp_MUR)), border = F, col = 'gray')
  lines(t, fmean_MUR, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  detach(age.income)
  ####################################################
}






if (fossil == 'yes') {
  ####################################################
  ######## Fossil data real application ###############
  ####################################################
  library(SemiPar) # semi parametric regression 
  data(fossil)
  attach(fossil)
  xtot <- age # input
  ytot <- strontium.ratio # total data
  ntot <- length(xtot) # nb of total data
  ntr <- floor(ntot * 0.8) # nb training data
  nte <- ntot - ntr # nb test data
  N1 <- 150 # size of the 1st subdomain
  p <- 30 # KLE truncation parameter
  M <- 10 # nb of subdomains
  nbsim <- 5000 # nb of simulation
  trial <- 25 # nb of replicates
  nu <- 3.5 # smoothness kernel parameter
  l <- l_est(nu, range(xtot), 0.5) # length-scale parameter
  tol <- 1e-12
  sigN <- sd(ytot) # standard deviation
  my_knots <- seq(from = min(xtot), to = max(xtot), length = (M * N1))
  ## split data
  timeLS <- rep(NA, trial) # run time Large scale
  timeMUR <- rep(NA, trial) # run time naive MUR
  for (Q in 1 : trial) {
    print(paste0('replicate = ', Q))
    set.seed(2 * Q)
    ind <- sample.int(ntot, ntr)
    xtr <- xtot[ind]
    ytr <- ytot[ind]
    xte <- xtot[-ind]
    yte <- ytot[-ind]
    A <- fctv(xtr, u = my_knots, M, N1)
    timeLS[Q] <- system.time(LS.KLE_MUR_v(nbsim, u = my_knots, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol))[3]
    timeMUR[Q] <- system.time(MUR(nbsim, u = my_knots, A, y = ytr, nu, l, sigN, tol))[3]
  }
  
  ## Illustration
  par(mar = c(3.1, 3.1, 1.9, 1.1)) # adapt margins
  post_samp_LS <- LS.KLE_MUR_v(nbsim, u = my_knots, A, y = ytr, N1, p, M, nu, l, sigN, tausq = 1, tol)
  post_samp_MUR <- MUR(nbsim, u = my_knots, A, y = ytr, nu, l, sigN, tol)
  t <- seq(from = min(xtot), to = max(xtot), length = 100)
  Y_LS <- fctv(t, u = my_knots, M, N1) %*% post_samp_LS
  Y_MUR <- fctv(t, u = my_knots, M, N1) %*% post_samp_MUR
  tmp_LS <- apply(Y_LS, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_LS <- tmp_LS[1, ]
  fmean_LS <- tmp_LS[2, ]
  f_upp_LS <- tmp_LS[3, ]
  tmp_MUR <- apply(Y_MUR, 1, quantile, probs = c(0.025, 0.5, 0.975))
  f_low_MUR <- tmp_MUR[1, ]
  fmean_MUR <- tmp_MUR[2, ]
  f_upp_MUR <- tmp_MUR[3, ]
  ## Large-scale KLE_MUR
  plot(xtr, ytr, pch = '*', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_LS, f_upp_LS, ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeLS), 1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'age', ylab = 'strontium ratio', line = 2)
  polygon(c(t, rev(t)), y = c(f_low_LS, rev(f_upp_LS)), border = F, col = 'gray')
  lines(t, fmean_LS, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  ## Naive Matheron's update rule
  plot(xtr, ytr, pch = '*', lwd = 2, lty = 1, col = 'black',
       ylim = range(f_low_MUR, f_upp_MUR, ytr), xlab = '', ylab = '')
  mtext(text =  paste("Average Time (s) = ", round(mean(timeMUR), 1)), side = 3, line = 0.8, cex = 1)
  mtext(text =  paste("Dimension N = ", M * N1), side = 3, line = 0.1, cex = 1)
  title(xlab = 'age', ylab = 'strontium ratio', line = 2)
  polygon(c(t, rev(t)),y = c(f_low_MUR, rev(f_upp_MUR)), border = F, col = 'gray')
  lines(t, fmean_MUR, type = 'l', lty = 2, lwd = 2)
  points(xtr, ytr, pch = '*')
  detach(fossil)
  ####################################################
}



## end
