#####################################################################
##### Function for MCMC samples using LS KLE algo & MUR #############
#####################################################################


source('all-fcts.R')
# require(Rcpp)



### Function for drawing posterior samples using LS and MUR without hyperparameter updates (nu & ell):

########### Nonparametric functions estimation ##############
LS.KLE.MUR.function <- function(y,x,N1,p,M,mcmc,brn,thin,nu,l,sig.in,xi.in,tau.in,
                             tau.fix,sig.fix,xi.fix,sseed,verbose,return.plot,tol){
  # y:Response variable; x: vector to form design matrix \Psi (n X N+1)
  # N1: the number of knots of 1st subdomain
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # nu.in, l.in, tau.in, sig.in, xi.in : initial values (supplied by user or the default values)
  # tau.fix, sig.fix, xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on nu,l,xi,tau,sig and fhat with posterior mean, 95% CI of fhat
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  n <- length(y)
  N <- M*N1
  y <- y[order(x)]
  x <- sort(x)
  delta <- (max(x)-min(x))/(N-1)
  my_knots <- seq(min(x),max(x),by=delta)
  X <- fctv(x,my_knots,M,N1)
  
  if (missing(return.plot))
    return.plot <- TRUE
  
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  
  if (missing(verbose))
    verbose <- TRUE
  
  if (missing(tol))
    tol <- 1e-8
  
  if (missing(l))
    l <- l_est(nu,c(my_knots[1],my_knots[length(my_knots)]),0.05)
  
  # prior covariance K:
  K <- kMat(my_knots,nu,l)
  # prior precision:
  K_inv <- tinv(K)
  
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  if (!missing(tau.fix)){
    tau.in <- tau.fix
    ## Hassan added two lines
    Gamma <- tau.in*K
    GX <- Gamma%*%t(X)
  }
  
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  
  if (missing(xi.fix) && missing(xi.in))
    xi.in <- mvrnorm(1,rep(0,N),K)
  
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  em <- mcmc+brn
  ef <- mcmc/thin
  xi_sam <- matrix(NA,N,ef)
  tau_sam <- rep(NA,ef)
  sig_sam <- rep(NA,ef)
  fhat_sam <- matrix(NA,n,ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em){
    # sampling Xi:
    y_tilde <- y 
    if (missing(xi.fix)){
      fprior <- LS.KLE(my_knots,N1,p,M,nu=nu,l=l,tau=tau,tol=tol,sseedLS=i)
      # as.vector(samp.WC(my_knots,nu_out,l_out,tau))
      if (missing(tau.fix)){
        Gamma <- tau*K
        GX <- Gamma%*%t(X)
        invXGX <- chol2inv(chol(X%*%GX+sig*diag(n)))
        xi_out <- fprior+GX%*%invXGX%*%(y-X%*%fprior)
      }
      else {
        invXGX <- chol2inv(chol(X%*%GX+sig*diag(n)))
        xi_out <- fprior+GX%*%invXGX%*%(y-X%*%fprior)
        # xi_out <- LS.KLE_MUR(my_knots,f=fprior,A=X,y,N1,p,M,nu,l,tausq=tau,sqrt(sig))
      }
    } else {
      xi_out <- xi_in
    }
    
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    set.seed(123)
    # sampling \sigma^2:
    if (missing(sig.fix))
      sig <- 1/rgamma(1,shape = n/2,rate = sum(y_star^2)/2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N/2, rate = (t(xi_out)%*%K_inv%*%xi_out)/2)
    
    # storing MCMC samples:
    if (i > brn && i%%thin == 0){
      xi_sam[,(i-brn)/thin] <- xi_out
      sig_sam[(i-brn)/thin] <- sig
      tau_sam[(i-brn)/thin] <- tau
      fhat_sam[,(i-brn)/thin] <- Xxi
    }
    
    if (i%%1000 == 0 && verbose){
      print(i)
    }
    
    # renewing the initial value:
    xi_in <- xi_out
  }; tm <- proc.time()-ptm
  
  
  fmean <- rowMeans(fhat_sam)
  qnt <- apply(fhat_sam,1,function(x) quantile(x,c(0.025,0.975),na.rm = TRUE))
  f_low <- qnt[1,]
  f_upp <- qnt[2,]
  ub <- max(f_low,f_upp,fmean,y)
  lb <- min(f_low,f_upp,fmean,y)
  
  if (return.plot){
    par(mfrow=c(1,1))
    par(mar=c(2.1,2.1,2.1,1.1)) # adapt margins
    plot(x,y,pch='*',lwd=2,lty=1,col='black',
         ylim=range(ub,lb),xlab='',ylab='')
    polygon(c(x,rev(x)),y=c(f_low, rev(f_upp)),border=F,col='gray')
    lines(x,fmean,type='l',lty=2,lwd=2)
    points(x,y,pch='*')
  }
  
  return(list("time"=tm,"xi_sam"=xi_sam,"sig_sam"=sig_sam,"tau_sam"=tau_sam,
              "fhat_sam"=fhat_sam,"fmean"=fmean,"f_low"=f_low,"f_upp"=f_upp))
}
###########





### Function for drawing posterior samples using LS  and MUR with hyperparameter updates (nu & ell):
########### Nonparametric functions estimation ##############
LS.KLE.MUR.hyp <- function(y,x,N1,p,M,mcmc,brn,thin,nu.in,l.in,sig.in,xi.in,tau.in,
                        tau.fix,sig.fix,xi.fix,sseed,verbose,return.plot,tol){
  # y:Response variable; x: vector to form design matrix \Psi (n X N+1)
  # N1: the number of knots of 1st subdomain
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # nu.in, l.in, tau.in, sig.in, xi.in : initial values (supplied by user or the default values)
  # tau.fix, sig.fix, xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  
  # OUTPUT: Posterior samples on nu,l,xi,tau,sig and fhat with posterior mean, 95% CI of fhat
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  n <- length(y)
  N <- M*N1
  y <- y[order(x)]
  x <- sort(x)
  delta <- (max(x)-min(x))/(N-1)
  my_knots <- seq(min(x),max(x),by=delta)
  X <- fctv(x,my_knots,M,N1)
  
  if (missing(return.plot))
    return.plot <- TRUE  
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(nu.in))
    nu.in <- 0.75
  if (missing(l.in))
    l.in <- l_est(nu.in,c(my_knots[1],my_knots[length(my_knots)]),0.05)
  
  # prior covariance K:
  # K <- kMat(my_knots,nu.in,l.in)
  # prior precision:
  # K_inv <- tinv(K)
  
  if (missing(mcmc))
    mcmc <- 5000
  if (missing(brn))
    brn <- 1000
  if (missing(thin))
    thin <- 1
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (!missing(xi.fix))
    xi.in <- xi.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  if (missing(xi.fix) && missing(xi.in))
    xi.in <- rep(0,N)
  
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  nu_in <- nu.in
  l_in <- l.in

  em <- mcmc+brn
  ef <- mcmc/thin
  xi_sam <- matrix(NA,N,ef)
  tau_sam <- rep(NA,ef)
  sig_sam <- rep(NA,ef)
  nu_sam <- rep(NA,ef)
  ell_sam <- rep(NA,ef)
  fhat_sam <- matrix(NA,n,ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em){
    # sampling from \nu and \ell
    MH.out <- nu.MH2(nu_in,l_in,tau,xi_in,my_knots,range.nu=c(0.5,2.5),range.l=c(0.1,1),sd.nu=0.05,sd.l=0.05,
                    seed=i)
    nu_out <- MH.out$nu
    l_out <- MH.out$l
    L_inv <- MH.out$L_inv
    
    # sampling Xi:
    y_tilde <- y 
    if (missing(xi.fix)){
      fprior <- LS.KLE(my_knots,N1,p,M,nu=nu_out,l=l_out,tau=tau,tol=tol,sseedLS=i)
      # fprior <- as.vector(samp.WC(my_knots,nu_out,l_out,tau,sseedWC=i)) ## wood and chen FFT
      Gamma <- tau*kMat(my_knots,nu_out,l_out)
      GX <- Gamma%*%t(X)
      invXGX <- chol2inv(chol(X%*%GX+sig*diag(n)))
      xi_out <- as.vector(fprior+GX%*%invXGX%*%(y-X%*%fprior))
      # xi_out <- LS.KLE_MUR(my_knots,f=fprior,A=X,y,N1,p,M,nu,l,tausq=tau,sqrt(sig))
    } else {
      xi_out <- xi_in
    }
    
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    set.seed(123)
    # sampling \sigma^2:
    if (missing(sig.fix))
      sig <- 1/rgamma(1,shape = n/2,rate = sum(y_star^2)/2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N/2, rate = ((sum((t(L_inv)%*%xi_out)^2))/2))
    
    
    # storing MCMC samples:
    if (i > brn && i%%thin == 0){
      xi_sam[,(i-brn)/thin] <- xi_out
      sig_sam[(i-brn)/thin] <- sig
      tau_sam[(i-brn)/thin] <- tau
      nu_sam[(i-brn)/thin] <- nu_out
      ell_sam[(i-brn)/thin] <- l_out
      fhat_sam[,(i-brn)/thin] <- Xxi
    }
    
    if (i%% 1000==0 && verbose){
      print(i)
    }
    
    # renewing the initial value:
    xi_in <- xi_out
    nu_in <- nu_out
    l_in <- l_out
    
  } 
  tm <- proc.time()-ptm
  
  
  fmean <- rowMeans(fhat_sam)
  qnt <- apply(fhat_sam,1,function(x) quantile(x,c(0.025,0.975),na.rm = TRUE))
  f_low <- qnt[1,]
  f_upp <- qnt[2,]
  ub <- max(f_low,f_upp,fmean,y)
  lb <- min(f_low,f_upp,fmean,y)
  
  if (return.plot){
    par(mfrow=c(1,1))
    par(mar=c(2.1,2.1,2.1,1.1)) # adapt margins
    plot(x,y,pch='*',lwd=2,lty=1,col='black',
         ylim=range(ub,lb),xlab='',ylab='')
    polygon(c(x,rev(x)),y=c(f_low, rev(f_upp)),border=F,col='gray')
    lines(x,fmean,type='l',lty=2,lwd=2)
    points(x,y,pch='*')
  }
  
  return(list("time"=tm,"xi_sam"=xi_sam,"sig_sam"=sig_sam,"tau_sam"=tau_sam,
              "nu_sam"=nu_sam,"ell_sam"=ell_sam,
              "fhat_sam"=fhat_sam,"fmean"=fmean,"f_low"=f_low,"f_upp"=f_upp))
}



## end
