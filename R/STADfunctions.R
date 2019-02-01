## ------------------------------------------------------------------- ##
##  R functions accompanying the file "STADcode.R".
##  Authors: Ugofilippo Basellini & Giancarlo Camarda
## ------------------------------------------------------------------- ##

## Function for smoothing mortality rates from deaths and exposures with 
## monotonic constraint. Age range can be greater than data, and monotonicity
## is needed for extrapolation
lmx_smooth <- function(age,y,e,w,ncoef){
  ## length of age range
  x <- age
  mx <- length(x)
  ## length of data
  my <- length(y)
  ## extend data if needed
  if (mx > my){
    yA <- c(y, rep(99, mx-my))
    eA <- c(e, rep(99, mx-my))
    wA <- c(w, rep(0, mx-my))
  }else{
    yA <- y
    eA <- e
    wA <- w
  }
  ## starting value: glm with weights
  options(warn=-1)
  fit0 <- glm(round(yA) ~ x, offset=log(eA),
    weights = wA, family=poisson()) 
  options(warn=0)
  # including monotonicity
  Dmon <- diff(diag(ncoef), diff=1)
  wmon <- rep(0, nrow(Dmon))
  Wmon <- diag(wmon)
  kappa <- 10^4
  eps <- 0.1
  epsvec <- rep(eps, ncoef - 1)
  ## selecting lambda
  lambdas <- 10^seq(-1,6,0.5)
  nl <- length(lambdas)
  BICs <- numeric(nl)
  COEFs <- matrix(0, ncoef, nl)
  ## difference matrix
  D <- diff(diag(ncoef), diff=2)
  tDD <- t(D)%*%D
  ## loop over the lambdas
  for(l in 1:nl){
    P <- lambdas[l] * tDD
    ## starting linear predictor
    eta <- fit0$coef[1] + fit0$coef[2]*x
    for(it in 1:100){
      ## monotonicity penalty
      Pmon0 <- t(Dmon) %*% Wmon %*% Dmon
      Pmon <- kappa * Pmon0
      vmon <- kappa * t(Dmon) %*% (wmon * epsvec)
      ## regression part
      mu <- exp(eta)*eA
      z <- wA * ((yA - mu)/mu + eta)
      z[is.na(z)] <- 0
      w <- c(wA*mu)
      BtWB <- t(B) %*% (w * B)
      BtWBpP <- BtWB + P + Pmon
      BtWz <- t(B) %*% (w * z)
      a <- solve(BtWBpP, BtWz + vmon)
      eta.old <- eta
      eta <- B %*% a
      ## update monotonicity weights
      diff.a <- diff(a) >= eps
      wmon <- rep(1, nrow(Dmon))
      wmon[diff.a] <- 0
      Wmon <- diag(wmon)
      ## convergence critera  
      deta <- max(abs((eta.old - eta)/eta.old))
      # cat(deta, "\n")
      if(it>=4 & deta<=10^-5) break   
    }
    ## deviance
    y1 <- yA
    mu1 <- mu
    y1[y1 == 0] <- 10^(-5)
    mu1[mu1 == 0] <- 10^(-5)
    dev <- 2 * sum(wA*(y1 * log(y1/mu1)),na.rm = T)
    ## ED
    H <- solve(BtWBpP, BtWB)
    h <- diag(H)
    ed <- sum(h)
    ## BIC
    BICs[l] <- dev + log(sum(wA)) * ed
    ## fitted coef
    COEFs[,l] <- a
  }
  ## evaluate mortality at finer grid ~ hazard function
  pmin <- which.min(BICs)
  coef.hat <- COEFs[,pmin]
  out <- list(coef=coef.hat)
  return(out)
}

## function for aligning fx given a shifting parameter
fx_shift <- function(age,fx,shift,ndx=25,deg=3){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(age)[1],4)
  ages.add <- 0
  ## augment basis
  if (shift > 0) ages.add = ceiling(shift + 5)
  if (shift < 0) ages.add = ceiling(abs(shift - 5))
  ## weights
  w <- c(rep(0,ages.add/delta),rep(1,ms),rep(0,ages.add/delta))
  ## augmented data
  xA <- c(rev(seq(from=xs[1]-delta, by=-delta,length=ages.add/delta)), xs, 
          seq(from=xs[ms]+delta, by=delta, length=ages.add/delta))
  fA <- c(rep(999, ages.add/delta), fx, rep(999, ages.add/delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  nb <- ndx+deg
  ## new B-splines bases
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## smoothing + extrapolation
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  P <- 10^-5 * tDD
  betasA <- solve(t(BA)%*%(w*BA)+P, t(BA)%*%(w*log(fA)))
  ## shifting
  ws.i <- xs - shift
  ## evaluating B on shifted x
  Bshift <- MortSmooth_bbase(ws.i, xl=xmin, xr=xmax, ndx, deg)
  ## shifted density
  fshift <- exp(Bshift %*% betasA)
  out <- fshift
  if (shift == 0) out <- fx
  return(out)
}

## function for obtaining coefficients of standard
coeff_stand <- function(age,fx,ages.add.l=30,ages.add.r=20,ndx=30,deg=3){
  ## Extrapolate age range to left and right + derive coefficients
  ## length of age range
  xs <- age
  delta <- round(diff(xs)[1],4)
  ms <- length(xs)
  ## weights
  w <- c(rep(0,ages.add.l/delta),rep(1,ms),rep(0,ages.add.r/delta))
  ## augmented data
  xA <- c(rev(seq(from=xs[1]-delta, by=-delta,length=ages.add.l/delta)), 
          xs, 
          seq(from=xs[ms]+delta, by=delta, length=ages.add.r/delta))
  fA <- c(rep(999, ages.add.l/delta), fx, 
          rep(999, ages.add.r/delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  nb <- ndx+deg
  ## B-splines bases standard expanded
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## simple smoothing + extrapolation
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  P <- 10^-5 * tDD
  ## finding betas for the standard distribution expanded
  betasA <- solve(t(BA)%*%(w*BA)+P, t(BA)%*%(w*log(fA)))
  ## effective dimension
  tBWB <- t(BA)%*%(w*BA)
  tBWBpP <- tBWB + P 
  H <- solve(tBWBpP, tBWB)
  h <- diag(H)
  ed <- sum(h)
  ## return
  out <- list(betasA=betasA,ed=ed)
  return(out)
}

## function to compute mx from dx
mx_from_dx <- function(dx,ax=NULL){
  ## dimension of dx
  m <- length(dx)
  ## template vectors
  lx <- Lx <- mx <- rep(NA,m)
  ## set the radix of life table
  lx[1] <- sum(dx,na.rm = T)
  ## compute l(x+1)=l(x)-d(x) 
  for (i in 2:m){
    lx[i] <- lx[i-1] - dx[i-1]
  }
  ## set ax = 1/2
  if (is.null(ax)) ax <- rep(1/2,m)
  ## compute Lx = l(x+1) + ax*dx
  Lx[-m] <- lx[-1]+ax[-m]*dx[-m]
  ## compute mx
  mx <- dx/Lx
  ## return mx value
  return(mx)
}

## Function to compute dx from mx 
dx_from_mx <- function(age,mx){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(xs)[1],4)
  ## identity matrix and C matrix
  I <- diag(ms)
  C <- lower.tri(I, diag = TRUE)
  C[C==1] <- -delta 
  ## compute dx
  dx <- mx * exp(C%*%mx)
  return(dx)
}

## optimization function for bL and bU
MLE_obj_FUN <- function(par,x,xA,
                        Mstand, shat, xlo, xup,coeff.stand, Dx, Ex, 
                        xmin, xmax, ndx, deg){
  ## starting b
  bL <- par[1]
  bU <- par[2]
  ## segment a linear transformation function
  ## below the mode
  aL <- Mstand - bL * (shat+Mstand)
  wL <- aL + bL*xlo
  ## above the mode
  aU <- Mstand - bU * (shat+Mstand)
  wU <- aU + bU*xup
  ## unique transformation function
  wb <- c(wL, wU)
  wb <- sort(wb)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xmin,xmax,ndx=ndx,deg=deg)
  ## transformed density
  dwb <- as.vector(exp(Bwb%*%coeff.stand))
  dwb <- dwb[xA%in%x]
  if (sum(dwb) > 1) dwb <- dwb/sum(dwb)
  ## hazard
  eta <- log(mx_from_dx(dx=dwb))
  mu <- exp(eta)
  ## minimise minus the Log-Likelihood (maximise the LL)
  Lij <- -sum(Dx * eta[1:length(Dx)]-Ex*mu[1:length(Ex)], na.rm = T)
  return(Lij)
}

## function for constructing a classic (& rather general) lifetable
## source: Gaincarlo Camarda webpage
## https://sites.google.com/site/carlogiovannicamarda/r-stuff/life-expectancy-confidence-interval
lifetable.mx <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

