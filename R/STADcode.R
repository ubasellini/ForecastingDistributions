## ------------------------------------------------------------------- ##
##  R code to estimate and forecast adult mortality using the 
##  STAD model described in: Basellini U. and Camarda C.G. (2019), 
##  "Modelling and forecasting adult age-at-death distributions", 
##  Population Studies
##  
##  Authors: Ugofilippo Basellini & Carlo Giovanni Camarda
##  
##  sessionInfo() details:
##  
##  R version 3.4.2 (2017-09-28)
##  Platform: x86_64-apple-darwin15.6.0 (64-bit)
##  Running under: macOS High Sierra 10.13.5
##  
##  locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
##  attached base packages:
##  splines  stats  graphics  grDevices  utils  datasets 
##  methods  base     
## 
##  other attached packages:
##  colorspace_1.3-2  vars_1.5-2  lmtest_0.9-35  urca_1.3-0  strucchange_1.5-1
##  sandwich_2.4-0  zoo_1.8-0  MASS_7.3-47  demography_1.20  forecast_8.3  
##  MortalitySmooth_2.3.4  lattice_0.20-35  svcm_0.1.2  Matrix_1.2-11
##
## ------------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library(MortalitySmooth)
library(demography)
library(vars)
library(colorspace)

## load STAD functions
source("STADfunctions.R")

## -- READING THE DATA -------------

## Read mortality data   
cou <- "SWE"   ## choose any country from HMD
username <- username   ## set your HMD credentials
password <- password   ## set your HMD credentials
FullData <- hmd.mx(cou, username, password)

## original ages
xo <- 30:110
mo <- length(xo)
## expanded ages
x <- 30:120
m <- length(x)
## expanded ages at finer grid
delta <- 0.1
xs <- seq(min(x), max(x), delta)
ms <- length(xs)
## years
years <- 1980:2014
n <- length(years)
## colors
coly <- rainbow_hcl(n)

## subset age and year-specific data
FittingData <- extract.years(FullData, years=years)
FittingData <- extract.ages(FittingData, ages=xo)

## final starting data: actual death counts and exposures
## (here you can change to male/total or use your own data)
E <- FittingData$pop$female
MX.act <- FittingData$rate$female
Z <- E*MX.act

## weigths to avoid issue with zero
WEI <- matrix(1,mo,n)
WEI[E==0] <- 0

## log death rates
lMX <- log(MX.act)

## B-splines parameters
xl <- min(x)
xr <- max(x)
xmin <- round(xl - 0.01 * (xr - xl),3)
xmax <- round(xr + 0.01 * (xr - xl),3)
ndx <- floor(m/3)
deg <- 3
nbx <- ndx+deg

## B-splines bases
B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
Bs <- MortSmooth_bbase(xs, xmin, xmax, ndx, deg)

## -- FITTING -------------

## smooth mortality for each y, forcing monotonicity at all ages
lMX.smooth <- matrix(NA, ms, n)
for(i in 1:n){
  lmx.smooth.coeff <- lmx_smooth(age=x,y=Z[,i],e=E[,i],w=WEI[,i],ncoef = nbx)$coef
  lMX.smooth[,i] <- Bs %*% lmx.smooth.coeff
}
matplot(xs, lMX.smooth, lty=1, t="l", col=coly,xlab="Age",
        ylab="log-mortality",main="Log-Mortality (Smooth)",cex.lab=1.25)

## compute density from smooth mortality rates
FXs <- matrix(0, nrow=ms, ncol=n)
for(i in 1:n){
  FXs[,i] <- dx_from_mx(age=xs,mx=exp(lMX.smooth[,i]))
}
matplot(xs, FXs, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Observed Distributions (Smooth)",cex.lab=1.25)

## compute modal age at death
M <- xs[apply(FXs, 2, which.max)]
plot(years, M, t="o", lwd=2, pch=16, main="Modal Age at Death")

## ESTIMATE THE STANDARD DISTRIBUTION
## 1. Alignment procedure: get shifting parameter
s <- M - M[1]

## derive aligned distributions
FXs.align <- matrix(0, nrow=ms, ncol=n)
for(i in 1:n){
  FXs.align[,i] <- fx_shift(age=xs,fx=FXs[,i],shift=-s[i],ndx = ndx,deg = deg)
}
matplot(xs, FXs.align, lty=1, t="l", col=coly, main="Aligned Distributions",
        xlab="Age",ylab="fx",cex.lab=1.25)

## 2. Standard = mean of the aligned density
FXallmean <- apply(FXs.align, 1, mean)
FXstand <- FXallmean
lines(xs, FXstand, lwd=3)
legend("topleft", c("Standard"),col=1, lwd=3,lty = 1,
       bg="white", pt.cex=1.2,bty="n",cex = 1.5)

## 3. Mode of the standard (extrapolate to left)
Mstand <- M[1]

## 4. Coefficients of the standard: need to augment x-axis
ages.add.l <- 40
ages.add.r <- 30
delta1 <- 1

## define new augmented age-axis 
xA <- c(rev(seq(from=x[1]-delta1, by=-delta1,length=ages.add.l/delta1)), 
        x, 
        seq(from=x[m]+delta1, by=delta1, length=ages.add.r/delta1))
mA <- length(xA)

## new B-splines parameters on augmented axis
xlA <- min(xA)
xrA <- max(xA)
xminA <- round(xlA - 0.01 * (xrA - xlA),3)
xmaxA <- round(xrA + 0.01 * (xrA - xlA),3)
ndxA <- floor((xA[mA]-xA[1])/3)

## augmented B-splines
BA <- MortSmooth_bbase(xA, xminA, xmaxA, ndx=ndxA, deg=deg)
nbxA <- ncol(BA)

## Derive coefficients of the standard
Standard <- coeff_stand(age=xs,fx=FXstand,ndx=ndxA,deg=deg,
                           ages.add.l=ages.add.l,ages.add.r=ages.add.r)
coeff_Stand <- Standard$betasA

## MLE ESTIMATION 
## Break point of the age axis for each year 
PLO <- PUP <-  list()
XLO <- XUP <- list()
for(i in 1:n){
  PLO[[i]] <- which(xA<=floor(M[i]))
  PUP[[i]] <- which(xA>floor(M[i]))
  XLO[[i]] <- xA[PLO[[i]]]
  XUP[[i]] <- xA[PUP[[i]]]
}

## empty vectors and matrices to store results 
PLOT=TRUE
ss <- bsLO <- bsUP <- numeric(n)
DXstad <- lMXstad <- matrix(NA, m, n)

## Estimation
for(i in 1:n){
  cat("fitting year", years[i], "\n")
  conv.stad <- FALSE
  ## starting values
  if (i==1){
  ## start from 1 if first year  
  start.value <- c(1,1)
  }else{
  ## start from previously estimated pars
  start.value <- c(bsLO[i-1],bsUP[i-1])
  }
  ## MLE
  opt <- optim(par=start.value, fn=MLE_obj_FUN,x=x,xA=xA,
               Mstand=Mstand, shat=s[i],
               xlo=XLO[[i]], xup=XUP[[i]],
               coeff.stand=coeff_Stand, Dx=Z[,i], Ex=E[,i],
               xmin=xminA, xmax=xmaxA,
               ndx=ndxA, deg=deg)
  if (opt$convergence != 0) break
  
  ## assign 
  shat <- s[i]
  bLhat <- opt$par[1]
  bUhat <- opt$par[2]

  ## compute dx:
  ## segment a linear transformation function
  ## below the mode
  aL <- Mstand - bLhat * (shat+Mstand)
  wL <- aL + bLhat*XLO[[i]]
  ## above the mode
  aU <- Mstand - bUhat * (shat+Mstand)
  wU <- aU + bUhat*XUP[[i]]
  ## unique transformation function
  wb <- c(wL, wU)
  wb <- sort(wb)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xminA,xmaxA,ndx=ndxA,deg=deg)
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff_Stand))
  fwb <- fwb[xA%in%x]
  ## check if density is greater than one
  if (sum(fwb) > 1) fwb <- (fwb)/sum(fwb)
  ## hazard
  eta <- log(mx_from_dx(fwb))
  ## save parameters, density and mx
  ss[i] <- shat
  bsLO[i] <- bLhat
  bsUP[i] <- bUhat
  DXstad[,i] <- fwb
  lMXstad[,i] <- eta
  ## plotting
  if (PLOT==TRUE){
    plot(xo,lMX[,i],pch=16,main=years[i],ylim = range(lMXstad,lMX,finite=T),
         xlim=range(x),xlab="Age")
    lines(xs,lMX.smooth[,i],col=2,lwd=2,lty=2)
    lines(x,lMXstad[,i],col=4,lwd=2,lty=1)
    legend("topleft",c("Observed","Smooth","STAD"),pch = c(16,NA,NA),
           lty = c(0,2,1),col = c(1,2,4),bty="n",cex = 1.3,lwd = 3)
  }
}

## plot estimated parameters
par(mfrow=c(1,3))
plot(years,ss,t="o",pch=16,main="s",cex.main=2,ylab="")
plot(years,bsLO,t="o",pch=16,ylim=range(bsLO,bsUP),cex.main=2,main=expression(b[L]),ylab="")
plot(years,bsUP,t="o",pch=16,ylim=range(bsLO,bsUP),cex.main=2,main=expression(b[U]),ylab="")
par(mfrow=c(1,1))

## actual and fitted life expectancy
e30.hat <- e30.act <- numeric(n)
for (i in 1:n){
  e30.hat[i] <- lifetable.mx(x=xo,mx=exp(lMXstad[1:mo,i]),sex="F")$ex[1]
  e30.act[i] <- lifetable.mx(x=xo,mx=exp(lMX[,i]),sex="F")$ex[1]
}
plot(years,e30.act,ylim=range(e30.act,e30.hat),pch=16,
     main=paste(cou,"- E30"),ylab="E30",cex.lab=1.25)
points(years,e30.hat,col=2,pch=4,lwd=2)
legend("topleft",c("Observed","STAD"),pch = c(16,4),
       lty = c(0,0),col = c(1,2),bty="n",cex = 1.3,lwd = 3)

## -- FORECASTING -------------

## uncertainty interval & seed for reproducibility
lev <- 80
lev.p <- lev/100
set.seed(2018) 

## forecast horizon
y.fore <- (years[n]+1):2040
n.fore <- length(y.fore)
coly.fore <- rainbow_hcl(n.fore)

## create ts dataframe for s, bL and bU
s.ts <- ts(ss, start = years[1])
df.par <- data.frame(bsLO=bsLO,bsUP=bsUP)
df.par.ts <- ts(df.par, start = years[1])

## Forecasting s: univariate best ARIMA
ds <- 0

## Testing for unit root in original series with augmented DF test
df.test.s <- ur.df(s.ts,lags=2,'trend')
if (df.test.s@teststat[1] > df.test.s@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1

## ts model for s
s.mod1 <- auto.arima(s.ts, d=ds,max.p=3,max.q=3,trace=TRUE)
pred.s <- forecast(s.mod1, h=n.fore,level=lev)
plot(pred.s)

## Forecasting bL and bU: VAR model 
var.typ <- "both"
var.p <- 1
var <- VAR(df.par.ts, p=var.p, type=var.typ)
pred.var <- forecast(var,h=n.fore,level=lev)
plot(pred.var)

## BOOTSTRAP FOR CONSTRUCTING PI 
n.simul <- 200    ## increase for smoother PI
coly.sim <- rainbow_hcl(n.simul)

## define three bootstapping matrices
bootstrap.s <- bootstrap.bsLO <-
  bootstrap.bsUP <- matrix(NA,nrow=n.fore,ncol=n.simul)

## bootsrap S
for(i in 1:n.simul){
  ## generate simulation with bootsrapping
  s.sim <- simulate(s.mod1, nsim=n.fore,
                          future=TRUE, bootstrap=TRUE)
  ## derive the bootsrap values
  bootstrap.s[,i] <- s.sim
}

## bootsrap BL and BU
## estimated residuals
res.bL <- resid(var)[,1]
res.bU <- resid(var)[,2]
for(i in 1:n.simul){
  ## var dimension
  n.var <- n - var.p
  ## sampling residuals
  res.bL.S <- sample(res.bL, n.var, replace=TRUE)
  res.bU.S <- sample(res.bU, n.var, replace=TRUE)
  ## fictitious responses
  bL.S <- c(fitted(var)[,1] + res.bL.S)
  bU.S <- c(fitted(var)[,2] + res.bU.S)
  ## re-fit the model
  df.par.S <- data.frame(bsLO=bL.S,bsUP=bU.S)
  df.par.ts.S <- ts(df.par.S, start = years[1+n.var])
  var.S <- VAR(df.par.ts.S, p=var.p, type=var.typ)
  ## forecast
  pred.var.S <- forecast(var.S,h=n.fore,level=lev)
  ## derive the bootsrap values
  bootstrap.bsLO[,i] <- pred.var.S$forecast$bsLO$mean
  bootstrap.bsUP[,i] <- pred.var.S$forecast$bsUP$mean
}

## plot bootstrap
matplot(y.fore,bootstrap.s,col="grey",t="l",lty=1,lwd=0.8)
matplot(y.fore,bootstrap.bsLO,col="grey",t="l",lty=1,lwd=0.8)
matplot(y.fore,bootstrap.bsUP,col="grey",t="l",lty=1,lwd=0.8)

## matrix and arrays to store bootsrap results
dwb.boot <- array(NA,c(m,n.simul,n.fore))
lmx.boot <- array(NA,c(m,n.simul,n.fore))
e30.mat <- matrix(NA,nrow=n.simul,ncol=n.fore)

## for each simulation and year, compute age-at-death distribution,
## hazard and life expectancy at age 30
for(j in 1:n.fore) {
  for(i in 1:n.simul){
    ## get values of s, bsLO, bsUP for each forecasted year
    s.boot <- bootstrap.s[j,i]
    bsLO.boot <- bootstrap.bsLO[j,i]
    bsUP.boot <- bootstrap.bsUP[j,i]
    ## compute new mode
    M.boot <- s.boot + Mstand
    ## divide forecast distribution in two pieces
    PLO[[n+j]] <- which(xA<=floor(M.boot))
    PUP[[n+j]] <- which(xA>floor(M.boot))
    XLO[[n+j]] <- xA[PLO[[n+j]]]
    XUP[[n+j]] <- xA[PUP[[n+j]]]

    ## segment a linear transformation function
    ## below the mode of the forecast year
    aL <- Mstand - bsLO.boot * (s.boot + Mstand)
    wL <- aL + bsLO.boot*XLO[[n+j]]
    ## above the mode of the forecast year
    aU <- Mstand - bsUP.boot * (s.boot + Mstand)
    wU <- aU + bsUP.boot*XUP[[n+j]]
    ## unique transformation function
    wbhat <- c(wL, wU)
    wbhat <- sort(wbhat)
    ## evaluating B on shifted x
    Bs.fore <- MortSmooth_bbase(wbhat, xl=xminA, xr=xmaxA, ndx=ndxA, deg)
    ## shifted density
    fwb.boot <- exp(Bs.fore %*% coeff_Stand)
    fwb.boot <- fwb.boot[xA%in%x]
    ## check if density is greater than one
    if (sum(fwb.boot) > 1) fwb.boot <- fwb.boot/sum(fwb.boot)
    ## save density, mx and e30
    dwb.boot[,i,j] <- fwb.boot
    eta.boot <- as.vector(log(mx_from_dx(fwb.boot)))
    lmx.boot[,i,j] <- eta.boot
    e30.mat[i,j] <- lifetable.mx(x=xo,mx=exp(eta.boot[1:mo]),sex="F")$ex[1]
  }
  cat("bootstrapping year", y.fore[j], "\n")
}

## derive median and PI of distribution
dwb.fore.MEAN <- dwb.fore.UPS <- dwb.fore.LOS <- matrix(NA,nrow=m,ncol=n.fore)
for(j in 1:n.fore){
  dwb.fore.MEAN[,j] <- apply(dwb.boot[,,j], 1, median)
  dwb.fore.UPS[,j] <- apply(dwb.boot[,,j], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  dwb.fore.LOS[,j] <- apply(dwb.boot[,,j], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}

## check that the transformation is sustained 
## (if something is strange here, augment ages.add.l and ages.add.r)
matplot(x,log(dwb.fore.MEAN),col=coly.fore,t="l",lty=1)
matplot(x,log(dwb.boot[,,n.fore]),col=coly.sim,t="l",lty=1)

## derive median and PI of mortality rates
lmx.fore.MEAN <- lmx.fore.UPS <- lmx.fore.LOS <- matrix(NA,nrow=m,ncol=n.fore)
for(j in 1:n.fore){
  lmx.fore.MEAN[,j] <- apply(lmx.boot[,,j], 1, median)
  lmx.fore.UPS[,j] <- apply(lmx.boot[,,j], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  lmx.fore.LOS[,j] <- apply(lmx.boot[,,j], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}

## compute E30 MEAN and PI
e30.fore.mean <- apply(e30.mat, 2, median)
e30.fore.up <- apply(e30.mat, 2, quantile, prob=1- (1-lev.p)/2,na.rm=T)
e30.fore.low <- apply(e30.mat, 2, quantile, prob=(1-lev.p)/2,na.rm=T)

## plot E30
my.col <- 4
my.colT <- adjustcolor(my.col,alpha.f = 0.3)
plot(years, e30.act, t="n",ylim=range(e30.act,e30.fore.up),
     xlim=range(years,y.fore),main=paste("e30 -",cou),xlab="Years", ylab="")
grid()
## STAD 
xx <- c(c(years[n], y.fore), rev(c(years[n], y.fore)))
yy <- c(c(e30.act[n],e30.fore.up), rev(c(e30.act[n],e30.fore.low)))
polygon(xx, yy, border = my.colT, col=my.colT)
lines(c(years[n], y.fore), c(e30.act[n],e30.fore.mean), col=my.col, lwd=3)
## actual & fitted
lines(years, e30.act, col=1, pch=16, cex=1.1, t="p", lwd=2)
points(years, e30.hat, col=2, pch=4, cex=1.1, t="p", lwd=2)
## legend
legend("topleft", c("Actual","Fitted",paste("STAD with",lev,"% PI")),
       col=c(1, 2,my.col), lwd=c(3,3,2),lty = c(NA,NA,1),
       pch=c(16,4,-1), bg="white", pt.cex=1.2,
       bty="n",cex = 1.5)

## END




