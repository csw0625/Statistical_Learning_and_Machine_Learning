## Non-Linear Model with non-constant variance
## Solving Non-constant Problem with Gamma function modeling


library(datasets)
library(ggplot2)

data=attenu ; dim(data)

x1 = data$mag ; x2 = data$dist ; y = data$accel


####initiate with linear regression since there is no information####
fit1=lm(accel~mag + dist,data=data) ; summary(fit1) ; par(mfrow=c(2,2)) ; plot(fit1)

par(mfrow=c(1,2)) ; plot(x1,y) ; plot(x2,y) 
dev.off()


##By seeing plot, residual is not symmetric to 0 -> Relationship of X and Y is not linear(mean function is not linear)
##If unfolding residual line to linear, It would be showing some patterns -> Non-constant variance 

##First of all, Need to solve mean function so that the residual have 0 mean 



####Using GAM to implement non-linear regression####
library(gam)
fit2 = gam(accel~s(mag,5) + s(dist,5), data=data) ; par(mfrow=c(1,2)) ; plot(fit2) ##Smoothing Spline=5

##By checking plot, It seems that degree of freedom 5 is too much for 'mag'.Need to check whether it is linear or non-linear
##'dist' is apparently seems to decrease exponentially. -> Use exponential

###model comparison for 'mag' (linear or not)
fit2_1=gam(accel~mag +s(dist,5),data=data)

anova(fit2,fit2_1) ## If the test is not significant, linear term will do. Else, Non-linear relation is appropiate 

###Test result tells Non-linear -> quadratic form would be good


### mag: Quadratic, dist:exponential

####Nonlinear Model with constant variance
#Y=accel, X1=mag, X2=dist
#Nonlinear parametric model: Y=Beta1 + Beta2*X1 + Beta3*X1^2 + Beta4*exp(-beta5*X2)
#Under Normality and constant Variance assumption

##Seeing mag, it drastically increases but slows down as mag becomes bigger at some point

f=function(beta,X){
  
  X1=X[,1] ; X2=X[,2]
  beta[1] + beta[2]*X1 + beta[3]*X1^2 + beta[4]*exp(-beta[5]*X2)   ## negative sign in exponential term is optional
}

##beta is 5*1 vector 

###objective function###
RSS = function(beta,Y,X) sum((Y-f(beta,X))^2)

###Gradient vector of the objective function###
grv = function(beta,Y,X){
  X1 = X[,1] ; X2 = X[,2]
  R = Y-f(beta,X)
  c(-2*sum(R),-2*sum(R*X1),-2*sum(R*X1^2),-2*sum(R*exp(-beta[5]*X2)),   ##derivatives for 5 parameters(If possible)
    2*beta[4]*sum(R*X2*exp(-beta[5]*X2)))
}

###optimization###                          
X=cbind(attenu$mag,attenu$dist) ; colnames(X)=c('mag','dist') # X matrix
Y=attenu$accel

ml1=optim(rep(0.1,5),RSS,gr=grv,method='BFGS',X=X,Y=Y,hessian = T) ; ml1 
## initial value 0.1 for 5 parameters Using BFGS(unconstrained optimization)
## BFGS uses approximate Hessian(2nd derivative) for step size (step size as inverse of Hessian)
## convergence 0 means successful optimization (1 means iteration had been reached->May not be full converged)

beta_hat=ml1$par ; beta_hat


Yhat=f(beta_hat,X) ###Fitted value, 

r = Y - Yhat ; par(mfrow=c(1,1)) ; plot(Yhat,r,ylim=c(-0.5,0.5)) ; lines(c(0,10),c(0,0),col='red') ###Residual plot

##megaphone shape appears ->Definitely Violating constant variance assumption



### Non-linear parametric model with non-constant variance modeling ###
library(matrixcalc) ##To check whether a variance matrix is singular or not

### objective function for mean function: Generalized least sqaure method (minimizing problem)
obj.mean = function(beta,X,Y,S) t(Y-f(beta,X)) %*% solve(S) %*% (Y-f(beta,X))  #S: Covariance matrix

###Gradient vector of the objective function
gr.mean = function(beta,Y,X,S){
  
  sigma2=diag(S)  ## Still Assuming independence of variance 
  X1 = X[,1] ; X2= X[,2]
  R = Y-f(beta,X)
  c(-2*sum(R/sigma2), -2*sum(R*X1/sigma2), -2*sum(R*X1^2/sigma2),
    -2*sum(R*exp(-beta[5]*X2)/sigma2),
    2*beta[4]*sum(R*X2*exp(-beta[5]*X2)/sigma2))
}

## Linear variance function: |r|=gam1+gam2*Yhat
## For linear variance function, we can consider absolute residuals instead of squared residual
## gam.hat = (Z^t W Z)^(-1) Z^t w |r|

##iteration procedure

beta.new=ml1$par    
#initial parameter under assumption that the new beta wouldn't be that much different with constant-variance beta
W = diag(rep(1,length(Y))) #initial value of weight matrix
mdif= 100000  


while(mdif > 0.000001){
  Yhat = f(beta.new,X)
  
  ##Variance function##
  r = Y - Yhat
  Z = cbind(1,Yhat)
  gam.hat = solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W %*% abs(r)
  sigma = Z %*% gam.hat
  S = diag(as.vector(sigma^2))
  
  if (is.non.singular.matrix(S)) W = solve(S)
  else W = solve(S+0.000000001)*diag(rep(1,nrow(S)))   ##To make matrix to Non-singular
  
  ml2 = optim(beta.new, obj.mean , gr=gr.mean, method='BFGS' ,Y=Y, X=X, S=S)
  beta.old = beta.new
  beta.new = ml2$par
  mdif = max(abs(beta.new - beta.old))  
  ## if the absolute value of difference of previous beta and beta after the iteration is smaller than stopping criterion,
  ## iteration stops
}

Yhat = f(beta.new,X) ; sigma = Z %*% gam.hat ; r = (Y-Yhat)/sigma

###new residual plot
plot(Yhat,r,ylim=c(-4,4)) ; lines(c(0,10),c(0,0),col='red')

