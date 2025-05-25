## Model with Time correlation
## Solving dependency of residuals with AR(1) Modeling
tsdat = read.table('tsdat.txt',header=T)

fit = lm(y ~ x, data=tsdat)
summary(fit)

par(mfrow=c(2,2))
plot(fit)

# Durbin-Watson test
# install.packages('lmtest')
library(lmtest)

dwtest(fit)

# Check ACF & PACF
# install.packages('astsa')
library(astsa)

# AR(p): ACF: Exponentially decreasing; PACF: Non-zero values at first p lags.
# MA(q): ACF: Non-zero values at first q lags; PACF: Exponentially decreasing.
# ARMA(p,q): ACF: Similar to ACF of AR(p); PACF: Similar to PACF of MA(q).

acf2(residuals(fit))

ar1 = sarima (residuals(fit), 1,0,0, no.constant=T)   #AR(1)
ar1$fit

# Method 1
# MLE: Multivariate normal distribution
X = cbind(1,tsdat$x)
Y = tsdat$y
n = length(Y)
S = diag(rep(1,n))    # initial covariance matrix

mdif = 1000000
beta.old = rep(100000,2)

while(mdif > 0.0000001)
{
  beta.new = as.vector(solve(t(X) %*% solve(S) %*% X) %*%t(X) %*% solve(S) %*% Y)
  r = as.vector(Y - (X %*% beta.new))
  ar1 = sarima (r, 1,0,0, no.constant=T, details=F)
  alpha = ar1$fit$coef
  sigma2 = ar1$fit$sigma2
  
  mdif = max(abs(beta.new - beta.old))
  beta.old = beta.new
  # Construct covariance matrix
  S = matrix(nrow=n,ncol=n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      if (i == j) S[i,j] = 1
      if (i != j) S[i,j] = alpha^(abs(i-j))
    }
  }
  S = (sigma2 / (1-alpha^2)) * S
}

round(beta.new,4)

# Method 2
# MLE: Product of conditional distribution (Approximation)
# Y_t | Y_t-1 ~ N(X_t*beta + alpha*epsilon_t-1, sigma^2)

fit = lm(y ~ x, data=tsdat)

Yt = tsdat$y[2:n]
Xt = tsdat$x[2:n]
et = residuals(fit)[1:(n-1)]
mdif = 10000
b.old = rep(0,3)

while(mdif > 0.0000001)
{
  fit.temp = lm(Yt ~ Xt + et)
  b.new = fit.temp$coefficient
  mdif = max(abs(b.new[1:2] - b.old[1:2]))
  
  et = (Y - X %*% b.new[1:2])[1:(n-1)]
  b.old = b.new
}

round(b.new,4)