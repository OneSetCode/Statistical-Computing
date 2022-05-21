#Problem 1 ---------------------------------------------------------------------
#compute the Laplace Approximation function
LaplaceApprox <- function(variable, response, data, mode)
{
  x = data[, variable];
  y = data[, response];
  l_ = LogLikelihood_(x, y, mode);
  result = log(2 * pi) + l_ - (1 / 2) * LogDet(-Hessian(x, y, mode));
  return(result)
}

#Problem 2 ---------------------------------------------------------------------
#compute the Metropolis-Hastings sampling algorithm
MHSample <- function(variable, response, data, mode, n)
{
  x = data[, variable]
  y = data[, response]
  nowBeta = mode
  cov = - solve(Hessian(x, y, mode))
  sample = matrix(0, n, 2)
  
  for(i in 1:n)
  {
    candBeta = mvrnorm(1, nowBeta, cov)
    prob = LogLikelihood_(x, y, candBeta) - LogLikelihood_(x, y, nowBeta)
    if(prob >= 1)
    {
      nowBeta = candBeta
    }
    else
    {
      if(log(runif(1, 0, 1)) <= prob)
      {
        nowBeta = candBeta
      }
    }
    sample[i, ] = nowBeta
  }
  return(sample)
}

#Compute the Metropolis-Hastings estimate
MHEstimate <- function(variable, response, data, mode, n)
{
  sample = MHSample(variable, response, data, mode, n)
  B0 = mean(sample[, 1])
  B1 = mean(sample[, 2])
  return(c(B0,B1))
}


#Functions needed --------------------------------------------------------------
#compute the Log-Likelihood star l* function
LogLikelihood_ <- function(x, y, beta)
{
  P = Pi(x, beta);
  LogLikelihood = sum(y * log(P)) + sum((1 - y) * log(1 - P));
  return(-log(2 * pi) - 1/2 * (sum(beta^2)) + LogLikelihood)
}

#compute the Pi function
Pi <- function(x, beta)
{
  regr = cbind(rep(1, length(x)), x) %*% beta;
  return(exp(regr) / (1 + exp(regr)));
}

#compute the log-det function
LogDet = function(m)
{
  return(sum(log(eigen(m)$values)));
}

#compute the Hessian matrix
Hessian = function(x, y, beta)
{
  P_ = Pi_(x, beta);
  hess = matrix(c(-sum(P_) - 1, -sum(P_ * x), -sum(P_ * x), -sum(P_ * x^2) - 1), 2, 2, byrow = TRUE);
  return(hess);
}

#compute the Pi_ function
Pi_ <- function(x, beta)
{
  regr = cbind(rep(1, length(x)), x) %*% beta;
  return(exp(regr) / (1 + exp(regr))^2);
}

#compute the mode with Newton-Raphson Algorithm
require(MASS)
NRMode <- function(variable, response, data)
{
  x = data[, variable]
  y = data[, response]
  nowBeta = matrix(c(0,0),2,1)
  
  while(1)
  {
    gradient = matrix(c(sum(y - Pi(x, nowBeta)) - nowBeta[1], sum((y - Pi(x, nowBeta)) * x) - nowBeta[2]), 2, 1)
    nextBeta = nowBeta - solve(Hessian(x, y, nowBeta)) %*% gradient
    
    if(abs(nextBeta[1]-nowBeta[1]) < 0.0001 & abs(nextBeta[2]-nowBeta[2]) < 0.0001)
    {
      break
    }
  }
  return(nowBeta)
}

