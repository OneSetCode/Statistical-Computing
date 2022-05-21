data = read.table("534binarydata.txt", header = FALSE);
library(rcdd);
#-------------------------------------------------------------------------------
getLogisticAIC <- function(y, vars, data)
{
  if(0 == length(vars))
  {
    out = glm(data[,y] ~ 1, family = binomial);
  }
  else
  {
    out = suppressWarnings(glm(data[,y] ~ as.matrix(data[,as.numeric(vars)]),
                               family = binomial));
  }
  AIC = out$deviance+2*(1+length(vars));
  converged = out$converged;
  return(ifelse(converged, AIC, NA));
}
#-------------------------------------------------------------------------------
isValidLogisticRCDD <- function(response,explanatory,data)
{
  if(0==length(explanatory))
  {
    return(TRUE);
  }
  logisticreg = suppressWarnings(glm(data[,response] ~ 
                                       as.matrix(data[,as.numeric(explanatory)]),
                                     family=binomial(link=logit),x=TRUE));
  tanv = logisticreg$x;
  tanv[data[,response] == 1, ] <- (-tanv[data[,response] == 1, ]);
  vrep = cbind(0, 0, tanv);
  lout = linearity(vrep, rep = "V");
  return(length(lout)==nrow(data));
}
#-------------------------------------------------------------------------------
# get the initial logistic model
initModel = function(response, data)
{
  vars = setdiff(1:ncol(data),response);
  valid = 1;
  while(valid)
  {
    initVars = sample(vars,sample(vars,1));
	
    if(isValidLogisticRCDD(response,initVars, data))
    {
      return(initVars);
    }
  }	
}
#-------------------------------------------------------------------------------
#get the valid neighbors
validNeighbors = function(nowVar, response, data)
{
  vars = setdiff(1:ncol(data),response)
  validNgbVars = list();
  for(i in vars)
  {
    if(i %in% nowVar)
      ngbVar = setdiff(nowVar, i)
    else
      ngbVar = c(nowVar,i)
    
    if(isValidLogisticRCDD(response,ngbVar,data))
      validNgbVars = append(validNgbVars,list(ngbVar))
  }
  return(validNgbVars)
}
#-------------------------------------------------------------------------------
MC3search = function(response, data, n)
{
  #get the initial model and AIC
  bestVar = nowVar = initModel(response, data)
  bestAIC = nowAIC = getLogisticAIC(response, nowVar, data)
  
  for(i in 1:n)
  {
    #get neighbors
    neigbors = validNeighbors(nowVar, response, data)
    #sample a model
    nextVar = unlist(neigbors[sample.int(length(neigbors),1)])
    #get the next neighbors
    nextVarNgb = validNeighbors(nextVar, response, data)
    #compute pA'
    nextAIC = getLogisticAIC(response,nextVar,data)
    pA_ = -nextAIC-log(length(nextVarNgb))
    #compute pA(r)
    pAr = -getLogisticAIC(response,nowVar,data)-log(length(neigbors))
    #compare
    if(pA_ > pAr)
    {
      nowVar = nextVar
      nowAIC = nextAIC
    }
    else
    {
      u = runif(n = 1, min = 0, max = 1)
      if(log(u) < pA_ - pAr)
      {
        nowVar = nextVar
        nowAIC = nextAIC
      }
    }
    #determine the best model so far
    if(nowAIC < bestAIC)
    {
      bestVar = nowVar
      bestAIC = nowAIC
    }
  }
  return(list(bestAIC = bestAIC, bestModel = sort(bestVar)))
}
#-------------------------------------------------------------------------------
set.seed(123);
out = list();
for(i in 1:10)
{
  out = append(out,MC3search(61, data, 25))
}
out
#See the results below.Through the results, I can observe that the models obtained 
#from this method can be quite different from each other, and we might need more 
#iterations to gain more precise model results.

# $bestAIC
# [1] 96.22832
# 
# $bestModel
# [1]  3  5  8 11 12 13 14 15 17 18 19 25 27 28 29 34 37 40 42 46 49 50 51 54 58 59
# 
# $bestAIC
# [1] 64.73364
# 
# $bestModel
# [1]  3  4 12 13 14 16 18 22 23 25 33 34 46 49 50 51 52 54 59
# 
# $bestAIC
# [1] 66.63069
# 
# $bestModel
# [1]  5  7 15 19 21 23 25 27 29 32 34 38 42 44 45 46 49 52 55
# 
# $bestAIC
# [1] 75.45252
# 
# $bestModel
# [1]  2  3  6  8 11 13 16 20 21 22 23 24 29 32 36 39 42 44 45 48 51 54 55 57 60
# 
# $bestAIC
# [1] 76.48669
# 
# $bestModel
# [1]  4  6 10 11 17 20 21 25 26 27 28 29 30 31 33 35 37 39 41 42 45 49 55 56 59
# 
# $bestAIC
# [1] 102.0004
# 
# $bestModel
# [1]  8 10 13 17 18 20 29 31 32 39 40 41 42 44 51 52 56 57
# 
# $bestAIC
# [1] 94.0068
# 
# $bestModel
# [1]  2  3  5  6  7 11 12 13 15 16 19 22 25 28 30 31 33 37 40 41 42 45 50 53 54 55 58
# 
# $bestAIC
# [1] 92.72992
# 
# $bestModel
# [1]  1  8 11 17 18 21 28 32 36 47 49 51 56 59
# 
# $bestAIC
# [1] 75.28514
# 
# $bestModel
# [1]  1  2  3  8  9 10 11 12 16 23 24 26 28 40 41 47 48 49 50 51 53 59 60
# 
# $bestAIC
# [1] 67.35596
# 
# $bestModel
# [1]  1  3  4  7  8 11 15 23 25 37 42 47 53 54

















