#data=read.table("534binarydata.txt");

#Question 1 ---------------------------------------------------------------
getLogisticAIC <- function(response,explanatory,data)
{
  #check if the regression has no explanatory variables
  if(0==length(explanatory))
  {
    #regression with no explanatory variables
    deviance = glm(data[,response] ~ 1,family=binomial(link=logit))$deviance;
  }
  else
  {
    #regression with at least one explanatory variable
    deviance = glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),family=binomial(link=logit))$deviance;
  }
  return(deviance+2*(1+length(explanatory)));
}


#Question 2 ---------------------------------------------------------------
forwardSearchAIC <- function(response,data,lastPredictor)
{
  
  #fit the first regression
  bestAIC = getLogisticAIC(response,c(),data);
  bestExp = c();
  #predictors that are not in the model yet
  remainExp = setdiff(1:lastPredictor,y=bestExp);
  
  while(length(remainExp)>0) 
  {
    #record if a new predictor will be added
    changed = 0
    #the index of the new-added predictor (if there is one)
    newInd = NULL
    for(i in remainExp) 
    {
      #new predictors list after adding the predictor i
      newExp = append(bestExp,i);
      if(getLogisticAIC(response,newExp,data) < bestAIC)
      {
        bestAIC = getLogisticAIC(response,newExp,data)
        changed = 1
        newInd = i
      }
    }
    if(changed == 1)
    {
      bestExp = append(bestExp,newInd);
      remainExp = setdiff(1:lastPredictor,y=bestExp);
    }
    else
    {
      #no new predictor added, return current predictor list and AIC
      return(list(aic=bestAIC,reg=bestExp)); 
    }
  }
  
  return(list(aic=bestAIC,reg=bestExp)); 
}



#Question 3 ---------------------------------------------------------------
backwardSearchAIC <- function(response,data,lastPredictor)
{
  #fit the first regression
  bestAIC = getLogisticAIC(response,1:lastPredictor,data);
  bestExp = 1:lastPredictor;
  
  while(length(bestExp)>0) 
  {
    #record if a new predictor will be deleted
    changed = 0
    #the index of the new-deleted predictor (if there is one)
    newInd = NULL
    for(i in bestExp) 
    {
      #new predictors list after deleting the predictor i
      newExp = bestExp[-which(bestExp==i)];
      if(getLogisticAIC(response,newExp,data) < bestAIC)
      {
        bestAIC = getLogisticAIC(response,newExp,data)
        changed = 1
        newInd = i
      }
    }
    if(changed == 1)
    {
      bestExp = bestExp[-which(bestExp==newInd)];    
    }
    else
    {
      #no new predictor deleted, return current predictor list and AIC
      return(list(aic=bestAIC,reg=bestExp)); 
    }
  }
  
  return(list(aic=bestAIC,reg=bestExp)); 
}



#Question 4 ---------------------------------------------------------------
#Answer:
#  No, the logistic regression models identified in Problems 2 and 3 are not the same.
#  Yes, they have the same AIC.
#  See the following code for repreating the forward and backward greedy procedures with BIC.
#  With BIC, we obtain the same two logistic regressions as that from AIC respectively, and with a same BIC 54.96934.
  

  getLogisticBIC <- function(response,explanatory,data)
  {
    #check if the regression has no explanatory variables
    if(0==length(explanatory))
    {
      #regression with no explanatory variables
      deviance = glm(data[,response] ~ 1,family=binomial(link=logit))$deviance;
    }
    else
    {
      #regression with at least one explanatory variable
      deviance = glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),family=binomial(link=logit))$deviance;
    }
    return(deviance+log(nrow(data))*(1+length(explanatory)));
  }
  
  
  forwardSearchBIC <- function(response,data,lastPredictor)
  {
    
    #fit the first regression
    bestBIC = getLogisticBIC(response,c(),data);
    bestExp = c();
    #predictors that are not in the model yet
    remainExp = setdiff(1:lastPredictor,y=bestExp);
    
    while(length(remainExp)>0) 
    {
      #record if a new predictor will be added
      changed = 0
      #the index of the new-added predictor (if there is one)
      newInd = NULL
      for(i in remainExp) 
      {
        #new predictors list after adding the predictor i
        newExp = append(bestExp,i);
        if(getLogisticBIC(response,newExp,data) < bestBIC)
        {
          bestBIC = getLogisticBIC(response,newExp,data)
          changed = 1
          newInd = i
        }
      }
      if(changed == 1)
      {
        bestExp = append(bestExp,newInd);
        remainExp = setdiff(1:lastPredictor,y=bestExp);
      }
      else
      {
        #no new predictor added, return current predictor list and BIC
        return(list(bic=bestBIC,reg=bestExp)); 
      }
    }
    
    return(list(bic=bestBIC,reg=bestExp)); 
  }
  
  
  
  backwardSearchBIC <- function(response,data,lastPredictor)
  {
    
    #fit the first regression
    bestBIC = getLogisticBIC(response,1:lastPredictor,data);
    bestExp = 1:lastPredictor;
    
    while(length(bestExp)>0) 
    {
      #record if a new predictor will be deleted
      changed = 0
      #the index of the new-deleted predictor (if there is one)
      newInd = NULL
      for(i in bestExp) 
      {
        #new predictors list after deleting the predictor i
        newExp = bestExp[-which(bestExp==i)];
        if(getLogisticBIC(response,newExp,data) < bestBIC)
        {
          bestBIC = getLogisticBIC(response,newExp,data)
          changed = 1
          newInd = i
        }
      }
      if(changed == 1)
      {
        bestExp = bestExp[-which(bestExp==newInd)];    
      }
      else
      {
        #no new predictor deleted, return current predictor list and BIC
        return(list(bic=bestBIC,reg=bestExp)); 
      }
    }
    
    return(list(bic=bestBIC,reg=bestExp)); 
  }
  