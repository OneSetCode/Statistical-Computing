#Question 1-------------------------
logdet <- function(R){
  det=prod(eigen(R)$values)
  return(log(det))
};

#Question 2--------------------------
data=read.table("erdata.txt");
logmarglik <- function(d,A){
  n=nrow(d)
  A_=length(A);
  D1=matrix(unlist(d[1]));
  Da=as.matrix(d[A]);
  Ma=diag(A_)+t(Da) %*% Da;
  logp=lgamma((n+A_+2)/2)-lgamma((A_+2)/2)-0.5*logdet(Ma)-
    ((n+A_+2)/2)*log(1+t(D1)%*%D1-t(D1)%*%Da%*%solve(Ma)%*%t(Da)%*%D1);
  return(logp)
};
logmarglik(data,c(2,5,10));

