predict.msda <- function(obj, x) {
  theta<-obj$theta
  mu<-obj$mu
  prior<-obj$prior
  mubar <- sweep(mu[, -1], 1, mu[, 1], "+")/2
  n <- nrow(x)
  p<-ncol(x)
  x.train<-obj$x
  y.train<-obj$y
  nclass <- length(prior)
  nlambda <- length(theta)
  pred<-matrix(0,n,nlambda)
  pred[1]<-which.max(prior)
  for(i in 1:nlambda){
    nz<-sum(theta[[i]][,1]!=0)
    if(nz==0){pred[i]<-which.max(prior)}else{
      xfit<-x.train%*%theta[[i]][,1:(min(nclass-1,nz))]
      suppressWarnings(l<-lda(xfit,y.train))
      pred[,i]<-predict(l,x%*%theta[[i]][,1:(min(nclass-1,nz))])$class}
  }    
  pred
} 
