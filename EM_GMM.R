EM_GMM <- function(X,K,max_it=50){
  # Initialization
  p = ncol(X)
  n = nrow(X)
  prop = rep(NA,K)
  mu = matrix(NA,K,p)
  Sigma = array(NA,dim = c(K,p,p))
  
  # Random init
  T = t(rmultinom(n,K,rep(1/K,K)))
  
  # The EM loop
  for (it in 1:max_it){       cat('.')
    # M step
    for (k in 1:K){
      nk = sum(T[,k])
      prop[k] = nk / n
      mu[k,] = 1/nk * colSums(T[,k] %*% matrix(1,1,p) * X)
      Ak = T[,k] * as.matrix(X - matrix(1,n,1)%*%mu[k,])
      Bk = as.matrix(X - matrix(1,n,1)%*%mu[k,])
      Sigma[k,,] = 1/nk * t(Ak) %*% Bk
    }
    # E step
    for (k in 1:K) T[,k] = prop[k] * dmvnorm(X,mu[k,],Sigma[k,,])
    T = T / rowSums(T) %*% matrix(1,1,K)
    # Visualization
    plot(X,pch=19,col=max.col(T))
    for(k in 1:K) points(mu[k,1],mu[k,2],pch="*",cex=5,col=k)
    Sys.sleep(0.25)
  }
  return(list(prop=prop,mu=mu,Sigma=Sigma,T=T))
}

#data(iris); X = iris[,2:3]
X = rbind(rmvnorm(100,mean = c(0,0)),
          rmvnorm(100,mean = c(-10,-10)),
          rmvnorm(100,mean = c(-10,10)))

out = EM_GMM(X,3)
