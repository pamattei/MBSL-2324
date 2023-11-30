library(MASS)
em <- function(X,K,nit=50,eps=1e-5){
	# Initialisation
	N = nrow(X)
	p = ncol(X)
	T = matrix(NA,N,K)
	ll = c(Inf)
	T = t(rmultinom(N,1,rep(1/K,K)))
	
	# Boucle
	for (i in 1:nit){
		cat('*')
		prms = em.mstep(X,K,T)
		resEstep = em.estep(X,K,prms)
    T = resEstep$T
    ll[i+1] = resEstep$ll
    # Visualization
    plot(X,pch=19,col=max.col(T))
    for(k in 1:K) points(prms$m[k,1],prms$m[k,2],pch="*",cex=5,col=k)
    Sys.sleep(0.05)
		if (abs(ll[i+1]-ll[i]) < N*eps) break
	}
	plot(res$ll,type='b')
	bic = em.bic(ll[i],K,p,N)
	cat('\n')
	list(T=T,cls=max.col(T),prms=prms,ll=ll,bic=bic)
}

em.estep <- function(X,K,prms){
	N = nrow(X); p = ncol(X)
	T = Q = matrix(NA,N,K)
	for (k in 1:K){
    Xk = as.matrix(X - matrix(1,N,1)%*%prms$m[k,])
		Q[,k] = log(prms$prop[k]) - p/2*log(2*pi) - 1/2*log(det(prms$S[k,,])) - 1/2*diag(Xk %*% ginv(prms$S[k,,]) %*% t(Xk))
	}
	ll = sum(log(rowSums(exp(Q-apply(Q,1,max))))+apply(Q,1,max))
	for (k in 1:K) {T[,k] = 1 / rowSums(exp(-(Q[,k]*matrix(1,N,K)-Q)))}
  list(T=T,ll=ll)
}

em.mstep <- function(X,K,T){
	p = ncol(X)
	N = nrow(X)
	prop = c()
	m = matrix(NA,K,p)
	S = array(NA,c(K,p,p))
	for (k in 1:K){
		prop[k] = sum(T[,k]) / N
		m[k,] = colSums(T[,k] * X) / sum(T[,k])
		Xmk = X - matrix(1,N,1) %*% m[k,]
		S[k,,] = (t(T[,k] * Xmk) %*% as.matrix(Xmk)) / sum(T[,k])
	}
	list(prop=prop,m=m,S=S)
}

em.bic <- function(ll,K,p,n){
  nu = K-1 + K*p + K*p*(p+1)/2
  bic = ll - nu/2*log(n) 
}

###################################################################
X = rbind(rmvnorm(100,mean = c(0,0)),
          rmvnorm(100,mean = c(-10,-10)),
          rmvnorm(100,mean = c(-10,10)))
res = em(X,3)