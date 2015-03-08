#simulate from factor model
#' @param n number of SNPs
#' @param d number of tissues
#' @param betasd the measurement error
#' @param esd the variance of the residual matrix E in Y=LF+E
factor_sim=function(n=1000,d=3,betasd=1,esd=0.1){

  temp=rep(list(c(0,1)),d)
  configs = expand.grid(temp) # all possible 2^d combinations
  F=as.matrix(configs)
  k = nrow(F) # number of factors

  possible_loadings = diag(k) #set possible loadings to be "sparse" (loaded on one factor each)
  z = sample(nrow(F),n,replace=TRUE) # randomly sample factor to be loaded on for each SNP
  L = possible_loadings[z,] 
  effects = c(rnorm(n/2,0,2),rnorm(n/2,0,0.1)) # a mixture of big and small effects
  L = L*effects #multiply each lambda by its effects
  E = matrix(rnorm(n*d,0,0.1),nrow=n)
  beta = L %*% F + E

  betahat = beta + rnorm(n*d,0,betasd)
  return(list(beta=beta,betahat=betahat,L=L,F=F,E=E))
}

test = factor_sim(1000,2)
#visualize the results
plot(test$beta[,1],test$beta[,2])
plot(test$betahat[,1],test$betahat[,2])

#do one with lower measurement error
test = factor_sim(1000,2,betasd=0.2)
plot(test$beta[,1],test$beta[,2])
plot(test$betahat[,1],test$betahat[,2])
