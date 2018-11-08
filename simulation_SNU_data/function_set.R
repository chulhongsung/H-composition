### Definition of loss, gradient, hessian

Loss <- function(beta, Y, X, lambda_1, lambda_2, nu_tmp, A, rho, gamma_tmp){
  L  <- -t(Y) %*% (X %*% beta) + colSums(log(1+exp(X%*%beta))) + 
    lambda_1*sum(abs(gamma_tmp[1:123])) + lambda_2*sum(abs(gamma_tmp[124:492])) +
    t(nu_tmp) %*% (A%*%beta - gamma_tmp) + (rho/2)*t(A%*%beta - gamma_tmp)%*%(A%*%beta - gamma_tmp)
  return(L)
}

gradient <- function(beta, X, Y)
{
  G <- t(X)%*%(-Y+(1/(1+exp(-X%*%beta))))
  return(t(G))
}

likelihood_hessian <- function(beta, X, n)
{
  H <- list()
  for( i in seq(n))
  {
    S <- as.numeric(1/(1+exp(-X[i,]%*%beta)))
    H[[i]] <- (X[i,]%*%t(X[i,])) * S * (1-S)
    
  }
  return(H)
}