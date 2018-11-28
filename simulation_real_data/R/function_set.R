### Definition of loss, gradient, hessian

setting_fun <- function(partition, i, newz, p, l)
{
  part <- partition[task_num,][i,]
  
  set.seed(part$seed)
  
  id <- sample(1:83, 60) 
  
  train.x <- newz[id,]
  test.x <- newz[-id,]
  
  train.y <- y[id]
  test.y <- y[-id]
  
  lambda_1 <- part$lambda_1
  
  lambda_2 <- part$lambda_2
  
  lambda_vec <- c(rep(lambda_1, p), rep(lambda_2, p*l))
  
  init_list <- list("train.x" = train.x, "test.x" = test.x, "train.y" = train.y, "test.y" = test.y, "lambda_1" = lambda_1, "lambda_2" = lambda_2, "lambda_vec" = lambda_vec) 
  
  return(init_list)
}

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