rm(list = ls())
gc(reset = TRUE)

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(quadprog)){ install.packages('quadprog')}; require(quadprog)

load('H-composition simulation.RData')

set.seed(2018)

p = 82; l = 3; n = 1000

### Loss function for beta

# Loss <- function(beta){
#   L  <- -t(Y) %*% (X %*% beta) + colSums(log(1+exp(X%*%beta)))
#   return(L)
# }

### Gradient
gradient <- function(beta)
{
  G <- t(X)%*%(-Y+(1/(1+exp(-X%*%beta))))
  return(t(G))
}


likelihood_hessian <- function(beta)
{
  H <- list()
  for( i in seq(n))
  {
    S <- as.numeric(1/(1+exp(-X[i,]%*%beta)))
    H[[i]] <- (X[i,]%*%t(X[i,])) * S * (1-S)
    
  }
  return(H)
}

### 
rho <- 0.5

lambda_1 <- 5

lambda_2 <- 5

lambda_vec = c(rep(lambda_1, p), rep(lambda_2, p*3))

### initial gamma, nu, u and beta matrix

gamma_tmp <- rnorm(p*4, 0, 1)

nu_tmp <- rnorm(p*4, 0, 1)

u_tmp <- nu_tmp / rho

beta_tmp <- rnorm(p, 0, 0.01)

### loop 
i <- 1

while( i <= 1000)
{
  d <- u_tmp - gamma_tmp
  
  j <- 1
  
  # while (TRUE)
  #   {
  
  Grad <- gradient(beta_tmp)
  
  H <- likelihood_hessian(beta_tmp)
  
  Hessian <- Reduce('+', H)
  
  Dmat <- rho*(t(A)%*%A) + Hessian
  
  Amat <- matrix(1, nrow = length(beta))
  
  dvec <- t(beta_tmp)%*%Hessian - Grad - rho*t(d) %*% A
  
  beta_new <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$solution
  
  #   criteria <- max(abs((rho*t(A)%*%A + Hessian)%*%beta_new + t(Grad) - Hessian%*%beta_tmp + rho*t(A)%*%d + 1)) 
  #   
  #   if( criteria <= 1e-6 )
  #   {
  #     cat(criteria, '\n',
  #         beta_new[1], '\n')
  #     beta_tmp <- beta_new
  #     break
  #   } else{
  # 
  #   j <- j + 1
  # 
  #   beta_tmp <- beta_new
  # 
  #     }
  # }
  
  beta_tmp <- beta_new
  
  d_tilde <- A %*% beta_tmp + u_tmp
  
  gamma_tmp = ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
  
  u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
  
  if( i %% 100 == 0)
  {
    cat( ' Epoch:: ', i, '\n', 
         'Gradient::', '', '\n',
         'By pylum', '\n', tapply(gamma_tmp[1:82], gl[[1]], sum), '\n', '\n',
         'By pylum & class', '\n', tapply(gamma_tmp[1:82], gl[[2]], sum), '\n', '\n',
         'By pylum & class & order', '\n', tapply(gamma_tmp[1:82], gl[[3]], sum) , '\n','\n',
         '======================================================================================', '\n','\n')
  }
  
  i <- i + 1
  
}

sum(beta_new)
plot(beta, beta_new[1:82])

