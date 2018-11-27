library(dplyr)
library(quadprog)

setwd("/home/hong/H-composition")

source("function_set.R")

load('simulation_setting_SNU_cv.RData')

test_simul <- function(file)
{
  
  lambda_1 <- seq(1e-3, 2, length.out = 20)[11]
  
  lambda_2 <- c(0, seq(1e-3, 1e-1, length.out = 20),seq(2*1e-1, 1, length.out = 5))[file]
  
  p <- 123; n <- 60; l <- 3; rho <- 0.5
  
  lambda_vec <- c(rep(lambda_1, p), rep(lambda_2, p*l))
  
  gamma_tmp <- rnorm(p*(l+1), 0, 1)
  
  nu_tmp <- rnorm(p*(l+1), 0, 1)
  
  u_tmp <- nu_tmp / rho
  
  beta_tmp <- rep(0, p)
    
  i <- 1
  
  while( i <= 10000)
  {
    d <- u_tmp - gamma_tmp
    
    j <- 1
    
    while (TRUE)
    {
      
      Grad <- gradient(beta_tmp, train.x, train.y)
      
      H <- likelihood_hessian(beta_tmp, train.x, n)
      
      Hessian <- Reduce('+', H)
      
      Dmat <- rho*(t(A)%*%A) + Hessian
      
      Amat <- matrix(1, nrow = length(beta_tmp))
      
      dvec <- t(beta_tmp)%*%Hessian - Grad - rho*t(d) %*% A
      
      beta_new <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$solution
      
      ### Lagrangian multiplier
      object <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$Lagrangian
      
      stationarity_beta <- max(abs((rho*t(A)%*%A + Hessian)%*%beta_new + t(Grad) - Hessian%*%beta_tmp + rho*t(A)%*%d - object))
      
      ### KKT condtion
      if( stationarity_beta <= 1e-6 ) ### Convergence
      {
        beta_tmp <- beta_new
        #cat( j, 'Convergence', '\n')
        break
      } else{
        
        j <- j + 1
        
        beta_tmp <- beta_new
        
        if( j >= 10 ){
          break
        }
        
        
      }
    }
    
    beta_tmp <- beta_new
    
    d_tilde <- A %*% beta_tmp + u_tmp
    
    gamma_tmp <- ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
    
    stationarity_gamma <- all(abs(rho*(gamma_tmp - d_tilde)) <= lambda_vec)
    
    u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
    
    i <- i + 1
    
  }
  
  prop <- test.x %*% beta_tmp
  
  yhat <- if_else(prop >= 0.5, 1, 0)
  
  error_rate <- mean(yhat != test.y)
  
  test1.result <- t(c(error_rate, lambda_1, lambda_2, t(beta_tmp)))
  
  test1.result <- as.data.frame(test1.result)
  
  test1.result

}
