rm(list = ls())
gc(reset = TRUE)

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(quadprog)){ install.packages('quadprog')}; require(quadprog)

load('H-composition simulation.RData')

### Loss function for beta

# Loss <- function(beta, X, Y){
#   L  <- -t(Y) %*% (X %*% beta) + colSums(log(1+exp(X%*%beta)))
#   return(L)
# }

### Gradient
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


### Simulation loop

p = 82; l = 3; n = 1000

result_table <- c()

for (lambda_1 in c(0,seq_len(5)))
{
  for ( lambda_2 in c(0,seq_len(20)))
  {
    rho <- 0.5
    
    lambda_vec <- c(rep(lambda_1, p), rep(lambda_2, p*3))
    
    set.seed(2018)
    
    gamma_tmp <- rnorm(p*4, 0, 1)
    
    nu_tmp <- rnorm(p*4, 0, 1)
    
    u_tmp <- nu_tmp / rho
    
    beta_tmp <- rep(0, p)
    
    i <- 1
    
    while( i <= 10000)
    {
      d <- u_tmp - gamma_tmp
      
      j <- 1
      
      while (TRUE)
      {
        
        Grad <- gradient(beta_tmp, X, Y)
        
        H <- likelihood_hessian(beta_tmp, X, n)
        
        Hessian <- Reduce('+', H)
        
        Dmat <- rho*(t(A)%*%A) + Hessian
        
        Amat <- matrix(1, nrow = length(beta))
        
        dvec <- t(beta_tmp)%*%Hessian - Grad - rho*t(d) %*% A
        
        beta_new <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$solution
        
        ### Lagrangian multiplier 
        object <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$Lagrangian
        
        stationarity_beta <- max(abs((rho*t(A)%*%A + Hessian)%*%beta_new + t(Grad) - Hessian%*%beta_tmp + rho*t(A)%*%d - object))
        
        ### KKT condtion
        if( stationarity_beta <= 1e-6 ) ### Convergence
        {
          beta_tmp <- beta_new
          # cat( j, 'Convergence', '\n')
          break
        } else{
          
          j <- j + 1
          
          beta_tmp <- beta_new
          
        }
      }
      
      beta_tmp <- beta_new
      
      d_tilde <- A %*% beta_tmp + u_tmp
      
      gamma_tmp <- ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
      
      stationarity_gamma <- all(abs(rho*(gamma_tmp - d_tilde)) <= lambda_vec)
      
      u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
      
      if( i %% 10000 == 0)
      {
        cat( ' Epoch:: ', i, '\n', 
             'Beta KKT condition Stationarity::', stationarity_beta, '\n',
             'Gamma KKT condition Stationarity::', stationarity_gamma, '\n', '\n',
             'By pylum', '\n', tapply(beta_tmp[1:82], gl[[1]], sum), '\n', '\n',
             'By pylum & class', '\n', tapply(beta_tmp[1:82], gl[[2]], sum), '\n', '\n',
             'By pylum & class & order', '\n', tapply(beta_tmp[1:82], gl[[3]], sum) , '\n','\n',
             '======================================================================================', '\n','\n')
      }
      
      i <- i + 1
      
    }
    result_table <- rbind(result_table, c(lambda_1, lambda_2, tapply(beta_tmp, gl[[1]], sum), tapply(beta_tmp, gl[[2]], sum), tapply(beta_tmp, gl[[3]], sum)))
  }
  
}

save(result_table, file = 'simulation_result.RData')

# Check 