args <- (commandArgs(TRUE))

start <- Sys.time()

setwd("~")
load('simulation_data.RData')

source("function_set.R")
require(dplyr)
require(quadprog)

today <- gsub('-', '', substr(Sys.time(), start = 1, stop = 10))

eval(parse(text = args[[1]]))

task_num <- seq((file_idx-1)*50+1, file_idx*50, length.out = 50)  

partition <- expand.grid(seed = 1:20, lambda_1 = seq(0, 2, length.out = 10), lambda_2 = c(0, seq(1e-3, 1e-1, length.out = 13), seq(1e-1, 1, length.out = 7)[2:7]))

eval(parse(text = paste0('test.result', file_idx, ' <- c( )')))

p <- 123; n <- 60; l <- 3; rho <- 0.5

for ( i in 1:50 )
{
  gamma_tmp <- rnorm(p*(l+1), 0, 1)
  
  nu_tmp <- rnorm(p*(l+1), 0, 1)
  
  u_tmp <- nu_tmp / rho
  
  beta_tmp <- rep(0, p)
  
  set <- settting_fun(partition, i, newz, p, l)
  
  train.x <- set$train.x
  
  test.x <- set$test.x
  
  train.y <- set$train.y
  
  test.y <- set$test.y
  
  lambda_vec <- set$lambda_vec 
  
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
        
        break
        
      } else{
        
        j <- j + 1
        
        beta_tmp <- beta_new
        
        if( j >= 5 ){
          break
        }
        
        
      }
    }
    
    beta_tmp <- beta_new
    
    d_tilde <- A %*% beta_tmp + u_tmp
    
    gamma_tmp <- ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
    
    u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
    
    i <- i + 1
    
  }
  
  prop <- test.x %*% beta_tmp
  
  yhat <- if_else(prop >= 0.5, 1, 0)
  
  error_rate <- mean(yhat != test.y)
  
  test.result <- t(c(error_rate, set$lambda_1, set$lambda_2, t(beta_tmp)))
  
  eval(parse(text = paste0('test.result', file_idx, ' <- rbind(test.result', file_idx, ', test.result)')))
  
}

eval(parse(text = paste0('test.result', file_idx, ' <- as.data.frame(test.result)')))

dir.create(paste0("~",today))

setwd(paste0("~",today))

eval(parse(text = paste0('save(test.result', file_idx, ',', paste0("file =","'", paste0('test.result',file_idx,'.RData',"')")))))

end <- Sys.time() 
  
end - start