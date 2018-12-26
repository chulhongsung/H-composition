start <- Sys.time()

setwd("/home/hong/H-composition")
load('simulation_data.RData')

source("function_set.R")

require(dplyr)
require(quadprog)

args <- (commandArgs(TRUE))

file_idx <- eval(parse(text=args[1]))

task_num <- seq((file_idx-1)*100+1, file_idx*100, length.out = 100)  

partition <- expand.grid(seed = 1:100, lambda_1 = seq(0, 2, length.out = 10), lambda_2 = c(0, seq(1e-3, 1e-1, length.out = 13), seq(1e-1, 1, length.out = 7)[2:7]))

eval(parse(text = paste0('cv.result', file_idx, ' <- c( )')))

p <- 123; n <- 48; l <- 3; rho <- 0.5

for (s in 1:100)
{
  gamma_tmp <- rnorm(p*(l+1), 0, 1)
  
  nu_tmp <- rnorm(p*(l+1), 0, 1)
  
  u_tmp <- nu_tmp / rho
  
  beta_tmp <- rep(0, p)
  
  set <- setting_cv_fun(partition, task_num, s, newz, p, l)
  
  seed <- set$seed
  
  cv.x <- set$cv.x
  
  cv.y <- set$cv.y
  
  lambda_1 <- set$lambda_1
  
  lambda_2 <- set$lambda_2
  
  lambda_vec <- set$lambda_vec 
  
  idx <- set$idx
  
  error_rate <- c(seed, lambda_1, lambda_2, rep(0, 5))
  
  for ( k in seq_len(5))
  {
    
    train.x <- cv.x[-idx[[k]],]
    
    train.y <- cv.y[-idx[[k]]]
    
    valid.x <- cv.x[idx[[k]],]
    
    valid.y <- cv.y[idx[[k]]]
    
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
    prop <- valid.x %*% beta_tmp
    
    yhat <- if_else(prop >= 0.5, 1, 0)
    
    error_rate[k+3] <- mean(yhat != valid.y)
    
  }
  
  eval(parse(text = paste0('cv.result', file_idx, ' <- rbind(cv.result', file_idx, ', error_rate)')))
}

eval(parse(text = paste0('cv.result', file_idx, ' <- as.data.frame(cv.result', file_idx, ')' )))

setwd('/home/hong/20181213cv')

eval(parse(text = paste0('save(cv.result', file_idx, ',', paste0("file =","'", paste0('cv.result',file_idx,'.RData',"')")))))

end <- Sys.time() 

end - start
