rm(list = ls())
gc(reset = T)

### required library
library(dplyr)
library(quadprog)

### set directory 
dir = '' 
setwd(dir)

source("function_set.R")

load('simulation_setting_SNU_cv.RData')

### 5-fold cross-validation

set.seed(1)

id <- rep(1:5, length.out = 60)

id <- sample(id)

lambda_1 <- 0.05

### file : lambda2 indexing 
file = 1

eval(parse(text = paste0('error_rate', file, ' <- c()')))

lambda_2 <- seq(1e-3, 1e-1, length.out = 20)[file]

### setting 
p <- 123; n <- 48; l <- 3; rho <- 0.5

lambda_vec = c(rep(lambda_1, p), rep(lambda_2, p*l))

gamma_tmp <- rnorm(p*(l+1), 0, 1)

nu_tmp <- rnorm(p*(l+1), 0, 1)

u_tmp <- nu_tmp / rho

beta_tmp <- rep(0, p)

### simulation loop
for ( k in seq_len(5))
{
  ### 5-fold cross-validation set
  train <- train.x[id != k, ]
  y.train <- train.y[id != k]
  
  valid <- train.x[id == k, ] 
  y.valid <- train.y[id == k]
  
  i <- 1
  
  while( i <= 10000)
  {
    d <- u_tmp - gamma_tmp
    
    j <- 1
    
    while (TRUE)
    {
      
      Grad <- gradient(beta_tmp, train, y.train)
      
      H <- likelihood_hessian(beta_tmp, train, n)
      
      Hessian <- Reduce('+', H)
      
      Dmat <- rho*(t(A)%*%A) + Hessian
      
      Amat <- matrix(1, nrow = length(beta_tmp))
      
      dvec <- t(beta_tmp)%*%Hessian - Grad - rho*t(d) %*% A
      
      beta_new <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$solution
      
      ### Lagrangian multiplier
      object <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$Lagrangian
      
      stationarity_beta <- max(abs((rho*t(A)%*%A + Hessian)%*%beta_new + t(Grad) - Hessian%*%beta_tmp + rho*t(A)%*%d - object))
      
      ### beta KKT condtion
      if( stationarity_beta <= 1e-6 ) ### Convergence
      {
        beta_tmp <- beta_new
        #cat( j, 'Convergence', '\n')
        break
      } else{
        
        j <- j + 1
        
        beta_tmp <- beta_new
        
        if( j >= 10 ){
          #cat(stationarity_beta, '\n','Not Convergence', '\n')
          break
        }
        
        
      }
    }
    
    beta_tmp <- beta_new
    
    d_tilde <- A %*% beta_tmp + u_tmp
    
    gamma_tmp = ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
    
    ### gamma KKT condition
    stationarity_gamma = all(abs(rho*(gamma_tmp - d_tilde)) <= lambda_vec)
    
    u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
    
    if( i %% 10000 == 0)
    {
      cat( ' Epoch:: ', i, '\n',
           'Beta KKT condition Stationarity::', stationarity_beta, '\n',
           'Gamma KKT condition Stationarity::', stationarity_gamma, '\n', '\n',
           'By pylum', '\n', tapply(beta_tmp[1:123], gl[[1]], sum), '\n', '\n',
           'By pylum & class', '\n', tapply(beta_tmp[1:123], gl[[2]], sum), '\n', '\n',
           'By pylum & class & order', '\n', tapply(beta_tmp[1:123], gl[[3]], sum) , '\n','\n',
           '======================================================================================', '\n','\n')
    }
    
    i <- i + 1
    
  }
  
  prop <- valid %*% beta_tmp
  yhat <- if_else(prop >= 0.5, 1, 0)
  
  ### error_rate
  error_rate <- mean(yhat != y.valid)
  
  eval(parse(text = paste0('error_rate', file, ' <- cbind(error_rate', file, ',error_rate)'))) 
}

eval(parse(text = paste0('error_rate', file, ' <- as.data.frame(cbind(lambda_1, lambda_2, error_rate', file,'))')))

### set save directory
dir2 <- ''
setwd(dir2)
eval(parse(text = paste0('save(error_rate', file, ',', paste0("file =","'", paste0('error_rate',file,'_20181107.RData',"')")))))


