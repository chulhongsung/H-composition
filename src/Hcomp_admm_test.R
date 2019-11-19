if(!require(quadprog)) install.packages("quadprog"); library(quadprog)
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)

args <- commandArgs(TRUE)

file_idx <- as.numeric(args[1])

task_num <- seq((file_idx-1)*100+1, file_idx*100, length.out = 100)  

partition <- expand.grid(seed = 1:20, lambda_1 = seq(0, 1, length.out = 10), lambda_2 = c(0, seq(1e-3, 1e-1, length.out = 13), seq(1e-1, 1, length.out = 7)[2:7]))

eval(parse(text = paste0('test.result', file_idx, ' <- c( )')))

setwd("/home/hong/H-composition/data")

load('bmi.RData')

lv <- lv %>% arrange(lv2,lv3,lv4,lv5,lv6)

gl <- lv %>% select(g1, g2, g3, g4)

m_matrix <- function(i, gl){
  k <- max(gl[[i]])
  m <- matrix(rep(0, k * length(gl[[i]])), ncol = k)
  for ( j in seq_len(k)){
    m[gl[[i]] == j,j] <- 1
  }
  return(m)
}

l = 3

M_matrix_list <- lapply(seq_len(l), function(x){m_matrix(x, gl)}) 

pi_matrix <- lapply(M_matrix_list, function(x){x%*%solve(t(x)%*%x)%*%t(x)})

A = rbind(diag(373), pi_matrix[[1]], pi_matrix[[2]], pi_matrix[[3]])

z = z[,as.character(lv$lv6)]

n = 114

p = 373

for(s in 1:100){
  
  part <- partition[task_num,][s,]
  
  seed <- part$seed
  
  set.seed(seed)
  
  idx <- sample(seq_len(n), 3/10*n)
  
  lambda_1 <- part$lambda_1
  
  lambda_2 <- part$lambda_2
  
  lambda_vec <- c(rep(lambda_1, p), rep(lambda_2, p*l))
  
  train.x <- z[-idx,]
  
  train.y <- y[-idx]
  
  test.x <- z[idx,]
  
  test.y <- y[idx]
  
  rho <- 50
  
  tmp_nu <- rnorm(nrow(A), 0, 0.01)
  
  tmp_gamma <- rnorm(nrow(A), 0, 0.01)
  
  tmp_nu <- rep(0, nrow(A))
  
  tmp_gamma <- rep(0, nrow(A))
  
  tmp_beta <- rep(0, ncol(A))
  
  tmp_u <- tmp_nu/rho
  
  ### Quadratic programming
  
  Dmat <- t(train.x)%*%train.x + rho*(t(A)%*%A)
  
  Amat <- matrix(1, nrow = length(tmp_beta))
  
  i <- 1
  
  while(i <= 10000){
    
    tmp_d <- tmp_u - tmp_gamma 
    
    dvec <- (-1)*(rho*t(tmp_d)%*%A - t(y)%*%train.x)
    
    QP_fit <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)
    
    tmp_beta <- QP_fit$solution
    
    QP_dual <- QP_fit$Lagrangian
    
    kkt1 <- all(abs(Dmat%*%tmp_beta - t(dvec)) - QP_dual < 1e-6)
    
    # sum(tmp_beta)
    
    d_tilde <- A %*% tmp_beta + tmp_u
    
    tmp_gamma <- if_else(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde)*(lambda_vec/rho), 0)
    
    # (abs(rho*tmp_gamma - rho*d_tilde) <= lambda_vec)
    
    kkt2 = near(abs(rho*tmp_gamma - rho*d_tilde)[!(abs(rho*tmp_gamma - rho*d_tilde) <= lambda_vec)], lambda_vec[!(abs(rho*tmp_gamma - rho*d_tilde) <= lambda_vec)])
    
    tmp_u <- tmp_u + (A %*% tmp_beta - tmp_gamma)
    
    kkt3 = max(abs(A %*% tmp_beta - tmp_gamma)) <= 1e-6
    
    kkt <- all(c(kkt1, kkt2, kkt3))
    
    if(kkt == TRUE) break 
    
  }
  
  if(i==10000) cat("Not satisfied kkt condition",'\n')
  
  mse = mean((test.x%*%tmp_beta - test.y)^2)
  
  result_table <- c(seed, lambda_1, lambda_2, mse, t(tmp_beta))
  
  eval(parse(text = paste0('test.result', file_idx, ' <- rbind(test.result', file_idx, ', result_table)')))
  
  cat("Task number", file_idx, ',', s, '% Completed','\n')
}

eval(parse(text = paste0('test.result', file_idx, ' <- as.data.frame(test.result', file_idx, ')' )))

setwd('/home/hong/H-composition/test_result')

eval(parse(text = paste0('save(test.result', file_idx, ',', paste0("file =","'", paste0('test.result',file_idx,'.RData',"')")))))