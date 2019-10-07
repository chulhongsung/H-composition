rm(list = ls())
gc(reset = T)

setwd("C:/Users/UOS/Desktop/H-composition/raw_data")

load('bmi.RData')

if(!require(quadprog)) install.packages("quadprog"); library(quadprog)
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)

lv <- lv %>% arrange(lv2,lv3,lv4,lv5,lv6)

gl <- lv %>% select(g1, g2, g3, g4)

dim(z) # 114 373

z = z[,as.character(lv$lv6)]

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

dim(A) # 1492 X 373

lambda1 <- 0.003

lambda2 <- 0.05

lambda_vec <- c(rep(lambda1, ncol(A)), rep(lambda2, l*ncol(A)))

### Initial 

rho <- 1

tmp_nu <- rnorm(nrow(A), 0, 0.01)

tmp_gamma <- rnorm(nrow(A), 0, 0.01)

tmp_nu <- rep(0, nrow(A))

tmp_gamma <- rep(0, nrow(A))

tmp_beta <- rep(0, ncol(A))

tmp_u <- tmp_nu/rho

### Quadratic programming

Dmat <- t(z)%*%z + rho*(t(A)%*%A)

Amat <- matrix(1, nrow = length(tmp_beta))

tmp_d <- tmp_u - tmp_gamma 

dvec <- (-1)*(rho*t(tmp_d)%*%A - t(y)%*%z)

QP_fit <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)

tmp_beta <- QP_fit$solution

# object <- QP_fit$Lagrangian

d_tilde <- A %*% tmp_beta + tmp_u

tmp_gamma <- if_else(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde)*(lambda_vec/rho), 0)

tmp_u <- tmp_u + (A %*% tmp_beta - tmp_gamma)

### Loss function

(1/2)*norm((y-z%*%tmp_beta),"2")^2 + sum(lambda_vec*abs(tmp_gamma)) +
  rho*t(tmp_u)%*%(A %*% tmp_beta - tmp_gamma) + (rho/2)*norm((A %*% tmp_beta - tmp_gamma), "2")^2

### Hierarchical constraints
tapply(tmp_beta, gl$g1, sum)

tapply(tmp_beta, gl$g2, sum)

tapply(tmp_beta, gl$g3, sum)

