rm(list = ls())
gc(reset = TRUE)

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(quadprog)){ install.packages('quadprog')}; require(quadprog)

load('H-composition3.RData')

set.seed(2018)

gl <- list()
gl[[1]] <- ug %>% group_by(pylum) %>% group_indices() ### k1 = 5
gl[[2]] <- ug %>% group_by(pylum, class) %>% group_indices() ### k2 = 11
gl[[3]] <- ug %>% group_by(pylum, class, order) %>% group_indices() ### k3 = 15
gl[[4]] <- ug %>% group_by(pylum, class, order, family) %>% group_indices() ### k4 = 38
gl[[5]] <- ug %>% group_by(pylum, class, order, family, genus) %>% group_indices() ### k5 = 87

p = 82; l = 3; n = 1000

#M_matrix
m_matrix <- function(i){
  k <- max(gl[[i]])
  m <- matrix(rep(0, k * length(gl[[i]])), ncol = k)
  for ( j in seq_len(k)){
    m[gl[[i]] == j,j] <- 1
  }
  return(m)
}

M_matrix_list <- lapply(seq_len(l), m_matrix ) 

#Pi matrix
pi_matrix <- lapply(M_matrix_list, function(x){ x %*% solve( t(x) %*% x ) %*% t(x)})

### Z
Z <- matrix(rnorm(p * n, mean = 0, sd = 1), nrow = n)

### real beta
beta <-  c()

beta_generator <- function(l)
{
  for ( i in seq_len(max(gl[[l]])))
  {
    if(sum(gl[[l]] == i) %% 2 == 0){
      beta_j <- rep(0, sum(gl[[l]] == i))
      beta <- c(beta, beta_j)
    } else{ 
      beta_j <- rep(c(1,-1), sum(gl[[l]] == i)/2 )
      beta <- c(beta, beta_j, 0)
    }
  }
  return(beta)
}

beta <- beta_generator(3)

### Y logistic distribution
prob <- exp(Z %*% beta)/(1 + exp(Z %*% beta))

Y <- rbinom(n = n, size = 1,prob = prob)

### X matrix
alpha <- c()

for ( i in seq_len(5))
{
  alpha[i] <- rexp(n = 1, rate = i)
  eval(parse(text = paste0('x', i, ' <- alpha[', i, '] * exp(Z[, gl[[1]] ==', i,'])')))
}

X <- cbind(x1,x2,x3,x4,x5)
X <- log(X)

### A matrix
I <- diag(1, nrow = p)

A <- rbind(I, pi_matrix[[1]], pi_matrix[[2]], pi_matrix[[3]])



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

set.seed(2018)

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

  while (TRUE)
    {

    Grad <- gradient(beta_tmp)

    H <- likelihood_hessian(beta_tmp)
    
    Hessian <- Reduce('+', H)
    
    Dmat <- rho*(t(A)%*%A) + Hessian
    
    Amat <- matrix(1, nrow = length(beta))
    
    dvec <- t(beta_tmp)%*%Hessian - Grad - rho*t(d) %*% A
    
    beta_new <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$solution
    
    ### Lagrangian multiplier
    object <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)$Lagrangian
      
    criteria <- max(abs((rho*t(A)%*%A + Hessian)%*%beta_new + t(Grad) - Hessian%*%beta_tmp + rho*t(A)%*%d - object))
    
    ### KKT condition
    if( criteria <= 1e-6 ) 
    { 
      ### Convergence
      beta_tmp <- beta_new
      break 
    } else{

    j <- j + 1

    beta_tmp <- beta_new

      }
  }
  
  beta_tmp <- beta_new
  
  d_tilde <- A %*% beta_tmp + u_tmp
  
  gamma_tmp = ifelse(abs(d_tilde) > lambda_vec/rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
  
  u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
  
  if( i %% 100 == 0)
  {
    cat( ' Epoch:: ', i, '\n', 
         'Gradient::', '', '\n',
         'By pylum', '\n', tapply(beta_tmp[1:82], gl[[1]], sum), '\n', '\n',
         'By pylum & class', '\n', tapply(beta_tmp[1:82], gl[[2]], sum), '\n', '\n',
         'By pylum & class & order', '\n', tapply(beta_tmp[1:82], gl[[3]], sum) , '\n','\n',
         '======================================================================================', '\n','\n')
  }
  
  i <- i + 1

}

sum(beta_new)
plot(beta, beta_new[1:82])
cbind(gamma_tmp[1:82], beta)

