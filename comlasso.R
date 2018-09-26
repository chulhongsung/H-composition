rm(list = ls())
gc(reset = TRUE)

#load(data)

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)

## data preprocessing

# g <-  bf_x.o[selv,]
# 
# ng <- g %>% select(genus) %>% group_by(genus) %>% tally() %>% filter(n ==1)
# 
# g <- g[!(g$genus %in% ng$genus), ]
# 
# g$g1 <- g %>% group_by(pylum) %>% group_indices() #5
# g$g2 <- g %>% group_by(pylum, class) %>% group_indices() #11
# g$g3 <- g %>% group_by(pylum, class, order) %>% group_indices() #16
# g$g4 <- g %>% group_by(pylum, class, order, family) %>% group_indices() #39
# g$g5 <- g %>% group_by(pylum, class, order, family, genus) %>% group_indices() #89
# 
# g <- g[g$g2 != 2,]
# g <- g[g$g4 != 32,]
# 
# ug <- g %>% distinct(pylum, class, order, family, genus)%>% arrange(pylum, class, order, family, genus)
# 
# ug$g1 <- ug %>% group_by(pylum) %>% group_indices() #5
# 
# ug <- ug[ug$g1 != 4,]
# 
# ug$g2 <- ug %>% group_by(pylum, class) %>% group_indices() #11
# 
# ug <- ug[ug$g2 != 7 & ug$g2 != 8,]
# 
# ug$g3 <- ug %>% group_by(pylum, class, order) %>% group_indices() #16
# 
# ug <- ug[ug$g3 != 2 & ug$g3 != 11,]
# 
# ug$g4 <- ug %>% group_by(pylum, class, order, family) %>% group_indices() #39
# 
# ug$g5 <- ug %>% group_by(pylum, class, order, family, genus) %>% group_indices() #87

set.seed(526)

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

### Loss function for beta

Loss <- function(beta){
  L  <- -t(Y) %*% (Z %*% beta) + colSums(log(1+exp(Z%*%beta))) + (rho/2)*(t(A%*%beta)%*%(A%*%beta) + 2*t(d)%*%A%*%beta)
  return(L)
}

### Newton's method under equality constrains

### Gradient
gradient <- function(beta, d_tmp)
{
  G <- -t(Y -(1/(1+exp(-Z%*%beta))))%*%Z + t(rho*(t(A)%*%(A%*%beta + d_tmp)))
  return(G)
}

### Hessian
likelihood_hessian <- function(beta)
{
  H <- list()
  for( i in seq(n))
  {
    S <- as.numeric(1/(1+exp(-Z[i,]%*%beta)))
    H[[i]] <- (Z[i,]%*%t(Z[i,])) * S * (1-S)
    
  }
  return(H)
}

### 
rho <- 0.5

lambda_1 <- 3

lambda_2 <- 3

lambda_vec = c(rep(lambda_1, p), rep(lambda_2, p*3))

### initial gamma, nu, u and beta matrix

gamma_tmp <- rnorm(p*4, 0, 1)

nu_tmp <- rnorm(p*4, 0, 1)

u_tmp <- nu_tmp / rho

beta_tmp <- rep(c(1,-1), 41)

### constraint
cons <- rep(1, 82)

### backtracking line search
bt_size <- 0.5
a <- 0.5

### loop 
i <- 1

while( i <= 500)
{
  d <- u_tmp - gamma_tmp
  
  G <- gradient(beta_tmp, d)
  
  H <- likelihood_hessian(beta_tmp)
  
  Hessian <- Reduce('+', H) + rho*t(A)%*%A
  
  ### Newton's method in linear equality constraints
  
  mat_tmp <- rbind(cbind(Hessian, c(Loss(beta_tmp))*cons), cbind(t(cons),0))
  
  v <- - solve(mat_tmp) %*% rbind(t(G), t(cons)%*%beta_tmp)
  
  ### backtracking line search
  t <- 1
  
  while (TRUE) {
    
    beta_candidate <- beta_tmp + t * v[1:82]
    
    if(Loss(beta_candidate)[1,1] <= Loss(beta_tmp) + a*t*G%*%v[1:82])
      {
        break
      } else{ 
        t <- bt_size*t
      }
  }
  
  beta_tmp <- beta_candidate
  
  d_tilde <- A %*% beta_tmp + u_tmp
  
  gamma_tmp = ifelse(abs(d_tilde) > lambda_vec /rho, d_tilde - sign(d_tilde) * (lambda_vec/rho), 0)
  
  u_tmp <- u_tmp + (A %*% beta_tmp - gamma_tmp)
  
  if( i %% 50 == 0)
  {
    cat( ' Epoch:: ', i, '\n', 
         'Gradient::', max(abs(G)), '\n',
         'Loss::', Loss(beta_tmp) + lambda_1*sum(abs(gamma_tmp[1:82])) + lambda_2*sum(abs(gamma_tmp[83:length(gamma_tmp)])) +
           t(rho*u_tmp) %*% (A%*%beta_tmp - gamma_tmp) + (rho/2)*t((A%*%beta_tmp - gamma_tmp))%*%(A%*%beta_tmp - gamma_tmp), '\n', '\n',
         'By pylum', '\n', tapply(beta_tmp[1:82], gl[[1]], sum), '\n', '\n',
         'By pylum & class', '\n', tapply(beta_tmp[1:82], gl[[2]], sum), '\n', '\n',
         'By pylum & class & order', '\n', tapply(beta_tmp[1:82], gl[[3]], sum), '\n','\n',
         '===================================================================================================================', '\n','\n')
  }
  
  i <- i + 1

}

plot(beta, beta_tmp[1:82])
