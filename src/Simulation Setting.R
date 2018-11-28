rm(list = ls())
gc(reset = T)

require(dplyr)

#load(data)
load('SNU.RData')

## data preprocessing

g <- bf_x.o[selv,]

ng <- g %>% dplyr::select(genus) %>% group_by(genus) %>% tally() %>% filter(n ==1)

g <- g[!(g$genus %in% ng$genus), ]

g$g1 <- g %>% group_by(pylum) %>% group_indices() #5
g$g2 <- g %>% group_by(pylum, class) %>% group_indices() #11
g$g3 <- g %>% group_by(pylum, class, order) %>% group_indices() #16
g$g4 <- g %>% group_by(pylum, class, order, family) %>% group_indices() #39
g$g5 <- g %>% group_by(pylum, class, order, family, genus) %>% group_indices() #89

g <- g[g$g2 != 2,]
g <- g[g$g4 != 32,]

ug <- g %>% distinct(pylum, class, order, family, genus)%>% arrange(pylum, class, order, family, genus)

ug$g1 <- ug %>% group_by(pylum) %>% group_indices()
ug <- ug[ug$g1 != 4,]

ug$g2 <- ug %>% group_by(pylum, class) %>% group_indices() 
ug <- ug[ug$g2 != 7 & ug$g2 != 8,]

ug$g3 <- ug %>% group_by(pylum, class, order) %>% group_indices() 
ug <- ug[ug$g3 != 2 & ug$g3 != 11,]

ug$g4 <- ug %>% group_by(pylum, class, order, family) %>% group_indices() 
ug$g5 <- ug %>% group_by(pylum, class, order, family, genus) %>% group_indices() 

gl <- list()
gl[[1]] <- ug %>% group_by(pylum) %>% group_indices() 
gl[[2]] <- ug %>% group_by(pylum, class) %>% group_indices() 
gl[[3]] <- ug %>% group_by(pylum, class, order) %>% group_indices()
gl[[4]] <- ug %>% group_by(pylum, class, order, family) %>% group_indices() 
gl[[5]] <- ug %>% group_by(pylum, class, order, family, genus) %>% group_indices() 

## Set simulation
set.seed(1)

p = 82; l = 3; n = 1000

#M_matrix
m_matrix <- function(i, gl){
  k <- max(gl[[i]])
  m <- matrix(rep(0, k * length(gl[[i]])), ncol = k)
  for ( j in seq_len(k)){
    m[gl[[i]] == j,j] <- 1
  }
  return(m)
}

M_matrix_list <- lapply(seq_len(l), function(x){m_matrix(x, gl)}) 

#Pi matrix
pi_matrix <- lapply(M_matrix_list, function(x){ x %*% solve( t(x) %*% x ) %*% t(x)})

### Z
Z <- matrix(rnorm(p * n, mean = 0, sd = 1), nrow = n)

### real beta


beta_generator <- function(l, gl)
{
  beta <-  c()
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

beta <- beta_generator(3, gl)

beta <- beta/10

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

save(gl, Z, X, Y, A, beta,  file = 'H-composition simulation.RData')




