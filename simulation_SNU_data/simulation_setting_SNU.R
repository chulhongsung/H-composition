rm(list = ls())
gc(reset = T)

setwd("C:/Users/UOS/Desktop/data/H-composition")
load('SNU.RData')

library(dplyr)

lv2 <- bf_x.o$pylum[selv]

lv3 <- bf_x.o$class[selv]

lv4 <- bf_x.o$order[selv]

lv5 <- bf_x.o$family[selv]

lv6 <- bf_x.o$genus[selv]

g <- cbind(lv2, lv3, lv4, lv5, lv6)

g <- data.frame(g)

g <- g %>% arrange(lv2, lv3, lv4, lv5, lv6)

### Odoribacter,  Oscillibacter, Parasutterella 
## g[g$lv6 == 'Odoribacter',] Bacteroidetes Bacteroidia Bacteroidales Porphyromonadaceae Odoribacter
## g[g$lv6 == 'Oscillibacter',] Firmicutes Clostridia Clostridiales Ruminococcaceae Oscillibacter
## g[g$lv6 == 'Parasutterella',] Proteobacteria Betaproteobacteria Burkholderiales Sutterellaceae Parasutterella

### lv6는 똑같은데 앞선 lv이 다른 경우를 같게 만드는 작업
g[g$lv6 == 'Odoribacter' & g$lv3 == 'Bacteroidetes',1:4] <- c('Bacteroidetes', 'Bacteroidia', 'Bacteroidales', 'Porphyromonadaceae')
g[g$lv6 == 'Oscillibacter' & g$lv3 == 'Bacilli', 1:4] <- data.frame(lv2 = c("Firmicutes", "Firmicutes"), lv3 = c('Clostridia','Clostridia'),
                                                                    lv4 = c('Clostridiales', 'Clostridiales'), lv5= c('Ruminococcaceae', 'Ruminococcaceae'))
g[g$lv6 == 'Parasutterella' & g$lv5 == 'Sutterellaceae', 1:4] <- c('Proteobacteria', 'Betaproteobacteria', 'Burkholderiales', 'Alcaligenaceae')

g <- g %>% distinct(lv2, lv3, lv4, lv5, lv6)

g$g2 <- g %>% group_by(lv2) %>% group_indices()
g$g3 <- g %>% group_by(lv2, lv3) %>% group_indices()
g$g4 <- g %>% group_by(lv2, lv3, lv4) %>% group_indices()
g$g5 <- g %>% group_by(lv2, lv3, lv4, lv5) %>% group_indices()
g$g6 <- g %>% group_by(lv2, lv3, lv4, lv5, lv6) %>% group_indices()

p = 123; l = 3; n = 83  

gl <- list()
gl[[1]] <- g$g2
gl[[2]] <- g$g3
gl[[3]] <- g$g4
gl[[4]] <- g$g5

m_matrix <- function(i, gl){
  k <- length(unique(gl[[i]]))
  m <- matrix(rep(0, k * length(gl[[i]])), ncol = k)
  for ( j in seq_len(k)){
    m[gl[[i]] == j,j] <- 1
  }
  return(m)
}

newz <- log(z[,as.character(g$lv6)])

M_matrix_list <- lapply(seq_len(l), function(x) m_matrix(x, gl) ) 

pi_matrix <- lapply(M_matrix_list, function(x){ x %*% solve( t(x) %*% x ) %*% t(x)})

I <- diag(1, nrow = p)

A <- rbind(I, pi_matrix[[1]], pi_matrix[[2]], pi_matrix[[3]])

### y = (-1, 1) -> (0, 1)
y <- if_else(y == 1, 1, 0)

### total(83), train(60), test set(23)
set.seed(1)

id <- sample(1:83, 60) 

train.x <- newz[id,]
test.x <- newz[-id,]

train.y <- y[id]
test.y <- y[-id]

save(gl, train.x, test.x, train.y, test.y , A, file = 'simulation_setting_SNU_cv.RData')