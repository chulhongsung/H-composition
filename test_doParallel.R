rm(list = ls())
gc(reset = T)

require(doParallel)

today <- gsub('-', '', substr(Sys.time(), start = 1, stop = 10))

numCores <- 26

cl <- makeCluster(numCores)

registerDoParallel(cores = cl)

setwd("/home/hong/H-composition")

source('test1_simulation.R')

eval(parse(text = paste0('result',today, "<- foreach(file = 1:25, .combine = 'rbind') %dopar% {test_simul(file = file)}" )))

setwd(paste0("/home/hong/simulation_result/seed1/",today))

eval(parse(text = paste0('save(result', today, ',', paste0("file =","'", paste0('result',today,'.RData',"')")))))
