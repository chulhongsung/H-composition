rm(list = ls())
gc(reset = T)

require(ggplot2)
require(reshape2)
require(gridExtra)
library(stringr)
library(dplyr)
library(kableExtra)

### set directory where cv simulation error_rate file downloaded
dir2 = '' 

setwd(dir2)

### find error_rate + 'file'.RData  
inx <- grep(list.files(), pattern = '^error')
for( i in inx ){ load(list.files()[i])}

b_matrix <- c()

inx2 <- grep(ls(), pattern = '^error')

### rbind error_rate + 'file'.RData 
for (i in inx2){b_matrix <- rbind(b_matrix, get(ls()[i]))}

error_table <- as.data.frame(b_matrix)

colnames(error_table) <- c('lambda1', 'lambda2', 'error_rate1', 'error_rate2', 'error_rate3', 'error_rate4', 'error_rate5')

### compute mean_error_rate
error_table <- error_table %>% mutate(mean_error = rowSums(error_table[,3:7])/5)

save(error_table, file = 'final_error_table.RData')
