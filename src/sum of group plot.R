rm(list = ls())
gc(reset = T)

require(stringr)
require(reshape2)
require(dplyr)
require(ggplot2)
require(gridExtra)

### set directory where simulation result table downloaded
dir = ''
setwd(dir)

### find simulation result table
inx <- grep(list.files(), pattern = '^simul')
for( i in inx ){ load(list.files()[i])}

b_matrix <- c()

inx2 <- grep(ls(), pattern = '^result')

for (i in inx2+1){b_matrix <- rbind(b_matrix, get(ls()[i]))}

final_result_table <- as.data.frame(b_matrix)

save(final_result_table, file = 'final_result_table.RData')

# load('final_result_table.RData')


### make df for sum of group plot 
plot_table <- melt(data = final_result_table[,1:44], id.vars = c('lambda1', 'lambda2'), value.name = 'sum')

colnames(plot_table)[3] <- 'group'

### fix a specific lambda1 
k = 1
lambda_1 = unique(plot_table$lambda1)[k]

### sum of group pylum plot 
sum_pylum <- ggplot(data = plot_table[plot_table$lambda1 == lambda_1 & plot_table$lambda2 != min(plot_table$lambda2) & grepl('^pylum', plot_table$group),]) + ### min(lambda2) too small making 1/lambda2 too large
  geom_line(aes(x = 1/lambda2, y = sum, col = group, linetype = group), size = 1) +
  geom_hline(yintercept = 0) + scale_x_continuous( breaks = 1/plot_table$lambda2[c(1:5,10,15,20)], labels = as.character(round(1/plot_table$lambda2, 2)[c(1:5,10,15,20)])) +
  theme_classic()

### sum of group class plot 
sum_class <- ggplot(data = plot_table[plot_table$lambda1 == lambda_1 & plot_table$lambda2 != min(plot_table$lambda2) & grepl('^class', plot_table$group),]) +
  geom_line(aes(x = 1/lambda2, y = sum, col = group, linetype = group), size = 1) +
  geom_hline(yintercept = 0) + scale_x_continuous( breaks = 1/plot_table$lambda2[c(1:5,10,15,20)], labels = as.character(round(1/plot_table$lambda2, 2)[c(1:5,10,15,20)])) +
  theme_classic()

### sum of group order plot 
sum_order <- ggplot(data = plot_table[plot_table$lambda1 == lambda_1 & plot_table$lambda2 != min(plot_table$lambda2) & grepl('^order', plot_table$group),]) +
  geom_line(aes(x = 1/lambda2, y = sum, col = group), size = 0.5) +
  geom_hline(yintercept = 0) + scale_x_continuous( breaks = 1/plot_table$lambda2[c(1:5,10,15,20)], labels = as.character(round(1/plot_table$lambda2, 2)[c(1:5,10,15,20)])) +
  theme_classic()

### plot
grid.arrange(sum_pylum, sum_class, sum_order, nrow = 3)