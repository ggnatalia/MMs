
library(ggplot2)
library(reshape2)

projectName = 'Freshwaters_master'
enviro = 'Aquatic.Freshwaters'

originFile = 'SpeciesperEnviro.tsv'
origin = read.table(originFile, header = T, row.names=1, stringsAsFactors = F, check.names = F)
env = origin[, enviro, drop = F]
env[env==0] = 0.005 # we add this to include a little probability of absence taxons


# Create new df
df = as.data.frame(matrix(0, nrow = nrow(env), ncol = 101))
rownames(df) = rownames(env)
colnames(df) = c('origin', paste0('R_', c(1:100)))
df$origin = env[,1]



for (i in 1:100){
    finalFile = paste0(projectName , i, '.pctax.tsv')
    final = t(read.table(finalFile , sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F))

    for (tx in rownames(df)){
        if (tx %in% rownames(final)){
            df[tx, paste0('R_', i)] = final[tx,'counts']
        } else { df[tx, paste0('R_', i)] = 0 }
    } 
}

dfnorm = as.data.frame(apply(df, 2, function(x) {x/sum(x)}*100))
dfnorm_sort = dfnorm[order(dfnorm$origin, decreasing = T),]



mean_of_replicas = as.data.frame(apply(dfnorm[, c(2:101)], 1, mean))
sd_of_replicas = as.data.frame(apply(dfnorm[, c(2:101)], 1, sd))

data2plot = cbind(dfnorm[, 'origin'], mean_of_replicas, sd_of_replicas, rownames(dfnorm))
colnames(data2plot) = c('origin', 'mean', 'sd', 'tax')





sink(paste0(projectName , '.', 'scatter_lm_taxonomy_subset_comparisons_100replicas.txt'))
print(summary(lm(data2plot$origin ~ data2plot$sd)))
sink()


p = ggplot(data2plot, aes(origin, mean)) + theme_classic()
#p = p + geom_point() + geom_text(aes(label=tax),hjust=0, vjust=0)
p = p + geom_point() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
p = p + geom_text(data = data2plot[data2plot$mean > 1,], aes(origin, mean, label = tax), hjust = 0, vjust = 0)
# Add regression line
p = p + geom_smooth(method = lm, formula = y ~ x)
p = p + labs(title = 'Subsetting taxa for an specific environment: \'Aquatic Freshwaters\'', x = 'Frequency of the different taxa in the reference of \'Aquatic Freshwaters\'', y = 'Mean of taxa abundance in 100 mocks of \'Aquatic Freshwaters\'')
p = p + theme(plot.title = element_text(hjust = 0.5))
p


ggsave(file = paste0(projectName , '.', 'scatter_lm_taxonomy_subset_comparisons2.svg'), plot = p, width = 10, height = 8)

dfnorm_sort$tax = rownames(dfnorm_sort)
data2plot2 = melt(dfnorm_sort, id.vars = c('origin', 'tax'))

p = ggplot(data2plot2, aes(x = as.factor(origin), y = value, label = tax, color = variable)) + theme_classic()
#p = p + geom_point() + geom_text(aes(label=tax),hjust=0, vjust=0)
p = p + geom_boxplot() 
p = p + labs(title = 'Subsetting taxa for an specific environment: \'Aquatic Freshwaters\'', x = 'Frequency of the different taxa in the reference of \'Aquatic Freshwaters\'', y = 'Mean of taxa abundance in 100 mocks of \'Aquatic Freshwaters\'')
p = p + theme(plot.title = element_text(hjust = 0.5))
p



plot()
