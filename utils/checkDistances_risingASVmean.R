# Check number os ASVs per species
library('reshape2')
library('ggplot2')
library('dplyr')
library('data.table')



data2plot = list()
d2plot_basic = as.data.frame(matrix(0, nrow = 100, ncol = 4))
colnames(d2plot_basic) = c('R','meanASVs', 'sdASVs', 'originalASVs')
projectName = 'Freshwaters2_master_asvsmean'


j = 1
for (asv in 1:10){
    for (i in 1:10){ ###### trial
    ASV = paste0('ASV_', asv)
    R = paste0('R_', i)
    finalFile = paste0(projectName , '_', asv, '.', i, '.distances.tsv')
    distances = read.table(finalFile, row.names = 1, header = T, stringsAsFactors = F,  sep = '\t', check.names = F)
    
    ind = which(upper.tri(distances, diag = TRUE), arr.ind = TRUE)
    M = cbind(ind, distances[ind])
    
    data = as.data.frame(M)
    data[,'row'] = rownames(distances)[unname(data[,'row'])]
    data[,'col'] = colnames(distances)[unname(data[,'col'])]
    colnames(data) = c('seq1','seq2','d')
    
    data$species1 = paste0(sapply(strsplit(data$seq1, '[.]'), "[", 1), '.', sapply(strsplit(data$seq1, '[.]'), "[", 2))
    data$species2 = paste0(sapply(strsplit(data$seq2, '[.]'), "[", 1), '.', sapply(strsplit(data$seq2, '[.]'), "[", 2))
    
    
    # Remove seq1==seq2
    
    data2 = data[!(data$seq1==data$seq2),]
    
    df = as.data.frame(matrix(0, nrow = length(levels(as.factor(data2$species1))), ncol = 6))
    colnames(df) = c('Sp', 'max_distance', 'min_distance', 'nASVsperSp', 'ASVs', 'OriginalASVs')
    rownames(df) = levels(as.factor(data2$species1))
    
    for (sp in levels(as.factor(data2$species1))){
        if (dim(data2[data2$species1 == sp & data2$species2 == sp,])[1] > 0){
            max_distance = max(data2[data2$species1 == sp & data2$species2 == sp,'d'])
            min_distance = min(data2[data2$species1 == sp & data2$species2 == sp,'d'])
            dumb = data2[data2$species1 == sp & data2$species2 == sp, c('seq1','seq2')]
            nASVsperSp = length(unique(c(dumb$seq1, dumb$seq2)))
            
            df[sp, 'Sp'] = sp
            df[sp,'max_distance'] = max_distance
            df[sp, 'min_distance'] = min_distance
            df[sp, 'nASVsperSp'] = as.numeric(nASVsperSp) #not count the original reference 
            df[sp, 'ASVs'] = paste(unique(c(dumb$seq1, dumb$seq2)), collapse = ',')
            df[sp, 'OriginalASVs'] = asv
        }
        else{df[sp,] = c(sp, 0, 0, 1,sp, asv)}
    }
    data2plot[[ASV]][[R]] = df
    d2plot_basic[j, ] = c(R, mean(as.numeric(df$nASVsperSp), na.rm=TRUE), sd(as.numeric(df$nASVsperSp), na.rm=TRUE), df$OriginalASVs)
    j=j+1
}
}




d2plot_basic$meanASVs = as.numeric(d2plot_basic$meanASVs)
d2plot_basic$sdASVs = as.numeric(d2plot_basic$sdASVs)
d2plot_basic$originalASVs = factor(d2plot_basic$originalASVs, levels = c(1:10))


sink(paste0(projectName , '.', 'average_ASVs.txt'))
print(summary(lm(as.numeric(originalASVs) ~ meanASVs, data = d2plot_basic)))
sink()

coefs = coef(lm(as.numeric(originalASVs) ~ meanASVs, data = d2plot_basic))


p = ggplot(d2plot_basic, aes(x = originalASVs, y = meanASVs))  + theme_light()
p = p + geom_boxplot() + scale_y_continuous(breaks = c(1:10))
p = p + labs(title = 'Average number of ASVs', x = 'Requested average number of ASVs', y = 'Average number of ASVs per species')
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + geom_abline(intercept = coefs[1], slope = coefs[2], color = "red")
p

ggsave(file = paste0(projectName , '.', 'boxplot_ASVsaccuracy.svg'), plot = p, width = 10, height = 8)
