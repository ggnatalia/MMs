# Check number os ASVs per species
library('reshape2')
library('ggplot2')
library('dplyr')
library('data.table')



data2plot = list()
d2plot_basic = as.data.frame(matrix(0, nrow = 100, ncol = 3))
colnames(d2plot_basic) = c('R','meanASVs', 'sdASVs')
projectName = 'Freshwaters_master'
ASVsmean = 2

for (i in 1:100){ ###### trial
    R = paste0('R_', i)
    finalFile = paste0(projectName , i, '.distances.tsv')
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
    
    df = as.data.frame(matrix(0, nrow = length(levels(as.factor(data2$species1))), ncol = 5))
    colnames(df) = c('Sp', 'max_distance', 'min_distance', 'nASVsperSp', 'ASVs')
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
        df[sp, 'nASVsperSp'] = as.numeric(nASVsperSp)  #not count the original reference 
        df[sp, 'ASVs'] = paste(unique(c(dumb$seq1, dumb$seq2)), collapse = ',')
        }
        else{df[sp,] = c(sp, 0, 0, 1,sp)}
    }
    data2plot[[R]] = df
    d2plot_basic[i, ] = c(R, mean(as.numeric(df$nASVsperSp), na.rm=TRUE), sd(as.numeric(df$nASVsperSp), na.rm=TRUE))
}


d2plot = rbindlist(data2plot)
d2plot = data.frame(lapply(data2plot, "length<-", max(lengths(data2plot))))
svg('Mean of ASVs')
boxplot(as.numeric(d2plot_basic$meanASVs))
dev.off()
#t.test(as.numeric(d2plot_basic$meanASVs), mu = 0, alternative = "one.sided")
#mean_of_ASVs = mean(df$nASVsperSp)

error = abs(ASVsmean - mean(as.numeric(d2plot_basic$meanASVs)))/ASVsmean
