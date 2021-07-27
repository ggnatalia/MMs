# Analyze entropy per position


data = read.table('../DB/Johnson2019Evaluation16S.tsv', sep = '\t', header=T, dec = ',')
regions = read.table('../DB/variableRegions_Ecoli_Brosius1978.tsv', sep = '\t', row.names = 1)
colnames(regions) = c('start', 'end')
#plot(x = data$Base_Position, y = data$Entropy, ylim=c(0,1), pch = 19, xlab = 'Position', ylab = 'Entropy')
#smoothScatter(x = data$Base_Position, y = data$Entropy, ylim=c(0,1), pch = 19, xlab = 'Position', ylab = 'Entropy')    

barplot(names = data$Base_Position, height = data$Entropy, ylim=c(0,0.5), pch = 19, xlab = 'Position', ylab = 'Entropy')
abline(v = regions$start, col = 'red')
abline(v = regions$end, col = 'blue')
