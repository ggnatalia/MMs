#!/usr/bin/Rscript


## Collect arguments
args = commandArgs(TRUE)
projectName = args[1]
projectPath = args[2]
mockName = args[3]
mockPath = args[4]
enviro = args[5]
rank = as.numeric(args[6]) + 1
N = as.numeric(args[7])


finalFile = paste0(projectPath , '/' , projectName , '.pctax.tsv')
originFile = '/home/natalia/Projects/natalia/opt/MMs/DB/SpeciesperEnviro.tsv'
print(finalFile)
print(originFile)


final = t(read.table(finalFile , sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F))
origin = read.table(originFile, header = T, row.names=1, stringsAsFactors = F, check.names = F)
env = origin[, enviro, drop = F]
env[env==0] = 0.005 # we add this to include a little probability of absence taxons




df = as.data.frame(matrix(0, nrow = nrow(env), ncol = 2))
rownames(df) = rownames(env)
colnames(df) = c('origin','final')


df$origin = env[,1]

for (tx in rownames(df)){
    if (tx %in% rownames(final)){
        df[tx, 'final'] = final[tx,'counts']
    } else { df[tx, 'final'] = 0 }
}

dfnorm = as.data.frame(apply(df, 2, function(x) {x/sum(x)}))
#dfnorm = round(dfnorm, 3)

svg(filename=paste0(projectPath , '/' , projectName , '.', 'scatter_lm_taxonomy_subset_comparisons.svg'),  width=10, height=8, pointsize = 12)
par(mfrow=c(2,2))
plot(dfnorm[, 'origin'] ~ dfnorm[,'final'], xlab = 'Origin', ylab = 'Final', pch = 16)
dev.off()


sink(paste0(projectPath , '/' , projectName , '.', 'scatter_lm_taxonomy_subset_comparisons.txt'))
print(summary(lm(dfnorm[, 'origin'] ~ dfnorm[,'final'])))
sink()










