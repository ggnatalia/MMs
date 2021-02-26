
# Script to check Shannon Index relation to lognormal according to S, H and number of read


library(ggplot2)
library(RColorBrewer)               #for brewer.pal()
library(vegan)
library(reshape2)
library(dplyr)

Hvalues = seq(1, 10, 2) #x # Shannon 0 makes no sense => 1 species
Svalues = seq(100, 5000, 500)  #y                 
reads = c(1000, 10000, 100000, 1000000)
reads = 10000

summary_fun = function(H, S, reads, meanlog = 0, replicas = 5) {
    if ((log(S) - H) < 0 ){
        print(S)
        print(H)
        return(c(NA, NA, NA, NA))
    }
    else{
        sdlog = sqrt(2*(log(S) - H))
        LN = as.data.frame(replicate(replicas, rlnorm(S, meanlog = meanlog, sdlog = sdlog))) # S = 1000, meanlog = 0, sdlog = 1, N =   100 replicates
        colnames(LN) = paste0('R_', c(1:ncol(LN)))
        rownames(LN) = paste0('ASV_', c(1:nrow(LN)))
        
        # Normalize counts to 1
        LN_norm = apply(LN, 2, function(x) x/sum(x))
        df = as.data.frame(matrix(0, nrow = nrow(LN_norm), ncol = ncol(LN_norm)))
        colnames(df) = colnames(LN_norm)
        rownames(df) = rownames(LN_norm)
        #counts = apply(LN_norm, 2, function(x) table(sample(rownames(LN_norm), r, replace = T, x))) # Useless when one replica don't have all the ASVs=> different length, not a df
        for (j in colnames(df)) {
            counts = table(sample(rownames(df), reads, replace = T, LN_norm[, j]))
            #print(j)
            for (value in names(counts)){
                df[value, j] = counts[value]
            }
        }
        mean_H = mean(apply(df, 2, function(x) diversity(x, index = 'shannon')))
        sd_H = mean(apply(df, 2, function(x) diversity(x, index = 'shannon')))
        mean_S = apply(df, 2, function(x) sum(x > 0))
        sd_S = apply(df, 2, function(x) sum(x > 0))
        return(c(mean_H, sd_H, mean_S, sd_S))
    }
}


data2plot = list()

for (r in reads){
    gg = expand.grid(H=Hvalues, S=Svalues)
    gg$H_real_mean = rep(NA,nrow(gg))
    gg$H_real_sd = rep(NA,nrow(gg))
    gg$S_real_mean = rep(NA,nrow(gg))
    gg$H_real_sd = rep(NA,nrow(gg))
    
    for (i in 1:nrow(gg)){
        results = summary_fun(gg$H[i], gg$S[i], r)
        gg$H_real_mean[i] = results[1]
        gg$H_real_sd[i] = results[2]
        gg$S_real_mean[i] = results[3]
        gg$S_real_sd[i] = results[4]
    }
    
    gg$Herror = (gg$H-gg$H_real_mean)/gg$H
    gg$Serror = (gg$S-gg$S_real_mean)/gg$S
    gg$reads = r
    data2plot[[as.character(r)]] = gg
    # SHANNON
    #brks <- cut(gg$Hexp_vs_Hreal, breaks= seq(0, 2, len=21))
    #brks <- gsub(","," - ",brks,fixed=TRUE)
    #gg$brks <- gsub("\\(|\\]","",brks)  # reformat guide labels
    #p = ggplot(gg,aes(H,S)) + geom_tile(aes(fill=brks))
    #p = p + scale_fill_manual("H/Hexp", values = (heat.colors(10)))
    #p = p + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    #p

    
    d2plot_H = gg %>% select(H, S, Herror) %>% acast(H~S, value.var = 'Herror')
    #svg(paste0('Herror_', r,'.svg'), width = 10, height = 10, pointsize = 12)
    filled.contour(d2plot_H,
                   color.palette = colorRampPalette(c('firebrick1','tan1', 'lemonchiffon','tan1', 'firebrick1')),
                   plot.title = title(main = 'Herror',  ylab = 'S', xlab = 'H'),
                   plot.axes = {
                       axis(1, at = seq(from=0, to = 1, length.out = nrow(d2plot_H)), labels = rownames(d2plot_H), las=2) #axis = 1 below
                       axis(2, at = seq(from=0, to = 1, length.out = ncol(d2plot_H)), labels = colnames(d2plot_H), cex.axis=0.9, las=2) #axis = 2 left
                       #axis(2, at =seq(from=0,to=1, length.out = 5), labels =c(4,3,2,1,0), cex.axis=0.9, las=2)
                       
                   }
    )
    #dev.off()
    
    d2plot_S = gg %>% select(H, S, Serror) %>% acast(H~S, value.var = 'Serror', fill=0)
    #svg(paste0('Serror_', r,'.svg'), width = 10, height = 10, pointsize = 12)
    filled.contour(d2plot_S, 
                   color.palette = colorRampPalette(c('firebrick1','tan1', 'lemonchiffon','tan1', 'firebrick1')),
                   plot.title = title(main = 'Serror',  ylab = 'S', xlab = 'H'),
                   plot.axes = {
                       axis(1, at = seq(from=0, to = 1, length.out = nrow(d2plot_S)), labels = rownames(d2plot_S), las=2) #axis = 1 below
                       axis(2, at = seq(from=0, to = 1, length.out = ncol(d2plot_S)), labels = colnames(d2plot_S), cex.axis=0.9, las=2) #axis = 2 left
                   #axis(2, at =seq(from=0,to=1, length.out = 5), labels =c(4,3,2,1,0), cex.axis=0.9, las=2)
                       
                   }
    )
    #dev.off()
}

















## Make vector of colors for values smaller than 0 (20 colors)
rc1 <- colorRampPalette(colors = c('firebrick1','tan1', 'lemonchiffon'), method = 'linear')

## Make vector of colors for values larger than 0 (180 colors)
rc2 <- colorRampPalette(colors = c('lemonchiffon','tan1', 'firebrick1'), method = 'linear')

## Combine the two color palettes
rampcols <- c(rc1, rc2)

mypal <- colorNumeric(palette = rampcols, domain = gg$Herror)

## If you want to preview the color range, run the following code
filled.contour(d2plot_H,
               color.palette = mypal,
               plot.title = title(main = 'Herror',  ylab = 'S', xlab = 'H'),
               plot.axes = {
                   axis(1, at = seq(from=0, to = 1, length.out = nrow(d2plot_H)), labels = rownames(d2plot_H), las=2) #axis = 1 below
                   axis(2, at = seq(from=0, to = 1, length.out = ncol(d2plot_H)), labels = colnames(d2plot_H), cex.axis=0.9, las=2) #axis = 2 left
                   #axis(2, at =seq(from=0,to=1, length.out = 5), labels =c(4,3,2,1,0), cex.axis=0.9, las=2)
                   
               }
)

ggmelt = melt(gg, id.var = c('H','S'), measure.vars = c('Herror'))
p = ggplot(ggmelt, aes(H, S)) + theme_light()
p = p + geom_raster(aes(fill = value), interpolate = T)
p = p + labs(title = 'Shannon reference values', x = 'H', y = 'S', fill = 'H relative error')
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + scale_fill_gradient2(low = 'firebrick', mid = 'lemonchiffon', high = 'tan1', na.value = 'white')
p
