
# Script to check Shannon Index relation to lognormal according to S, H and number of read


library(ggplot2)
library(RColorBrewer)               #for brewer.pal()
library(vegan)
library(reshape2)
library(dplyr)

Hvalues = seq(1, 10, 0.1) #x # Shannon 0 makes no sense => 1 species
Svalues = seq(100, 5000, 100)  #y                 
reads = c(1000, 10000, 100000, 1000000)

summary_fun = function(H, S, reads, meanlog = 0, replicas = 100) {
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
        sd_H = sd(apply(df, 2, function(x) diversity(x, index = 'shannon')))
        mean_S = apply(df, 2, function(x) sum(x > 0))
        #sd_S = apply(df, 2, function(x) sum(x > 0)) # revisar
        sd_S = NA
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

    
    #d2plot_H = gg %>% select(H, S, Herror) %>% acast(H~S, value.var = 'Herror')
    #svg(paste0('Herror_', r,'.svg'), width = 10, height = 10, pointsize = 12)
    
    #Hggmelt = melt(gg, id.var = c('H','S'), measure.vars = c('Herror'))
    #p = ggplot(Hggmelt, aes(H, S)) + theme_light()
    #p = p + geom_raster(aes(fill = value), interpolate = T)
    #p = p + labs(title = 'Shannon reference values', x = 'H', y = 'S', fill = 'H relative error')
    #p = p + theme(plot.title = element_text(hjust = 0.5))
    #p = p + scale_fill_gradient2(low = 'firebrick', mid = 'yellowgreen', high = 'firebrick', na.value = 'white', breaks = seq(0,1,0.1))
    #p

    #Sggmelt = melt(gg, id.var = c('H','S'), measure.vars = c('Serror'))
    #p = ggplot(Sggmelt, aes(H, S)) + theme_light()
    #p = p + geom_raster(aes(fill = value), interpolate = T)
    #p = p + labs(title = 'Shannon reference values', x = 'H', y = 'S', fill = 'S relative error')
    #p = p + theme(plot.title = element_text(hjust = 0.5))
    #p = p + scale_fill_gradient2(low = 'firebrick', mid = 'yellowgreen', high = 'firebrick', na.value = 'white')
    #p
    
    #ggmelt = melt(gg, id.var = c('H','S'), measure.vars = c('Herror', 'Serror'))
    #p = ggplot(ggmelt, aes(H, S)) + theme_light() + facet_grid(.~variable)
    #p = p + geom_raster(aes(fill = value), interpolate = T)
    #p = p + labs(title = 'Shannon reference values', x = 'H', y = 'S', fill = 'S relative error')
    #p = p + theme(plot.title = element_text(hjust = 0.5))
    #p = p + scale_fill_gradient2(low = 'firebrick', mid = 'yellowgreen', high = 'firebrick', na.value = 'white')
    #p
}


dd = do.call(rbind.data.frame, data2plot)
dd$Herror = (dd$H-dd$H_real_mean)/dd$H
dd$Serror = (dd$S-dd$S_real_mean)/dd$S
ggmelt = melt(dd, id.var = c('H','S', 'reads'), measure.vars = c('Herror', 'Serror'))

p = ggplot(ggmelt, aes(H, S)) + theme_light() + facet_grid(reads~variable)
p = p + geom_raster(aes(fill = value), interpolate = T)
p = p + labs(title = 'Shannon reference values', x = 'H', y = 'S', fill = 'Relative error')
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + scale_fill_gradient2(low = 'firebrick', mid = 'yellowgreen', high = 'firebrick', na.value = 'white')
p

ggsave(file = paste0(projectName , '.', 'relativeError.svg'), plot = p, width = 10, height = 8)



# La idea de usar errores relativos es que precisamente iguala las diferencias cuando se miran a una escala pequeña, por ejemplo H, 1(pedido) y 3 (obs) o con números más altos (100 S o 50S)

# Calculate number of reads

rr = dd



gg2 = expand.grid(H = Hvalues, S = Svalues)
gg2$ValidH = rep(NA, nrow(gg2))
gg2$ValidS = rep(NA, nrow(gg2))
rownames(gg2) = paste0(gg2$H, '_', gg2$S)


for (r in sort(reads)){
    print(r)
    srr = rr[rr$reads == r, ]
    # Fix a relative error to be considered acceptable: Ex: 0.10
    for (i in c(1:nrow(srr))) {
        print(r)
        #print(srr[i,])
        #print(gg2[paste0(srr$H[i], '_', srr$S[i]),])
        if (is.na(srr$Herror[i])){
            print(1)
            next()
        } else if(gg2[paste0(srr$H[i], '_', srr$S[i]), 'ValidH'] != 0){
            print(2)
            next()
        } else {
            print(3)
            if (abs(srr$Herror[i]) < 0.1){
                print('fff')
                gg2[paste0(srr$H[i], '_', srr$S[i]), 'ValidH'] = r # save minimun number of reads to satisfy this criteria
            }else if (abs(srr$Herror[i]) > 0.1){
                print('gg')
                gg2[paste0(srr$H[i], '_', srr$S[i]), 'ValidH'] = 0
            }else{next()}
        }
        print(gg2[paste0(srr$H[i], '_', srr$S[i]), 'ValidH'])
    }
}
    

    ggplot(gg2, aes(x = H, y = S, fill = ValidH)) + geom_tile()






