
library(vegan)
library(ggplot2)

Species = c(100, 500, 1000, 5000, 10000)
Shannon = c(1, 2, 3, 4, 5, 6, 7, 8, 9 )
Reads = c(1000, 10000, 100000, 1000000)
meanlogs = c(0) 
replicas = 100
#S = 500
#H = 3




    data = list()
    data_to_plot_S = list()

    for (S in Species){
        l_H = list()
        data_to_plot_H = list()
        for (H in Shannon){
            if ((log(S) - H) <0 ){
                print(S)
                print(H)
                next
            }
            sdlog = sqrt(2*(log(S) - H))
            
            l_reads = list()
            data_to_plot_reads = list()
            for (r in Reads){
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
                    counts = table(sample(rownames(df), r, replace = T, LN_norm[, j]))
                    #print(j)
                    for (value in names(counts)){
                        df[value, j] = counts[value]
                    }
                }
                
                
                d2plot = as.data.frame(apply(df, 2, function(x) diversity(x, index = 'shannon')), drop = T)
                colnames(d2plot) = 'H'
                d2plot$S = apply(df, 2, function(x) sum(x > 0))
                d2plot$expected_H = H
                d2plot$expected_S = S
                d2plot$reads = as.character(r)
                df$Reads = as.character(r)
                #data_to_plot[[as.character(r)]] = d2plot
                data_to_plot_reads[[as.character(r)]] = d2plot
                l_reads[[as.character(r)]] = df
            }
            df_reads = do.call(rbind, l_reads)
            df_reads$Shannon = as.character(H)
            l_H[[as.character(H)]] = df_reads
            
            d2plot_reads = do.call(rbind, data_to_plot_reads)
            d2plot_reads$Shannon = as.character(H)
            data_to_plot_H[[as.character(H)]] = d2plot_reads
        }
    df_H = do.call(rbind, l_H)
    df_H$S = S
    data[[as.character(S)]] = df_H
    
    d2plot_H = do.call(rbind, data_to_plot_H)
    d2plot_H$Species = as.character(S)
    data_to_plot_S[[as.character(S)]] = d2plot_H
    
    }




    df2plot = do.call(rbind, data_to_plot_S)


    df2plot$diffH = (df2plot$expected_H - df2plot$H)/df2plot$expected_H
    df2plot$reads2 = as.numeric(df2plot$reads)
    df2plot$reads = factor(df2plot$reads, levels = sort(as.numeric(unique(df2plot$reads))))
    
    
    df2plot$Species = factor(df2plot$Species, levels = sort(as.numeric(unique(df2plot$Species))))
    df2plot$Shannon = factor(df2plot$Shannon, levels = sort(as.numeric(unique(df2plot$Shannon))))
    
    
    p = ggplot(df2plot, aes(x = (reads2), y = diffH, group = reads2)) + theme_light()
    p = p + geom_boxplot(color = 'darkolivegreen4') + geom_hline(yintercept = 0, linetype="dashed", color = "darkgreen", size=0.5)
    p = p + facet_grid(Species~Shannon)
    #p = p + geom_point(mapping = aes(x = reads, y = S /1000 - 5), colour = 'blue') 
    #p = p + scale_y_continuous(name = "expected H - H", limits = c(-5, 5), sec.axis = sec_axis(~ (. + 5 )*1000, name = "S"))
    
    p = p + geom_boxplot(mapping = aes(x = (reads2), y = (log10(S))), colour = 'dodgerblue')  
    p = p + geom_point(mapping = aes(x = (reads2), y = (log10(expected_S)), group = reads2), colour = 'dodgerblue4', stat = "identity",   position = "identity") 
    
    p = p + scale_y_continuous(name = "H relative error", limits = c(-5, 5), sec.axis = sec_axis(~ (.+0), name = "log10(S)"))
    p = p + theme(plot.title = element_text(hjust = 0.5), axis.title.y.left = element_text(color = "darkolivegreen4"), axis.title.y.right = element_text(color = "dodgerblue4"))
    p = p + ggtitle(paste0('meanlog = ', meanlog))
    #p
    ggsave(file=paste0('meanlog_', meanlog, '.svg'), plot= p, width=16, height=10)



