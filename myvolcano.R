myvolcano=function(df,gene.plot=10,logfc.cutoff=0,p.cutoff=0.01){
  df$color = '#BEBEBE'
  df$group = 'Other'
  df$label = rownames(df)
  df$gene = rownames(df)
  up = df %>% filter(., avg_logFC > logfc.cutoff & (p_val < p.cutoff))
  down = df %>% filter(.,avg_logFC < -logfc.cutoff & (p_val < p.cutoff))
  df[up$gene,'color'] = '#FF0000'
  df[up$gene,'group'] = 'up'
  df[down$gene,'color'] = '#0000FF'
  df[down$gene,'group'] = 'down'
  df$p_val = -log10(df$p_val)
  df = df[order(df$group,-df$p_val),]
  top = df %>% group_by(group) %>% top_n(n = gene.plot, wt = p_val)
  top = top %>% filter(.,group !='Other')
  df[!df$gene %in% (top$gene),'label'] = ''
  
  png(file = paste('volcano_',gene.plot,'.png'),width = 650,height = 500)
  p <- ggplot(df, aes(avg_logFC, p_val,color = group)) + geom_point()+
    labs(x='logFC',y='-Log10 P')
  #labs(x='logFC',y='p-value')
  p <- p + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.line = element_line(colour = 'black',size = 0.8),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = 'grey',size=0.6,linetype = 2),
          axis.ticks = element_line(colour = 'black',size = 0.8),
          legend.title = element_blank(),legend.text = element_text(size = 14),
          legend.position = 'top',legend.key.size = unit(1.5, 'lines'))
  p = p + scale_color_manual(values = setNames(df$color, df$group))
  
  #p = p + geom_label_repel(aes(label = label,colour = 'red'),box.padding = 0.35,point.padding = 0.5,segment.color = 'grey50')
  p = p + geom_text_repel(aes(label = df$label),color= 'black',size = 5, box.padding = 0.35,point.padding = 0.5,segment.color = 'grey70')
  #p = p + geom_hline(yintercept=0, size = 0.8, color = "grey")
  #p = p + geom_vline(xintercept=0, size = 0.8, color = "grey")
  #p = p + geom_abline(intercept = 0, slope = 1, color="grey", 
  #                    size=0.8)
  #p = p + guides(colour = guide_legend(override.aes = list(size = 5),nrow=3))
  p = p + theme(legend.key=element_rect(fill=NA))
  print(p)
  dev.off()
  return(top$gene)
}
