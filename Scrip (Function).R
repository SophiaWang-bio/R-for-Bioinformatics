#PCA function
make_pc1_pc2 = function(colour_groups, e_data)
{
  # scale data
  e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  
  # run PCA
  xx = prcomp(t(e_data_scaled))
  pca_coordinates = data.frame(xx$x)
  
  # get % variation
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
  y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
  
  # plot  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = colour_groups)) +
    geom_point(size = 3) +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) +
    theme_bw()
  
  return(ggp)
}

# save plot
save_plot= function(plot, plot.path, plot.height, plot.width)
{
  png(plot.path, width = plot.width, height = plot.height)
  print(plot)
  dev.off()
}

#get DE and plot
do_DE = function(e_data, de, ss_data)
{
  # setup a list object, to store all our results
  de_object = list("data" = list(), "plots" = list())
  
  
  #create master
  master = merge(e_data, de, by.x = 0, by.y = 0)
  
  # re_dataove nas
  master = na.omit(master)
  colnames(master)[1] = "SYMBOL"
  row.names(master) = master[,"SYMBOL"]

  # sort by p
  master = master[order(master[,"p"]),]
  
  # create new columns
  master$mean = rowMeans(master[,row.names(ss_data)])
  master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1)
  master$mlog10p = -log10(master$p.adj)
  
  # create sig_genes
  master_sig = subset(master, sig==TRUE)
  sig_genes = row.names(master_sig)
  
  # create e_data symbols scaled
  e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  
  # create e_data sig
  e_data_symbols_sig = e_data[sig_genes,]
  e_data_scaled_sig = e_data_scaled[sig_genes,]
  
  
  # parse up and down
  master_up = subset(master,p.adj < 0.05 & log2fold > 1)
  master_down = subset(master,p.adj < 0.05 & log2fold < -1)
  master_up_top5 = master_up[1:5,]
  master_down_top5 = master_down[1:5,]
  
  ## add the tables to the list ####
  de_object$data$e_data = e_data 
  de_object$data$ss_data = ss_data
  de_object$data$de = de
  de_object$data$sig_genes = sig_genes
  de_object$data$master = master
  de_object$data$master_sig = master_sig
  de_object$data$e_data_symbols_sig = e_data_symbols_sig
  de_object$data$e_data_scaled = e_data_scaled
  de_object$data$e_data_scaled_sig = e_data_scaled_sig
  de_object$data$master_up = master_up
  de_object$data$master_up_top5 = master_up_top5
  de_object$data$master_down = master_down
  de_object$data$master_down_top5 = master_down_top5
  
  
  ####-- theme --####
  
  my_theme = theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.background = element_blank(), 
    legend.background = element_rect(fill="transparent", colour=NA),
    legend.key = element_rect(fill="transparent", colour=NA),
    plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
    title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
    axis.text.y = element_text(size = 11, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
    axis.text.x = element_text(size = 11, margin = margin(t = 5),hjust=0.5,vjust=1, family="Arial", face="bold",colour="black"), 
    axis.title.y = element_text(size = 12, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="Arial", face="bold"),
    axis.title.x = element_text(size = 12, margin = margin(t = 10),hjust=0.5,vjust=1, family="Arial", face="bold"),
    legend.text=element_text(size=12, family="Arial", face="bold"),
    legend.title=element_blank(), 
    legend.key.size=unit(1,"line"),
    plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
    strip.text.x = element_text(size = 12, family="Arial", face="bold", vjust=1),
    panel.spacing = unit(1, "lines")
  )
  
  
  
  ####-- Volcano With Legend --####
  
  ggp.volcano = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
    geom_point(aes(colour = "a")) +
    geom_point(data=master_up, aes(colour = "b")) +
    geom_point(data=master_down, aes(colour = "c")) +
    geom_text_repel(data=master_up_top5, aes(label=SYMBOL, colour = "b"), vjust=-0.5, show.legend=FALSE) + 
    geom_text_repel(data=master_down_top5, aes(label=SYMBOL, colour = "c"), vjust=1.5, show.legend=FALSE) + 
    scale_color_manual(values=c("black", "red", "blue"), labels = c("No change", "Up", "Down"), name = "") + 
    geom_vline(xintercept=-1,linetype="dashed") +
    geom_vline(xintercept=1,linetype="dashed") +
    geom_hline(yintercept=-log10(0.05),linetype="dashed") +
    my_theme +
    labs(title = "Volcano Plot", x= "Log2 fold change", y= "-Log10 p")
  
  ##-- MA With Legend --####
  
  # make plot
  ggp.ma = ggplot(master, aes(x=log10(mean), y=log2fold)) + 
    
    # adds the dots
    geom_point(aes(colour = "a")) +
    geom_point(data=master_up, aes(colour = "b")) +
    geom_point(data=master_down, aes(colour = "c")) +
    
    # choose colours and legend names
    scale_color_manual(values=c("black", "red", "blue"), labels = c("No change", "Up", "Down"), name = "") + 
    
    # adds the fancy lines
    geom_hline(yintercept=1,linetype="dashed") +
    geom_hline(yintercept=-1,linetype="dashed") +
    
    # adds the theme and axis titles
    my_theme +
    labs(title = "MA", x= "Log10 mean express_dataion", y= "Log2 fold change")
  
  
  
  ##-- PCA --####
  
  # run PCA
  xx = prcomp(t(e_data_scaled))
  pca_coordinates = data.frame(xx$x)
  
  # get % variation 
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  
  x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
  y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
  
  
  # plot
  ggp.pca = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss_data$sample_group)) + 
    geom_point() + 
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) + 
    my_theme
  
  
  
  ####-- Express_dataion Density --####
  
  ## gets a groups column - so we can colour the plot
  e_data_groups = e_data
  for (index in 1:ncol(e_data))
  {
    group = ss_data[index,"sample_group"]
    e_data_groups[,index] = group
  }
  
  # finally melt the new table
  e_data_groups = melt(e_data_groups, id.vars = NULL)
  
  
  ## PLOT
  
  # melts e_data symbols
  e_data_symbols.m = melt(e_data)
  
  # adds the group column
  e_data_symbols.m$groups = e_data_groups$value
  
  # plot
  ggp.density = ggplot(e_data_symbols.m, aes(x = log10(value+0.01), fill = groups)) + 
    geom_density(alpha = 0.75) + 
    facet_wrap(~variable, ncol=3) + 
    my_theme + 
    theme(strip.background = element_rect(fill="transparent", size=0), legend.position = "none") + 
    labs(x = "Express_dataion (log10)", y = "Density")
  
  
  
  ####-- Facet Boxplot --####
  
  # get the gene list
  candidate_genes = row.names(master[1:10,])
  candidate_genes = c(candidate_genes)
  
  # create the table
  gene_data = e_data_scaled[candidate_genes, ]
  gene_data = data.frame(t(gene_data))
  gene_data$groups = ss_data$sample_group
  
  # melt
  gene_data.m = melt(gene_data, id.vars = "groups")
  
  # makes the plot
  ggp.top10.boxplot = ggplot(gene_data.m, aes(x=groups, y=value, fill = groups)) + 
    geom_boxplot()  + 
    my_theme + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~variable, ncol=5) + 
    labs(x = "", y = "Express_dataion (z-score)")
  
  
  
  ####-- Heatmap --#### 
  
  # makes a matrix
  hm.matrix = as.matrix(e_data_scaled_sig)
  
  # does the y clustering
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  
  # does the y clustering
  x.dist = Dist(t(hm.matrix), method="spearman")
  x.cluster = hclust(x.dist, method="average")
  x.dd = as.dendrogram(x.cluster)
  x.dd.reorder = reorder(x.dd,0,FUN="average")
  x.order = order.dendrogram(x.dd.reorder)
  
  # include the y and x orders here
  hm.matrix_clustered = hm.matrix[y.order,x.order]
  
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  # colour palette
  colours = c("blue","black","yellow")
  palette = colorRampPalette(colours)(100)
  
  # plot
  ggp.heatmap = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = palette) + 
    labs(x="", y="") + 
    theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'), axis.text.y = element_blank(), axis.ticks=element_blank())
  
  
  ## Add plots to results list ####
  de_object$plots$ggp.density = ggp.density
  de_object$plots$ggp.volcano = ggp.volcano
  de_object$plots$ggp.ma = ggp.ma
  de_object$plots$ggp.pca = ggp.pca
  de_object$plots$ggp.top10.boxplot = ggp.top10.boxplot
  de_object$plots$ggp.heatmap = ggp.heatmap
  
  
  ## return the results
  return(de_object)
}

#heat map
make_heatmap = function(e_data, genes)
{
  
  # libraries
  library(amap)
  library(reshape2)
  library(ggplot2)
  
  # parses
  e_data_scaled_candidates = na.omit(data.frame(t(scale(t(e_data[genes,])))))
  hm.matrix = as.matrix(e_data_scaled_candidates)
  
  # does the y clustering
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  # plot
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = colorRampPalette(c("blue","black","yellow"))(100)) + 
    labs(x="", y="") + 
    theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'), axis.text.x = element_blank(), axis.ticks=element_blank())
  
  return(ggp)
}

#boxchart
make_boxplot = function(e_data, candidate_genes, g_data)
{
  
  # libraries
  library(reshape2)
  library(ggplot2)
  
  # create the table
  gene_data = na.omit(data.frame(t(scale(t(e_data[candidate_genes,])))))
  
  gene_data = data.frame(t(gene_data))
  gene_data$groups = g_data
  
  # melt
  gene_data.m = melt(gene_data, id.vars = "groups")
  
  # makes the plot
  ggp = ggplot(gene_data.m, aes(x=variable, y=value, fill = groups)) + 
    geom_boxplot()  + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(ggp)
  
}

#pathway
do_pathway = function(organism_db, genes,e_data, g_data)
{
  # libraries
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(STRINGdb)
  
  # converts from ensembl Symbols to Entrez
  sig_genes_entrez = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism_db)
  
  # gets the enrichment
  pathway_data = enrichGO(gene = sig_genes_entrez$ENTREZID,OrgDb = organism_db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  
  # extract the enrichment table from go enrich 
  gene_sets = pathway_data$geneID
  description = pathway_data$Description
  p.adj = pathway_data$p.adjust
  ora_results = data.frame(cbind(gene_sets, description, p.adj))
  
  # bar and dotplot
  ggp.barplot = barplot(pathway_data, showCategory=10)
  ggp.dotplot = dotplot(pathway_data, showCategory=10)
  #list
  pathway_results = list("tables" = list(), "plots" = list()) 
  # store tables and plots
  pathway_results$tables$pathway_data = pathway_data
  pathway_results$tables$ora_results = ora_results
  pathway_results$plots$ggp.barplot = ggp.barplot 
  pathway_results$plots$ggp.dotplot = ggp.dotplot 
  ## Take the top 10 ontologies, and make the boxplot & heatmap,
  for (row_index in 1:10)
  {
    # get the genes
    enriched_gene_set = as.character(ora_results[row_index,1])
    candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
    
    # make the plots
    ggp.boxplot = make_boxplot(e_data, candidate_genes, g_data)
    
    ## store in a list
    ontology_result = list()
    ontology_result$candidate_genes = candidate_genes
    ontology_result$ggp.boxplot = ggp.boxplot
    
    # put list into results
    pathway_results$plots[[paste0("ontology",row_index)]] = ontology_result
  }
  
  # return the results
  return(pathway_results)
}

#function of metagene
make_metagene=function(e_data,genes,g_data)
{
  signature_1_em = e_data[genes,]
  # get the metagene for each cluster
  signature_1_metagene = data.frame(colMeans(signature_1_em))
  names(signature_1_metagene) = "meta_expression"
  signature_1_metagene$group = g_data
  ggp =ggplot(signature_1_metagene, aes(x=group, y=meta_expression, fill=group)) + geom_boxplot()
  return(ggp)
}

#plot_signature
plot_signature = function(organism_db,genes,e_data, g_data)
{
  ggp.heatmap = make_heatmap(e_data,genes)
  pathway_results = do_pathway(organism_db,genes,e_data, g_data)
  ggp.metagene = make_metagene(e_data,genes,g_data)
  results = list("pathway" = list(),"plots" = list())
  #pathway_results = list("tables" = list(), "plots" = list()) 
  # store tables and plots
  results$pathway$pathway_results = pathway_results
  results$plots$ggp.heatmap = ggp.heatmap
  results$plots$ggp.metagene = ggp.metagene 
  return(results)
}

#rug
rug_function=function(g_data)
{
  rug_colours = c("#F8766D","lightgreen","#619CFF" )
  # rug for discrete variable
  rug_data = as.matrix(as.numeric(factor(g_data)))
  rug_data = melt(rug_data)
  #plot
  ggp = ggplot(rug_data , aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_gradientn(colours = rug_colours)
  return(ggp)
}
  




