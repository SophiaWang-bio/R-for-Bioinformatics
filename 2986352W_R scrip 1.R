##run the scrip function.R first##
#library####

library(ggplot2)     # For creating elegant data visualizations
library(reshape2)    # For reshaping data
library(sva)         # For PCA plot
library(ggrepel)     # For adding non-overlapping text labels to ggplot2 plots
library(amap)        # For heatmap
library(eulerr)      # For Venn plot
library(survival)    # For survival analysis
library(survminer)   # For survival analysis
library(dplyr)       # For survival analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)

#load and organize data####
setwd("C:/Users/user/Desktop/UoG/semester 2/clinic data/1Final assignment omics/data") 
DE_GB_1_vs_HC = read.table("DE_GB_1_vs_HC.csv", header = TRUE, row.names = 1, sep = "\t") 
DE_GB_2_vs_GB_1 = read.table("DE_GB_2_vs_GB_1.csv", header = TRUE, row.names = 1, sep = "\t") 
DE_GB_2_vs_HC = read.table("DE_GB_2_vs_HC.csv", header = TRUE, row.names = 1, sep = "\t") 
ss = read.table("ss.csv", header = TRUE, sep = "\t") 
ss_patient = read.table("ss_per_patient.csv", header = TRUE, sep = "\t") 
row.names(ss) = ss$sample
em =  read.table("em.csv", header = TRUE, sep = "\t")
row.names(em) = em$ID
em = em[,-1]

#reorder the em table into [HC, GB_1, GB_2]
ss_hc = subset(ss, sample_group == "HC" ) 
ss_GB1 = subset(ss, sample_group == "GB_1" )
ss_GB2 = subset(ss, sample_group == "GB_2" )
ss_reorder = rbind(ss_hc,ss_GB1,ss_GB2)
order_names = row.names(ss_reorder)
em_reorder = em[,order_names]

#rug####
ggprug = rug_function(ss_reorder$sample_group)
save_plot(ggprug,"rug.png", 500,500 )

#survival plot####
surv_object = Surv(time = ss_patient$days_survived, event = ss_patient$censored)
fit1 = survfit(surv_object ~ sample_group, data = ss_patient)
ggp = ggsurvplot(fit1, data = ss_patient, pval = TRUE)
save_plot(ggp,"survival plot.png", 500,500 )

#A vs B Analysis (GB1 vs HC)####

#use DE function to get PCA plot, volcano plot, ma plot, heatmap, boxplot
results_GB1_HC = do_DE(em, DE_GB_1_vs_HC, ss)
#get plots from results_GB1_HC
pca_GB1_HC = results_GB1_HC$plots$ggp.pca #PCA
save_plot(pca_GB1_HC,"pca_GB1_HC_GB2.png", 300,300 ) #save plot
vol_GB1_HC = results_GB1_HC$plots$ggp.volcano #volcano
save_plot(vol_GB1_HC,"vol_GB1_HC.png", 500,500 )  #save plot
ma_GB1_HC = results_GB1_HC$plots$ggp.ma  #MA plot
save_plot(ma_GB1_HC,"ma_GB1_HC.png", 500,500 )  #save plot
results_GB1_HC$plots$ggp.heatmap    # heatmap
boxplottop10_GB1_HC = results_GB1_HC$plots$ggp.top10.boxplot #boxpolot of top 10 significant genes
save_plot(boxplottop10_GB1_HC,"boxplottop10_GB1_HC111.png", 500,500 ) #save plot

#top 5 signifcant genes that upregulate
sig_genes_up = results_GB1_HC$data$master_up_top5$SYMBOL #get the top 5 gene names
#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_up = make_heatmap(em_reorder, sig_genes_up) #heatmap
ggp_boxplot_up = make_boxplot(em, sig_genes_up, ss$sample_group) #boxplot
save_plot(ggp_boxplot_up,"boxplot_up5_GB1_HC.png", 500,500 ) #save plot

#top 5 significant genes that downregulate
sig_genes_down = results_GB1_HC$data$master_down_top5$SYMBOL #get the top 5 gene names
#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_down = make_heatmap(em_reorder, sig_genes_down) #heatmap
ggp_boxplot_down = make_boxplot(em, sig_genes_down, ss$sample_group) #boxplot
save_plot(ggp_boxplot_down,"boxplot_down5_GB1_HC.png", 500,500 ) #save plot

#A vs B Analysis (GB2 vs HC)####

#use DE function to get PCA plot, volcano plot, ma plot, heatmap, boxplot
results_GB2_HC = do_DE(em, DE_GB_2_vs_HC, ss)

#get plots from results_GB2_HC
vol_GB2_HC = results_GB2_HC$plots$ggp.volcano # volcano plot
save_plot(vol_GB2_HC,"vol_GB2_HC.png", 500,500 ) #save plot
ma_GB2_HC = results_GB2_HC$plots$ggp.ma #ma plot
save_plot(ma_GB2_HC,"ma_GB2_HC.png", 500,500 ) #save plot
results_GB2_HC$plots$ggp.heatmap  # heatmap
boxplottop10_GB2_HC = results_GB2_HC$plots$ggp.top10.boxplot  #boxpolot of top 10 significant genes
save_plot(boxplottop10_GB2_HC,"boxplottop10_GB2_HC.png", 500,500 ) #save plot

#top 5 signifcant genes that upregulate
sig_genes_up = results_GB2_HC$data$master_up_top5$SYMBOL #get the top 5 gene names

#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_up = make_heatmap(em_reorder, sig_genes_up) #heatmap
ggp_boxplot_up = make_boxplot(em, sig_genes_up, ss$sample_group) #boxplot
save_plot(ggp_boxplot_up,"boxplot_up5_GB2_HC.png", 500,500 ) #save plot

#top 5 signifcant genes that downregulate
sig_genes_down = results_GB2_HC$data$master_down_top5$SYMBOL #get the top 5 gene names

#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_down = make_heatmap(em_reorder, sig_genes_down) #heatmap
ggp_boxplot_down = make_boxplot(em, sig_genes_down, ss$sample_group) #boxplot
save_plot(ggp_boxplot_down,"boxplot_down5_GB2_HC.png", 500,500 ) #save plot

#do de (GB1 vs GB2)####

#use DE function to get PCA plot, volcano plot, ma plot, heatmap, boxplot
results_GB1_GB2 = do_DE(em, DE_GB_2_vs_GB_1, ss)

#get plots from results_GB1_GB2
vol_GB2_GB_1 = results_GB1_GB2$plots$ggp.volcano #volcano plot
save_plot(vol_GB2_GB_1,"vol_GB2_GB1.png", 500,500 ) #save plot
ma_GB2_GB_1 = results_GB1_GB2$plots$ggp.ma #ma plot
save_plot(ma_GB2_GB_1,"ma_GB2_GB_1.png", 500,500 ) #save plot
results_GB1_GB2$plots$ggp.heatmap #heatmap
boxplottop10_GB2_GB_1 = results_GB1_GB2$plots$ggp.top10.boxplot   #boxpolot of top 10 significant genes
save_plot(boxplottop10_GB2_GB_1,"boxplottop10_GB2_GB_1.png", 500,500 ) #save plot

#top 5 signifcant genes that upregulate
sig_genes_up = results_GB1_GB2$data$master_up_top5$SYMBOL #get the top 5 gene names

#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_up = make_heatmap(em_reorder, sig_genes_up) #heatmap
ggp_boxplot_up = make_boxplot(em, sig_genes_up, ss$sample_group) #boxplot
save_plot(ggp_boxplot_up,"boxplot_up5_GB2_GB_1.png", 500,500 ) #save plot

#top 5 signifcant genes that downregulate
sig_genes_down = results_GB1_GB2$data$master_down_top5$SYMBOL #get the top5 gene names

#use make_heatmap and make_boxplot function to make plots
ggp_heatmap_down = make_heatmap(em_reorder, sig_genes_down) #heatmap 
ggp_boxplot_down = make_boxplot(em, sig_genes_down, ss$sample_group) #boxplot
save_plot(ggp_boxplot_down,"boxplot_down5_GB2_GB_1.png", 500,500 ) #save plot

#pathway GB1 VS H####

#get significant genes name
sig_genes_GB1_h = results_GB1_HC$data$sig_genes

#use do pathway function
pathway_GB1_h = do_pathway(org.Hs.eg.db,sig_genes_GB1_h,em, ss$sample_group)

ggp = pathway_GB1_h$plots$ggp.barplot #get the barplot from pathway_GB1_h
save_plot(ggp,"pathway_GB1_h_barplot.png", 500,500 ) #save plot
ggp1 = pathway_GB1_h$plots$ontology1$ggp.boxplot #get the boxplot of ontology1 from pathway_GB1_h
save_plot(ggp1,"GB_H_ontology1.png", 1000,1000 )  #save plot

#pathway GB1 VS GB2 by ####

#get significant genes name
sig_genes_GB1_2 = results_GB1_GB2$data$sig_genes

#use do pathway function
pathway_results = do_pathway(org.Hs.eg.db, sig_genes_GB1_2, em, ss$sample_group)

ggp = pathway_results$plots$ggp.barplot #get the barplot from pathway_GB1_2
save_plot(ggp,"pathway_GB1_GB2_barplot.png", 500,500 ) #save plot
ggp1 = pathway_results$plots$ontology1$ggp.boxplot #get the boxplot of ontology1 from pathway_GB1_h
save_plot(ggp1,"GB_GB2_ontology1.png", 1000,1000 )  #save plot

#pathway h VS GB2####

#get significant genes name
sig_genes_GB2_HC = results_GB2_HC$data$sig_genes 

#use do pathway function
pathway_h_gb2 = do_pathway(org.Hs.eg.db,sig_genes_GB2_HC,em, ss$sample_group)
ggp = pathway_h_gb2$plots$ggp.barplot #get the barplot from pathway_GB2_hc
save_plot(ggp,"pathway_GB2_h_barplot.png", 500,500 ) #save plot
ggp1 = pathway_h_gb2$plots$ontology1$ggp.boxplot   #get the boxplot from pathway_GB2_hc
save_plot(ggp1,"GB2_H_ontology1.png", 1000,1000 )  #top 10 ontologies

#[A&B]n####

#make new column to label significant or not
DE_GB_1_vs_HC$sig = as.factor(DE_GB_1_vs_HC$p.adj < 0.05 & abs(DE_GB_1_vs_HC$log2fold) > 1)
DE_GB_2_vs_GB_1$sig = as.factor(DE_GB_2_vs_GB_1$p.adj < 0.05 & abs(DE_GB_2_vs_GB_1$log2fold) > 1)
DE_GB_2_vs_HC$sig = as.factor(DE_GB_2_vs_HC$p.adj < 0.05 & abs(DE_GB_2_vs_HC$log2fold) > 1)

#gb1 gb2 vs gb1 hc (1) ####

#Venn plot of gb1 gb2 vs gb1 hc

#make master 1 and organized
master1 = merge(DE_GB_2_vs_GB_1, DE_GB_1_vs_HC, by = 0, suffixes = c(".1_v_2",".1_v_h"))
colnames(master1)[1] = "SYMBOL"
row.names(master1) = master1[,"SYMBOL"]

# get names of significant genes
sig_1_v_2 = row.names(subset(master1, sig.1_v_2 == TRUE ))
sig_1_v_h = row.names(subset(master1, sig.1_v_h == TRUE ))

# Create Venn plot between gb1 gb2 vs gb1 hc
venn_data = list("GB1_v_GB2" = sig_1_v_2, "GB1_v_Health" = sig_1_v_h)
venn_1 = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE)
save_plot(venn_1,"venn_1.png", 500,500 ) #save plot

# get the gene names form different district of Venn plot
sig_1_v_2_only = row.names(subset(master1, sig.1_v_2 == TRUE & sig_1_v_h == FALSE ))
sig_1_v_h_only = row.names(subset(master1, sig.1_v_2 == FALSE & sig.1_v_h == TRUE ))
sig_both = row.names(subset(master1, sig.1_v_2 == TRUE & sig.1_v_h == TRUE )) #both significant genes

#make heatmap of both significant genes
ggp_heatmap_12vs1h = make_heatmap(em_reorder, sig_both)
save_plot(ggp_heatmap_12vs1h,"heatmap_12vs1h.png", 500,500 ) #save plot

#p value calculation
group1 = nrow(subset(master1, sig.1_v_2 == TRUE )) 
group2 = nrow(subset(master1, sig.1_v_h == TRUE ))
overlap = nrow(subset(master1, sig.1_v_2 == TRUE & sig.1_v_h == TRUE ))
total = nrow(master1)
p = phyper(overlap-1, group2, total-group2, group1,lower.tail =  FALSE) 

#fold - fold plot
ggp_fold = ggplot(master1, aes(x = log2fold.1_v_2, y = log2fold.1_v_h)) + geom_point()
save_plot(ggp_fold,"gb1 gb2 vs gb1 hc fold.png", 500,500 ) #save plot

#correlation
cor.test(master1$log2fold.1_v_2, master1$log2fold.1_v_h)

# Differential expression signatures

# get the genes that go up then back down
signature_1 = row.names(subset(master1, (sig.1_v_2 == TRUE & log2fold.1_v_2 > 1) & (sig.1_v_h == TRUE & log2fold.1_v_h < -1) ))
# get the genes that go down then back up
signature_2 = row.names(subset(master1, (sig.1_v_2 == TRUE & log2fold.1_v_2 < -1) & (sig.1_v_h == TRUE & log2fold.1_v_h > 1) ))

# use plot_signature function
gb1_2_vs_gb1_h_sig1 = plot_signature(org.Hs.eg.db, signature_1, em_reorder, ss$sample_group)

#get the plots from gb1_2_vs_gb1_h_sig1 
sig1_heatmap = gb1_2_vs_gb1_h_sig1$plots$ggp.heatmap #heatmap
save_plot(sig1_heatmap,"sig1_heatmap.png", 500,500 ) #save plot
sig1_metagene = gb1_2_vs_gb1_h_sig1$plots$ggp.metagene #metagene
save_plot(sig1_metagene,"sig1_metagene.png", 500,500 ) #save plot
sig1_pathway = gb1_2_vs_gb1_h_sig1$pathway$pathway_results$plots$ggp.barplot #barplot of pathway
save_plot(sig1_pathway,"sig1_pathway.png", 500,500 ) #save plot
sig1_ontology1 = gb1_2_vs_gb1_h_sig1$pathway$pathway_results$plots$ontology1 #boxplot of ontology1
save_plot(sig1_ontology1,"sig1_ontology1.png", 500,500 ) #save plot

# use plot_signature function
gb1_2_vs_gb1_h_sig2 = plot_signature(org.Hs.eg.db, signature_2, em_reorder, ss$sample_group)

#get the plot from gb1_2_vs_gb1_h_sig2
sig2_heatmap = gb1_2_vs_gb1_h_sig2$plots$ggp.heatmap #heatmap
save_plot(sig2_heatmap,"sig2_heatmap.png", 500,500 ) #save plot
sig2_metagene = gb1_2_vs_gb1_h_sig2$plots$ggp.metagene #metagene
save_plot(sig2_metagene,"sig2_metagene.png", 500,500 ) #save plot
sig2_pathway = gb1_2_vs_gb1_h_sig2$pathway$pathway_results$plots$ggp.barplot #barplot of pathway
save_plot(sig2_pathway,"sig2_pathway.png", 500,500 ) #save plot
sig2_ontology = gb1_2_vs_gb1_h_sig2$pathway$pathway_results$plots$ontology1 #boxplot of ontology1
save_plot(sig2_ontology,"sig2_ontology1.png", 300,300 ) #save plot

#gb1 gb2 vs gb2 hc (2) ####

#Venn plot of gb1 gb2 vs gb2 hc

#make master 2 and organized
master2 = merge(DE_GB_2_vs_GB_1, DE_GB_2_vs_HC, by = 0, suffixes = c(".1_v_2",".2_v_h"))
colnames(master2)[1] = "SYMBOL"
row.names(master2) = master2[,"SYMBOL"]

# get names of significant genes
sig_1_v_2 = row.names(subset(master2, sig.1_v_2 == TRUE ))
sig_2_v_h = row.names(subset(master2, sig.2_v_h == TRUE ))

# Create Venn plot between gb1 gb2 vs gb1 hc
venn_data = list("GB1_v_GB2" = sig_1_v_2, "GB2_v_Health" = sig_2_v_h)
venn_2 = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE)
save_plot(venn_2,"venn_2.png", 500,500 ) #save plot

# get the gene names form different district of Venn plot
sig_1_v_2_only = row.names(subset(master2, sig.1_v_2 == TRUE & sig.2_v_h == FALSE ))
sig_2_v_h_only = row.names(subset(master2, sig.1_v_2 == FALSE & sig.2_v_h == TRUE ))
sig_both = row.names(subset(master2, sig.1_v_2 == TRUE & sig.2_v_h == TRUE )) #both significant genes

#make heatmap of both significant genes
heatmap_12vs2h = make_heatmap(em_reorder.s, sig_both)
save_plot(heatmap_12vs2h,"heatmap_12vs2h.png", 500,500 ) #save plot

#p value calculation
group1 = nrow(subset(master2, sig.1_v_2 == TRUE ))
group2 = nrow(subset(master2, sig.2_v_h == TRUE ))
overlap = nrow(subset(master2, sig.1_v_2 == TRUE & sig.2_v_h == TRUE ))
total = nrow(master2)
p_2 = phyper(overlap-1, group2, total-group2, group1,lower.tail =  FALSE) #2.816184e-17

#fold - fold plot
ggp_fold = ggplot(master2, aes(x = log2fold.1_v_2, y = log2fold.2_v_h)) + geom_point()
save_plot(ggp_fold2,"gb1 gb2 vs gb2 hc fold.png", 500,500 ) #save plot

#correlation
cor.test(master2$log2fold.1_v_2, master2$log2fold.2_v_h)

# Differential expression signatures

# get the genes that no change then go down #ONLY ONE GENE (EFEMP1)
signature_1_2 = row.names(subset(master2, (sig.1_v_2 ==  TRUE & log2fold.2_v_h < -1) & (sig.2_v_h == TRUE & log2fold.2_v_h < -1) ))
# get the genes that no change then go up
signature_2_2 = row.names(subset(master2, (sig.1_v_2 == TRUE & log2fold.2_v_h > 1) & (sig.2_v_h == TRUE & log2fold.2_v_h > 1) ))

# use plot_signature function
gb1_2_vs_gb2_h_sig2 = plot_signature(org.Hs.eg.db, signature_2_2, em_reorder, ss$sample_group)

#get the plots from gb1_2_vs_gb2_h_sig2 
sig2_heatmap = gb1_2_vs_gb2_h_sig2$plots$ggp.heatmap #heatmap
save_plot(sig2_heatmap,"sig2_2h_heatmap.png", 500,500 ) #save plot
sig2_metagene = gb1_2_vs_gb2_h_sig2$plots$ggp.metagene  #metagene
save_plot(sig2_metagene,"sig2_2h_metagene.png", 500,500 ) #save plot
sig2_pathway = gb1_2_vs_gb2_h_sig2$pathway$pathway_results$plots$ggp.barplot #barplot of pathway
save_plot(sig2_pathway,"sig2_2h_pathway.png", 500,500 ) #save plot
sig2_ontology1 = gb1_2_vs_gb2_h_sig2$pathway$pathway_results$plots$ontology1 #boxplot of ontology1
save_plot(sig2_ontology1,"sig2_2h_ontology1.png", 300,300 ) #save plot


