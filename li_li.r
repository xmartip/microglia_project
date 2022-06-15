#This paper had no GSE nor supplementary information the last time I checked it,
#but I have some SeuratObjects.
#The 1st one (fig1b.RData) is scRNA-seq on microglia from different regions of 
#developing human brain and early fate decision of microglia progenitors
#The 2nd one (fig4b.RData) is Fetal microglia with regional specification exit 
#from embryonic relative resting state
#I guess it can be seen as if the two datasets are set in a pseudo-chronological
#order, but I have to ask anyway.

##Libraries + WD
require(dplyr)
require(patchwork)
require(Seurat)
require(ggplot2)
setwd("/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/li_li_2022")

##Load objects
load("fig1b.RData")
mg_genes8_22 <- readRDS("heatmap_order_genes8_22_hard.rds")
mg_genes1_7 <- readRDS("heatmap_order_genes1_7_hard.rds")
ps_table <- readRDS("ps_table.rds")
deg_clusters <- readRDS("marker_genes_li_clusters.rds")

##Create metadata column for minimally decent clustering
#The current idents (currrently not a metadata column)
Microglia.Brain.Region$cluster <- Idents(Microglia.Brain.Region)
#The combined state+cluster column
Microglia.Brain.Region$state_cluster <- paste(Microglia.Brain.Region$state, Microglia.Brain.Region$cluster)
#Assign this last column as Idents
Idents(Microglia.Brain.Region) <- Microglia.Brain.Region$state_cluster

##Translate the symbols into IDs and store it as a character (just in case)
id_list1_7 <- as.character(mg_genes1_7)
id_list8_22 <- as.character(mg_genes8_22)

##Store the DotPlot fuction data into an object for further manipulation
preplot1_7 <- DotPlot(Microglia.Brain.Region, features = id_list1_7, cluster.idents = TRUE)
preplot8_22 <- DotPlot(Microglia.Brain.Region, features = id_list8_22, cluster.idents = TRUE)


##Add the PS column
preplot1_7$data$PS <- ps_table$PS[match(preplot1_7$data$features.plot,
                                        ps_table$Name)]
preplot8_22$data$PS <- ps_table$PS[match(preplot8_22$data$features.plot,
                                         ps_table$Name)]

##Find Gene markers (done and hushed)
#deg_clusters <- FindAllMarkers(Microglia.Brain.Region, min.pct = 0.5, logfc_threshold = 0.5, test.use = "wilcox")
#deg_clusters <- subset(deg_clusters, avg_log2FC>0.5 & p_val_adj<=0.05)
#saveRDS(deg_clusters, file = 'marker_genes_li_clusters.rds')

##Add cluster column to preplot$data because deg$cluster is not the same as preplot$id
##this time.
preplot1_7$data$cluster <- Microglia.Brain.Region$cluster[match(preplot1_7$data$id, Microglia.Brain.Region$state_cluster)]
preplot8_22$data$cluster <- Microglia.Brain.Region$cluster[match(preplot8_22$data$id, Microglia.Brain.Region$state_cluster)]

##Add significance
for (i in 1:nrow(deg_clusters)) {
  gene <- deg_clusters$gene[i]
  cluster <- deg_clusters$cluster[i]
  preplot1_7$data$deg[preplot1_7$data$features.plot == gene & preplot1_7$data$cluster == cluster] <- "Yes"
  preplot8_22$data$deg[preplot8_22$data$features.plot == gene & preplot8_22$data$cluster == cluster] <- "Yes"
}
#For the ggplot code to work there mustn't be any NA (I don't know why)
preplot1_7$data$deg[is.na(preplot1_7$data$deg)] <- "No"
preplot8_22$data$deg[is.na(preplot8_22$data$deg)] <- "No"

##Set the order of the genes
order_genes1_7 <- mg_genes1_7
order_genes8_22 <- mg_genes8_22

#Set the <0.05 % expr values to NA.
preplot1_7$data$pct.exp[preplot1_7$data$pct.exp < 5] <- NA
preplot8_22$data$pct.exp[preplot8_22$data$pct.exp < 5] <- NA

##Order clusters in this order: 
order_clusters <- c("Origin C2", "Origin C3", "Cycling C6", "Cycling C9", "Cycling C10",
                    "Intermediate C11", "Intermediate C12", "Intermediate C13",
                    "Intermediate C17", "Immune_spec C18", "Immune_spec C19", "Immune_spec C20", "Neuron C4", 
                    "Neuron C5", "Neuron C7", "Neuron C8", "Neuron C14", "Neuron C15", "Neuron C16")

##Actual plot 1-7
my_plot1_7 <- ggplot(preplot1_7$data, aes(x=factor(id, levels = order_clusters),
                                          y=factor(features.plot, levels = order_genes1_7))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  geom_point(data = preplot1_7$data[preplot1_7$data$deg=="Yes", ], aes(x=factor(id, levels = order_clusters),
                                                                       y=factor(features.plot, levels = order_genes1_7)),
             shape = "*", size=4, color="black") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="#237A57", high = "#ff6a00") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")) + 
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE) +
  labs(color="Avg. scaled exp.", size="Exp. %")
my_plot1_7

##Actual plot 8-22
my_plot8_22 <- ggplot(preplot8_22$data, aes(x=factor(id, levels = order_clusters),
                                            y=factor(features.plot, levels = order_genes8_22))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  geom_point(data = preplot8_22$data[preplot8_22$data$deg=="Yes", ], aes(x=factor(id, levels = order_clusters),
                                                                         y=factor(features.plot, levels = order_genes8_22)),
             shape = "*", size=4, color="black") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="#237A57", high = "#ff6a00") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")) + 
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE) +
  labs(color="Avg. scaled exp.", size="Exp. %")
my_plot8_22

##Save plots
ggsave(my_plot1_7, height = 9, width = 6,
       filename = "mg_genes1_7_li.pdf",
       dpi = 700)
print("mg_genes1_7_li.pdf saved")
ggsave(my_plot8_22, height = 6, width = 6,
       filename = "mg_genes8_22_li.pdf",
       dpi = 700)
print("mg_genes8_22_li.pdf saved")
