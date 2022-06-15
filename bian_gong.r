#Paper with 2 GEOs: GSE133345 (STRT-seq) and GSE137010 (scRNA-seq)

##########################DEVELOPING MG SPECIFIC GENES##########################

#Using Supplementary Table 2 (S2). Cluster Mac_4 corresponds to developing MG.
#S2 corresponds to STRT-seq data

##Libraries + WD
require(dplyr)
require(patchwork)
require(ggplot2)
require(readxl)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/bian_gong_2020/GSE133345')

##Load objects
bian_mg_markers <- readRDS('bian_mg_markers.rds')
trevino_mg_markers <- readRDS('mg_markers.rds')
ps_table <- readRDS('ps_table.rds')

##Prepare the bian_mg_markers object to make it more comfortable to work with.
#(already done and saved, so it's hushed)
#bian_mg_markers <- read_excel('S2.xlsx',sheet = 2)
#bian_mg_markers <- bian_mg_markers %>% filter(cluster=="Mac_4" & avg_logFC>=1 &
#                                                p_val_adj<0.001)
#bian_mg_markers$...1<-NULL
#saveRDS(bian_mg_markers, file = 'bian_mg_markers.rds')

##Prepare the trevino_mg_markers_dataset
trevino_mg_markers <- trevino_mg_markers %>% filter(avg_log2FC>=1 &
                                                      p_val_adj<0.001)

##Join the 2 dataframes by gene symbol
human_mg_markers <- left_join(trevino_mg_markers,bian_mg_markers, 
                                by=c("gene_symbol" = "gene"))
human_mg_markers <- na.omit(human_mg_markers)

##Add PS
human_mg_markers$PS <- ps_table$PS[match(human_mg_markers$gene_symbol,ps_table$Name)]

##Better column names
human_mg_markers$cluster <- NULL
column_names <- c("pval_brain","avg_log2FC_brain","pct.1_brain", "pct.2_brain",
                  "adj_pval_brain","gene_symbol","pval_macro","avg_logFC_macro",
                  "pct.1_macro", "pct.2_macro","adj_pval_macro","PS")
colnames(human_mg_markers) <- column_names
#Change order of column with dplyr's relocate
human_mg_markers <- human_mg_markers %>% relocate(gene_symbol, .before = pval_brain)

############################PREPARING SEURAT OBJECT############################

##Libraries + WD
require(Seurat)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/bian_gong_2020/GSE133345')

##Load Objects
counts <- read.delim(file = "GSE133345_Quality_controled_UMI_data_of_all_1231_embryonic_cells.txt", sep = " ")
metadata <- read.delim(file = "GSE133345_Annotations_of_all_1231_embryonic_cells_updated_0620.txt", sep = " ")

##Create Seurat Object
bian_seurat <- CreateSeuratObject(counts=counts,
                                  project="bian",
                                  meta.data=metadata,
                                  min.cells = 0) # Genes need to be expressed in at least 0 cells

##Normalization
bian_seurat <- NormalizeData(bian_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

##Variable Features (just wanna check)
bian_seurat <- FindVariableFeatures(object = bian_seurat)
length(VariableFeatures(bian_seurat))
##Save Seurat Object (for cluster, up until normalization)
saveRDS(bian_seurat, file = "bian_seurat.rds")

##Scale Data (of variable features, namely bian_seurat_toy.rds)
bian_seurat <- ScaleData(bian_seurat, vars.to.regress = c('nUMI','nGene'))
##Save Seurat Object
saveRDS(bian_seurat, file = "bian_seurat_toy.rds")

###########################SCALE ALL GENES (CLUSTER)############################

##Packages, options and WD
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(Seurat)
options(stringsAsFactors = F)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD OK")

##Load objects
bian_seurat<-readRDS('bian_seurat.rds')
print("Objects OK")

##Scale Data
all.genes <- rownames(bian_seurat)
print("Scaling all genes...")
bian_seurat <- ScaleData(object = bian_seurat, features = all.genes, vars.to.regress = c('nUMI','nGene'))

#Save the modified objects
saveRDS(bian_seurat, file = 'bian_seurat.rds')
print("R job done")

###########################DEVELOPMENT STAGES##################################

##Libraries + WD
require(Seurat)
require(dplyr)
require(patchwork)
require(ggplot2)
require(readxl)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/bian_gong_2020/GSE133345')

##Libraries + WD (for cluster)
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Seurat)
require(patchwork)
require(ggplot2)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD loaded correctly")

##Load Objects
mg_genes1_7 <- readRDS('heatmap_order_genes1_7_hard.rds')
mg_genes8_22 <- readRDS('heatmap_order_genes8_31_hard.rds')
ps_table <- readRDS('ps_table.rds')
bian_seurat <- readRDS('bian_seurat_toy.rds')
#deg_clusters <- read_excel('S2.xlsx', sheet = 2)
#saveRDS(deg_clusters, file = "marker_genes_bian_clusters.rds")
deg_clusters <- readRDS("marker_genes_bian_clusters.rds")

##Set clusters as indentations
Idents(bian_seurat) <- bian_seurat@meta.data$cluster

##Translate the symbols into IDs and store it as a character (just in case)
id_list1_7 <- as.character(mg_genes1_7)
id_list8_22 <- as.character(mg_genes8_22)

##Store the DotPlot fuction data into an object for further manipulation
preplot1_7 <- DotPlot(bian_seurat, features = id_list1_7, cluster.idents = TRUE)
preplot8_22 <- DotPlot(bian_seurat, features = id_list8_22, cluster.idents = TRUE)


##Add the PS column
preplot1_7$data$PS <- ps_table$PS[match(preplot1_7$data$features.plot,
                                        ps_table$Name)]
preplot8_22$data$PS <- ps_table$PS[match(preplot8_22$data$features.plot,
                                         ps_table$Name)]

##Add significance
for (i in 1:nrow(deg_clusters)) {
  gene <- deg_clusters$gene[i]
  cluster <- deg_clusters$cluster[i]
  preplot1_7$data$deg[preplot1_7$data$features.plot == gene & preplot1_7$data$id == cluster] <- "Yes"
  preplot8_22$data$deg[preplot8_22$data$features.plot == gene & preplot8_22$data$id == cluster] <- "Yes"
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
order_clusters <- c("YSMP","ErP","MkP",
                    "GMP","Myeloblast","Monocyte", "Mac_1", "Mac_2",
                    "Mac_3", "Mac_4", "HSPC", "CD7loP", "CD7hiP",
                    "ILC", "Mast cell")

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
       filename = "mg_genes1_7_bian.pdf",
       dpi = 700)
print("mg_genes1_7_bian.pdf saved")
ggsave(my_plot8_22, height = 6, width = 6,
       filename = "mg_genes8_22_bian.pdf",
       dpi = 700)
print("mg_genes8_22_bian.pdf saved")
print("R job done")
