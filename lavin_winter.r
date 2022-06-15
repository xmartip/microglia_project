#GEO: GSE63340

########PS ENRICHMENT ANALYSIS FOR DIFFERENTIALLY EXPRESSED MACROPHAGE GENES####

##Libraries + WD
require(dplyr)
require(patchwork)
require(ggplot2)
require(readxl)
require(biomaRt)
require(reshape2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE63340')

##Load objects
macrophages <- read_excel("S1_adapted.xlsx")
ps_table <- readRDS('ps_table.rds')

##Obtain genes in "human version"
#(Translation from Genbank Nucleotide Accession to Gene Symbol done in:
#https://biodbnet.abcc.ncifcrf.gov/db/db2db.php). Results saved in:
#genbank_to_symbols.tsv. It may not be totally accurate because of the version,
#however, I should find a tool that lets you include it as a parameter.
mouse_symbols <- read.delim("genbank_to_symbols_deg.tsv")
#Substitute the gene symbols for the ones that are actually useful
macrophages$Mouse <- mouse_symbols$Gene.Symbol[match(macrophages$Genes,
                                                       mouse_symbols$GenBank.Nucleotide.Accession)]

#Man it just doesn't work, I have found that bioDBnet has a tool that does the
#same, I'm gonna try it. It seems to work, and maybe even better than biomaRt,
#because it has the reference to what you have searched.
human_symbols <- read.delim("genbank_to_human_deg.tsv")

#Add the human gene symbols as a new column
macrophages$Human <- human_symbols$Gene.Symbol[match(macrophages$Genes,
                                                     human_symbols$GenBank.Nucleotide.Accession)]

##Assign PS to the genes. Note, for some reason empty names get assigned PS31,
##but luckily it doesn't matter in this case
macrophages$PS <- ps_table$PS[match(macrophages$Human,
                                    ps_table$Name)]
##Drop rows with no PS. This should't be necessary, but if I don't apply the next
##one (the PS filtering) returns the NA rows for some reason
macrophages <- macrophages[!is.na(macrophages$PS),]

##Drop the genes from PS1-7 (before chordata) and from PS23-31 (no longer shared
##with humans). Mouse and human differ in the Order (Rodentia vs Primates), but
##have the same Class (Mammalia), Subclass (Theria) and Infraclass (Eutheria),
##according to https://www.itis.gov Victor's table has extra strata: Boroeutheria
##and Euarchontoglires, which (at least according to Wikipedia) are also shared
##between them. This means that they share genes up to PS22.
human_mouse_macro_genes <- macrophages[macrophages$PS > 7 &
                                         macrophages$PS < 23,]

#As a curiosity, I want to know which genes have a human equivalent but are, in
#theory, primate-forward specific. This means there is an error either in the
#PS classification, or in the refseq ID --> symbol translation.
weird_cases <- macrophages[macrophages$PS > 22 &
                             macrophages$PS < 31,]

##Selected PS enrichment of DEG for every type of Macrophage

#Create empty lists
contingency_tables <- list()
contingency_names <- c()
fisher_tests <- list()
#Set a counter to assign the tables and test results to a list
counter = 0
#Navigate through cluster
for (i in 1:11) {
  #Navigate through PS
  for (j in 8:22) {
    counter <- counter + 1
    #Present in the PS and Cluster
    a <- nrow(human_mouse_macro_genes[human_mouse_macro_genes$`11 Clusters`==i &
                                        human_mouse_macro_genes$PS==j,])
    #Present in the PS but not in the Cluster
    b <- nrow(human_mouse_macro_genes[!human_mouse_macro_genes$`11 Clusters`==i &
                                        human_mouse_macro_genes$PS==j,])
    #Present in the Cluster but not in the PS
    c <- nrow(human_mouse_macro_genes[human_mouse_macro_genes$`11 Clusters`==i &
                                        !human_mouse_macro_genes$PS==j,])
    #Not present in the PS nor the Cluster
    d <- nrow(human_mouse_macro_genes[!human_mouse_macro_genes$`11 Clusters`==i &
                                        !human_mouse_macro_genes$PS==j,])
    #Create the contingency table (in the form of a matrix but then turn it into
    #a dataframe for comfort)
    test_matrix <- as.data.frame(matrix(data = c(a,b,c,d), ncol = 2, byrow = TRUE))
    #Assign rownames and column names in cse the table wants to be consulted manually
    rownames(test_matrix) <- c(paste0("PS",j),paste0("not_PS",j))
    colnames(test_matrix) <- c(paste0("C",i),paste0("not_C",i))
    #Assign the table to the contingency table list
    contingency_tables[[counter]] <- test_matrix
    #Store the name in a vector to assign the names later all at the same time
    #(only way I know how to do this)
    contingency_names <- c(contingency_names,paste0("PS",j,"_C",i))
    #Perform the Fisher exact test and store it into the Fisher list.
    fisher <- as.list(fisher.test(test_matrix))
    fisher_tests[[counter]] <- fisher
  }
}
#Assign the names of the elements of each list
names(contingency_tables) <- contingency_names
names(fisher_tests) <- contingency_names

##Correct the p-values (FDR adjustment)

#First I should create a dataframe with the p-values
rownames_vector <- c("PS8","PS9","PS10","PS11","PS12","PS13","PS14","PS15","PS16",
                 "PS17","PS18","PS19","PS20","PS21","PS22")
colnames_vector <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
fisher_pvalues <- as.data.frame(matrix(0, ncol = 0, nrow = length((rownames_vector))))
#Navigate through the clusters (11)
for (i in 1:11) {
  #Empty vector for column values
  column_values <- c()
  #Navigate through the PS
  for (j in 1:length(rownames_vector)) {
    #Set an index for the correct extraction of the p-values of the Fisher list
    index <- length(rownames_vector)*(i-1)+j
    #Fill the vector with the correct p-values
    column_values <- c(column_values, fisher_tests[[index]][["p.value"]])
  }
  #Fill the columns of the dataframe accordingly
  fisher_pvalues[i] <- column_values
}
#Assing rownames and column names
rownames(fisher_pvalues) <- rownames_vector
colnames(fisher_pvalues) <- colnames_vector

#Assign rownames to a new column
fisher_pvalues$`PS` <- rownames_vector

#Melt the dataframes so that all the p-values are in one column
molten_fisher_pvalues <- melt(fisher_pvalues, id.vars = "PS", variable.name = 'Cluster')

#Correct the p-values FDR-style
molten_fisher_pvalues$adj_pvalue <- p.adjust(molten_fisher_pvalues$value,
                                    method = "fdr")

##Plot the results

#minus-log-transform the p-values
molten_fisher_pvalues$log_pvalue <- -log10(molten_fisher_pvalues$value)
molten_fisher_pvalues$log_adj_pvalue <- -log10(molten_fisher_pvalues$adj_pvalue)

#Rename "value" column for "mean_spec"
names(molten_fisher_pvalues)[names(molten_fisher_pvalues) == "value"] <- "pvalue"

#Add a new column with categorical data (colors) depending on the p-values. Uses dplyr.
  #Green for p_adj OK
  #Blue for p_adj not OK but p OK
  #Red for p not OK
molten_fisher_pvalues <- molten_fisher_pvalues %>% mutate(Color =
                                                    case_when(adj_pvalue <= 0.05 ~ "p_adj<=0.05", 
                                                              adj_pvalue > 0.05 & pvalue <= 0.05  ~ "p_adj>0.05 & p<=0.05",))

#Actual dot plot
p1 <- ggplot(molten_fisher_pvalues[!is.na(molten_fisher_pvalues$Color),], aes(x = factor(Cluster, levels = colnames_vector),
                                                                      y = factor(PS, levels = rownames_vector))) + 
  geom_point(aes(size = log_pvalue, color = Color)) + theme_light() + 
  scale_color_manual(values = c("gold","black")) + scale_y_discrete(drop=FALSE, position = "right") +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_line(color = "#201b1b", size = 0.2, linetype = 1),
        legend.position = "left", axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = unit(c(0.2,0.05,0.2,0.2), "cm")) +
  labs(color = "Significance", size = "-log10(p-value)")
p1
ggsave(p1, height = 6, width = 5, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/fisher_dotplot.pdf",
       dpi = 400)

##Add barplot with number of genes per PS, and per Cluster probably, (a.k.a
##the marginals)

#Create the genes per PS dataframe
genes_per_ps <- c()
for (i in 8:22) {
  genes_per_ps <- c(genes_per_ps, nrow(human_mouse_macro_genes[human_mouse_macro_genes$PS==i,]))
}
ps_names <- c("PS8","PS9","PS10","PS11","PS12","PS13","PS14","PS15","PS16",
              "PS17","PS18","PS19","PS20","PS21","PS22")
ps_barplot <- data.frame(ps_names,genes_per_ps)
ps_barplot$log <- log2(genes_per_ps)

#Create the genes per Cluster dataframe
genes_per_cluster <- c()
for (i in 1:11) {
  genes_per_cluster <- c(genes_per_cluster, nrow(human_mouse_macro_genes[human_mouse_macro_genes$`11 Clusters`==i,]))
}
cluster_names <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
cluster_barplot <- data.frame(cluster_names,genes_per_cluster)
cluster_barplot$log <- log2(genes_per_cluster)
#In this case, for reasons the order must be inverse
cluster_order <- c("C11","C10","C9","C8","C7","C6","C5","C4","C3","C2","C1")

#Actual genes per PS barplot
p2 <- ggplot(ps_barplot, aes(x=factor(ps_names, levels=ps_names), y=genes_per_ps)) + 
  geom_bar(position = position_dodge(), stat="identity", width = 0.5, show.legend = FALSE,
           fill = "#37ccff") +
  coord_flip() +
  geom_text(aes(label=genes_per_ps), hjust = 1.1, size=2.9)+theme_light() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))
p2
#Save the plot
ggsave(p2, height = 6, width = 3, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/fisher_barplot_ps.pdf",
       dpi = 400)

#Actual genes per Cluster barplot
p3 <- ggplot(cluster_barplot, aes(x=factor(cluster_names, levels=cluster_order), y=genes_per_cluster)) + 
  geom_bar(position = position_dodge(), stat="identity", width = 0.5, show.legend = FALSE,
           fill = "#37ccff") +
  coord_flip() + scale_y_reverse() + scale_x_discrete(name = "", position = "top") +
  geom_text(aes(label=genes_per_cluster), hjust = -0.2, size=2.9)+theme_light() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"))
p3
#Save the plot
ggsave(p3, height = 5, width = 3, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/fisher_barplot_cluster.pdf",
       dpi = 400)
####################ADD HEATMAP TO TRUE MG GENES DOTPLOT#######################

#Dotplot available in TRUE MG SPECIFIC GENES section from trevino_muller.r script

##Libraries + WD
require(dplyr)
require(patchwork)
require(ggplot2)
require(readxl)
require(reshape2)
require(viridis)
require(RColorBrewer)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE63340')

##Load objects
macrophages <- read_excel("GSE63340_OUT_merged4.norm.xlsx")
ps_table <- readRDS('ps_table.rds')
mg_genes <- readRDS('heatmap_order_genes8_31_soft.rds')

##Obtain genes in "human version"
#(Translation from Genbank Nucleotide Accession to Gene Symbol done in:
#https://biodbnet.abcc.ncifcrf.gov/db/db2db.php). Results saved in:
#genbank_to_symbols.tsv. It may not be totally accurate because of the version,
#however, I should find a tool that lets you include it as a parameter.
mouse_symbols <- read.delim("genbank_to_symbols_full.tsv")
#Substitute the gene symbols for the ones that are actually useful
macrophages$Mouse <- mouse_symbols$Gene.Symbol[match(macrophages$Genes,
                                                     mouse_symbols$GenBank.Nucleotide.Accession)]

#Add the human gene symbols as a new column (also done manually via BioDBnet)
human_symbols <- read.delim("genbank_to_human_full.tsv")

#Add the human gene symbols as a new column
macrophages$Human <- human_symbols$Gene.Symbol[match(macrophages$Genes,
                                                     human_symbols$GenBank.Nucleotide.Accession)]

##Assign PS to the genes
macrophages$PS <- ps_table$PS[match(macrophages$Human,
                                    ps_table$Name)]

##Drop rows with no PS. This shouldn't be necessary, but if I don't apply the next
##one (the PS filtering) returns the NA rows for some reason
macrophages <- macrophages[!is.na(macrophages$PS),]

#Select genes from PS8-22 (the ones between PS of interest and shared with Humans)
human_mouse_macro_genes <- macrophages[macrophages$PS > 7 &
                                         macrophages$PS < 23,]

#Select also for PS1-7
human_mouse_macro_genes1_7 <- macrophages[macrophages$PS <= 7,]

#Genes with either wrong PS assignation or wrong translation due to versions.
weird_cases <- macrophages[macrophages$PS > 22 &
                             macrophages$PS < 31,]

###For PS8-22

##Heatmap of true MG genes

#Filtering the mg_macro_genes dataframe by mg_genes
mg_macro_genes <- human_mouse_macro_genes[match(mg_genes, human_mouse_macro_genes$Human),]
mg_macro_genes <- mg_macro_genes[!is.na(mg_macro_genes$Human),]

#Cleaning the table so that the heatmap is easier (probably not necessary)
heatmap_df <- mg_macro_genes
heatmap_df <- heatmap_df[,-c(1,2,21,23)]
#This names are in order, so it's more comfortable to change the colnames
heatmap_columns <- c("Neutrophile 1","Neutrophile 2","Monocyte 1","Monocyte 2",
                     "Microglia 1","Microglia 2","Peritoneum 1","Peritoneum 2",
                     "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                     "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                     "Lung 1","Lung 2","Symbols")
colnames(heatmap_df) <- heatmap_columns

#I think I have to melt the table
molten_heatmap_df <- melt(heatmap_df, id.vars = "Symbols", variable.name = 'Macrophage')
#Add the PS here (if you add it before, the melt step removes it)
molten_heatmap_df$PS <- mg_macro_genes$PS[match(molten_heatmap_df$Symbols,mg_macro_genes$Human)]
#Use the log10, otherwise some gene takes all the glory. Set the -inf values of
#the 0-expression genes to 0, otherwise the plot returns a horrible brown color
molten_heatmap_df$logExpression <- (log10(molten_heatmap_df$value))
#Set the -inf values of the 0-expression genes to 0, otherwise the plot returns 
#a horrible brown color
molten_heatmap_df <- molten_heatmap_df %>% 
  mutate(logExpression = if_else(logExpression < 0, 0, logExpression))

#Set the macrophages in a better order
order_macro <- c("Microglia 1","Microglia 2","Lung 1","Lung 2",
                   "Peritoneum 1","Peritoneum 2",
                   "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                   "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                   "Neutrophile 1","Neutrophile 2",
                   "Monocyte 1","Monocyte 2")

#Missing genes (add manually). Without this I think I won't be able to match the
#heatmap with the dotplot. I basically need them because of the PS separation,
#because MNDA is not in any other present PS and as such it doesn't appear.
missing_genes <- subset(mg_genes, !(mg_genes %in% unique(molten_heatmap_df$Symbols)))

#Add them to molten_heatmap_df for every Macrophage, add the PS, but add an NA to
#logExpression so that it appears white. This allows me to drop the scale_y_discrete(drop = FALSE)
#condition of the ggplot, which prevented the facet_grid to free_scale_y properly
for (gene in missing_genes) {
  for (macrophage in order_macro) {
    missing_ps <- ps_table$PS[match(gene,ps_table$Name)]
    molten_heatmap_df <- molten_heatmap_df %>% add_row(Symbols = gene, Macrophage = macrophage,
                                  PS = missing_ps)
  }
}

#Actual heatmap plot
p1 <- ggplot(data = molten_heatmap_df, aes(x = factor(Macrophage, levels = order_macro),
                       y = factor(Symbols, levels = mg_genes),
                       fill = logExpression)) + 
  geom_tile() + theme_light() + #scale_y_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0),"cm"),
        axis.text.y = element_blank(), strip.background = element_blank()) +
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE) +
  labs(fill = "Log. Norm.\nExpression") + 
  scale_fill_viridis(discrete=FALSE) + ggtitle("Heatmap Microglia Specific Genes PS8-22")
p1
##Save objects
ggsave(p1, height = 6, width = 5,
       filename = "heatmap_mg_macro_genes8_31_soft.pdf",
       dpi = 700)

###For PS1-7

##Change the list
mg_genes <- readRDS('heatmap_order_genes1_7_soft.rds')
##Heatmap of true MG genes

#Filtering the mg_macro_genes dataframe by mg_genes
mg_macro_genes <- human_mouse_macro_genes1_7[match(mg_genes, human_mouse_macro_genes1_7$Human),]
mg_macro_genes <- mg_macro_genes[!is.na(mg_macro_genes$Human),]

#Cleaning the table so that the heatmap is easier (probably not necessary)
heatmap_df <- mg_macro_genes
heatmap_df <- heatmap_df[,-c(1,2,21,23)]
#This names are in order, so it's more comfortable to change the colnames
heatmap_columns <- c("Neutrophile 1","Neutrophile 2","Monocyte 1","Monocyte 2",
                     "Microglia 1","Microglia 2","Peritoneum 1","Peritoneum 2",
                     "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                     "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                     "Lung 1","Lung 2","Symbols")
colnames(heatmap_df) <- heatmap_columns

#I think I have to melt the table
molten_heatmap_df <- melt(heatmap_df, id.vars = "Symbols", variable.name = 'Macrophage')
#Add the PS here (if you add it before, the melt step removes it)
molten_heatmap_df$PS <- mg_macro_genes$PS[match(molten_heatmap_df$Symbols,mg_macro_genes$Human)]
#Use the log10, otherwise some gene takes all the glory. Set the -inf values of
#the 0-expression genes to 0, otherwise the plot returns a horrible brown color
molten_heatmap_df$logExpression <- (log10(molten_heatmap_df$value))
#Set the -inf values of the 0-expression genes to 0, otherwise the plot returns 
#a horrible brown color
molten_heatmap_df <- molten_heatmap_df %>% 
  mutate(logExpression = if_else(logExpression < 0, 0, logExpression))

#Set the macrophages in a better order
order_macro <- c("Microglia 1","Microglia 2","Lung 1","Lung 2",
                 "Peritoneum 1","Peritoneum 2",
                 "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                 "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                 "Neutrophile 1","Neutrophile 2",
                 "Monocyte 1","Monocyte 2")

#Missing genes (add manually). Without this I think I won't be able to match the
#heatmap with the dotplot. I basically need them because of the PS separation,
#because MNDA is not in any other present PS and as such it doesn't appear.
missing_genes <- subset(mg_genes, !(mg_genes %in% unique(molten_heatmap_df$Symbols)))

#Add them to molten_heatmap_df for every Macrophage, add the PS, but add an NA to
#logExpression so that it appears white. This allows me to drop the scale_y_discrete(drop = FALSE)
#condition of the ggplot, which prevented the facet_grid to free_scale_y properly
for (gene in missing_genes) {
  for (macrophage in order_macro) {
    missing_ps <- ps_table$PS[match(gene,ps_table$Name)]
    molten_heatmap_df <- molten_heatmap_df %>% add_row(Symbols = gene, Macrophage = macrophage,
                                                       PS = missing_ps)
  }
}

#Actual heatmap plot
p1 <- ggplot(data = molten_heatmap_df, aes(x = factor(Macrophage, levels = order_macro),
                                           y = factor(Symbols, levels = mg_genes),
                                           fill = logExpression)) + 
  geom_tile() + theme_light() + #scale_y_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0),"cm"),
        axis.text.y = element_blank(), strip.background = element_blank()) +
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE) +
  labs(fill = "Log. Norm.\nExpression") + 
  scale_fill_viridis(discrete=FALSE) + ggtitle("Heatmap Microglia Specific Genes PS8-22")
p1
##Save objects
ggsave(p1, height = 9, width = 5,
       filename = "heatmap_mg_macro_genes1_7_soft.pdf",
       dpi = 700)

#########################EMERCENCY BAR PLOT###################################

##Objects + WD
require(dplyr)
require(patchwork)
require(ggplot2)
require(readxl)
require(reshape2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE63340')

##Load Objects
macrophages <- read_excel("GSE63340_OUT_merged4.norm.xlsx")

##Select 3 genes: One that is clearly expressed in MG (TREM2), one that is more or less
##expressed in all macrophages (I'm going to choose TYROBP), and FOXP2. Done
##almost manually because because the other methods didn't work
barplot_genes <- macrophages[c(19262,15330,19369),]
barplot_genes <- macrophages[c(5365,15330,19369),]
#Put the symbols in mouse-way
barplot_genes$Symbols <- c("Tmem119", "Tyrobp", "Foxp2")
barplot_genes$Symbols <- c("Trem2", "Tyrobp", "Foxp2")
#Remove useless info right now
barplot_genes$Genes <- NULL

##Prepare the plot

#Make the columns name prettier
colnames(barplot_genes) <- c("Symbols","Neutrophile 1","Neutrophile 2","Monocyte 1","Monocyte 2",
                     "Microglia 1","Microglia 2","Peritoneum 1","Peritoneum 2",
                     "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                     "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                     "Lung 1","Lung 2")
#Melt the plot so that I can log-transform the expression values
molten_barplot_genes <- melt(barplot_genes, id.vars = "Symbols", variable.name = 'Macrophage')
molten_barplot_genes$Log <- log10(molten_barplot_genes$value)
#Set the -inf values to 0
molten_barplot_genes <- molten_barplot_genes %>% 
  mutate(Log = if_else(Log < 0, 0, Log))

#Add Fill column to color the replicates with the same color
molten_barplot_genes <- molten_barplot_genes %>% mutate("Macro. Tissue" =
                                                    case_when(Macrophage == "Neutrophile 1" ~ "Neutrophile", 
                                                              Macrophage == "Neutrophile 2" ~ "Neutrophile",
                                                              Macrophage == "Monocyte 1" ~ "Monocyte",
                                                              Macrophage == "Monocyte 2" ~ "Monocyte",
                                                              Macrophage == "Microglia 1" ~ "Microglia",
                                                              Macrophage == "Microglia 2" ~ "Microglia",
                                                              Macrophage == "Peritoneum 1" ~ "Peritoneum",
                                                              Macrophage == "Peritoneum 2" ~ "Peritoneum",
                                                              Macrophage == "Spleen 1" ~ "Spleen",
                                                              Macrophage == "Spleen 2" ~ "Spleen",
                                                              Macrophage == "Kupffer 1" ~ "Kupffer",
                                                              Macrophage == "Kupffer 2" ~ "Kupffer",
                                                              Macrophage == "L. Intest. 1" ~ "L. Intest.",
                                                              Macrophage == "L. Intest. 2" ~ "L. Intest.",
                                                              Macrophage == "S. Intest. 1" ~ "S. Intest.",
                                                              Macrophage == "S. Intest. 2" ~ "S. Intest.",
                                                              Macrophage == "Lung 1" ~ "Lung",
                                                              Macrophage == "Lung 2" ~ "Lung"))
#Specify the order in which I want the mcarophages
order_macro <- c("Microglia 1","Microglia 2","Lung 1","Lung 2",
                 "Peritoneum 1","Peritoneum 2",
                 "Spleen 1","Spleen 2","Kupffer 1","Kupffer 2",
                 "L. Intest. 1","L. Intest. 2","S. Intest. 1","S. Intest. 2",
                 "Neutrophile 1","Neutrophile 2",
                 "Monocyte 1","Monocyte 2")

##Actual bar plot
my_barplot <- ggplot(molten_barplot_genes, aes(x=factor(Macrophage, levels=order_macro), y=Log, fill=`Macro. Tissue`)) + 
  geom_bar(stat="identity") + facet_wrap(~Symbols, scales="free_y", ncol=1) +
  theme_light() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Log. Norm. Expression (ESAT software units)") + 
  scale_fill_discrete(breaks=c('Microglia', 'Lung', 'Peritoneum', 'Spleen', 'Kupffer',
                               'L. Intest.', 'S. Intest.', 'Neutrophile', 'Monocyte'))

my_barplot
##Save the plot
ggsave(my_barplot, height = 4, width = 4,
       filename = "mouse_genes_piled_barplot_2.pdf",
       dpi = 700)
