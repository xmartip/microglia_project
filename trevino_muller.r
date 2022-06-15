#Actual utility script from the Dataset obtained from Cell article Chromatin and 
#gene-regulatory dynamics of the developing human cerebral cortex at single-cell
#resolution, From Trevino A., Müller F. et al.

#Objective of the script is to find marker genes from Microglia cluster, translate
#them into their actual name, and assign them a specificity score for Chordata vs.
#non-Chordata. For the specificity score, a signature matrix must be made beforehand.

#Update: The specificity scores and non-parametric statistical tests were not a
#good approach, so we decided to finish the EWCE tutorial workflow and apply a
#bootstrap enrichment test. The results seem better and are more logical.

########################PREPARING OBJECTS######################################

#Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(EWCE)
library(Matrix)

#Working directory
setwd('/home/onecellardoor/Dropbox/imim/data/GSE162170')

#Load objects
#dhcc<-readRDS('dhcc_v2.rds')
#scaled_data<-readRDS('scaled_data_v2.rds')
annotation_df<-readRDS('annotation_df.rds')
tr_table_gencode_v27<-readRDS('tr_table_gencode_v27.rds')
counts<-readRDS('counts.rds') #From the dhcc_V2, it's a dgCMatrix
#subset_counts <- counts[1000:3000,] <-- Small set for experimenting
#saveRDS(subset_counts, file = '2000_counts.rds')
#counts <- readRDS('2000_counts.rds')
#counts<-readRDS('counts_df.rds') <-- Directly from the .tsv file, I can't load it.
scT<-readRDS('sct_2000_counts.rds')
corrected_counts<-readRDS('corrected_2000_counts.rds')
spec_matrix <- readRDS('spec_matrix.rds')
spec_matrix_normed <- readRDS('spec_matrix_normed.rds')


#After generate_celltype_data run (either here or in the cluster)
celltypedata <- load(file='CellTypeData_dhcc.rda') #From the 2000 genes subset
celltypedata <- load(file='CellTypeData_ewce_all_genes.rda')

#Okay, so this is the hard part. I will try to use the EWCE or edgeR package from
#Bioconductor to get the signature matrix. The manual can be accessed via:
browseVignettes("EWCE")
edgeRUsersGuide() #Done, it downloads a PDF file that I saves in "utilities" ->
# "manuals".

#I am going to try EWCE first. Installation:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EWCE")

#Load objects
dhcc<-readRDS('dhcc_v2.rds')
scaled_data<-readRDS('scaled_data_v2.rds')
annotation_df<-readRDS('annotation_df.rds')
tr_table_gencode_v27<-readRDS('tr_table_gencode_v27.rds')
counts<-readRDS('counts.rds') #From the dhcc_V2, it's a dgCMatrix
subset_counts <- counts[1000:3000,] #Small set for experimenting
saveRDS(subset_counts, file = '2000_counts.rds')
counts <- readRDS('2000_counts.rds')
counts<-readRDS('counts_df.rds') #Directly from the .tsv file, I can't load it.
scT<-readRDS('sct_2000_counts.rds')
corrected_counts<-readRDS('corrected_2000_counts.rds')
spec_matrix <- readRDS('spec_matrix.rds')
spec_matrix_normed <- readRDS('spec_matrix_normed.rds')

#If you wish to take the @scale.data values out of Seurat (to import into another
#package), you can access them using the GetAssayData function (or just take the
#object@scale.data matrix).
scaled_data <- GetAssayData(object=dhcc, slot="scale.data")
saveRDS(scaled_data, file='scaled_data_v2.rds')

#Subsetting the metadata that I want with "select" function from package dplyr:
annotation_df <- select(dhcc@meta.data, Cell.ID, Name)
saveRDS(annotation_df, file='annotation_df.rds')

###########################SPECIFICITY MATRIX##################################

#Load libraries
library(dplyr)
library(patchwork)
library(EWCE)
library(sctransform)
library(Matrix)

#Working directory
setwd('/home/onecellardoor/Dropbox/imim/data/GSE162170')

##Load objects

annotation_df<-readRDS('annotation_df.rds')
tr_table_gencode_v27<-readRDS('tr_table_gencode_v27.rds')
counts<-readRDS('counts.rds')
corrected_counts_normed <- readRDS('corrected_counts_normed.rds')
corrected_counts <- readRDS('sct_corrected.rds')
spec_matrix <- readRDS('spec_matrix.rds')
spec_matrix_normed <- readRDS('spec_matrix_normed.rds')
print("Objects OK")

##Actual job

#Normalizig and scaling via scTransform, as opposed to logNormalization + scaling
#n_genes set to NULL to include all genes (default is 2000)
print("The number of genes won't match because there are genes that are not expressed in ANY cell")
scT = sctransform::vst(counts, n_genes = NULL, min_cells = 1, return_cell_attr = TRUE)
saveRDS(scT, file='sct.rds')
print("sctransform OK and saved as sct.rds")
#Correcting the scT
corrected_counts = correct_counts(scT, counts)
saveRDS(corrected_counts, file='sct_corrected.rds')
print("corrected_counts OK and saves as corrected_counts.rds")

#Transposition step
corrected_counts<-readRDS('sct_corrected.rds')
corrected_counts_normed =
  Matrix::t(Matrix::t(corrected_counts)*(1/Matrix::colSums(corrected_counts)))
saveRDS(corrected_counts_normed, file = 'corrected_counts_normed.rds')

#generate_celltype_data saves an rda workspace that is an object, specify where
#with savePath option

annotLevels = list(cell_type=annotation_df$Name)
spec_matrix_normed <- generate_celltype_data(exp=corrected_counts_normed, annotLevels=annotLevels,
                      groupName="ewce_all_genes_normed", savePath='/users/genomics/xmarti/data')
spec_matrix <- generate_celltype_data(exp=corrected_counts, annotLevels=annotLevels,
               groupName="ewce_all_genes", savePath='/users/genomics/xmarti/data')



#load the workspace a.k.a. an object called "ctd"
celltypedata <- load(file='CellTypeData_ewce_all_genes_normed.rda')
spec_matrix_normed <- ctd[[1]][["specificity"]]
saveRDS(spec_matrix_normed, file='spec_matrix_normed.rds')
print("spec_matrix_normed OK")

celltypedata <- load(file='CellTypeData_ewce_all_genes.rda')
my_spec_matrix <- ctd[[1]][["specificity"]]
saveRDS(spec_matrix, file='spec_matrix.rds')
print("spec_matrix OK")

#Translate de ENSEMBL IDs into Gene Symbols
ensembl_ids <- rownames(spec_matrix_normed)
gene_symbols <- tr_table_gencode_v27$Gene.name[match(ensembl_ids,
                tr_table_gencode_v27$`Gene.stable.ID`)]
rownames(spec_matrix_normed) <- gene_symbols
spec_matrix_normed <- as.data.frame(spec_matrix_normed)
mg_spec_normed <- select(spec_matrix_normed, MG)
saveRDS(mg_spec_normed, file = 'mg_spec_normed.rds')

ensembl_ids <- rownames(spec_matrix)
gene_symbols <- tr_table_gencode_v27$Gene.name[match(ensembl_ids,
                tr_table_gencode_v27$`Gene.stable.ID`)]
rownames(spec_matrix) <- gene_symbols
spec_matrix <- as.data.frame(spec_matrix)
mg_spec <- select(spec_matrix, MG)
saveRDS(mg_spec, file = 'mg_spec.rds')

######################STATISTICAL TEST (NOT USED)###############################

#I have the norm-scaled specificity matrix. As of now, it seems a reasonable way
#to perform the across-phylostrata microglia-specific gene enrichment analysis.

#The legend for the phylostrata table is as follows:
# PS 1 means a protein has sequence homology with (at least one species of) bacteria
# 8 means a protein is shared with Chordates
# 18 - Mammals
# 23 - Primates
# 31 is human-specific: no homologs in other species

## Libraries + WD

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

## Object preparation
#USE READ.DELIM!!!!!!!!!!!!! READ.TABLE HAS A WARNING AND TRIMMED 2/3 OF THE
#GENES
ps_table <- read.delim(file = 'Hs_Ens89+2102_PS_seq_etc_hg38.txt', sep = '\t', header = TRUE)
saveRDS(ps_table, file = 'ps_table.rds')

## Load objects

ps_table<-readRDS('ps_table.rds')
mg_spec <- readRDS('mg_spec.rds')
mg_spec_normed <- readRDS('mg_spec_normed.rds')
#mg_spec_with_ps <- readRDS('mg_spec_with_ps.rds')
mg_spec_normed_with_ps <- readRDS('mg_spec_normed_with_ps.rds')

## Adding the PS to the MG genes
#mg_spec$Name <- rownames(mg_spec)
#mg_spec$PS <- ps_table$PS[match(mg_spec$Name, ps_table$Name)]

mg_spec_normed$Name <- rownames(mg_spec_normed)
mg_spec_normed$PS <- ps_table$PS[match(mg_spec_normed$Name, ps_table$Name)]

## Getting just the names with a PS value assigned (dplyr used)
#mg_spec_with_ps <- dplyr::filter(mg_spec, !is.na(PS) | PS != "")
#saveRDS(mg_spec_with_ps, file = 'mg_spec_with_ps.rds')

mg_spec_normed_with_ps <- dplyr::filter(mg_spec_normed, !is.na(PS) | PS != "")
saveRDS(mg_spec_normed_with_ps, file = 'mg_spec_normed_with_ps.rds')

#Separating the sets into phylostrata

for(i in 1:31) {
  #assign allows you to associate the iterator value with the object
  assign(paste0("ps_", i), 
         subset(mg_spec_normed_with_ps, PS == i, select=c(MG, Name, PS)))
}

#Generating a data frame with the number of genes per each phylostratum. Done
#manually because I have ended faster than doing it the fancy way
n_genes <- data.frame(PS=c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8",
                                 "PS9","PS10","PS11","PS12","PS13","PS14","PS15","PS16","PS17",
                                 "PS18","PS19","PS20","PS21","PS22","PS23","PS24","PS25","PS26",
                                 "PS27","PS28","PS29","PS30","PS31"), N_GENES=c(9334,3897,566,
                                  905,542,565,68,119,0,48,713,0,209,57,1,56,139,39,76,163,44,6,6,1,30,
                                  17,2,0,1,0,27))

#Save table because it is useful and may be used in my thesis report. I checked
#if it saves correctly and is able to be read later. It does... in R, not it Excel
#(maybe because it has rownames I don't know)
saveRDS(n_genes, file = 'ngenes_per_ps.rds')
write.table(n_genes, 
            file = '/home/onecellardoor/Dropbox/imim/records/tables/gse162170_ngenes_per_ps.tsv',
            sep = "\t")

#Plotting the density function of the gene specificity of each phylostratum
for(i in 1:31) {
  #assign allows you to associate the iterator value with the object
  assign(paste0("density_plot_", i), 
         ggplot(paste0("ps_", i), aes(x=MG)) +   geom_density() + 
           geom_label(x=0.75, y=8, label=paste0("Obs = ",nrow(ps_16))))
}

plot1<-ggplot(ps_1to7, aes(x=MG)) + geom_density() + ggtitle("Density plot of PS 1-7") +
  xlab("Specificity") + ylab("Density") + 
  geom_label(x=0.50, y=20, label=paste0("N_Genes = ",nrow(ps_1to7)))

plot2<-ggplot(ps_8to11, aes(x=MG)) + geom_density() + ggtitle("Density plot of PS 8-11") +
  xlab("Specificity") + ylab("Density") + 
  geom_label(x=0.50, y=20, label=paste0("N_Genes = ",nrow(ps_8to11)))
grid.arrange(plot1, plot2, ncol=2)

plot3<-ggplot(ps_8, aes(x=MG)) + geom_density() + ggtitle("Density plot of PS 8") +
  xlab("Specificity") + ylab("Density") + 
  geom_label(x=0.50, y=20, label=paste0("N_Genes = ",nrow(ps_8)))
grid.arrange(plot1, plot3, ncol=2)

#Non-chordates (PS1 - PS7)
ps_1to7 <- subset(mg_spec_normed_with_ps, PS < 8, select=c(MG, Name, PS))
#Provisional "chordates" (PS8 - PS11)
ps_8to11 <- subset(mg_spec_normed_with_ps, PS >= 8 & PS < 12, select=c(MG, Name, PS))

## Applying the K-S test of chordates vs "non chordates"
ks.test(ps_1to7$MG, ps_8$MG, alternative = "two.sided")
wilcox.test(ps_1to7$MG, ps_8$MG, alternative = "two.sided")

############AUTOMATED FOR ALL CELL TYPES AND PS (NOT USED)######################

## I already have the specificity matrix with the gene symbols, so I just have to
## do the significance test and the plots for every cell type and PS.

#Load libraries + Set WD
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

#Load objects
spec_matrix <- readRDS('spec_matrix_normed.rds')
ps_table <- readRDS('ps_table.rds')

# Adding the PS to the spec_matrix genes when possible
spec_matrix$Gene_Name <- rownames(spec_matrix)
spec_matrix$PS <- ps_table$PS[match(spec_matrix$Gene_Name, ps_table$Name)]

#Select only those genes with PS assigned. With this step, we now have a matrix
#with 8267 genes, whereas the original spec_matrix had 28009 genes = almost 20000
#gene dropout (no PS available). UPDATE: THE PS TABLE WAS WRONG, NOW IT IS
#CORRECT AND THE RESULTING FILTERED MATRIX HAS 17631 GENES.
spec_matrix_with_ps <- dplyr::filter(spec_matrix, !is.na(PS) | PS != "")
#saveRDS(spec_matrix_with_ps, file = "spec_matrix_with_ps.rds")

#Generate a dataframe for each cell type. I will store them in a list so I can
#iterate over it
cell_type_specs_with_ps <- list()
for (i in colnames(spec_matrix_with_ps)[1:23]){
  print(colnames(spec_matrix_with_ps[i]))
  cell_type_specs_with_ps[[i]] <- select(spec_matrix_with_ps, i, Gene_Name, PS)
}

#Separate each cell type gene set into PS. Here I do a list of lists, and fill
#the latter lists with dataframes via the subset() function.
individual_ps_per_cell_type <- list()
for (i in names(cell_type_specs_with_ps)) {
  individual_ps_per_cell_type[[i]] <- list()
  #PS 1-7 together
  individual_ps_per_cell_type[[i]][[1]] <- subset(cell_type_specs_with_ps[[i]], 
                                                  PS < 8, 
                                                  select=c(colnames(cell_type_specs_with_ps[[i]])[1], 
                                                           "Gene_Name", "PS"))
 
  #PS 8-11 together
  individual_ps_per_cell_type[[i]][[2]] <- subset(cell_type_specs_with_ps[[i]], 
                                                  PS >= 8 & PS < 12, 
                                                  select=c(colnames(cell_type_specs_with_ps[[i]])[1], 
                                                           "Gene_Name", "PS"))
  #Create a vector with the list names to assign them all at the same time
  #(it's the only way that works, 1 on 1 names() and assign() function do not)
  list_names <- c(paste0(names(individual_ps_per_cell_type[i])," PS1-7"),
                  paste0(names(individual_ps_per_cell_type[i])," PS8-11"))
  #PS 8-31 individual
  for(j in 1:max(spec_matrix_with_ps$PS)) {
    individual_ps_per_cell_type[[i]][[j]] <- subset(cell_type_specs_with_ps[[i]], 
                                                    PS == j, 
                                                    select=c(colnames(cell_type_specs_with_ps[[i]])[1], 
                                                             "Gene_Name", "PS"))
  }
  #When doing the previous loop, the list automatically created empty slots for
  #positions 3-7, and the vector of names and the objects of the list must have the
  #same number, so here I must create "useless" names for the rest of the dataframes
  #to have the correct ones.
  for(j in 1:max(spec_matrix_with_ps$PS)) {
    list_names <- c(list_names, paste0(names(individual_ps_per_cell_type[i])," PS",j))
  }
  #Now we assign the names of all the dataframes in the list
  names(individual_ps_per_cell_type[[i]]) <- list_names
}

##Eliminate empty list elements and empty dataframes (they exist but have 0 rows)

#Reason: ks.text and wilcox.test return an error and cause the program to stop
#running if at least one of the two input data vectors is empty.

#Removing the NULL elements of the lists
for (i in names(individual_ps_per_cell_type)) {
  individual_ps_per_cell_type[[i]] = individual_ps_per_cell_type[[i]][-which(sapply(individual_ps_per_cell_type[[i]], is.null))]
}

#Removing the dataframes with 0 rows == 0 genes
for (i in names(individual_ps_per_cell_type)) {
  individual_ps_per_cell_type[[i]] = individual_ps_per_cell_type[[i]][sapply(individual_ps_per_cell_type[[i]], nrow)>0]
}

##Create a dataframe with the mean specificity per cell type per PS

#Define row names (but add them in the end because it cannot be done when it's
#empty). Done manually, otherwise this step was very complex.
rownames_mean_df <- c("PS1-7","PS8-11","PS8","PS10","PS11","PS13","PS14","PS16","PS17","PS18",
                 "PS19","PS20","PS21","PS22","PS23","PS24","PS25","PS26","PS27",
                 "PS29","PS31")
mean_spec_df <- as.data.frame(matrix(0, ncol = 0, nrow = length((rownames_mean_df))))

#Fill the dataframe
for (i in names(individual_ps_per_cell_type)) {
  #Create empty vectors where the data will be stored to transfer to dataframe
  mean_values <- c()
 
  for (j in 1:length(individual_ps_per_cell_type[[i]])) {
    #Apply the tests for PS1-7 against PS8-11 and each PS with at least one gene
    #for each cell type
    my_mean <- mean(individual_ps_per_cell_type[[i]][[j]][[1]])

    #Concatenate each p-value and D obtained to the storage vector
    mean_values <- c(mean_values, my_mean)
    }
  #Add new columns with the name of the cell type with their corresponding data
  #values
  mean_spec_df[i] <- mean_values
  }

#Assign the row names to the statistical tests dataframe.
rownames(mean_spec_df) <- rownames_mean_df

##Plot the mean specificities

#Create a new column being the rownames (for the melt() step)
mean_spec_df$`X_Axis` <- rownames_mean_df
#melt() is a reshape2 function that basically converts a dataframe into another
#dataframe with less columns but more rows, the data being conserved
molten_mean_spec_df <- melt(mean_spec_df, id.vars = "X_Axis", variable.name = 'Cell_type')
#Actual plot (not working as intended for now)
ggplot(molten_mean_spec_df, aes(x = factor(X_Axis, levels=rownames_mean_df), y = value, group = Cell_type)) + 
  geom_line(aes(colour = Cell_type)) + theme_bw() +
  labs(y= "Mean Specificity", x = "Phylostrata", group = "Cell type") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##Save p-values from K-S and Wilcoxon tests in a dataframe

#Define row names (but add them in the end because it cannot be done when it's
#empty). Done manually, otherwise this step was very complex.
rownames_df <- c("PS8-11","PS8","PS10","PS11","PS13","PS14","PS16","PS17","PS18",
                 "PS19","PS20","PS21","PS22","PS23","PS24","PS25","PS26","PS27",
                 "PS29","PS31")
#Create dataframes where the statistics will go
ks_pvalues <- as.data.frame(matrix(0, ncol = 0, nrow = length((rownames_df))))
ks_d <- as.data.frame(matrix(0, ncol = 0, nrow = length((rownames_df))))
wilcox_pvalues <- as.data.frame(matrix(0, ncol = 0, nrow = length((rownames_df))))


#Fill the dataframes
for (i in names(individual_ps_per_cell_type)) {
  #Create empty vectors where the data will be stored to transfer to dataframe
  pvalues_ks <- c()
  d <- c()
  pvalues_wilcox <- c()
  
  for (j in 2:length(individual_ps_per_cell_type[[i]])) {
    #Apply the tests for PS1-7 against PS8-11 and each PS with at least one gene
    #for each cell type
    ks <- ks.test(individual_ps_per_cell_type[[i]][[1]][[1]], 
            individual_ps_per_cell_type[[i]][[j]][[1]], 
            alternative = "two.sided")
    wilcox <- wilcox.test(individual_ps_per_cell_type[[i]][[1]][[1]], 
                individual_ps_per_cell_type[[i]][[j]][[1]], 
                alternative = "two.sided")
    #Concatenate each p-value and D obtained to the storage vector
    pvalues_ks <- c(pvalues_ks, ks[["p.value"]])
    d <- c(d, ks[["statistic"]][["D"]])
    pvalues_wilcox <- c(pvalues_wilcox, wilcox[["p.value"]])
  }
  #Add new columns with the name of the cell type with their corresponding data
  #values
  ks_pvalues[i] <- pvalues_ks
  ks_d[i] <- d
  wilcox_pvalues[i] <- pvalues_wilcox
}

#Assign the row names to the statistical tests dataframes.
rownames(ks_pvalues) <- rownames_df
rownames(ks_d) <- rownames_df
rownames(wilcox_pvalues) <- rownames_df

##Correct the p-values (FDR method)

#Assign rownames to a new column
ks_pvalues$`PS` <- rownames_df
wilcox_pvalues$`PS` <- rownames_df
#Melt the dataframes so that all the p-values are in one column
molten_ks_pvalues <- melt(ks_pvalues, id.vars = "PS", variable.name = 'Cell_type')
molten_wilcox_pvalues <- melt(wilcox_pvalues, id.vars = "PS", variable.name = 'Cell_type')
#Correct the p-values FDR-style
molten_ks_pvalues$value <- p.adjust(molten_ks_pvalues$value,
                                    method = "fdr")
molten_wilcox_pvalues$value <- p.adjust(molten_wilcox_pvalues$value,
                                    method = "fdr")
#Reconstruct the table with dcast() (it's like the opposite of melt() in reshape2)
#It return a dataframe but it has no rownames and the PS order is altered (it is
#now ordered alphanumerically)
ks_pvalues_corrected <- dcast(molten_ks_pvalues, PS ~ Cell_type, value.var="value")
wilcox_pvalues_corrected <- dcast(molten_wilcox_pvalues, PS ~ Cell_type, value.var="value")

##Plot the p-values

#This one should be a dot plot of -log10(p-values) with the
#color degradation #corresponding to the -log10(p-value) and the size of the plot 
#being proportional to the mean specificity (this is provisional, it may not be
#great nor very useful)

#Create new molten dataframes for the -log10(pvalues)
molten_log_ks <- molten_ks_pvalues
molten_log_wilcox <- molten_wilcox_pvalues

#minus-log-transform the p-values
molten_log_ks$log_pvalue <- -log10(molten_log_ks$value)
molten_log_ks$value <- NULL
molten_log_wilcox$log_pvalue <- -log10(molten_log_wilcox$value)
molten_log_wilcox$value <- NULL

#Combine the p-values molten dataframes with their mean_spec value. I suppose this
#is necessary for the dot plot to perform correctly.
molten_log_ks <- left_join(molten_log_ks, molten_mean_spec_df, 
                                     by = c("PS" = "X_Axis", "Cell_type"))
molten_log_wilcox <- left_join(molten_log_wilcox, molten_mean_spec_df, 
                           by = c("PS" = "X_Axis", "Cell_type"))

#Rename "value" column for "mean_spec"
names(molten_log_ks)[names(molten_log_ks) == "value"] <- "mean_spec"
names(molten_log_wilcox)[names(molten_log_wilcox) == "value"] <- "mean_spec"

#Actual plot
order <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
           "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
           "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")

ggplot(molten_log_ks, aes(x = factor(Cell_type, levels = order), y = factor(PS, levels = rownames_df))) + 
  geom_point(aes(size = mean_spec, color = log_pvalue)) + theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low = "cadetblue", high = "tomato") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0,0.2), "cm")) +
  labs(color = "-log10(p-value)", size = "Mean Spec.") + ggtitle("K-S")

#######################BOOTSTRAP CLUSTER##################################

#I already have the specificity matrix with the gene symbols, but I need to
#finish the EWCE workflow I suppose, because with just the specificity matrix
#it's not working great. So I basically have to adapt the parameters that Skene
#et al use in his paper into mine, and it would be as follows:
#The sct_data is the ctd (this step should be the same)
#The background gene set is all the expressed genes in my data (ignore mouse orthologs)
#The set of genes (hits) is a vector containing the MG specific genes in PS8
#Reps should be 10000, but I may do 100 at first for speed while testing
#annotLevel is the only one that I have

##Libraries + wd
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(reshape2)
require(Matrix)
require(patchwork)
require(RNOmni, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
require(EWCE, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
setwd('/users/genomics/xmarti/data')
print("Libraries and WD OK")

##Load objects
#corrected_counts_normed <- readRDS('corrected_counts_normed.rds')
#annotation_df<-readRDS('annotation_df.rds')
ps_table<-readRDS('ps_table.rds')
#tr_table_gencode_v27<-readRDS('tr_table_gencode_v27.rds')
load('CellTypeData_ewce_all_genes_normed_translated.rda')
print("Objects OK")

##Prepare necessary objects

#Translating the rownames of corrected_counts_normed matrix
#ids <- rownames(corrected_counts_normed)
#gene_names <- tr_table_gencode_v27$Gene.name[match(ids, tr_table_gencode_v27$Gene.stable.ID)]
#rownames(corrected_counts_normed) <- gene_names

#I need to repeat the generate_celltype_data step
#annotLevels = list(cell_type=annotation_df$Name)
#spec_matrix_normed <- generate_celltype_data(exp=corrected_counts_normed, 
#                                             annotLevels=annotLevels,
#                                             groupName="ewce_all_genes_normed_translated", 
#                                             savePath='/users/genomics/xmarti/data')

#Background set
background <- unique(rownames(ctd[[1]][["specificity"]]))
ps_with_more_than_4_genes <- c(1,2,3,4,5,6,7,8,10,11,13,14,16,17,18,19,20,21,22,
                               23,25,26,31)
#Create list to store results
bootstrap_enrichment_analysis <- list()
for (i in ps_with_more_than_4_genes) {
  
  #Hit set
  ps <- subset(ps_table, PS == i)
  hits <- ps$Name[match(background, ps$Name)]
  hits <- hits[!is.na(hits) & hits != ""]
  
  #Reps and annotLevel
  reps = 100000
  level = 1
  
  ##Actual job
  full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                           reps=reps,annotLevel=level, 
                                           genelistSpecies = "human",
                                           sctSpecies = "human")
  
  ##Save object into the list
  bootstrap_enrichment_analysis[[i]] <- full_results
}
#Hit set for PS 8-11
ps <- subset(ps_table, PS > 7 | PS < 12)
hits <- ps$Name[match(background, ps$Name)]
hits <- hits[!is.na(hits) & hits != ""]

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         reps=reps,annotLevel=level,
                                         genelistSpecies = "human",
                                         sctSpecies = "human")

##Save object into the list
bootstrap_enrichment_analysis8_11 <- full_results

##Save objects of interest
saveRDS(bootstrap_enrichment_analysis, file ='bootstrap_enrichment_analysis.rds')
print("bootstrap_enrichment_analysis.rds saved")
saveRDS(bootstrap_enrichment_analysis8_11, file ='bootstrap_enrichment_analysis8_11.rds')
print("bootstrap_enrichment_analysis8_11.rds saved")
print("R job done")

#############################BOOTSTRAP NO CLUSTER#############################

##Libraries + WD
require(dplyr)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(reshape2)
require(Matrix)
require(patchwork)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Load objects
bootstrap_enrichment_analysis <- readRDS('bootstrap_enrichment_analysis.rds')
#bootstrap_enrichment_analysis8_11 <- readRDS('bootstrap_enrichment_analysis8_11.rds')

##Add the extra PS 8-11 analysis (not done because it is redundant)
#bootstrap_enrichment_analysis8_11 <- list(bootstrap_enrichment_analysis8_11)
#bootstrap_enrichment_analysis <- append(bootstrap_enrichment_analysis, bootstrap_enrichment_analysis8_11)

#Create vector of names (could have been done automatically pretty easily I think)
list_names <- c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8","PS9","PS10","PS11",
                "PS12","PS13","PS14","PS15","PS16","PS17","PS18","PS19","PS20","PS21","PS22",
                "PS23","PS24","PS25","PS26","PS27","PS28","PS29","PS30","PS31")
#Assign the names
names(bootstrap_enrichment_analysis) <- list_names

#Eliminate empty lists
bootstrap_enrichment_analysis = bootstrap_enrichment_analysis[-which(sapply(bootstrap_enrichment_analysis, is.null))]

## Plot the results. Dot plots with p-value as color and fold-change as size

#Add the PS to the results tables (as a column equal to the list name)
for (i in names(bootstrap_enrichment_analysis)) {
  bootstrap_enrichment_analysis[[i]][["results"]]$PS <- i
}

#Create a dataframe with all the p-values and fold-changes per cell type per
#phylostratum (First create separate data frames with just what you want, then
#use rbind)
bootstrap_pval_fold <- list()
for (i in names(bootstrap_enrichment_analysis)) {
  bootstrap_pval_fold[[i]] <- select(bootstrap_enrichment_analysis[[i]][["results"]], CellType, p, fold_change, PS)
}
molten_dotplot_df <- as.data.frame(matrix(0, ncol = length(bootstrap_pval_fold[[1]]), nrow = 0))
#First rbind manual
molten_dotplot_df <- rbind(bootstrap_pval_fold[[1]],  bootstrap_pval_fold[[2]])
#Rest in a for loop
for (i in 3:length(bootstrap_pval_fold)) {
  molten_dotplot_df <- rbind(molten_dotplot_df,  bootstrap_pval_fold[[i]])
}
#Adjust the p-values. This may be tricky because they're not independent, they
#come "by group"
molten_dotplot_df$p_adj <- p.adjust(molten_dotplot_df$p,
                                    method = "fdr")

#Apply log2(fold_change)
molten_dotplot_df$log2_fold_change <- log(molten_dotplot_df$fold_change, base = 2)
#Add a new column with categorical data (colors) depending on the p-values. Uses dplyr.
#Green for p_adj OK
#Blue for p_adj not OK but p OK
#Red for p not OK
molten_dotplot_df <- molten_dotplot_df %>% mutate(Color =
                                                    case_when(p_adj <= 0.05 ~ "p_adj<=0.05", 
                                                              p_adj > 0.05 & p <= 0.05  ~ "p_adj>0.05 & p<=0.05",))

#Set the order of the cell types and PS
order_ct <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
              "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
              "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")
order_ps <- c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8","PS10","PS11",
              "PS13","PS14","PS16","PS17","PS18","PS19","PS20","PS21",
              "PS25","PS26","PS31")

#Actual dot plot
p1 <- ggplot(molten_dotplot_df[!is.na(molten_dotplot_df$Color),], aes(x = factor(CellType, levels = order_ct),
                                                                      y = factor(PS, levels = list_names))) + 
  geom_point(aes(size = log2_fold_change, color = Color)) + theme_light() + 
  scale_color_manual(values = c("gold","#c1c1c1")) + scale_y_discrete(drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_line(color = "#201b1b", size = 0.2, linetype = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")) +
  labs(color = "Significance", size = "log2(Fold Change)")
p1
#Save the plot
ggsave(p1, height = 6, width = 6, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/ewce_dotplot.pdf",
       dpi = 400)


##Bar plot to the left of it to show how many genes each PS has.

#Create the necessary object to know this. It is a bit of a mess so I will just
#do this here and save the object (it must be done through the workflow of 
#AUTOMATED FOR ALL CELL TYPES AND PS until you reach individual_ps_per_cell_type)
#genes_per_ps <- c()
#for (i in 1:length(individual_ps_per_cell_type[["CGE IN"]])){
#  genes_per_ps <- c(genes_per_ps, length(rownames(individual_ps_per_cell_type[["CGE IN"]][[i]])))
#}
#saveRDS(genes_per_ps, file='genes_per_ps.rds')

#Load the objects
genes_per_ps <- readRDS('ngenes_per_ps.rds')
ps_names <- c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8","PS9","PS10","PS11",
              "PS12","PS13","PS14","PS15","PS16","PS17","PS18","PS19","PS20","PS21","PS22",
              "PS23","PS24","PS25","PS26","PS27","PS28","PS29","PS30","PS31")
order <- ps_names

#Create the data frame for the barplot
genes_per_ps <- genes_per_ps %>% mutate(Color =
                                      case_when(genes_per_ps$N_GENES >= 5 ~ "OK", 
                                                genes_per_ps$N_GENES < 5  ~ "not_OK",))
genes_per_ps$log <- log2(genes_per_ps$N_GENES+1)
##Actual bar plot
p2 <- ggplot(genes_per_ps, aes(x=factor(PS, levels=order), y=log, fill=Color)) + 
  geom_bar(position = position_dodge(), stat="identity", width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("OK" = "#37ccff", "not_OK" = "tomato")) +
  coord_flip() + scale_y_reverse() + scale_x_discrete(name = "", position = "top") +
  geom_text(aes(label=N_GENES), hjust = -0.2, size=2.9)+theme_light() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"))
p2
#Save the plot
ggsave(p2, height = 6, width = 3, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/ewce_barplot.pdf",
       dpi = 400)

##Combine plots (result not good enough, I suppose I will have to adjust manually
#via Inkscape)
plot_grid(p2, p1, labels = c('A', 'B'))

########################BOOTSTRAP GROUPS PS CLUSTER############################

##Libraries + WD
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Matrix)
require(patchwork)
require(RNOmni, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
require(EWCE, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
setwd('/users/genomics/xmarti/data')
print("Libraries and WD OK")

##Load objects
load('CellTypeData_ewce_all_genes_normed_translated.rda')
ps_table<-readRDS('ps_table.rds')
print("Objects OK")

##Prepare necessary objects

#Background set
background <- unique(rownames(ctd[[1]][["specificity"]]))
#Create PS groups (done manually, I think there is no other way)
n_groups <- c(1:5)
ps1_7 <- subset(ps_table, PS < 8)
ps8_11 <- subset(ps_table, PS > 7 & PS < 12)
ps12_17 <- subset(ps_table, PS > 11 & PS < 18)
ps18_22 <- subset(ps_table, PS > 17 & PS < 23)
ps23_31 <- subset(ps_table, PS > 22)
#Make a list to iterate
ps_groups <- list(ps1_7 = ps1_7,
                  ps8_11 = ps8_11,
                  ps12_17 = ps12_17,
                  ps18_22 = ps18_22,
                  ps23_31 = ps23_31)
#Create empty list to store results
bootstrap_groups <- list()
#Actual job
for (i in names(ps_groups)) {
  #Hit set
  hits <- ps_groups[[i]]$Name[match(background, ps_groups[[i]]$Name)]
  hits <- hits[!is.na(hits) & hits != ""]

  #Reps and annotLevel
  reps = 100000
  level = 1
  
  ##Actual job
  full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                           reps=reps,annotLevel=level, 
                                           genelistSpecies = "human",
                                           sctSpecies = "human")
  
  ##Save object into the list
  bootstrap_groups[[i]] <- full_results
}

##Save results
saveRDS(bootstrap_groups, file ='bootstrap_groups.rds')
print("bootstrap_groups.rds saved")
print("R job done")

########################BOOTSTRAP GROUPS NO CLUSTER############################

##Libraries + WD
require(dplyr)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(reshape2)
require(Matrix)
require(patchwork)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Load objects
bootstrap_groups <- readRDS('bootstrap_groups.rds')

## Plot the results. Dot plots with p-value as color and fold-change as size

#Add the PS to the results tables (as a column equal to the list name, done manually)
bootstrap_groups[[1]][["results"]]$PS <- "PS1-7"
bootstrap_groups[[2]][["results"]]$PS <- "PS8-11"
bootstrap_groups[[3]][["results"]]$PS <- "PS12-17"
bootstrap_groups[[4]][["results"]]$PS <- "PS18-22"
bootstrap_groups[[5]][["results"]]$PS <- "PS23-31"


#Create a dataframe with all the p-values and fold-changes per cell type per
#phylostratum (First create separate data frames with just what you want, then
#use rbind)
bootstrap_pval_fold <- list()
for (i in names(bootstrap_groups)) {
  bootstrap_pval_fold[[i]] <- select(bootstrap_groups[[i]][["results"]], CellType, p, fold_change, PS)
}
molten_dotplot_df <- as.data.frame(matrix(0, ncol = length(bootstrap_pval_fold[[1]]), nrow = 0))
#First rbind manual
molten_dotplot_df <- rbind(bootstrap_pval_fold[[1]],  bootstrap_pval_fold[[2]])
#Rest in a for loop
for (i in 3:length(bootstrap_pval_fold)) {
  molten_dotplot_df <- rbind(molten_dotplot_df,  bootstrap_pval_fold[[i]])
}
#Adjust the p-values. This may be tricky because they're not independent, they
#come "by group"
molten_dotplot_df$p_adj <- p.adjust(molten_dotplot_df$p,
                                    method = "fdr")

#Apply log2(fold_change)
molten_dotplot_df$log2_fold_change <- log(molten_dotplot_df$fold_change, base = 2)
#Add a new column with categorical data (colors) depending on the p-values. Uses dplyr.
#Green for p_adj OK
#Blue for p_adj not OK but p OK
#Red for p not OK
molten_dotplot_df <- molten_dotplot_df %>% mutate(Color =
                                                    case_when(p_adj <= 0.05 ~ "p_adj<=0.05", 
                                                              p_adj > 0.05 & p <= 0.05  ~ "p_adj>0.05 & p<=0.05",))

#Set the order of the cell types and PS
order_ct <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
              "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
              "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")
order_ps <- c("PS1-7","PS8-11","PS12-17","PS18-22","PS23-31")

#Actual dot plot
p1 <- ggplot(molten_dotplot_df[!is.na(molten_dotplot_df$Color),], aes(x = factor(CellType, levels = order_ct),
                                                                      y = factor(PS, levels = order_ps))) + 
  geom_point(aes(size = log2_fold_change, color = Color)) + theme_light() + 
  scale_color_manual(values = c("gold","#c1c1c1")) + scale_y_discrete(drop=FALSE) +
  scale_x_discrete(drop=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_line(color = "grey", size = 0.2, linetype = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")) +
  labs(color = "Significance", size = "log2(Fold Change)")
p1
#Save the plot
ggsave(p1, height = 2, width = 7, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/ewce_dotplot_groups.pdf",
       dpi = 400)


##Bar plot to the left of it to show how many genes each PS has.

#Load the objects
genes_per_ps <- readRDS('ngenes_per_ps.rds')

#Adjust the genes_per_ps
genes_per_ps1_7 <- sum(genes_per_ps$N_GENES[1:7])
genes_per_ps8_11 <- sum(genes_per_ps$N_GENES[8:11])
genes_per_ps12_17 <- sum(genes_per_ps$N_GENES[12:17])
genes_per_ps18_22 <- sum(genes_per_ps$N_GENES[18:22])
genes_per_ps23_31 <- sum(genes_per_ps$N_GENES[23:31])
genes_per_ps <- c(genes_per_ps1_7,genes_per_ps8_11,genes_per_ps12_17,
                  genes_per_ps18_22,genes_per_ps23_31)

ps_names <- c("PS1-7","PS8-11","PS12-17","PS18-22","PS23-31")
order <- ps_names

#Create the data frame for the barplot
barplot_df <- data.frame(genes_per_ps,ps_names)
barplot_df <- barplot_df %>% mutate(Color =
                                      case_when(genes_per_ps >= 5 ~ "OK", 
                                                genes_per_ps < 5  ~ "not_OK",))
barplot_df$log <- log2(barplot_df$genes_per_ps+1)
##Actual bar plot
p2 <- ggplot(barplot_df, aes(x=factor(ps_names, levels=order), y=log, fill=Color)) + 
  geom_bar(position = position_dodge(), stat="identity", width = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("OK" = "#37ccff", "not_OK" = "tomato")) +
  coord_flip() + scale_y_reverse() + scale_x_discrete(name = "", position = "top") +
  geom_text(aes(label=genes_per_ps), hjust = -0.2, size = 2.5) + theme_light() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"))
p2
#Save the plot
ggsave(p2, height = 2, width = 3, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/ewce_barplot_groups.pdf",
       dpi = 400)

######################BOOTSTRAP CELL-TYPE GROUPS CLUSTER#######################
##Libraries + wd
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Matrix)
require(patchwork)
require(RNOmni, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
require(EWCE, lib.loc = '/users/genomics/xmarti/packages/R/4.1')
setwd('/users/genomics/xmarti/data')
print("Libraries and WD OK")

##Load objects
corrected_counts_normed <- readRDS('corrected_counts_normed.rds')
annotation_df<-readRDS('annotation_df.rds')
ps_table<-readRDS('ps_table.rds')
tr_table<-readRDS('gse162170_tr_table.rds')
print("Objects OK")

##Prepare necessary objects

#Translating the rownames of corrected_counts_normed matrix
ids <- rownames(corrected_counts_normed)
gene_names <- tr_table$Symbol[match(ids, tr_table$ID)]
rownames(corrected_counts_normed) <- gene_names

#Adding clustered cell types to the annotation dataframe (excitatory neurons and
#interneurons separated)
annotation_df <- annotation_df %>% mutate("Clustered_CT" =
                                            case_when(Name == "Cyc. Prog." ~ "Prog.",
                                                      Name == "Early RG" ~ "Prog.",
                                                      Name == "Late RG" ~ "Prog.",
                                                      Name == "tRG" ~ "Prog.",
                                                      Name == "nIPC" ~ "Prog.",
                                                      Name == "SP" ~ "EN",
                                                      Name == "GluN1" ~ "EN",
                                                      Name == "GluN2" ~ "EN",
                                                      Name == "GluN3" ~ "EN",
                                                      Name == "GluN4" ~ "EN",
                                                      Name == "GluN5" ~ "EN",
                                                      Name == "GluN6" ~ "EN",
                                                      Name == "GluN7" ~ "EN",
                                                      Name == "GluN8" ~ "EN",
                                                      Name == "MGE IN" ~ "IN",
                                                      Name == "CGE IN" ~ "IN",
                                                      Name == "mGPC" ~ "mGPC",
                                                      Name == "OPC/Oligo" ~ "OPC/Oligo",
                                                      Name == "MG" ~ "MG",
                                                      Name == "Peric." ~ "Peric.",
                                                      Name == "EC" ~ "EC",
                                                      Name == "RBC" ~ "RBC",
                                                      Name == "VLMC" ~ "VLMC"))

#Adding clustered cell types to the annotation dataframe (excitatory neurons and
#interneurons together)
annotation_df <- annotation_df %>% mutate("Clustered_CT_2" =
                                            case_when(Name == "Cyc. Prog." ~ "Prog.",
                                                      Name == "Early RG" ~ "Prog.",
                                                      Name == "Late RG" ~ "Prog.",
                                                      Name == "tRG" ~ "Prog.",
                                                      Name == "nIPC" ~ "Prog.",
                                                      Name == "SP" ~ "Neuron",
                                                      Name == "GluN1" ~ "Neuron",
                                                      Name == "GluN2" ~ "Neuron",
                                                      Name == "GluN3" ~ "Neuron",
                                                      Name == "GluN4" ~ "Neuron",
                                                      Name == "GluN5" ~ "Neuron",
                                                      Name == "GluN6" ~ "Neuron",
                                                      Name == "GluN7" ~ "Neuron",
                                                      Name == "GluN8" ~ "Neuron",
                                                      Name == "MGE IN" ~ "Neuron",
                                                      Name == "CGE IN" ~ "Neuron",
                                                      Name == "mGPC" ~ "mGPC",
                                                      Name == "OPC/Oligo" ~ "OPC/Oligo",
                                                      Name == "MG" ~ "MG",
                                                      Name == "Peric." ~ "Peric.",
                                                      Name == "EC" ~ "EC",
                                                      Name == "RBC" ~ "RBC",
                                                      Name == "VLMC" ~ "VLMC"))

#I need to repeat the generate_celltype_data step
annotLevels = list(cell_type=annotation_df$Clustered_CT_2)
spec_matrix_normed <- generate_celltype_data(exp=corrected_counts_normed, 
                                             annotLevels=annotLevels,
                                             groupName="ewce_all_genes_normed_translated_clustered_ct", 
                                             savePath='/users/genomics/xmarti/data')

#Load the EWCE data now
load('CellTypeData_ewce_all_genes_normed_translated_clustered_ct.rda')

#Background set
background <- unique(rownames(ctd[[1]][["specificity"]]))
ps_with_more_than_4_genes <- c(1,2,3,4,5,6,7,8,10,11,13,14,16,17,18,19,20,21,22,
                               23,25,26,31)
#Create list to store results
bootstrap_enrichment_analysis_clustered_ct <- list()
for (i in ps_with_more_than_4_genes) {
  
  #Hit set
  ps <- subset(ps_table, PS == i)
  hits <- ps$Name[match(background, ps$Name)]
  hits <- hits[!is.na(hits) & hits != ""]
  
  #Reps and annotLevel
  reps = 100000
  level = 1 #This doesn't change because I always introduce a list with just 1 level 
  
  ##Actual job
  full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                           reps=reps,annotLevel=level, 
                                           genelistSpecies = "human",
                                           sctSpecies = "human")
  
  ##Save object into the list
  bootstrap_enrichment_analysis_clustered_ct[[i]] <- full_results
}

##Save objects of interest
saveRDS(bootstrap_enrichment_analysis_clustered_ct, 
        file ='bootstrap_enrichment_analysis_clustered_ct_2.rds')
print("bootstrap_enrichment_analysis_clustered_ct_2.rds saved")
print("R job done")

####################BOOTSTRAP CELL-TYPE GROUPS NO CLUSTER######################
##Libraries + WD
require(dplyr)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(reshape2)
require(Matrix)
require(patchwork)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Load objects
bootstrap_enrichment_analysis <- readRDS('bootstrap_enrichment_analysis_clustered_ct_2.rds')

#Create vector of names (could have been done automatically pretty easily I think)
list_names <- c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8","PS9","PS10","PS11",
                "PS12","PS13","PS14","PS15","PS16","PS17","PS18","PS19","PS20","PS21","PS22",
                "PS23","PS24","PS25","PS26","PS27","PS28","PS29","PS30","PS31")
#Assign the names
names(bootstrap_enrichment_analysis) <- list_names

#Eliminate empty lists
bootstrap_enrichment_analysis = bootstrap_enrichment_analysis[-which(sapply(bootstrap_enrichment_analysis, is.null))]

## Plot the results. Dot plots with p-value as color and fold-change as size

#Add the PS to the results tables (as a column equal to the list name)
for (i in names(bootstrap_enrichment_analysis)) {
  bootstrap_enrichment_analysis[[i]][["results"]]$PS <- i
}

#Create a dataframe with all the p-values and fold-changes per cell type per
#phylostratum (First create separate data frames with just what you want, then
#use rbind)
bootstrap_pval_fold <- list()
for (i in names(bootstrap_enrichment_analysis)) {
  bootstrap_pval_fold[[i]] <- select(bootstrap_enrichment_analysis[[i]][["results"]], CellType, p, fold_change, PS)
}
molten_dotplot_df <- as.data.frame(matrix(0, ncol = length(bootstrap_pval_fold[[1]]), nrow = 0))
#First rbind manual
molten_dotplot_df <- rbind(bootstrap_pval_fold[[1]],  bootstrap_pval_fold[[2]])
#Rest in a for loop
for (i in 3:length(bootstrap_pval_fold)) {
  molten_dotplot_df <- rbind(molten_dotplot_df,  bootstrap_pval_fold[[i]])
}
#Adjust the p-values.
molten_dotplot_df$p_adj <- p.adjust(molten_dotplot_df$p,
                                    method = "fdr")

#Apply log2(fold_change)
molten_dotplot_df$log2_fold_change <- log(molten_dotplot_df$fold_change, base = 2)
#Add a new column with categorical data (colors) depending on the p-values. Uses dplyr.
#Green for p_adj OK
#Blue for p_adj not OK but p OK
#Red for p not OK
molten_dotplot_df <- molten_dotplot_df %>% mutate(Color =
                                                    case_when(p_adj <= 0.05 ~ "p_adj<=0.05", 
                                                              p_adj > 0.05 & p <= 0.05  ~ "p_adj>0.05 & p<=0.05",))

#Change "Neuron" name for "EN" (excitatory neuron)
#molten_dotplot_df["CellType"][molten_dotplot_df["CellType"] == "Neuron"] <- "EN"
#Set the order of the cell types and PS
order_ct <- c("Prog.", "EN", "IN",  "mGPC", "OPC/Oligo", "MG","Peric.", "EC",
              "RBC", "VLMC")
order_ct <- c("Prog.", "Neuron",  "mGPC", "OPC/Oligo", "MG","Peric.", "EC",
              "RBC", "VLMC")
order_ps <- c("PS1","PS2","PS3","PS4","PS5","PS6","PS7","PS8","PS10","PS11",
              "PS13","PS14","PS16","PS17","PS18","PS19","PS20","PS21",
              "PS25","PS26","PS31")

#Actual dot plot
p1 <- ggplot(molten_dotplot_df[!is.na(molten_dotplot_df$Color),], aes(x = factor(CellType, levels = order_ct),
                                                                      y = factor(PS, levels = list_names))) + 
  geom_point(aes(size = log2_fold_change, color = Color)) + theme_light() + 
  scale_color_manual(values = c("gold","#c1c1c1")) + scale_y_discrete(drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_line(color = "#201b1b", size = 0.2, linetype = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")) +
  labs(color = "Significance", size = "log2(Fold Change)")
p1
#Save the plot
ggsave(p1, height = 6, width = 5, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/ewce_dotplot_clustered_ct_2.pdf",
       dpi = 700)


##Bar plot to the left of it to show how many genes each PS has. This is already
##done in the BOOTSTRAP NO CLUSTER section (the bar plot is exactly the same)

##############################MG MARKER GENES##################################

##Libraries + WD
require(dplyr)
require(patchwork)
#require(Seurat)
require(ggplot2)
#require(readxl)
require(biomaRt)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Libraries + WD (for cluster)
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Seurat)
require(patchwork)
setwd('/users/genomics/xmarti/data')

##Load objects
#dhcc <- readRDS('dhcc_v4.rds')
mg_markers <- readRDS("mg_markers.rds")
#tr_table <- readRDS("tr_table_gencode_v27.rds")
mg_spec_normed_with_ps <- readRDS("mg_spec_normed_with_ps.rds")
#mg_not_macro <- readRDS("mg_not_macro.rds")
human_genes_mg_not_macro <- readRDS("human_genes_mg_not_macro.rds")

##Find MG markers
#mg_markers <- FindMarkers(dhcc, ident.1 = "MG", min.pct = 0.25)

##Translate genes of mg_markers (already done and saved)
#ensembl_ids <- rownames(mg_markers)
#gene_symbols <- tr_table$Gene.name[match(ensembl_ids, tr_table$`Gene.stable.ID`)]
#rownames(mg_markers) <- gene_symbols
#mg_markers$gene_symbol <- gene_symbols

##Cross the dataframe of specificty and PS with the one with Markers (not very
##cleanly done)
#mg_markers$gene_symbol <- rownames(mg_markers) (done and saved)
mg_markers_spec_ps <- left_join(mg_markers,mg_spec_normed_with_ps, 
                                by=c("gene_symbol" = "Name"))
#Remove rows with at least 1 NA value
mg_markers_spec_ps <- na.omit(mg_markers_spec_ps)

##Violin/box plot of distribution of logFC vs. chunks of specificity of 0.10.

#Add columns with chunks of specificity (dplyr with case)
mg_markers_spec_ps <- mg_markers_spec_ps %>% 
  mutate(sp_chunks =
           case_when(MG <= 0.10 ~ "0-0.10", 
                     MG > 0.10  & MG <= 0.20 ~ "0.10-0.20",
                     MG > 0.20  & MG <= 0.30 ~ "0.20-0.30",
                     MG > 0.30  & MG <= 0.40 ~ "0.30-0.40",
                     MG > 0.40  & MG <= 0.50 ~ "0.40-0.50",
                     MG > 0.50  & MG <= 0.60 ~ "0.50-0.60",
                     MG > 0.60  & MG <= 0.70 ~ "0.60-0.70",
                     MG > 0.70  & MG <= 0.80 ~ "0.70-0.80",
                     MG > 0.80  & MG <= 0.90 ~ "0.80-0.90",
                     MG > 0.90 ~ "0.90-1"))

#Make sure that new column is a factor (this step is recommended in a tutorial)
mg_markers_spec_ps$sp_chunks <- as.factor(mg_markers_spec_ps$sp_chunks)

#Helper function to get the number (counts) of each group
get_box_stats <- function(y, upper_limit = quantile(y,0.95)*1.1) {
  return(data.frame(y = 0.95 * upper_limit, label = length(y)))
}
#Actual violin + box plot
violin_plot <- ggplot(mg_markers_spec_ps, aes(x=sp_chunks, y=avg_log2FC, fill=sp_chunks), show.legend = FALSE) + 
  geom_violin(trim = FALSE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() + 
  labs(x = "Specificity", y = "Average log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = -0.2, vjust = 1)
violin_plot
#Reasonable plot, but I would like the colors to be displayed according to the 
#number of genes for each chunk, with a gradient a like the rainbow: red for very
#few, blue/violet for a lot. Or just add the number of observations, that also
#works.

#Save the plot
ggsave(violin_plot, height = 6, width = 8, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/fc_sp_violin.pdf",
       dpi = 400)

##Loading the S4 table from Geirsdottir et al, with genes specific for MG but not
##for macrophages
#mg_not_macro <- read_excel('/home/onecellardoor/Dropbox/imim/data/geirsdottir_david_2019/s4.xlsx')
#saveRDS(mg_not_macro,'mg_not_macro.rds')

##Translating the mouse symbols of mm10/ensemble v95 into human gencode v27.

#Gabriel has a function in a script of him, maybe it works (or at least it should
#be adaptable). It doesn't work, and I really don't know how to do this one.
convertMouseGeneList <- function(x) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = x , mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,])
  
  #Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
#human_genes_mg_not_macro <- convertMouseGeneList(mg_not_macro$`Symbol mouse`)

##Cross the human genes of MG and not Macrophages vs the MG Markers and the PS
##table.

#vs MG Markers
true_mg_genes <- left_join(mg_markers,human_genes_mg_not_macro,
                           by=c("gene_symbol" = "HGNC.symbol"))
true_mg_genes <- na.omit(true_mg_genes) #5 (almost 6) good genes!

true_mg_genes_with_ps <- left_join(mg_markers_spec_ps,human_genes_mg_not_macro,
                                   by=c("gene_symbol" = "HGNC.symbol"))
true_mg_genes_with_ps <- na.omit(true_mg_genes_with_ps) #All have PS, but below
#PS8 :/ Maybe they adopt new functions in different gene networks. 

##Save objects
#saveRDS(mg_markers, file = "mg_markers.rds")
#saveRDS(human_genes_mg_not_macro, file = "human_genes_mg_not_macro.rds")

########################GENES VS CELL TYPE DOT PLOTS###########################

##Libraries + WD
require(dplyr)
require(patchwork)
require(Seurat)
require(ggplot2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Libraries + WD (for cluster)
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Seurat)
require(patchwork)
require(ggplot2)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD loaded correctly")

##Load objects
tr_table <- readRDS('gse162170_tr_table.rds')
mg_spec_normed_with_ps <- readRDS("mg_spec_normed_with_ps.rds")
print("Loading the big fella, this takes a while:")
dhcc <- readRDS('dhcc_v4.rds')
print("Objects loaded correctly")

##Create lists with found genes for every PS
#Create empty list for iteration
gene_names_per_ps <- list()
#Fill the list with dataframes of just 1 PS
for (i in 1:max(mg_spec_normed_with_ps$PS)) {
  gene_names_per_ps[[i]] <- mg_spec_normed_with_ps[mg_spec_normed_with_ps$PS == i 
                                                   & mg_spec_normed_with_ps$MG >= 0.1,]
}
#Assign names for easier titles and stuff later on
names(gene_names_per_ps) <- c(1:max(mg_spec_normed_with_ps$PS))
#Eliminate empty list objects
gene_names_per_ps <- gene_names_per_ps[sapply(gene_names_per_ps, function(x) dim(x)[1]) > 0]
#Translate the names into IDs for correct matching in the SeuratObject
tr_gene_names_per_ps <- gene_names_per_ps
for (i in 1:length(gene_names_per_ps)) {
  tr_gene_names_per_ps[[i]]$Name <- tr_table$ID[match(tr_gene_names_per_ps[[i]]$Name, tr_table$Symbol)]
}

##Order cell types in this order: 
order <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
           "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
           "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")

##Plotting loop
for (i in 8:length(tr_gene_names_per_ps)) {
  #Store the gene names as a character
  id_list <- as.character(tr_gene_names_per_ps[[i]]$Name)
  #Store the DotPlot fuction data into an object for further manipulation
  preplot <- DotPlot(dhcc, features = id_list, cluster.idents = TRUE)
  #Retranslate the IDs into symbols and store them as a factor (necessary)
  gene_list <- as.factor(tr_table$Symbol[match((preplot$data$features.plot),
                                                  tr_table$ID)])
  #Actual plot
  my_plot <- ggplot(preplot$data, aes(x=factor(id, levels = order), y=gene_list)) + 
    geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
    theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_color_gradient(low="cadetblue", high = "tomato") +  
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0,0.2),"cm")) +
    labs(color="Avg. scaled exp.", size="Exp. %") + 
    ggtitle(paste0("PS",names(tr_gene_names_per_ps[i])))
  #Save objects
  ggsave(my_plot, filename = paste0("PS",names(tr_gene_names_per_ps[i]),".pdf"))
}
print("Done!")

##Investigate the IDs not found because of a "poor translation" (already solved)
#Manually annotate the IDs that returned an error
#lost_IDs <- c("ENSG00000274145", "ENSG00000206485", "ENSG00000278718", "ENSG00000278186", 
#              "ENSG00000276982","ENSG00000278634", "ENSG00000277398",
#              "ENSG00000206478", "ENSG00000275539", "ENSG00000274914", 
#              "ENSG00000277092","ENSG00000276042", "ENSG00000276858",
#              "ENSG00000223833", "ENSG00000274110", "ENSG00000236697",
#              "ENSG00000223465", "ENSG00000274454")
#Translate into gene symbol
#lost_genes <- tr_table$Gene.name[match(lost_IDs, tr_table$Gene.stable.ID)]
#Find their PS
#lost_PS <- mg_spec_normed_with_ps$PS[match(lost_genes, mg_spec_normed_with_ps$Name)]                                 

##Gabriel suggested the idea that I just retranslated the rownames of the SeuratObject
##and stored these IDs with the symbols together for a perfect translation afterwards
#ID <- rownames(dhcc)
#Symbol <- tr_table$Gene.name[match(gse162170_ensembl_id, tr_table$Gene.stable.ID)]
#gse162170_tr_table <- data.frame(ID, Symbol)
#saveRDS(gse162170_tr_table, file = "gse162170_tr_table.rds")

##Split the PS11 into 3 (otherwise it's a mess) and plot the parts
splits <- list(tr_gene_names_per_ps[["11"]][1:43,]$Name,tr_gene_names_per_ps[["11"]][44:86,]$Name,
            tr_gene_names_per_ps[["11"]][87:130,]$Name)
counter<-0
for (i in splits) {
  counter <- counter + 1
  id_list <- as.character(i)
  #Store the DotPlot fuction data into an object for further manipulation
  preplot <- DotPlot(dhcc, features = id_list, cluster.idents = TRUE)
  #Retranslate the IDs into symbols and store them as a factor (necessary)
  gene_list <- as.factor(tr_table$Symbol[match((preplot$data$features.plot),
                                               tr_table$ID)])
  #Actual plot
  my_plot <- ggplot(preplot$data, aes(x=factor(id, levels = order), y=gene_list)) + 
    geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
    theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_color_gradient(low="cadetblue", high = "tomato") +  
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0,0.2),"cm")) +
    labs(color="Avg. scaled exp.", size="Exp. %") + 
    ggtitle(paste0("PS11_",counter))
  #Save objects
  ggsave(my_plot, filename = paste0("PS11_",counter,".pdf"))
}

############################TRUE MG SPECIFIC GENES#############################
##Libraries + WD
require(dplyr)
require(patchwork)
require(Seurat)
require(ggplot2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE162170')

##Libraries + WD (for cluster)
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Seurat)
require(patchwork)
require(ggplot2)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD loaded correctly")

##Load objects
spec_matrix_with_ps <- readRDS('spec_matrix_with_ps.rds')
tr_table <- readRDS('gse162170_tr_table.rds')
mg_markers <- readRDS('mg_markers.rds')
print("Loading the big fella, this takes a while:")
dhcc <- readRDS('dhcc_v4.rds')
print("Objects loaded correctly")

##Select genes with average log2FC >= 1 from the mg_markers obtained with
##FindMarkers function of Seurat

selected_mg_markers <- mg_markers[mg_markers$avg_log2FC>=1,]

##Filter by specificity. MG > 0.1, rest < 0.05. Returns weird stuff if only filtered
#with this method (see true_mg_spec_genes_dotplot.pdf). HARD FILTERING (SPEC)
true_mg_spec_genes <- spec_matrix_with_ps[spec_matrix_with_ps$MG >= 0.1
                                              & spec_matrix_with_ps$`CGE IN` <= 0.05
                                              & spec_matrix_with_ps$`Cyc. Prog.` <= 0.05
                                              & spec_matrix_with_ps$`Early RG` <= 0.05
                                              & spec_matrix_with_ps$EC <= 0.05
                                              & spec_matrix_with_ps$GluN1 <= 0.05
                                              & spec_matrix_with_ps$GluN2 <= 0.05
                                              & spec_matrix_with_ps$GluN3 <= 0.05
                                              & spec_matrix_with_ps$GluN4 <= 0.05
                                              & spec_matrix_with_ps$GluN5 <= 0.05
                                              & spec_matrix_with_ps$GluN6 <= 0.05
                                              & spec_matrix_with_ps$GluN7 <= 0.05
                                              & spec_matrix_with_ps$GluN8 <= 0.05
                                              & spec_matrix_with_ps$`Late RG` <= 0.05
                                              & spec_matrix_with_ps$`MGE IN` <= 0.05
                                              & spec_matrix_with_ps$mGPC <= 0.05
                                              & spec_matrix_with_ps$nIPC <= 0.05
                                              & spec_matrix_with_ps$`OPC/Oligo` <= 0.05
                                              & spec_matrix_with_ps$Peric. <= 0.05
                                              & spec_matrix_with_ps$RBC <= 0.05
                                              & spec_matrix_with_ps$SP <= 0.05
                                              & spec_matrix_with_ps$tRG <= 0.05
                                              & spec_matrix_with_ps$VLMC <= 0.05
                                              & spec_matrix_with_ps$PS > 7,]

true_mg_spec_genes1_7 <- spec_matrix_with_ps[spec_matrix_with_ps$MG >= 0.1
                                          & spec_matrix_with_ps$`CGE IN` <= 0.05
                                          & spec_matrix_with_ps$`Cyc. Prog.` <= 0.05
                                          & spec_matrix_with_ps$`Early RG` <= 0.05
                                          & spec_matrix_with_ps$EC <= 0.05
                                          & spec_matrix_with_ps$GluN1 <= 0.05
                                          & spec_matrix_with_ps$GluN2 <= 0.05
                                          & spec_matrix_with_ps$GluN3 <= 0.05
                                          & spec_matrix_with_ps$GluN4 <= 0.05
                                          & spec_matrix_with_ps$GluN5 <= 0.05
                                          & spec_matrix_with_ps$GluN6 <= 0.05
                                          & spec_matrix_with_ps$GluN7 <= 0.05
                                          & spec_matrix_with_ps$GluN8 <= 0.05
                                          & spec_matrix_with_ps$`Late RG` <= 0.05
                                          & spec_matrix_with_ps$`MGE IN` <= 0.05
                                          & spec_matrix_with_ps$mGPC <= 0.05
                                          & spec_matrix_with_ps$nIPC <= 0.05
                                          & spec_matrix_with_ps$`OPC/Oligo` <= 0.05
                                          & spec_matrix_with_ps$Peric. <= 0.05
                                          & spec_matrix_with_ps$RBC <= 0.05
                                          & spec_matrix_with_ps$SP <= 0.05
                                          & spec_matrix_with_ps$tRG <= 0.05
                                          & spec_matrix_with_ps$VLMC <= 0.05
                                          & spec_matrix_with_ps$PS <= 7,]

##Filter by spec, but just the one of MG. SOFT FILTERING (SPEC)
#true_mg_spec_genes <- spec_matrix_with_ps[spec_matrix_with_ps$MG >= 0.1
#                                     & spec_matrix_with_ps$PS > 7,]
#true_mg_spec_genes1_7 <- spec_matrix_with_ps[spec_matrix_with_ps$MG >= 0.1
#                                     & spec_matrix_with_ps$PS <= 7,]

####Combine mg markers with average log2FC >= 1 to the ones with spec>10 and other
## cell types < 0.05
true_mg_spec_genes <- left_join(true_mg_spec_genes,selected_mg_markers,
                                   by=c("Gene_Name"="gene_symbol"))
true_mg_spec_genes1_7 <- left_join(true_mg_spec_genes1_7,selected_mg_markers,
                                by=c("Gene_Name"="gene_symbol"))
true_mg_spec_genes <- true_mg_spec_genes[!is.na(true_mg_spec_genes$avg_log2FC),]
true_mg_spec_genes1_7 <- true_mg_spec_genes1_7[!is.na(true_mg_spec_genes1_7$avg_log2FC),]

##Sort the selected genes by alphabetic order and then PS (in this order).
##This is to group them into PS and an intuitive alphabetical order in the plot
true_mg_spec_genes <- true_mg_spec_genes[order(true_mg_spec_genes$Gene_Name, decreasing = TRUE),]
true_mg_spec_genes1_7 <- true_mg_spec_genes1_7[order(true_mg_spec_genes1_7$Gene_Name, decreasing = TRUE),]
true_mg_spec_genes <-true_mg_spec_genes[order(true_mg_spec_genes$PS),]
true_mg_spec_genes1_7 <-true_mg_spec_genes1_7[order(true_mg_spec_genes1_7$PS),]

###PS8-31

##Translate the symbols into IDs and store it as a character (just in case)
id_list <- as.character(tr_table$ID[match(true_mg_spec_genes$Gene_Name, tr_table$Symbol)])

##Store the DotPlot fuction data into an object for further manipulation
preplot <- DotPlot(dhcc, features = id_list, cluster.idents = TRUE)

##Add the PS as a new column of preplot$data

#The Genes come in ID format here so I have to translate them into symbols first
preplot$data$features.plot <- as.factor(tr_table$Symbol[match(preplot$data$features.plot,
                                             tr_table$ID)])
#Add the PS column
preplot$data$PS <- spec_matrix_with_ps$PS[match(preplot$data$features.plot,
                                                spec_matrix_with_ps$Gene_Name)]

#Create vector with cell types that are not MG
cell_types_not_mg <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
                       "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
                       "CGE IN", "mGPC", "OPC/Oligo", "Peric.", "EC", "RBC", "VLMC")

#Worst nested for loop in history. Navigates the molten data frame of the ggplot
#via ID and cell type and filters genes not very specific to MG.
for (gene in unique(preplot$data$features.plot)) {
  #Creates a counter to know if the gene is worthy or not (must pass all the tests)
  counter <- 0
  for (cell_type in cell_types_not_mg) {
    #Is the MG avg. expression 2x higher than that of all the others?
    difference <- preplot$data$avg.exp.scaled[preplot$data$features.plot==gene &
                                                preplot$data$id=="MG"] -
      preplot$data$avg.exp.scaled[preplot$data$features.plot==gene&preplot$data$id==cell_type]
    #Is the MG % of cell expression 3x higher than that of all the others?
    pct <- preplot$data$pct.exp[preplot$data$features.plot==gene &
                                         preplot$data$id=="MG"] /
      preplot$data$pct.exp[preplot$data$features.plot==gene&preplot$data$id==cell_type]
    #If both criteria are met, plus min. % exp of 20, add 1 to the counter
    if ((difference>2) & (pct>3) & 
      (preplot$data$pct.exp[preplot$data$features.plot==gene & preplot$data$id=="MG"] >= 20)) {
      counter <- counter + 1
    }
  }
  #If the criteria are met for all the cell types, the gene is worthy. Otherwise
  #all the rows with that gene are deleted.
  if (counter<length(cell_types_not_mg)) {
    preplot$data<-preplot$data[!(preplot$data$features.plot==gene),]
  }
}

#Set the order of the genes. This just defactorizes the features.plot and sets the
#order as it is in the true_mg_genes table, which is already the desired one.
order_genes <- as.character(unique(preplot$data$features.plot))
#Also save this for the heatmap (they share y axis)
saveRDS(order_genes, file = "heatmap_order_genes8_31_hard.rds")

#Set the <0.05 % expr values to NA.
preplot$data$pct.exp[preplot$data$pct.exp < 5] <- NA

##Order cell types in this order: 
order_celltype <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
           "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
           "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")

##Actual plot
my_plot <- ggplot(preplot$data, aes(x=factor(id, levels = order_celltype),
                                    y=factor(features.plot, levels = order_genes))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="cadetblue", high = "tomato") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0,0.2,0.2),"cm"),
        legend.position = "left", axis.text.y=element_text(hjust=-0.5)) + 
  scale_y_discrete(position = "right") + 
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE, switch = "y") +
  labs(color="Avg. scaled exp.", size="Exp. %") + 
  ggtitle("Microglia Specific Genes PS8-31")
my_plot
#Checked it, success this time (with the difficult filtering) and the facet_grid

##Save objects
ggsave(my_plot, height = 6, width = 6,
       filename = "true_mg_avg_pct_exp_genes_dotplot8_31_hard.pdf",
       dpi = 700)

###PS1-7

##Translate the symbols into IDs and store it as a character (just in case)
id_list <- as.character(tr_table$ID[match(true_mg_spec_genes1_7$Gene_Name, tr_table$Symbol)])

##Store the DotPlot fuction data into an object for further manipulation
preplot <- DotPlot(dhcc, features = id_list, cluster.idents = TRUE)

##Add the PS as a new column of preplot$data

#The Genes come in ID format here so I have to translate them into symbols first
preplot$data$features.plot <- as.factor(tr_table$Symbol[match(preplot$data$features.plot,
                                                              tr_table$ID)])
#Add the PS column
preplot$data$PS <- spec_matrix_with_ps$PS[match(preplot$data$features.plot,
                                                spec_matrix_with_ps$Gene_Name)]

#Create vector with cell types that are not MG
cell_types_not_mg <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
                       "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
                       "CGE IN", "mGPC", "OPC/Oligo", "Peric.", "EC", "RBC", "VLMC")

#Worst nested for loop in history. Navigates the molten data frame of the ggplot
#via ID and cell type and filters genes not very specific to MG.
for (gene in unique(preplot$data$features.plot)) {
  #Creates a counter to know if the gene is worthy or not (must pass all the tests)
  counter <- 0
  for (cell_type in cell_types_not_mg) {
    #Is the MG avg. expression 2x higher than that of all the others?
    difference <- preplot$data$avg.exp.scaled[preplot$data$features.plot==gene &
                                                preplot$data$id=="MG"] -
      preplot$data$avg.exp.scaled[preplot$data$features.plot==gene&preplot$data$id==cell_type]
    #Is the MG % of cell expression 3x higher than that of all the others?
    pct <- preplot$data$pct.exp[preplot$data$features.plot==gene &
                                  preplot$data$id=="MG"] /
      preplot$data$pct.exp[preplot$data$features.plot==gene&preplot$data$id==cell_type]
    #If both criteria are met, plus min. % exp of 20, add 1 to the counter
    if ((difference>2) & (pct>3) & 
        (preplot$data$pct.exp[preplot$data$features.plot==gene & preplot$data$id=="MG"] >= 20)) {
      counter <- counter + 1
    }
  }
  #If the criteria are met for all the cell types, the gene is worthy. Otherwise
  #all the rows with that gene are deleted.
  if (counter<length(cell_types_not_mg)) {
    preplot$data<-preplot$data[!(preplot$data$features.plot==gene),]
  }
}

#Set the order of the genes. This just defactorizes the features.plot and sets the
#order as it is in the true_mg_genes table, which is already the desired one.
order_genes <- as.character(unique(preplot$data$features.plot))
#Also save this for the heatmap (they share y axis)
saveRDS(order_genes, file = "heatmap_order_genes1_7_hard.rds")

#Set the <0.05 % expr values to NA.
preplot$data$pct.exp[preplot$data$pct.exp < 5] <- NA

##Order cell types in this order: 
order_celltype <- c("Cyc. Prog.", "Early RG", "Late RG", "tRG", "nIPC", "SP", "GluN1", 
                    "GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8", "MGE IN", 
                    "CGE IN", "mGPC", "OPC/Oligo", "MG", "Peric.", "EC", "RBC", "VLMC")

##Actual plot
my_plot <- ggplot(preplot$data, aes(x=factor(id, levels = order_celltype),
                                    y=factor(features.plot, levels = order_genes))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="cadetblue", high = "tomato") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0,0.2,0.2),"cm"),
        legend.position = "left", axis.text.y=element_text(hjust=-0.5)) + 
  scale_y_discrete(position = "right") + 
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE, switch = "y") +
  labs(color="Avg. scaled exp.", size="Exp. %") + 
  ggtitle("Microglia Specific Genes PS1_7")
my_plot
#Checked it, success this time (with the difficult filtering) and the facet_grid

##Save objects
ggsave(my_plot, height = 9, width = 6,
       filename = "true_mg_avg_pct_exp_genes_dotplot1_7_hard.pdf",
       dpi = 700)

print("Done!")

#Addition of heatmap in lavin_winter.r script (corresponds to whole new dataset)
