##Data in GSE141862

##########################SEURAT OBJECT (CLUSTER)###############################

##Creation of merged counts and Seurat object extracted from code of the authors

# Load libraries and set basic options
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(Seurat)
options(stringsAsFactors = F)
ct_target_dir <-"/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE141862/"
ct_target_dir <- "/users/genomics/xmarti/data/GSE141862/"

#####-----------------------------------------------------------------------------------------------#
## Functions ##
#####-----------------------------------------------------------------------------------------------#

# Make some quick QC graphs based on seurat QC
quickPoolQC <- function(countmatrix,pool,superpool) {
  # Create directory
  dir.create("qc/quick_qc_unfiltered", showWarnings = F)
  dir.create(paste0("qc/quick_qc_unfiltered/",superpool), showWarnings = F)
  # Create SingleCellExp object for easy QC
  seuset <- CreateSeuratObject(counts=countmatrix, project=paste0(pool))
  # Save graphs
  pdf(paste0("qc/quick_qc_unfiltered/",superpool,"/",pool,"_qc.pdf"), useDingbats=F)
  # Plot features by counts
  print(FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  # Plot total number of counts (library size)
  hist(seuset@meta.data$nCount_RNA, breaks=100,
       xlab="Library size (nCount_RNA)", ylab="Number of cells",
       main=paste("Median total nCount_RNA:", median(seuset@meta.data$nCount_RNA), "#Cells: ", ncol(seuset), "in", seuset@project.name)) #total_counts = total counts of cell, i.e. library size
  abline(v=median(seuset@meta.data$nCount_RNA), col="red")
  # Plot total number of features above detection limit (= genes)
  hist(seuset@meta.data$nFeature_RNA, breaks=100,
       xlab="# of features (nFeature_RNA)", ylab="Number of cells",
       main=paste("Median total nFeature_RNA:", median(seuset@meta.data$nFeature_RNA), "#Cells: ", ncol(seuset), "in", seuset@project.name)) #total_features_by_counts = number of features above detection limit (default=0)
  abline(v=median(seuset@meta.data$nFeature_RNA), col="red")
  # Violin plot for nCount_RNA and nFeature_RNA
  print(VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
  dev.off()
}

# Read all raw countfiles in one folder + adjust colnames
loadCountFiles <- function(ctdir,pool_file,runqc=F) {
  #Get pool and superpool
  pool <- substr(pool_file, 12, 29)
  superpool <- poolmeta[poolmeta$pool_id==pool,"superpool"]
  # Loading pool ...
  print(paste("Loading pool:",pool))
  # Read counts
  currentDf <- read.table(paste0(ctdir,pool_file), sep = "\t", header = T, row.names=1)
  # Change colnames to pool + cellbarcode. So if several pools are merged later, we always know where the cell came from.
  names(currentDf) <- paste(pool,colnames(currentDf), sep="_") 
  # Run quick QC on pool
  if(runqc==T) { quickPoolQC(currentDf,pool,superpool) }
  # Return pool dataframe
  return(currentDf)
}

# Get meta data for all cells of each donor
getMeta <- function(countDataDf) {
  # Get all cell_ids
  cell_ids <- colnames(countDataDf)
  # Make dataframe with cell_ids as rownames
  meta_df <- data.frame(row.names = cell_ids)
  # Extract pool_id from cell_id
  pool <- substr(cell_ids, 1, 14)
  # Add meta to each cell
  meta_df <- cbind(meta_df, poolmeta[match(pool,poolmeta$pool_id),])
  # Return meta dataframe
  return(meta_df)
}

#####-----------------------------------------------------------------------------------------------#
##Load countfiles + merge ##
#####-----------------------------------------------------------------------------------------------#

# Metadata and count files are available at NCBI GEO under accession number GSE141862.
# Infos can be found in the manuscript.

# Load target file and create info variables
poolmeta <- read.csv(paste0(ct_target_dir,"target.csv"))

# Define count directory
ctdir = paste0(ct_target_dir,"countfiles/")
# Make list of all countfiles in this directory (including subdirectories)
ct_file_list <- list.files(path=ctdir, pattern="*completeCounts.txt", recursive = TRUE)

# Does a merged count file exist already?
#all_cts = read.table("merged_counts.txt", header=T, check.names = F)

# The following code allows continueing of merging data in case script stopped (memory issues)
if(exists("all_cts")) {
  # Get all pool names that already exist in all_cts
  done = unique(substr(names(all_cts),1,18))
  # Remove already loaded pool from ct_file_list
  ct_file_list = ct_file_list[!grepl(paste(done,collapse="|"), ct_file_list)]
} else {
  # Otherwise create new all_cts dataframe to begin merging
  all_cts = data.frame()
}

# Load all pools, run quick QC, and merge to big dataframe. This will take very long.
# Recommended to run with at least 64gb of memory.
for(pool_file in ct_file_list) {
  # Merge pools in one dataframe
  all_cts <- base::merge(all_cts, loadCountFiles(ctdir,pool_file, runqc=F), by="row.names", all=T)
  # Genes that were not present in one pool but in the other are marked with NA, replace NA by 0
  all_cts <- replace(all_cts, is.na(all_cts), 0)
  # Move ensembl back to rownames and delete row.names column
  rownames(all_cts) <- all_cts$Row.names
  all_cts <- all_cts[-1]
  gc()
}

# Export raw merged counts
write.table(all_cts, paste0(ctdir,"merged_counts.tsv"), sep = "\t")

# Load gene annotation
#gene_annotation <- read.csv(paste0(ct_target_dir,"gene_annotation.csv"), header=T)
# Replace ensembl rownames by gene rownames
#rownames(all_cts) <- gene_annotation[match(rownames(all_cts), gene_annotation$gene_id),"gene_name"]

#####-----------------------------------------------------------------------------------------------#
## Create seurat object ##
#####-----------------------------------------------------------------------------------------------#

# Create combined seurat object
seu_all_raw <- CreateSeuratObject(counts=all_cts,
                              project="kracht",
                              meta.data=getMeta(all_cts),
                              min.cells = 0) # Genes need to be expressed in at least 0 cells
seu_all_raw$age[seu_all_raw$age=="9"] <- "09"
seu_all_raw$age_donor <- paste0(seu_all_raw$age,"_",seu_all_raw$donor)

## Save Seurat Object
saveRDS(seu_all_raw, file = "kracht_seurat.rds")

## Save Seurat Object (not done because I forgot to do it in the first part (it
## was just to save it), done in a separate script from the merged counts and the
## getMeta function

# Load libraries and set basic options
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(Seurat)
options(stringsAsFactors = F)
setwd('/users/genomics/xmarti/data/GSE141862')
# Load objects
poolmeta <- read.csv(paste0("target.csv"))
all_cts <- read.table('merged_counts.txt',sep = " ", header = TRUE)

# Functions
# Get meta data for all cells of each donor
getMeta <- function(countDataDf) {
  # Get all cell_ids
  cell_ids <- colnames(countDataDf)
  # Make dataframe with cell_ids as rownames
  meta_df <- data.frame(row.names = cell_ids)
  # Extract pool_id from cell_id
  pool <- substr(cell_ids, 1, 14)
  # Add meta to each cell
  meta_df <- cbind(meta_df, poolmeta[match(pool,poolmeta$pool_id),])
  # Return meta dataframe
  return(meta_df)
}

# Create Seurat Object
seu_all_raw <- CreateSeuratObject(counts=all_cts,
                                  project="kracht",
                                  meta.data=getMeta(all_cts),
                                  min.cells = 0) # Genes need to be expressed in at least 0 cells
seu_all_raw$age[seu_all_raw$age=="9"] <- "09"
seu_all_raw$age_donor <- paste0(seu_all_raw$age,"_",seu_all_raw$donor)

saveRDS(seu_all_raw, file = "kracht_seurat.rds")

## Fix the merge. Ultimately not done, decided to re-run previous script and save
## The seurat object
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE141862')
a <- data.frame(c(0,1,2),row.names=1)
b <- data.frame(c(3,4,5,6),rownames=1)
c <- data.frame()
c <- merge(c,a,by="row.names", all=T)
c <- merge(c,b,by="row.names", all=T)
currentDf <- read.table("countfiles/GSM4214857_103591-001-000-001_S1_completeCounts.txt", sep = "\t", header = T, row.names=1)
currentDf2 <- read.table("countfiles/GSM4214858_103591-001-000-002_S2_completeCounts.txt", sep = "\t", header = T, row.names=1)
merged <- merge(currentDf,currentDf2,by="row.names", all=T)
write.table(merged, paste0("merged.txt"), sep = "\t")

## Add metadata
#Read the object and the metadata file
kracht_seurat <- readRDS('kracht_seurat.rds')
kracht_metadata <- read.delim('metadata.csv', sep = ",")
#Create column for matching from the rownames
kracht_seurat@meta.data$cell_pool_id <- rownames(kracht_seurat@meta.data)
#Manual implementation of the metadata (the automated version was a pain in the ass)
kracht_seurat@meta.data$superpool <- kracht_metadata$superpool[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$pool_id <- kracht_metadata$pool_id[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$pool <- kracht_metadata$pool[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$donor <- kracht_metadata$donor[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$age <- kracht_metadata$age[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$celltype <- kracht_metadata$celltype[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$sex <- kracht_metadata$sex[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$percent.mt <- kracht_metadata$percent.mt[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$percent.ribo <- kracht_metadata$percent.ribo[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$S.Score <- kracht_metadata$S.Score[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$G2M.Score <- kracht_metadata$G2M.Score[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$Phase <- kracht_metadata$Phase[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$seurat_clusters <- kracht_metadata$seurat_clusters[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
kracht_seurat@meta.data$age_exact <- kracht_metadata$age_exact[match(kracht_seurat@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]

#Add cluster names
require(dplyr)
kracht_seurat@meta.data <- kracht_seurat@meta.data %>% mutate("Clustered_CT" =
                                                                case_when(seurat_clusters == "1" ~ "Prog.",
                                                                          seurat_clusters == "2" ~ "Prog.",
                                                                          seurat_clusters == "3" ~ "Prog.",
                                                                          seurat_clusters == "4" ~ "Prog.",
                                                                          seurat_clusters == "5" ~ "Prog.",
                                                                          seurat_clusters == "6" ~ "Prog.",
                                                                          seurat_clusters == "7" ~ "Prog.",
                                                                          seurat_clusters == "8" ~ "Prog.",
                                                                          seurat_clusters == "9" ~ "Prog.",
                                                                          seurat_clusters == "10" ~ "Prog.",
                                                                          seurat_clusters == "11" ~ "Prog.",
                                                                          seurat_clusters == "12" ~ "Prog.",
                                                                          seurat_clusters == "13" ~ "Prog.",
                                                                          seurat_clusters == "14" ~ "Prog.",
                                                                          seurat_clusters == "15" ~ "Prog.",
                                                                          seurat_clusters == "16" ~ "Prog."))

#Assign identations
Idents(kracht_seurat) <- kracht_seurat@meta.data$seurat_clusters

#Save SeuratObject
saveRDS(kracht_seurat, file = 'kracht_seurat.rds')

##############################QC & FILTERING###################################

##Code from the original paper's authors

# Load libraries, objects and set basic options
library(Seurat)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE141862')
seu_all_raw <- readRDS('kracht_seurat.rds')
options(stringsAsFactors = F)

#####-----------------------------------------------------------------------------------------------#
## Functions ##
#####-----------------------------------------------------------------------------------------------#

# Seurat quick basic QC including MT and Ribo graphs
runSeuratQC <- function(seuset, title) {
  # Plot features by counts
  print(FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "celltype"))
  # Plot total number of counts (library size)
  hist(seuset@meta.data$nCount_RNA, breaks=100,
       xlab="Library size (nCount_RNA)", ylab="Number of cells",
       main=paste("Median total nCount_RNA:", median(seuset@meta.data$nCount_RNA), "in", seuset@project.name)) #total_counts = total counts of cell, i.e. library size
  abline(v=median(seuset@meta.data$nCount_RNA), col="red")
  # Plot total number of features above detection limit (= genes)
  hist(seuset@meta.data$nFeature_RNA, breaks=100,
       xlab="# of features (nFeature_RNA)", ylab="Number of cells",
       main=paste("Median total nFeature_RNA:", median(seuset@meta.data$nFeature_RNA), "in", seuset@project.name)) #total_features_by_counts = number of features above detection limit (default=0)
  abline(v=median(seuset@meta.data$nFeature_RNA), col="red")
  # Violin plot for nCount_RNA and nFeature_RNA
  print(VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "celltype")  + theme(legend.position = "none"))
  print(VlnPlot(object = seuset, features = c("nFeature_RNA"), group.by="donor") + theme(legend.position = "none"))
  print(VlnPlot(object = seuset, features = c("nFeature_RNA"), group.by="superpool") + theme(legend.position = "none"))
  # Make tables of nFeatures, nCount, percent.mt, percent.ribo
  per_donor = data.frame(aggregate(cbind(nFeature_RNA,nCount_RNA, percent.mt, percent.ribo) ~ donor, seuset@meta.data, FUN=median))
  per_donor$nCells = sapply(per_donor$donor, function(y) nrow(seuset@meta.data[seuset@meta.data$donor==y,]))
  write.table(per_donor, paste0("qc/donor_qc_medians_",title,".csv"), col.names = NA, sep=",")
  
  per_superpool = data.frame(aggregate(cbind(nFeature_RNA,nCount_RNA, percent.mt, percent.ribo) ~ superpool, seuset@meta.data, FUN=median))
  per_superpool$nCells = sapply(per_superpool$superpool, function(y) nrow(seuset@meta.data[seuset@meta.data$superpool==y,]))
  write.table(per_superpool, paste0("qc/superpool_qc_medians_",title,".csv"), col.names = NA, sep=",")
  
  per_pool_id = data.frame(aggregate(cbind(nFeature_RNA,nCount_RNA, percent.mt, percent.ribo) ~ pool_id, seuset@meta.data, FUN=median))
  per_pool_id$nCells = sapply(per_pool_id$pool_id, function(y) nrow(seuset@meta.data[seuset@meta.data$pool_id==y,]))
  write.table(per_pool_id, paste0("qc/pool_id_qc_medians_",title,".csv"), col.names = NA, sep=",")
}

# Visualise cell filtering technique
visualiseCellFilter <- function(seuset) {
  tot_cts <- log10(seuset@meta.data$nCount_RNA)
  tot_cts_median <- median(tot_cts)
  tot_cts_mad <- mad(tot_cts)
  print(paste("Median total counts:",10^tot_cts_median))
  print(paste("Mad total counts:",10^tot_cts_mad))
  # Set thresholds
  cell_thresholds <- list(low = tot_cts_median-3*tot_cts_mad,
                          high = tot_cts_median+3*tot_cts_mad)
  # Visualise using histogram
  hist(tot_cts,col="grey80",breaks = 100,
       xlab="log10(nGene_RNA)",  main=paste("# of genes in ", seuset@project.name),
       ylab="Number of cells"
  )
  abline(v = cell_thresholds$low, col = "red")
  abline(v = cell_thresholds$high, col = "red")
  # Return thresholds
  return(cell_thresholds)
}

#####-----------------------------------------------------------------------------------------------#
## MT, RPS, RPL gene labelling and QC ##
#####-----------------------------------------------------------------------------------------------#

# Label MT and ribo genes
seu_all_raw[["percent.mt"]] <- PercentageFeatureSet(object=seu_all_raw,pattern="^MT-")
seu_all_raw[["percent.ribo"]] <- PercentageFeatureSet(object=seu_all_raw,pattern="^RPS|^RPL")

# Export graph: MT counts by total features counts
pdf("MT-Ribo_gene_content.pdf", useDingbats = F)
FeatureScatter(object = seu_all_raw, feature1 = "percent.mt", feature2 = "percent.ribo", group.by = "celltype")
print(VlnPlot(object = seu_all_raw, features = c("percent.mt", "percent.ribo"), ncol = 2, group.by="celltype"))
dev.off()

#####-----------------------------------------------------------------------------------------------#
## Run before-filtering QC ##
#####-----------------------------------------------------------------------------------------------#

# Run standard QC before filtering.
pdf("general_qc_unfiltered.pdf", useDingbats = F)
runSeuratQC(seu_all_raw,"unfiltered")
dev.off()

#####-----------------------------------------------------------------------------------------------#
## Cell filtering (Based on MT content and low/high counts) ##
#####-----------------------------------------------------------------------------------------------#

# Find thresholds for low/high count cells + plot graph
# Function will return low and high treshold based on median absolute deviation
pdf("cell_filtering_approach.pdf", useDingbats = F)
cell_thresholds <- visualiseCellFilter(seu_all_raw)
dev.off()

# MT gene expression threshold in %
MT_threshold <- 10

# Do the actual filtering
seu_all_filtered <- subset(x = seu_all_raw, subset = nCount_RNA > 10^cell_thresholds$low &
                             nCount_RNA < 10^cell_thresholds$high &
                             percent.mt < MT_threshold)

#####-----------------------------------------------------------------------------------------------#
## Export overview and QC ##
#####-----------------------------------------------------------------------------------------------#

# QC
pdf("qc/general_qc_filtered.pdf", useDingbats = F)
runSeuratQC(seu_all_filtered, "filtered")
dev.off()

# Create file with all filtering info
write.csv(data.frame(cells_before = nrow(seu_all_raw@meta.data),
                     cells_after = nrow(seu_all_filtered@meta.data),
                     below_threshold = nrow(subset(seu_all_raw, subset= nCount_RNA < 10^cell_thresholds$low)@meta.data),
                     above_threshold = nrow(subset(seu_all_raw, subset= nCount_RNA > 10^cell_thresholds$high)@meta.data),
                     mt_filter = nrow(subset(seu_all_raw, subset= percent.mt > MT_threshold)@meta.data)),
          "qc/overview_filtering.csv")

#####-----------------------------------------------------------------------------------------------#
## Pool QC scatterplots ##
#####-----------------------------------------------------------------------------------------------#

# Read QC tables
poolqc = read.csv('qc/pool_id_qc_medians_unfiltered.csv')
poolqc_filt = read.csv('qc/pool_id_qc_medians_filtered.csv')

# Plot scatterplots nFeature and nCounts
pdf("qc/general_qc_all.pdf", useDingbats = F)
plot(poolqc$nFeature_RNA, poolqc$nCount_RNA, main='Unfiltered')
plot(poolqc_filt$nFeature_RNA, poolqc_filt$nCount_RNA, main='Filtered')
dev.off()

#Whatever whatever, there is one more cell in metadata than it should, I remove
#it manually. The rest seemed to work fine except for the fact that the percent.mt
#and percent.ribo were reset to 0 for some reason (it's probably in the code)
test <- subset(x = seu_all_filtered, subset = seurat_clusters >= 1)

#Assigning the missing values because of OCD and removing useless columns. This
#helped me confirm that I kept the same cells.
test@meta.data$percent.mt <- kracht_metadata$percent.mt[match(test@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
test@meta.data$percent.ribo <- kracht_metadata$percent.ribo[match(test@meta.data$cell_pool_id, kracht_metadata$pool_cell_id)]
test@meta.data$age_donor <- NULL
#Overwriting the SeuratObject, I don't think it's gonna be a problem.
saveRDS(test, file = 'kracht_seurat.rds')
##################NORMALIZATION & SCALING (CLUSTER)############################

#I suppose I will copy the author's parameters.

##Packages, options and WD
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(Seurat)
options(stringsAsFactors = F)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD OK")

##Load objects
kracht_seurat<-readRDS('kracht_seurat.rds')
print("Objects OK")

##Normalization
kracht_seurat <- NormalizeData(kracht_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

##Variable Features (just wanna check)
kracht_seurat <- FindVariableFeatures(object = kracht_seurat, selection.method = "mean.var.plot")
length(VariableFeatures(kracht_seurat))

##Scale Data
all.genes <- rownames(kracht_seurat)
print("Scaling all genes...")
kracht_seurat <- ScaleData(object = kracht_seurat, vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA"))

#Save the modified objects
saveRDS(kracht_seurat, file = 'kracht_seurat_toy.rds')
print("R job done")

## FROM HERE ON NOW, KRACHT_SEURAT.RDS IS THE HEAVY FILE WITH ALL GENES SCALED,
## AND KRACHT_SEURAT_TOY.RDS IS THE VERSION TO TRY STUFF IN THIS LAPTOP.

#########################DEVELOPMENTAL STAGE (CLUSTER)##########################
##Libraries + WD (for cluster)
.libPaths('/users/genomics/xmarti/packages/R/4.1')
require(dplyr)
require(Seurat)
require(patchwork)
require(ggplot2)
setwd('/users/genomics/xmarti/data')
print("Libraries and WD loaded correctly")

##Load objects
tr_table <- readRDS('gse141862_tr_table.rds')
mg_genes1_7 <- readRDS('heatmap_order_genes1_7_hard.rds')
mg_genes8_22 <- readRDS('heatmap_order_genes8_31_hard.rds')
ps_table <- readRDS('ps_table.rds')
deg_clusters <- readRDS("marker_genes_kracht_clusters.rds")
print("Loading the big fella:")
kracht_seurat <- readRDS('kracht_seurat_toy.rds')
print("Objects loaded correctly")

##Translate the symbols into IDs and store it as a character (just in case)
id_list1_7 <- as.character(tr_table$ID[match(mg_genes1_7,tr_table$Symbol)])
id_list8_22 <- as.character(tr_table$ID[match(mg_genes8_22,tr_table$Symbol)])

##Store the DotPlot fuction data into an object for further manipulation
preplot1_7 <- DotPlot(kracht_seurat, features = id_list1_7, cluster.idents = TRUE)
preplot8_22 <- DotPlot(kracht_seurat, features = id_list8_22, cluster.idents = TRUE)

##Add the PS as a new column of preplot$data

#The Genes come in ID format here so I have to translate them into symbols first
preplot1_7$data$features.plot <- as.factor(tr_table$Symbol[match(preplot1_7$data$features.plot,
                                                              tr_table$ID)])
preplot8_22$data$features.plot <- as.factor(tr_table$Symbol[match(preplot8_22$data$features.plot,
                                                                 tr_table$ID)])
#Add the PS column
preplot1_7$data$PS <- ps_table$PS[match(preplot1_7$data$features.plot,
                                        ps_table$Name)]
preplot8_22$data$PS <- ps_table$PS[match(preplot8_22$data$features.plot,
                                        ps_table$Name)]

#Add cluster names for (a slightly better) clarity
preplot1_7$data <- preplot1_7$data %>% mutate(Cluster =
                                                case_when(id == 1 ~ "GW15-17 (1)",
                                                          id == 2 ~ "GW11-12 (1)",
                                                          id == 3 ~ "GW11-12 (2)",
                                                          id == 4 ~ "GW15-17 (2)",
                                                          id == 5 ~ "IEG",
                                                          id == 6 ~ "Cell Cycle",
                                                          id == 7 ~ "GW9-10 (1)",
                                                          id == 8 ~ "GW9-10 (2)",
                                                          id == 9 ~ "Myeloid Cells (1)",
                                                          id == 10 ~ "Myeloid Cells (2)",
                                                          id == 11 ~ "MRPL23-E",
                                                          id == 12 ~ "PARP4-E",
                                                          id == 13 ~ "MTX1-E",
                                                          id == 14 ~ "HB-E",
                                                          id == 15 ~ "ZP3-E",
                                                          id == 16 ~ "NAMPT-E"))
preplot8_22$data <- preplot8_22$data %>% mutate(Cluster =
                                                case_when(id == 1 ~ "GW15-17 (1)",
                                                          id == 2 ~ "GW11-12 (1)",
                                                          id == 3 ~ "GW11-12 (2)",
                                                          id == 4 ~ "GW15-17 (2)",
                                                          id == 5 ~ "IEG",
                                                          id == 6 ~ "Cell Cycle",
                                                          id == 7 ~ "GW9-10 (1)",
                                                          id == 8 ~ "GW9-10 (2)",
                                                          id == 9 ~ "Myeloid Cells (1)",
                                                          id == 10 ~ "Myeloid Cells (2)",
                                                          id == 11 ~ "MRPL23-E",
                                                          id == 12 ~ "PARP4-E",
                                                          id == 13 ~ "MTX1-E",
                                                          id == 14 ~ "HB-E",
                                                          id == 15 ~ "ZP3-E",
                                                          id == 16 ~ "NAMPT-E"))

##Add significance (experimental)
#require(readxl)
#deg_clusters <- read_excel('kracht_st_adapted.xlsx', sheet = 3)
#saveRDS(deg_clusters, file = "marker_genes_kracht_clusters.rds")

#Try with loop. Funny, I have iterate by position and not gene because there are
#genes that are markers for more than one cluster

for (i in 1:nrow(deg_clusters)) {
  gene <- deg_clusters$gene[i]
  cluster <- deg_clusters$cluster[i]
  preplot1_7$data$deg[preplot1_7$data$features.plot == gene & preplot1_7$data$id == cluster] <- "Yes"
  preplot8_22$data$deg[preplot8_22$data$features.plot == gene & preplot8_22$data$id == cluster] <- "Yes"
}
#For the ggplot code to work there mustn't be any NA (I don't know why)
preplot1_7$data$deg[is.na(preplot1_7$data$deg)] <- "No"
preplot8_22$data$deg[is.na(preplot8_22$data$deg)] <- "No"

#Set the order of the genes
order_genes1_7 <- mg_genes1_7
order_genes8_22 <- mg_genes8_22

#Set the <0.05 % expr values to NA.
preplot1_7$data$pct.exp[preplot1_7$data$pct.exp < 5] <- NA
preplot8_22$data$pct.exp[preplot8_22$data$pct.exp < 5] <- NA

##Order clusters in this order: 
order_clusters <- c("GW9-10 (1)","GW9-10 (2)","GW11-12 (1)",
                    "GW11-12 (2)","GW15-17 (1)","GW15-17 (2)", "IEG", "Cell Cycle",
                    "MRPL23-E", "PARP4-E", "MTX1-E", "HB-E", "ZP3-E", "NAMPT-E",
                    "Myeloid Cells (1)", "Myeloid Cells (2)")

##Actual plot 1-7
my_plot1_7 <- ggplot(preplot1_7$data, aes(x=factor(Cluster, levels = order_clusters),
                                    y=factor(features.plot, levels = order_genes1_7))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  geom_point(data = preplot1_7$data[preplot1_7$data$deg=="Yes", ], aes(x=factor(Cluster, levels = order_clusters),
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
my_plot8_22 <- ggplot(preplot8_22$data, aes(x=factor(Cluster, levels = order_clusters),
                                    y=factor(features.plot, levels = order_genes8_22))) + 
  geom_point(aes(size=pct.exp, color=avg.exp.scaled)) +
  geom_point(data = preplot8_22$data[preplot8_22$data$deg=="Yes", ], aes(x=factor(Cluster, levels = order_clusters),
                                                                         y=factor(features.plot, levels = order_genes8_22)),
             shape = "*", size=4, color="black") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="#237A57", high = "#ff6a00") +  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")) + 
  facet_grid(PS~., scales="free_y", space="free_y", as.table = FALSE) +
  labs(color="Avg. scaled exp.", size="Exp. %")
my_plot8_22

##Save objects
ggsave(my_plot1_7, height = 9, width = 6,
       filename = "mg_genes_dev_stage1_7.pdf",
       dpi = 700)
print("mg_genes_dev_stage1_7.pdf saved")
ggsave(my_plot8_22, height = 6, width = 6,
       filename = "mg_genes_dev_stage8_22.pdf",
       dpi = 700)
print("mg_genes_dev_stage8_22.pdf saved")
print("R job done")

########################DEVELOPMENTAL STAGE (not used)##########################

#Gene version is GrCh38, but this looks really old, August 2014 when the paper is
#from July 2020. The authors don't specify the subversion, so I'm just gonna apply
#the translation table from the Trevino paper and hope for the best.

##Libraries + WD
require(dplyr)
require(readxl)
require(Seurat)
require(reshape2)
setwd('/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/GSE141862')

##Load objects
tr_table <- readRDS('gse141862_tr_table.rds')
genes_stages <- read_excel('kracht_st_adapted.xlsx', sheet = 5)
mg_genes1_7 <- readRDS('heatmap_order_genes1_7_hard.rds')
mg_genes8_22 <- readRDS('heatmap_order_genes8_31_hard.rds')
deg_clusters <- read_excel('kracht_st_adapted.xlsx', sheet = 3)
kracht_seurat <- readRDS('kracht_seurat_toy.rds')
ps_table <- readRDS('ps_table.rds')

###Stages (pseudotime, S5 Table)

#I don't have the values of expression nor the cells that 
#correspond to each stage. Not very epic.

##Filtering by the genes names of the selected genes. Doesn't work properly
mg_genes1_7_stages <- subset(genes_stages, gene %in% mg_genes1_7)
mg_genes8_22_stages <- subset(genes_stages, gene %in% mg_genes8_22)

##Dropping empty rows
mg_genes1_7_stages <- mg_genes1_7_stages[!is.na(mg_genes1_7_stages$gene),]
mg_genes8_22_stages <- mg_genes8_22_stages[!is.na(mg_genes8_22_stages$gene),]

###DEG (S3 Table). Has avg_logFC, so heatmap could be applied. Only available
###for DEG genes, and not all of them.

#Filtering by the genes names of the selected genes. Doesn't work properly
deg_genes1_7 <- subset(deg_clusters, gene %in% mg_genes1_7)
deg_genes8_22 <- subset(deg_clusters, gene %in% mg_genes8_22)

#Dropping empty rows
deg_genes1_7 <- deg_genes1_7[!is.na(deg_genes1_7$gene),]
deg_genes8_22 <- deg_genes8_22[!is.na(deg_genes8_22$gene),]

###Expression by cluster. I must QC_norm_scale the data, but once this is done
###the expression for every gene present in the raw data is available (CLUSTER)

##Translating the ENSEMBL IDs into symbols via obtaining a dataset-specific
##translation table (done, saved and hushed)

#ID <- rownames(kracht_seurat)
#Symbol <- tr_table$Gene.name[match(gene_ids,tr_table$Gene.stable.ID)]
#gse141862_tr_table <- data.frame(ID, Symbol)
#saveRDS(gse141862_tr_table, file = "gse141862_tr_table.rds")

##Obtaining the ENSEMBL IDs of the dataset of the gene symbols. They all have
##matching ENSEMBLS IDs, so we're good here.
genes1_7_ids <- tr_table$ID[match(mg_genes1_7,tr_table$Symbol)]
genes8_22_ids <- tr_table$ID[match(mg_genes8_22,tr_table$Symbol)]

##Obtaining average expression through Seurat Function "AverageExpression"
avg_exp_genes1_7 <- AverageExpression(object = kracht_seurat, features = genes1_7_ids, group.by = "ident")
avg_exp_genes8_22 <- AverageExpression(object = kracht_seurat, features = genes8_22_ids, group.by = "ident")
#The previous function returns a matrix in a list, which is not great. I need 
#to create a dataframe out of this. Turns out it was pretty easy.
avg_exp_genes1_7 <- as.data.frame(avg_exp_genes1_7[[1]])
avg_exp_genes8_22 <- as.data.frame(avg_exp_genes8_22[[1]])
#Reorder the columns (pending, but not necessary)

#Applying DotPlot function
DotPlot(object = kracht_seurat, features = genes1_7_ids, cluster.idents = TRUE)
DotPlot(object = kracht_seurat, features = genes8_22_ids, cluster.idents = TRUE)

##Adding the gene symbols and PS as columns to the dataframes
#Symbols
avg_exp_genes1_7$Name <- tr_table$Symbol[match(rownames(avg_exp_genes1_7),tr_table$ID)]
avg_exp_genes8_22$Name <- tr_table$Symbol[match(rownames(avg_exp_genes8_22),tr_table$ID)]
#PS
avg_exp_genes1_7$PS <- ps_table$PS[match(avg_exp_genes1_7$Name,ps_table$Name)]
avg_exp_genes8_22$PS <- ps_table$PS[match(avg_exp_genes8_22$Name,ps_table$Name)]
