#This script has no dataset per se, the macaque and mouse tables were obtained
#from querying BioMart. 

##Libraries + WD
require(dplyr)
require(ggplot2)
require(ggrepel)
require(patchwork)
setwd("/home/onecellardoor/Documents/bioinformatics/practicum/imim/data/orthology")

##Load objects
macaque <- readRDS(file = "macaque_ortholegues_unique.rds")
mouse <- readRDS(file = "mouse_ortholegues.rds")
chicken <- readRDS(file = "chicken_ortholegues.rds")
mg_genes8_22 <- readRDS("heatmap_order_genes8_22_hard.rds")
mg_genes1_7 <- readRDS("heatmap_order_genes1_7_hard.rds")
ps_table <- readRDS("ps_table.rds")

##Preprocess the tables (done and hushed)
#macaque <- read.delim(file = "macaque_ortholegues_unique.tsv", sep = "\t")
#mouse <- read.delim(file = "mouse_ortholegues.tsv", sep = "\t")
#macaque <- macaque[complete.cases(macaque), ]
#mouse <- mouse[complete.cases(mouse), ]
#chicken <- read.delim(file = "chicken_ortholegues.tsv", sep = "\t")
#chicken <- chicken[complete.cases(chicken), ]
#saveRDS(macaque, file = "macaque_ortholegues_unique.rds")
#saveRDS(mouse, file = "mouse_ortholegues.rds")
#saveRDS(chicken, file = "chicken_ortholegues.rds")

##Add PS
macaque$PS <- ps_table$PS[match(macaque$Gene.name,ps_table$Name)]
mouse$PS <- ps_table$PS[match(mouse$Gene.name,ps_table$Name)]
chicken$PS <- ps_table$PS[match(chicken$Gene.name,ps_table$Name)]
#Fix the no-name = PS31 thing
macaque$PS[macaque$Gene.name==""] <- NA
mouse$PS[mouse$Gene.name==""] <- NA
chicken$PS[chicken$Gene.name==""] <- NA
#Remove rows the PS of which correspond to NA values (maybe counterproductive)
macaque <- macaque[!is.na(macaque$PS),]
mouse <- mouse[!is.na(mouse$PS),]
chicken <- chicken[!is.na(chicken$PS),]

##Subset specific MG genes. I don't know if this is actually useful
mg_macaque1_7 <- macaque[macaque$Gene.name %in% mg_genes1_7,]
mg_macaque8_22 <- macaque[macaque$Gene.name %in% mg_genes8_22,]
mg_mouse1_7 <- mouse[mouse$Gene.name %in% mg_genes1_7,]
mg_mouse8_22 <- mouse[mouse$Gene.name %in% mg_genes8_22,]
mg_chicken1_7 <- chicken[chicken$Gene.name %in% mg_genes1_7,]
mg_chicken8_22 <- chicken[chicken$Gene.name %in% mg_genes8_22,]

##Exploratory data violin plots of % of query target identical to non-human species gene

#Make sure that PS column is a factor (this step is recommended in a tutorial)
mouse$PS <- as.factor(mouse$PS)
macaque$PS <- as.factor(macaque$PS)
chicken$PS <- as.factor(chicken$PS)
#Helper function to get the number (counts) of each group
get_box_stats <- function(y, upper_limit = quantile(y,0.95)*1.1) {
  return(data.frame(y = 0.95 * upper_limit, label = length(y)))
}
#Actual violin + box plot for the mouse data
violin_plot_mouse <- ggplot(mouse[mouse$Mouse.orthology.confidence..0.low..1.high. == 1,], aes(x=PS, y=X.id..query.gene.identical.to.target.Mouse.gene, fill=PS), show.legend = FALSE) + 
  geom_violin(trim = TRUE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() + 
  labs(x = "Phylostratum", y = "% of identity to mouse gene") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = -0.4, vjust =0.5)
violin_plot_mouse
ggsave(violin_plot_mouse, height = 6, width = 10, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/violin_plots/mouse_ortholegues_violin.pdf",
       dpi = 700)

#Actual violin + box plot for the macaque data
violin_plot_macaque <- ggplot(macaque[macaque$Macaque.orthology.confidence..0.low..1.high. == 1,], aes(x=PS, y=X.id..query.gene.identical.to.target.Macaque.gene, fill=PS), show.legend = FALSE) + 
  geom_violin(trim = TRUE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() + 
  labs(x = "Phylostratum", y = "% of identity to macaque gene") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text")
violin_plot_macaque
ggsave(violin_plot_macaque, height = 6, width = 10, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/violin_plots/macaque_ortholegues_violin.pdf",
       dpi = 700)

#Actual violin + box plot for the chicken data
violin_plot_chicken <- ggplot(chicken[chicken$Chicken.orthology.confidence..0.low..1.high. == 1,], aes(x=PS, y=X.id..query.gene.identical.to.target.Chicken.gene, fill=PS), show.legend = FALSE) + 
  geom_violin(trim = TRUE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() + 
  labs(x = "Phylostratum", y = "% of identity to chicken gene") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = -0.3, vjust =0.8)
violin_plot_chicken
ggsave(violin_plot_chicken, height = 6, width = 10, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/violin_plots/chicken_ortholegues_violin.pdf",
       dpi = 700)


##################################PS8-17#######################################

##Statistical testing

#Subsetting the datasets for just the "big" chordata cluster (PS8-17)
mouse$PS <- as.integer(mouse$PS)
macaque$PS <- as.integer(macaque$PS)
chicken$PS <- as.integer(chicken$PS)
mouse8_17 <- subset(mouse, PS > 7 & PS < 18)
macaque8_17 <- subset(macaque, PS > 7 & PS < 18)
chicken8_17 <- subset(chicken, PS > 7 & PS < 18) #Check maximum shared PS with chicken, it's not 17 probably
mg_mouse8_17 <- subset(mg_mouse8_22, PS > 7 & PS < 18)
mg_macaque8_17 <- subset(mg_macaque8_22, PS > 7 & PS < 18)
mg_chicken8_17 <- subset(mg_chicken8_22, PS > 7 & PS < 18)

mouse8_17_not_mg <- mouse8_17[!mouse8_17$Gene.name %in% mg_genes8_22,]
macaque8_17_not_mg <- macaque8_17[!macaque8_17$Gene.name %in% mg_genes8_22,]
chicken8_17_not_mg <- chicken8_17[!chicken8_17$Gene.name %in% mg_genes8_22,]

#Applying Shapiro-Wilk test to see if the distributions are normal-lile
shapiro.test(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene) #Normal (probs by chance)
shapiro.test(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse8_17$Mouse.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene) #Not normal
shapiro.test(mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse8_17_not_mg$Mouse.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene) #Not normal
shapiro.test(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque8_17$Macaque.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene) #Not normal
shapiro.test(macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque8_17_not_mg$Macaque.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene) #Not normal
shapiro.test(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken8_17$Chicken.homology.type=='ortholog_one2one']) #Normal
shapiro.test(chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene) #Not normal
shapiro.test(chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken8_17_not_mg$Chicken.homology.type=='ortholog_one2one']) #Not normal

#Applyting Fligner-Killeen test to assess if the variances are the same (non parametric test)
fligner_mouse <- fligner.test(x = list(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene,
                                       mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene)) #Significant :/
fligner_mouse_one2one <- fligner.test(x = list(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse8_17$Mouse.homology.type=='ortholog_one2one'],
                                               mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse8_17_not_mg$Mouse.homology.type=='ortholog_one2one'])) #Non-significant
fligner_macaque <- fligner.test(x = list(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene,
                                         macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene)) #Non-significant
fligner_macaque_one2one <- fligner.test(x = list(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque8_17$Macaque.homology.type=='ortholog_one2one'],
                                                 macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque8_17_not_mg$Macaque.homology.type=='ortholog_one2one'])) #Non-significant
fligner_chicken <- fligner.test(x = list(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene,
                                         chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene)) #Significant :/
fligner_chicken_one2one <- fligner.test(x = list(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken8_17$Chicken.homology.type=='ortholog_one2one'],
                                                 chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken8_17_not_mg$Chicken.homology.type=='ortholog_one2one'])) #Significant :/

#Applying Mann-Whitney Wilcoxon (non parametric)to see whether there is a significant
#difference in evolution rate between mg_enriched chordate genes and the rest of chordate
#genes (big version)
wilcox_mouse <- wilcox.test(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene,
                            mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene,
                            alternative = "two.sided") #Not significant
wilcox_mouse_one2one <- wilcox.test(mg_mouse8_17$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse8_17$Mouse.homology.type=='ortholog_one2one'],
                                    mouse8_17_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse8_17_not_mg$Mouse.homology.type=='ortholog_one2one'],
                                    alternative = "two.sided") #Significant
wilcox_macaque <- wilcox.test(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene,
                              macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene,
                              alternative = "two.sided") #Not significant
wilcox_macaque_one2one <- wilcox.test(mg_macaque8_17$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque8_17$Macaque.homology.type=='ortholog_one2one'],
                                      macaque8_17_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque8_17_not_mg$Macaque.homology.type=='ortholog_one2one'],
                                      alternative = "two.sided") #Not significant
wilcox_chicken <- wilcox.test(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene,
                              chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene,
                              alternative = "two.sided") #Not significant
wilcox_chicken_one2one <- wilcox.test(mg_chicken8_17$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken8_17$Chicken.homology.type=='ortholog_one2one'],
                                      chicken8_17_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken8_17_not_mg$Chicken.homology.type=='ortholog_one2one'],
                                      alternative = "two.sided") #Significant


##Create new dataframe for violin plot

#Subset elements of interest from the dataframe
mouse_violin_df <- mouse8_17 %>% select(Gene.name, X.id..query.gene.identical.to.target.Mouse.gene, Mouse.homology.type)
colnames(mouse_violin_df) <- c("Gene","Pct_id", "Type")
macaque_violin_df <- macaque8_17 %>% select(Gene.name, X.id..query.gene.identical.to.target.Macaque.gene, Macaque.homology.type)
colnames(macaque_violin_df) <- c("Gene","Pct_id", "Type")
chicken_violin_df <- chicken8_17 %>% select(Gene.name, X.id..query.gene.identical.to.target.Chicken.gene, Chicken.homology.type)
colnames(chicken_violin_df) <- c("Gene","Pct_id", "Type")

#Add column corresponding to the non-human species
mouse_violin_df$Species <- "Mouse"
macaque_violin_df$Species <- "Macaque"
chicken_violin_df$Species <- "Chicken"

#Add column corresponding to the color of the dots
for (i in 1:length(mg_genes8_22)) {
  mouse_violin_df$Color[mouse_violin_df$Gene == mg_genes8_22[i]] <- "Red"
  macaque_violin_df$Color[macaque_violin_df$Gene == mg_genes8_22[i]] <- "Red"
  chicken_violin_df$Color[chicken_violin_df$Gene == mg_genes8_22[i]] <- "Red"
}
#Change the rest of the values (NA) of the color column to black
#mouse_violin_df$Color[is.na(mouse_violin_df$Color)] <- "Black"
#macaque_violin_df$Color[is.na(macaque_violin_df$Color)] <- "Black"

#Stack the two dataframes into a single one (appending them)
violin_df <- rbind(mouse_violin_df, macaque_violin_df, chicken_violin_df)

#Put the species as a factor
violin_df$Species <- as.factor(violin_df$Species)

#Create a new dataframe with just the one2one ortholegues
violin_one2one_df <- subset(violin_df, Type == "ortholog_one2one")

#Create new columns "All" and "Fill" in each dataframe
violin_df$All <- violin_df$Species
violin_df$Fill <- "Cian"
violin_one2one_df$All <- paste0(violin_one2one_df$Species, " 1 to 1")
violin_one2one_df$Fill <- "Gold"

#Stack the 2 dataframes together
violin_df <- rbind(violin_df, violin_one2one_df)

#Put the "species"All" column as a factor
violin_df$All <- as.factor(violin_df$All)

#Helper function to get the count of observations per factor
get_box_stats <- function(y, upper_limit = quantile(y,0.95)*1.1) {
  return(data.frame(y = 0.95 * upper_limit, label = length(y)))
}

#Add pvalues
pvalues <- c(wilcox_macaque$p.value, wilcox_macaque_one2one$p.value,
             wilcox_mouse$p.value, wilcox_mouse_one2one$p.value,
             wilcox_chicken$p.value, wilcox_chicken_one2one$p.value)
pval_df <- data.frame(pval = round(pvalues,4))
pval_df$All <- c("Macaque", "Macaque 1 to 1", "Mouse", "Mouse 1 to 1",
                 "Chicken", "Chicken 1 to 1")
for (i in pval_df$All) {
  violin_df$pval[violin_df$All == i] <- pval_df$pval[pval_df$All==i]
} 
pval_df$y = -0.05

adj_pvalues <- p.adjust(pvalues, method = "fdr")

#Set the same position a priori for both the points and the labels via saving the
#same seed
pos <- position_jitter(width = 0.3, seed = 3)

order_species <- c("Macaque", "Macaque 1 to 1", "Mouse", "Mouse 1 to 1", "Chicken",
                   "Chicken 1 to 1")

#Actual plot
violin_plot <- ggplot(data = violin_df, aes(x=factor(All, levels = order_species), y=Pct_id), show.legend = FALSE) + 
  geom_violin(trim = TRUE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() +
  geom_point(data = subset(violin_df, Color=='Red'), position=pos, 
             pch=21, aes(fill=Color), show.legend = F) +
  scale_fill_manual(values = c("Red" = "red", "Cian" = "cyan", "Gold" = "yellow")) +
  labs(x = "Species compared to human", y = "Human gene identity to compared species gene (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = -0.3) +
 geom_text(data = pval_df, aes(x=factor(All, levels = order_species), y=y, label = paste0("p-val=",pval)), size = 3.5) +
  geom_text_repel(data = subset(violin_df, Color=='Red'), aes(label=Gene[Color=='Red']), 
                  position = pos, max.overlaps = 30, size = 2.5)
violin_plot

#Save the plot
ggsave(violin_plot, height = 6, width = 10, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/violin_plots/mg_genes_8_17_violin.pdf",
       dpi = 700)

###################################PS1-7#######################################

##Statistical testing

#Subsetting the datasets for just the "big" chordata cluster (PS8-17)
mouse$PS <- as.integer(mouse$PS)
macaque$PS <- as.integer(macaque$PS)
chicken$PS <- as.integer(chicken$PS)
mouse1_7 <- subset(mouse, PS > 0 & PS < 8)
macaque1_7 <- subset(macaque, PS > 0 & PS < 8)
chicken1_7 <- subset(chicken, PS > 0 & PS < 8)

mouse1_7_not_mg <- mouse1_7[!mouse1_7$Gene.name %in% mg_genes1_7,]
macaque1_7_not_mg <- macaque1_7[!macaque1_7$Gene.name %in% mg_genes1_7,]
chicken1_7_not_mg <- chicken1_7[!chicken1_7$Gene.name %in% mg_genes1_7,]

#Applying Shapiro-Wilk test to see if the distributions are normal-like
shapiro.test(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene) #Not normal
shapiro.test(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse1_7$Mouse.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene) #Sample too big
shapiro.test(mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse1_7_not_mg$Mouse.homology.type=='ortholog_one2one']) #Sample too big
shapiro.test(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene) #Not normal
shapiro.test(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque1_7$Macaque.homology.type=='ortholog_one2one']) #Not normal
shapiro.test(macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene) #Sample too big
shapiro.test(macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque1_7_not_mg$Macaque.homology.type=='ortholog_one2one']) #Sample too big
shapiro.test(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene) #Not normal
shapiro.test(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken1_7$Chicken.homology.type=='ortholog_one2one']) #Normal
shapiro.test(chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene) #Sample too big
shapiro.test(chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken1_7_not_mg$Chicken.homology.type=='ortholog_one2one']) #Sample too big

#Applyting Fligner-Killeen test to assess if the variances are the same (non parametric test)
fligner_mouse <- fligner.test(x = list(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene,
                                       mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene)) #Non-significant
fligner_mouse_one2one <- fligner.test(x = list(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse1_7$Mouse.homology.type=='ortholog_one2one'],
                                               mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse1_7_not_mg$Mouse.homology.type=='ortholog_one2one'])) #Non-significant
fligner_macaque <- fligner.test(x = list(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene,
                                         macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene)) #Non-significant
fligner_macaque_one2one <- fligner.test(x = list(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque1_7$Macaque.homology.type=='ortholog_one2one'],
                                                 macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque1_7_not_mg$Macaque.homology.type=='ortholog_one2one'])) #Non-significant
fligner_chicken <- fligner.test(x = list(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene,
                                         chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene)) #Non-significant
fligner_chicken_one2one <- fligner.test(x = list(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken1_7$Chicken.homology.type=='ortholog_one2one'],
                                                 chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken1_7_not_mg$Chicken.homology.type=='ortholog_one2one'])) #Significant

#Applying Mann-Whitney Wilcoxon (non parametric)to see whether there is a significant
#difference in evolution rate between mg_enriched chordate genes and the rest of chordate
#genes (big version)
wilcox_mouse <- wilcox.test(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene,
                            mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene,
                            alternative = "two.sided") #Significant
wilcox_mouse_one2one <- wilcox.test(mg_mouse1_7$X.id..query.gene.identical.to.target.Mouse.gene[mg_mouse1_7$Mouse.homology.type=='ortholog_one2one'],
                                    mouse1_7_not_mg$X.id..query.gene.identical.to.target.Mouse.gene[mouse1_7_not_mg$Mouse.homology.type=='ortholog_one2one'],
                            alternative = "two.sided") #Significant
wilcox_macaque <- wilcox.test(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene,
                              macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene,
                              alternative = "two.sided") #Significant
wilcox_macaque_one2one <- wilcox.test(mg_macaque1_7$X.id..query.gene.identical.to.target.Macaque.gene[mg_macaque1_7$Macaque.homology.type=='ortholog_one2one'],
                                      macaque1_7_not_mg$X.id..query.gene.identical.to.target.Macaque.gene[macaque1_7_not_mg$Macaque.homology.type=='ortholog_one2one'],
                              alternative = "two.sided") #Significant
wilcox_chicken <- wilcox.test(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene,
                              chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene,
                              alternative = "two.sided") #Significant
wilcox_chicken_one2one <- wilcox.test(mg_chicken1_7$X.id..query.gene.identical.to.target.Chicken.gene[mg_chicken1_7$Chicken.homology.type=='ortholog_one2one'],
                                      chicken1_7_not_mg$X.id..query.gene.identical.to.target.Chicken.gene[chicken1_7_not_mg$Chicken.homology.type=='ortholog_one2one'],
                              alternative = "two.sided") #Significant

##Create new dataframe for violin plot

#Subset elements of interest from the dataframe
mouse_violin_df <- mouse1_7 %>% select(Gene.name, X.id..query.gene.identical.to.target.Mouse.gene, Mouse.homology.type)
colnames(mouse_violin_df) <- c("Gene","Pct_id","Type")
macaque_violin_df <- macaque1_7 %>% select(Gene.name, X.id..query.gene.identical.to.target.Macaque.gene, Macaque.homology.type)
colnames(macaque_violin_df) <- c("Gene","Pct_id","Type")
chicken_violin_df <- chicken1_7 %>% select(Gene.name, X.id..query.gene.identical.to.target.Chicken.gene, Chicken.homology.type)
colnames(chicken_violin_df) <- c("Gene","Pct_id","Type")

#Add column corresponding to the non-human species
mouse_violin_df$Species <- "Mouse"
macaque_violin_df$Species <- "Macaque"
chicken_violin_df$Species <- "Chicken"

#Add column corresponding to the color of the dots
for (i in 1:length(mg_genes1_7)) {
  mouse_violin_df$Color[mouse_violin_df$Gene == mg_genes1_7[i]] <- "Red"
  macaque_violin_df$Color[macaque_violin_df$Gene == mg_genes1_7[i]] <- "Red"
  chicken_violin_df$Color[chicken_violin_df$Gene == mg_genes1_7[i]] <- "Red"
}
#Change the rest of the values (NA) of the color column to black
#mouse_violin_df$Color[is.na(mouse_violin_df$Color)] <- "Black"
#macaque_violin_df$Color[is.na(macaque_violin_df$Color)] <- "Black"

#Stack the two dataframes into a single one (appending them)
violin_df <- rbind(mouse_violin_df, macaque_violin_df, chicken_violin_df)

#Put the species as a factor
violin_df$Species <- as.factor(violin_df$Species)

#Create a new dataframe with just the one2one ortholegues
violin_one2one_df <- subset(violin_df, Type == "ortholog_one2one")

#Create new columns "All" and "Fill" in each dataframe
violin_df$All <- violin_df$Species
violin_df$Fill <- "Cian"
violin_one2one_df$All <- paste0(violin_one2one_df$Species, " 1 to 1")
violin_one2one_df$Fill <- "Gold"

#Stack the 2 dataframes together
violin_df <- rbind(violin_df, violin_one2one_df)

#Put the "species"All" column as a factor
violin_df$All <- as.factor(violin_df$All)

#Helper function to get the count of observations per factor
get_box_stats <- function(y, upper_limit = quantile(y,0.95)*1.1) {
  return(data.frame(y = 0.95 * upper_limit, label = length(y)))
}

#Add pvalues
pvalues <- c(wilcox_macaque$p.value, wilcox_macaque_one2one$p.value,
             wilcox_mouse$p.value, wilcox_mouse_one2one$p.value,
             wilcox_chicken$p.value, wilcox_chicken_one2one$p.value)
pval_df <- data.frame(pval = round(pvalues,4))
pval_df$All <- c("Macaque", "Macaque 1 to 1", "Mouse", "Mouse 1 to 1",
                 "Chicken", "Chicken 1 to 1")
for (i in pval_df$All) {
  violin_df$pval[violin_df$All == i] <- pval_df$pval[pval_df$All==i]
} 
pval_df$y = -0.05

adj_pvalues <- p.adjust(pvalues, method = "fdr")

#Set the same position a priori for both the points and the labels via saving the
#same seed
pos <- position_jitter(width = 0.3, seed = 3)

order_species <- c("Macaque", "Macaque 1 to 1", "Mouse", "Mouse 1 to 1", "Chicken",
                   "Chicken 1 to 1")

#Actual plot
violin_plot <- ggplot(data = violin_df, aes(x=factor(All, levels = order_species), y=Pct_id), show.legend = FALSE) + 
  geom_violin(trim = TRUE) + geom_boxplot(width=0.1, outlier.shape = NA) + theme_bw() +
  geom_point(data = subset(violin_df, Color=='Red'), position=pos, 
             pch=21, aes(fill=Color), show.legend = F) +
  scale_fill_manual(values = c("Red" = "red", "Cian" = "cyan", "Gold" = "yellow")) +
  labs(x = "Species compared to human", y = "Human gene identity to compared species gene (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = -0.3) +
  geom_text(data = pval_df, aes(x=factor(All, levels = order_species), y=y, label = paste0("p-val=",pval)), size = 3.5)
  #geom_text_repel(data = subset(violin_df, Color=='Red'), aes(label=Gene[Color=='Red']), 
  #                position = pos, max.overlaps = 30, size = 2.5)
violin_plot

#Save the plot
ggsave(violin_plot, height = 6, width = 10, 
       filename = "/home/onecellardoor/Dropbox/imim/records/plots/violin_plots/mg_genes_1_7_violin.pdf",
       dpi = 700)
