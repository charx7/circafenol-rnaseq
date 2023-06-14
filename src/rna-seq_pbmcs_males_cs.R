library(readxl)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(gridExtra)
library(ggsci)
library(pheatmap)
library(reshape2)

# read the data
excel_path = "/home/rstudio/scripts/data/gene_count_matrix.xlsx"

count_data <- read_excel(
  excel_path,
  sheet = "count_data"
)

# convert tibble to dataframe to rename rows
count_data <- as.data.frame(count_data)

# select only subject columns
count_data_matrix <- as.data.frame(count_data[,-1])
gene_names <- as.character(count_data[, 1])
count_data_final <- data.frame(count_data_matrix, row.names=gene_names)
# restore column names because it appends an X char
colnames(count_data_final) <- colnames(count_data_matrix)

# remove to save memory
rm(count_data_matrix)
rm(gene_names)
rm(count_data)

# read col data
col_data <- read_excel(
  excel_path,
  sheet = "col_data"
)

col_data <- as.data.frame(col_data)
col_data <- data.frame(as.data.frame(as.data.frame(col_data[,-1])), row.names=col_data[,1])

# Add library size, which is the total amount of gene reads per sample
col_data$libSize <- colSums(count_data_final)

# Make sure that the colnames of your count table are matching with the row names of your colData
all(rownames(col_data) %in% colnames(count_data_final))

# get counts and col data by treatment
# chronodisruption -> to test for this we remove L11-GSPE
chronodisruption_sample_names <- rownames(col_data)[which(col_data$treatment != "L11-GSPE")]
chronodisruption_count_table <- count_data_final[, chronodisruption_sample_names]
chronodisruption_col_data <- col_data[chronodisruption_sample_names,]


# stimulated dataset -> GSPE vs VH on chronodisruption
stimulus_sample_names <- rownames(col_data)[which(col_data$treatment != "control")]
stimulus_count_table <- count_data_final[, stimulus_sample_names]
stimulus_col_data <- col_data[stimulus_sample_names,]

# chronodisruption GSPE vs control
chrono_gspe_sample_names <- rownames(col_data)[which(col_data$treatment != "L11-VH")]
chrono_gspe_count_table <-count_data_final[, chrono_gspe_sample_names]
chrono_gspe_col_data <- col_data[chrono_gspe_sample_names,]

# get deseq datasets
dds_chronodisruption <- DESeqDataSetFromMatrix(
  countData = chronodisruption_count_table,
  colData = chronodisruption_col_data,
  design = ~ treatment
)

dds_stimulus <- DESeqDataSetFromMatrix(
  countData = stimulus_count_table,
  colData = stimulus_col_data,
  design = ~ treatment
)

dds_chrono_gspe <- DESeqDataSetFromMatrix(
  countData = chrono_gspe_count_table,
  colData = chrono_gspe_col_data,
  design = ~ treatment
) 

# filter low count genes
keeps <- rowSums(counts(dds_chronodisruption) >= 10)
dds_chronodisruption <- dds_chronodisruption[keeps,]

keeps <- rowSums(counts(dds_stimulus) >= 10)
dds_stimulus <- dds_stimulus[keeps,]

keeps <- rowSums(counts(dds_chrono_gspe) >= 10)
dds_chrono_gspe <- dds_chrono_gspe[keeps,]

# Variance stabilizing transformations
vst_chronodisruption <- assay(vst(dds_chronodisruption, blind=FALSE))
vst_stimulus <- assay(vst(dds_stimulus, blind=FALSE))
vst_chrono_gspe <- assay(vst(dds_chrono_gspe, blind=FALSE))

## PCA
# To calculate the components by sample we need to transpose our matrix of normalized gene expression 
pcData_chrono <- prcomp(t(vst_chronodisruption))
pcVar_chrono <- summary(pcData_chrono)

pcData_stim <- prcomp(t(vst_stimulus))
pcVar_stim <- summary(pcData_stim)

pcData_chrono_gspe <- prcomp(t(vst_chrono_gspe))
pcVar_chrono_gspe <- summary(pcData_chrono_gspe)

# By getting the summary() of a prcomp object (in this case pcData) we can also obtain the total amount of variance explained by each of the components.
pcVar_chrono$importance
pcVar_stim$importance
pcVar_chrono_gspe$importance

# We can then extract the variance explained by components 1 and 2. 
varPC1_chrono <- pcVar_chrono$importance[2,1]
varPC2_chrono <- pcVar_chrono$importance[2,2]

varPC1_stim <- pcVar_stim$importance[2,1]
varPC2_stim <- pcVar_stim$importance[2,2]

varPC1_chrono_gspe <- pcVar_chrono_gspe$importance[2,1]
varPC2_chrono_gspe <- pcVar_chrono_gspe$importance[2,2]

pcPlotData_chrono <- data.frame(pcData_chrono$x[,1:4], col_data[rownames(pcData_chrono$x),])
pcPlotData_stim <- data.frame(pcData_stim$x[,1:4], col_data[rownames(pcData_stim$x),])
pcPlotData_chrono_gspe <- data.frame(pcData_chrono_gspe$x[,1:4], col_data[rownames(pcData_chrono_gspe$x),])


pcaPlot_chronodisruption <- ggplot(pcPlotData_chrono, aes(x=PC1 , y=PC2 , color=treatment))+
  geom_jitter(alpha=0.6)+
  #facet_grid(~stimulation, scales = "free")+
  xlab(paste0("PC1 explained variance = ", varPC1_chrono*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2_chrono*100, "%"))+
  scale_color_aaas()+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 8))

pcaPlot_chronodisruption

pcaPlot_stim <- ggplot(pcPlotData_stim, aes(x=PC1 , y=PC2 , color=treatment))+
  geom_jitter(alpha=0.6)+
  #facet_grid(~stimulation, scales = "free")+
  xlab(paste0("PC1 explained variance = ", varPC1_stim*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2_stim*100, "%"))+
  scale_color_aaas()+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 8))

pcaPlot_stim

pcaPlot_chrono_gspe <- ggplot(pcPlotData_chrono_gspe, aes(x=PC1 , y=PC2 , color=treatment))+
  geom_jitter(alpha=0.6)+
  #facet_grid(~stimulation, scales = "free")+
  xlab(paste0("PC1 explained variance = ", varPC1_chrono_gspe*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2_chrono_gspe*100, "%"))+
  scale_color_aaas()+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 8))

pcaPlot_chrono_gspe

# Unsupervised clustering

# Again we need to transpose our matrix to then calculate the distance between each of the samples.
sampleDists_chrono <- dist(t(vst_chronodisruption))
sampleDistMatrix_chrono <- as.matrix(sampleDists_chrono)

sampleDists_stim <- dist(t(vst_stimulus))
sampleDistMatrix_stim <- as.matrix(sampleDists_stim)

sampledists_chrono_gspe <- dist(t(vst_chrono_gspe))
sampleDistMatrix_chrono_gspe <- as.matrix(sampledists_chrono_gspe)


# By using brewer.pal() we can generate a palette of colors, for more colors check (http://colorbrewer2.org/)
colors <- colorRampPalette(brewer.pal(9, "GnBu"))(255)

pheatmap(
  sampleDistMatrix_chrono,
  main = "Chronodisrupted vs Control samples",
  show_colnames = FALSE,
  annotation = chronodisruption_col_data[,c("treatment", "type")],
  clustering_distance_rows=sampleDists_chrono,
  clustering_distance_cols=sampleDists_chrono,
  col=colors
)

pheatmap(
  sampleDistMatrix_stim,
  main = "L11GSPE vs L11VH",
  show_colnames = FALSE,
  annotation = stimulus_col_data[,c("treatment", "type")],
  clustering_distance_rows=sampleDists_stim,
  clustering_distance_cols=sampleDists_stim,
  col=colors
)

pheatmap(
  sampleDistMatrix_chrono_gspe,
  main = "Control vs L11GSPE",
  show_colnames = FALSE,
  annotation = chrono_gspe_col_data[,c("treatment", "type")],
  clustering_distance_rows=sampledists_chrono_gspe,
  clustering_distance_cols=sampledists_chrono_gspe,
  col=colors
)

## Differential expression analysis

# To properly compare control versus chronodisruption, we need to define that in our colData objects, 
# therefore we need to set the colData$type as a factor, and the firs level should be control 
chronodisruption_col_data$treatment <- factor(chronodisruption_col_data$treatment, levels= c("control", "L11-VH"))
stimulus_col_data$treatment <- factor(stimulated_col_data$treatment, levels= c("L11-VH", "L11-GSPE"))

# We now generate the new DESeq objects and run the differential expression analysis.
dds_chronodisruption <- DESeqDataSetFromMatrix(
  countData = chronodisruption_count_table,
  colData = chronodisruption_col_data,
  design = ~ treatment
)

dds_stim <- DESeqDataSetFromMatrix(
  countData = stimulus_count_table,
  colData = stimulus_col_data,
  design = ~ treatment
)


dds_chronodisruption <- DESeq(dds_chronodisruption)
res_chronodisruption <- results(dds_chronodisruption)
head(res_chronodisruption)

dds_stim <- DESeq(dds_stim)
res_stim <- results(dds_stim)
head(res_stim)

# Given the conditions we declared for differential expression analysis we can subset our list of 
de_chronodisruption <- res_chronodisruption[which(res_chronodisruption$padj <= 0.05),]
head(de_chronodisruption)

de_stim <- res_stim[which(res_stim$padj <= 0.05),]
head(de_stim)

# how many differentially expressed genes we got
nrow(de_chronodisruption)
nrow(de_stim)  # there are no DE genes in L11-VH vs L11-GPSE

# Vizualization Volcano plot
pData <- as.data.frame(res_chronodisruption[which(!is.na(res_chronodisruption$padj)),])
chronodisruption_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(aes(color= padj <= 0.05))+
  geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
  geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
  scale_color_d3()+
  ggtitle("Chronodisruption Volcano Plot")+
  theme_bw()

chronodisruption_Volcano

# sort
de_chronodisruption[order(de_chronodisruption$padj, decreasing = FALSE), ]

# plot disperssion estimates
plotDispEsts(dds_chronodisruption)

# The write.csv function will generate an Excel "friendly" file.
chrono_fileName <- "./scripts/results/diffExpGenes_chronodisruption.csv"
write.csv(de_chronodisruption, chrono_fileName)

# prettier volcano plot
pData$top10label <- NA
pData$top10label[order(pData$padj)[1:10]] <- rownames(pData)[order(pData$padj)[1:10]]

chronodisruption_Volcano <- ggplot(pData, aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(aes(color= padj <= 0.05))+
  geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
  geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
  scale_color_d3()+
  ggtitle("Chronodisruption Volcano Plot")+
  geom_text_repel(aes(label=top10label))+ ##add the lables in the top 10 
  theme_bw()+
  theme(legend.position = "bottom")

chronodisruption_Volcano

# Heatmap of differentially expressed genes
vst_chronodisruption <- assay(vst(dds_chronodisruption))
deGenes <- vst_chronodisruption[rownames(de_chronodisruption),]
pheatmap(deGenes, scale = "row", show_rownames = FALSE, main = "Differentially expressed genes chronodisruption")

# Boxplot of differentially expressed (top10) genes
top10Genes <- rownames(de_chronodisruption)[order(de_chronodisruption$padj)[1:10]]

# Using melt to literaly melt a wide data.frame into a long data.frame
pData <- melt(vst_chronodisruption[top10Genes,])
# Add your sample information. 
pData <- cbind(pData, chronodisruption_col_data[as.character(pData$Var2),])

top10_Plot <- ggplot(pData, aes(x= type, y= value))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(alpha=0.6)+
  facet_grid(~Var1, scale="free")+
  ylab("VST expression values of DE chronodisruption vs control genes")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

print(top10_Plot)

# Pathway analysis
chrono <- res_chronodisruption
summary(chrono)

# genes with NA for log2FoldChange
head(chrono[is.na(chrono$log2FoldChange), ])


# How many such genes?
table(chrono[is.na(chrono$log2FoldChange), ]$baseMean == 0)

# Data cleaning, remove unexpressed genes (baseMean = 0)
chrono <- chrono[!chrono$baseMean == 0, ]

library(purrr)  # functional map

split_gene_string <- function(gene_name) {
  # implicit return
  unlist(strsplit(unlist(strsplit(gene_name, "-"))[2], "\\|"))[1]
}	

# need to dynamically split the string and parse to get the HGNC gene symbols and
# then convert to ENTREZ IDs

