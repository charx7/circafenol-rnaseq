library(readxl)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(gridExtra)

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

# get deseq datasets
dds_chronodisruption <- DESeqDataSetFromMatrix(
  countData = chronodisruption_count_table,
  colData = chronodisruption_col_data,
  design = ~ treatment
)

# Variance stabilizing transformations
vst_chronodisruption <- assay(vst(dds_chronodisruption, blind=FALSE))

## PCA
# To calculate the components by sample we need to transpose our matrix of normalized gene expression 
pcData <- prcomp(t(vst_chronodisruption))
pcVar <- summary(pcData)

# By getting the summary() of a prcomp object (in this case pcData) we can also obtain the total amount of variance explained by each of the components.
pcVar$importance

# We can then extract the variance explained by components 1 and 2. 
varPC1 <- pcVar$importance[2,1]
varPC2 <- pcVar$importance[2,2]

pcPlotData <- data.frame(pcData$x[,1:4], col_data[rownames(pcData$x),])

pcaPlot_chronodisruption <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=treatment))+
  geom_jitter(alpha=0.6)+
  #facet_grid(~stimulation, scales = "free")+
  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
  scale_color_aaas()+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 8))

pcaPlot_chronodisruption

# TODO: repeat the analysis with gspe vs vh

# Unsupervised clustering

# Again we need to transpose our matrix to then calculate the distance between each of the samples.
sampleDists <- dist(t(vst_chronodisruption))
sampleDistMatrix <- as.matrix(sampleDists)

# By using brewer.pal() we can generate a palette of colors, for more colors check (http://colorbrewer2.org/)
colors <- colorRampPalette(brewer.pal(9, "GnBu"))(255)


pheatmap(
  sampleDistMatrix,
  main = "Chronodisrupted vs Control samples",
  show_colnames = FALSE,
  annotation = chronodisruption_col_data[,c("treatment", "type")],
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors
)


## Differential expression analysis

# To properly compare control versus celiac, we need to define that in ourt colData objects, 
# therefore we need to set the colData$type as a factor, and the firs level should be control 
chronodisruption_col_data$treatment <- factor(chronodisruption_col_data$treatment, levels= c("control", "L11-VH"))


# We now generate the new DESeq objects and run the differential expression analysis.
dds_chronodisruption <- DESeqDataSetFromMatrix(
  countData = chronodisruption_count_table,
  colData = chronodisruption_col_data,
  design = ~ treatment
)

dds_chronodisruption <- DESeq(dds_chronodisruption)
res_chronodisruption <- results(dds_chronodisruption)
head(res_chronodisruption)

# Given the conditions we declared for differential expression analysis we can subset our list of 
de_chronodisruption <- res_chronodisruption [which(res_chronodisruption$padj <= 0.05),]

# how many differentially expressed genes we got
nrow(de_chronodisruption)

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

