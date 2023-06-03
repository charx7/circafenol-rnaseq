install.packages("ggplot2") # for plotting
install.packages("ggsci")  # scientific journal color palettes
install.packages("pheatmap")  # heatmap

# Commands pakcages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")  # for differential expression analysis
