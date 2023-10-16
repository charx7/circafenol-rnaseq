# Commands pakcages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")  # for differential expression analysis

# other packages
install.packages("ggplot2") # for plotting
install.packages("ggsci")  # scientific journal color palettes
install.packages("pheatmap")  # heatmap
install.packages("ggsci")  # scientific style plots
install.packages("reshape2")  # dataframe transformations (melf func)
install.packages("purrr")  # map functionallity
install.packages("ggrepel")  # additions to ggplot
install.packages("VennDiagram")  # Venn Diagram
install.packages("gridExtra")
install.packages("readxl")  # to read excel files
