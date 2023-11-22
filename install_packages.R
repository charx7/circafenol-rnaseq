# Commands pakcages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")  # for differential expression analysis
BiocManager::install("KEGGREST")
#BiocManager::install("org.Hs.eg.db")  # Human gene annotation package
BiocManager::install("org.Rn.eg.db")  # rat gene annotation package
BiocManager::install("clusterProfiler") # Package for enrichment tests
BiocManager::install("topGO")  # Gene Ontology vizualization
BiocManager::install("Rgraphviz")  # graph viz

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
