library(readxl)
library(DESeq2)

# read the data
excel_path = "/home/rstudio/scripts/data/gene_keys.xlsx"
gene_names <- read_excel(
  excel_path,
  sheet = "Sheet1"
)

excel_path = "/home/rstudio/scripts/data/gene_count_matrix.xlsx"

count_data <- read_excel(
  excel_path,
  sheet = "count_data"
  
)
merged_df <- merge(gene_names, count_data, by = "gene_id")

columns <- names(count_data)
columns <- columns[2:22]

names <- c("gene name" , columns)
merged_df <- merged_df[names]
