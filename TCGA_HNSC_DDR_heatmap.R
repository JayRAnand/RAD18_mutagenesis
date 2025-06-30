library(here)
library(tidyverse)
library(TCGAbiolinks)
library(FirebrowseR)
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(readxl)
library(openxlsx)
library(AnnotationDbi)
library(clusterProfiler)
library(GOSemSim)
library(biomaRt)
library(pheatmap)
library(org.Hs.eg.db)
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(circlize)  # For color management
library(viridis)
library(RColorBrewer)
library(scales)

#using TCGA-HNSC sample downloaded using TCGAbiolinks
#assayNames(tcga_hnsc_data)

tcga_hnsc_fpkm_uq <- assay(tcga_hnsc_data, 'fpkm_uq_unstrand')
tcga_hnsc_fpkm_uq[1:10, 1:10]
# Save as a tab-separated .txt file
write.table(tcga_hnsc_fpkm_uq, file = "tcga_heatmap/tcga_hnsc_fpkm_uq.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# Save as a comma-separated .csv file
write.table(tcga_hnsc_fpkm_uq, file = "tcga_heatmap/tcga_hnsc_fpkm_uq.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)



#read file
hnsc_fpkm_uq <- read.table("tcga_heatmap/tcga_hnsc_fpkm_uq.txt", header = TRUE, sep = "\t", row.names = 1)
# Remove everything after and including the dot in the row names and assign back to matrix 
hnsc_fpkm_uq <- as.matrix(hnsc_fpkm_uq)
# Standardize hnsc_barcodes by replacing periods with hyphens
colnames(hnsc_fpkm_uq) <- gsub("\\.", "-", colnames(hnsc_fpkm_uq))
str(hnsc_fpkm_uq)
#ensemble ids
rownames(hnsc_fpkm_uq) <- sub("\\..*", "", rownames(hnsc_fpkm_uq))
hnsc_fpkm_uq[1:10,1:10]


#keep
#DDR
#remove duplicate FANCL gene
ddr_genes <- ddr_genes %>%
  group_by(Gene_Name) %>%
  filter(!(Gene_Name == "FANCL" & row_number() > 1)) %>%
  ungroup()

# Get the Ensembl Gene IDs from ddr_genes
ddr_ids_to_keep <- ddr_genes$Ensembl_Gene_ID
# Get the row names of hnsc_fpkm_uq
hnsc_ensembl_ids <- rownames(hnsc_fpkm_uq)
head(hnsc_ensembl_ids)
# Find the common Ensembl Gene IDs
common_ddr_ensembl_ids <- intersect(hnsc_ensembl_ids, ddr_ids_to_keep)
str(common_ddr_ensembl_ids)
# Filter hnsc_fpkm_uq to keep only the rows with common Ensembl Gene IDs
filtered_hnsc_fpkm_uq <- hnsc_fpkm_uq[common_ddr_ensembl_ids, ]
# Verify the result
str(filtered_hnsc_fpkm_uq)
#filter normal, HPV positive and HPV negative barcodes using coldata_filtered
# Get the barcodes from coldata_filtered
barcodes_to_keep <- coldata_filtered$barcode # Replace 'barcode' with the actual column name
str(barcodes_to_keep)
# Get the column names of filtered_hnsc_fpkm_uq
hnsc_barcodes <- colnames(filtered_hnsc_fpkm_uq)
str(hnsc_barcodes)
# Find the common barcodes
common_barcodes <- intersect(hnsc_barcodes, barcodes_to_keep)
str(common_barcodes)
# Filter filtered_hnsc_fpkm_uq to keep only the columns with common barcodes
filtered_hnsc_fpkm_uq_final <- filtered_hnsc_fpkm_uq[, common_barcodes]
# Verify the result
str(filtered_hnsc_fpkm_uq_final)
# Shorten column names to 15 digits
colnames(filtered_hnsc_fpkm_uq_final) <- substr(colnames(filtered_hnsc_fpkm_uq_final), 1, 15)
# Verify the changes
head(filtered_hnsc_fpkm_uq_final)
# Convert to data frame
filtered_hnsc_fpkm_uq_final_df <- as.data.frame(filtered_hnsc_fpkm_uq_final)
# Convert ENSEMBL IDs to gene symbols
gene_symbols_hnsc_ddr <- mapIds(org.Hs.eg.db,
                                keys = rownames(filtered_hnsc_fpkm_uq_final_df),
                                keytype = 'ENSEMBL',
                                column = 'SYMBOL')
# Add gene symbols as a column
filtered_hnsc_fpkm_uq_final_df$gene_symbols <- gene_symbols_hnsc_ddr
# Verify the changes
filtered_hnsc_fpkm_uq_final_df[1:5, 1:6] # Show gene symbols column
str(filtered_hnsc_fpkm_uq_final_df)
# Set rownames to be the values in the gene_symbols column
rownames(filtered_hnsc_fpkm_uq_final_df) <- filtered_hnsc_fpkm_uq_final_df$gene_symbols
# Remove the gene_symbols column from the data frame
filtered_hnsc_fpkm_uq_final_df$gene_symbols <- NULL
# Calculate Z-scores for each row (gene)
hnsc_ddr_z_score <- t(apply(filtered_hnsc_fpkm_uq_final_df, 1, function(x) (x - mean(x)) / sd(x)))
# Check the result
str(hnsc_ddr_z_score)
# Save the hnsc_ddr_z_score matrix as a .txt file
write.table(hnsc_ddr_z_score, file = "tcga_heatmap/hnsc_ddr_z_score.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)




#keep
#MAGEs
#load mage genes and ids
mages <- read.xlsx("mages_human.xlsx")
head(mages)
#filter DDR genes
# Get the Ensembl Gene IDs from ddr_genes
mage_ids_to_keep <- mages$ensemble_id
# Find the common Ensembl Gene IDs
common_mage_ensembl_ids <- intersect(hnsc_ensembl_ids, mage_ids_to_keep)
str(common_mage_ensembl_ids)
# Filter hnsc_fpkm_uq to keep only the rows with common Ensembl Gene IDs
filtered_hnsc_fpkm_uq_mages <- hnsc_fpkm_uq[common_mage_ensembl_ids, ]
# Verify the result
str(filtered_hnsc_fpkm_uq_mages)
#filter normal, HPV positive and HPV negative barcodes using coldata_filtered
#str(barcodes_to_keep)
# Get the column names of filtered_hnsc_fpkm_uq_mages
hnsc_barcodes_mages <- colnames(filtered_hnsc_fpkm_uq_mages)
str(hnsc_barcodes_mages)
# Find the common barcodes
common_barcodes_mages <- intersect(hnsc_barcodes_mages, barcodes_to_keep)
str(common_barcodes_mages)
# Filter filtered_hnsc_fpkm_uq to keep only the columns with common barcodes
filtered_hnsc_fpkm_uq_final_mages <- filtered_hnsc_fpkm_uq_mages[, common_barcodes_mages]
# Verify the result
str(filtered_hnsc_fpkm_uq_final_mages)
# Shorten column names to 15 digits
colnames(filtered_hnsc_fpkm_uq_final_mages) <- substr(colnames(filtered_hnsc_fpkm_uq_final_mages), 1, 15)
# Verify the changes
head(filtered_hnsc_fpkm_uq_final_mages)
# Convert to data frame
filtered_hnsc_fpkm_uq_final_mages_df <- as.data.frame(filtered_hnsc_fpkm_uq_final_mages)
# Convert ENSEMBL IDs to gene symbols
gene_symbols_hnsc_mages <- mapIds(org.Hs.eg.db,
                                keys = rownames(filtered_hnsc_fpkm_uq_final_mages_df),
                                keytype = 'ENSEMBL',
                                column = 'SYMBOL')
# Add gene symbols as a column
filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols <- gene_symbols_hnsc_mages
# Verify the changes
filtered_hnsc_fpkm_uq_final_mages_df[1:5, 1:6] # Show gene symbols column
cbind(rownames(filtered_hnsc_fpkm_uq_final_mages_df), filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols)
# Create a lookup table for gene_name based on ensemble_id from mages
ensemble_lookup <- mages[, c("ensemble_id_human", "gene_name")]
# Replace NA values in gene_symbols with corresponding gene_name based on matching ensemble_id
for (i in 1:nrow(filtered_hnsc_fpkm_uq_final_mages_df)) {
  if (is.na(filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols[i])) {
    # Get the ensemble_id from the row name (which is the row name itself in filtered_hnsc_fpkm_uq_final_mages_df)
    ensemble_id <- rownames(filtered_hnsc_fpkm_uq_final_mages_df)[i]
    
    # Find the corresponding gene_name from the mages dataframe using ensemble_id
    matching_gene_name <- ensemble_lookup$gene_name[ensemble_lookup$ensemble_id == ensemble_id]
    
    # Replace the NA value in gene_symbols with the corresponding gene_name
    if (length(matching_gene_name) > 0) {
      filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols[i] <- matching_gene_name
    } else {
      # Optional: print message if no match is found (useful for debugging)
      message("No match found for ensemble_id: ", ensemble_id)
    }
  }
}

# Check the result
cbind(rownames(filtered_hnsc_fpkm_uq_final_mages_df), filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols)
# Set rownames to be the values in the gene_symbols column
rownames(filtered_hnsc_fpkm_uq_final_mages_df) <- filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols
# Remove the gene_symbols column from the data frame
filtered_hnsc_fpkm_uq_final_mages_df$gene_symbols <- NULL
# Calculate Z-scores for each row (gene)
hnsc_mage_z_score <- t(apply(filtered_hnsc_fpkm_uq_final_mages_df, 1, function(x) (x - mean(x)) / sd(x)))
# Check the result
str(hnsc_mage_z_score)
# Save the hnsc_ddr_z_score matrix as a .txt file
write.table(hnsc_mage_z_score, file = "tcga_heatmap/hnsc_mage_z_score.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)





#keep
#TLSgenes
#load mage genes and ids
tlsgenes <- read.xlsx("TLSgenes.xlsx")
str(tlsgenes)
#filter DDR genes
# Get the Ensembl Gene IDs from tlsgenes
tls_ids_to_keep <- tlsgenes$ensembl_id_human
str(tls_ids_to_keep)
# Find the common Ensembl Gene IDs
common_tls_ensembl_ids <- intersect(hnsc_ensembl_ids, tls_ids_to_keep)
str(common_tls_ensembl_ids)
# Filter hnsc_fpkm_uq to keep only the rows with common Ensembl Gene IDs
filtered_hnsc_fpkm_uq_tls <- hnsc_fpkm_uq[common_tls_ensembl_ids, ]
# Verify the result
str(filtered_hnsc_fpkm_uq_tls)
#filter normal, HPV positive and HPV negative barcodes using coldata_filtered
#str(barcodes_to_keep)
# Get the column names of filtered_hnsc_fpkm_uq_mages
hnsc_barcodes_tls <- colnames(filtered_hnsc_fpkm_uq_tls)
str(hnsc_barcodes_tls)
# Find the common barcodes
common_barcodes_tls <- intersect(hnsc_barcodes_tls, barcodes_to_keep)
str(common_barcodes_mages)
# Filter filtered_hnsc_fpkm_uq to keep only the columns with common barcodes
filtered_hnsc_fpkm_uq_final_tls <- filtered_hnsc_fpkm_uq_tls [, common_barcodes_tls]
# Verify the result
str(filtered_hnsc_fpkm_uq_final_tls)
# Shorten column names to 15 digits
colnames(filtered_hnsc_fpkm_uq_final_tls) <- substr(colnames(filtered_hnsc_fpkm_uq_final_tls), 1, 15)
# Verify the changes
head(filtered_hnsc_fpkm_uq_final_tls)
# Convert to data frame
filtered_hnsc_fpkm_uq_final_tls_df <- as.data.frame(filtered_hnsc_fpkm_uq_final_tls)
# Convert ENSEMBL IDs to gene symbols
gene_symbols_hnsc_tls <- mapIds(org.Hs.eg.db,
                                  keys = rownames(filtered_hnsc_fpkm_uq_final_tls_df),
                                  keytype = 'ENSEMBL',
                                  column = 'SYMBOL')
# Add gene symbols as a column
filtered_hnsc_fpkm_uq_final_tls_df$gene_symbols <- gene_symbols_hnsc_tls
# Verify the changes
filtered_hnsc_fpkm_uq_final_tls_df[1:5, 1:6] # Show gene symbols column
cbind(rownames(filtered_hnsc_fpkm_uq_final_tls_df), filtered_hnsc_fpkm_uq_final_tls_df$gene_symbols)
filtered_hnsc_fpkm_uq_final_tls_df$gene_symbols
# Set rownames to be the values in the gene_symbols column
rownames(filtered_hnsc_fpkm_uq_final_tls_df) <- filtered_hnsc_fpkm_uq_final_tls_df$gene_symbols
# Remove the gene_symbols column from the data frame
filtered_hnsc_fpkm_uq_final_tls_df$gene_symbols <- NULL
# Calculate Z-scores for each row (gene)
hnsc_tls_z_score <- t(apply(filtered_hnsc_fpkm_uq_final_tls_df, 1, function(x) (x - mean(x)) / sd(x)))
# Check the result
str(hnsc_tls_z_score)
# Save the hnsc_ddr_z_score matrix as a .txt file
write.table(hnsc_tls_z_score, file = "tcga_heatmap/hnsc_tls_z_score.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)










#get TP53 mutation data
query_maf <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_maf)
maf <- GDCprepare(query_maf)
# Convert to maftools object
maf <- read.maf(maf)

# TP53 mutations
tp53_mutations <- maf@data %>%
  as_tibble() %>%
  filter(Hugo_Symbol == "TP53") %>%
  dplyr::select(Tumor_Sample_Barcode, Variant_Classification)

# Add Tumor_Sample_Barcode_short column and remove duplicates
tp53_mutations_filtered <- tp53_mutations %>%
  mutate(Tumor_Sample_Barcode_short = substr(Tumor_Sample_Barcode, 1, 15)) %>%
  distinct(Tumor_Sample_Barcode_short, .keep_all = TRUE)

# Verify the changes
head(tp53_mutations_filtered)
# 1. Shorten Barcodes (If Necessary)
shorten_barcode <- function(barcode) {
  return(substr(barcode, 1, 15))
}
coldata_filtered$barcode_short <- sapply(coldata_filtered$barcode, shorten_barcode)
tp53_mutations_filtered$Tumor_Sample_Barcode_short <- sapply(tp53_mutations_filtered$Tumor_Sample_Barcode, shorten_barcode)
# 2. Filter tp53_mutations for matches in coldata_filtered
matched_tp53_mutations <- tp53_mutations_filtered %>%
  filter(Tumor_Sample_Barcode_short %in% coldata_filtered$barcode_short)
# 3. Join the DataFrames
coldata_filtered_heatmap <- coldata_filtered %>%
  left_join(matched_tp53_mutations, by = c("barcode_short" = "Tumor_Sample_Barcode_short"))
coldata_filtered_heatmap$Tumor_Sample_Barcode <- NULL
#add new column for TP53.m
coldata_filtered_heatmap <- coldata_filtered_heatmap %>%
  mutate(TP53.m = ifelse(is.na(Variant_Classification), "no", "yes"))



# PIK3CA mutations
pik3ca_mutations <- maf@data %>%
  as_tibble() %>%
  filter(Hugo_Symbol == "PIK3CA") %>%
  dplyr::select(Tumor_Sample_Barcode, Variant_Classification)
# Add Tumor_Sample_Barcode_short column
pik3ca_mutations_filtered <- pik3ca_mutations %>%
  mutate(Tumor_Sample_Barcode_short = substr(Tumor_Sample_Barcode, 1, 15)) %>%
  distinct(Tumor_Sample_Barcode_short, .keep_all = TRUE)
# 1. Shorten Barcodes (If Necessary)
shorten_barcode <- function(barcode) {
  return(substr(barcode, 1, 15))
}
pik3ca_mutations_filtered$Tumor_Sample_Barcode_short <- sapply(pik3ca_mutations_filtered$Tumor_Sample_Barcode, shorten_barcode)
# 2. Filter tp53_mutations for matches in coldata_filtered
matched_pik3ca_mutations <- pik3ca_mutations_filtered %>%
  filter(Tumor_Sample_Barcode_short %in% coldata_filtered$barcode_short)
head(matched_pik3ca_mutations)
matched_pik3ca_mutations <- matched_pik3ca_mutations %>%
  rename(PIK3CA_variant = Variant_Classification)
# Verify the change
head(matched_pik3ca_mutations)
# 3. Join the DataFrames
coldata_filtered_heatmap <- coldata_filtered_heatmap %>%
  left_join(matched_pik3ca_mutations, by = c("barcode_short" = "Tumor_Sample_Barcode_short"))
coldata_filtered_heatmap$Tumor_Sample_Barcode <- NULL
#add new column for TP53.m
coldata_filtered_heatmap <- coldata_filtered_heatmap %>%
  mutate(PIK3CA.m = ifelse(is.na(PIK3CA_variant), "no", "yes"))






#delete
#copy number data
# Query for CNV data, filtering for both primary tumor and solid tissue normal
query_cnv_normal_tumor <- GDCquery(
  project = 'TCGA-HNSC',
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_cnv_normal_tumor)
cnv_normal_tumor <- GDCprepare(query_cnv_normal_tumor)
# 1. Extract the Copy Number Matrix
cnv_matrix <- assay(cnv_normal_tumor, "copy_number")
rownames(cnv_matrix) <- sub("\\..*", "", rownames(cnv_matrix))
cnv_matrix[1:10,1:10]
# 2. Filter for Genes of Interest
genes_of_interest_names <- c("CDKN2A", "CCND1", "PIK3CA")
genes_of_interest <- c("ENSG00000147889", "ENSG00000110092", "ENSG00000121879")
gene_filter <- rownames(cnv_matrix) %in% genes_of_interest
cnv_matrix_genes <- cnv_matrix[gene_filter, ]
cnv_matrix_genes[1:3, 1:3]
# Adjusted shortening function
shorten_barcode <- function(barcode) {
  return(substr(barcode, 1, 15))
}
# Shorten common_barcodes
common_barcodes_short <- sapply(common_barcodes, shorten_barcode)
# Shorten colnames(cnv_matrix_genes)
colnames_cnv_short <- sapply(colnames(cnv_matrix_genes), shorten_barcode)
# Filter for Samples of Interest
sample_filter <- colnames_cnv_short %in% common_barcodes_short
cnv_matrix_genes_samples <- cnv_matrix_genes[, sample_filter]
# Check the results
head(cnv_matrix_genes_samples)
# 4. Transpose the Matrix (Optional)
cnv_matrix_final <- t(cnv_matrix_genes_samples)
# Now cnv_matrix_final contains the copy number data for the genes and samples of interest.
head(cnv_matrix_final)
# Change column names
colnames(cnv_matrix_final) <- genes_of_interest_names
# Shorten row names of cnv_matrix_final
rownames(cnv_matrix_final) <- substr(rownames(cnv_matrix_final), 1, 15)
head(cnv_matrix_final)
# Convert cnv_matrix_final to a data frame with barcode as a column
cnv_df <- as.data.frame(cnv_matrix_final) %>%
  mutate(barcode_short = rownames(cnv_matrix_final))
# Filter cnv_df for barcodes present in coldata_filtered_heatmap$barcode_short
cnv_filtered <- cnv_df %>%
  filter(barcode_short %in% coldata_filtered_heatmap$barcode_short)
# Join cnv_filtered with coldata_filtered_heatmap
coldata_cnv <- coldata_filtered_heatmap %>%
  left_join(cnv_filtered, by = "barcode_short")
# Verify the joined data
head(coldata_cnv$PIK3CA)
#number of rows in coldata_cnv that has sample_type value of 'Primary Tumor'
#and have NA in either of these columns CDKN2A, CCND1, and PIK3CA
nrow(coldata_cnv %>% filter(sample_type == "Primary Tumor" & (is.na(CDKN2A) | is.na(CCND1) | is.na(PIK3CA))))
#58 samples
#Assign values of 1 if values is more than 2, 0 if 2, -1 if less than 2, and 3 if NA
#in columns "CDKN2A", "CCND1", "PIK3CA" 
#to new columns "CDKN2A_CN", "CCND1_CN", "PIK3CA_CN" for coldata_cnv
coldata_cnv <- coldata_cnv %>%
  mutate(
    CDKN2A_CN = case_when(
      is.na(CDKN2A) ~ 3, #NA
      CDKN2A > 2 ~ 1, #amplification
      CDKN2A == 2 ~ 0, #neutral
      CDKN2A < 2 ~ -1, #deletion
      TRUE ~ NA_real_  # Handle any unexpected cases
    ),
    CCND1_CN = case_when(
      is.na(CCND1) ~ 3, #NA
      CCND1 > 2 ~ 1, #amplification
      CCND1 == 2 ~ 0, #neutral
      CCND1 < 2 ~ -1, #deletion
      TRUE ~ NA_real_ 
    ),
    PIK3CA_CN = case_when(
      is.na(PIK3CA) ~ 3, #NA
      PIK3CA > 2 ~ 1, #amplification
      PIK3CA == 2 ~ 0, #neutral
      PIK3CA < 2 ~ -1, #deletion
      TRUE ~ NA_real_
    )
  )
# Verify the changes
head(coldata_cnv[, c("barcode_short", "HPV_status", "TP53.m", "PIK3CA.m", "CDKN2A_CN", "PIK3CA_CN", "CCND1_CN")])
# Select the desired columns
annotations_hm <- coldata_cnv[, c("barcode_short", "HPV_status","TP53.m", "PIK3CA.m", "CDKN2A_CN", "PIK3CA_CN", "CCND1_CN")]
head(annotations_hm)
# Save annotations_hm as a tab-separated .txt file
write.table(annotations_hm, file = "tcga_heatmap/annotations_hm.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#delete
# Define the TCGA project/disease of interest (HNSC in this case)
# Download GISTIC thresholded data
gistic_data <- getGistic(disease = "HNSC", type = "thresholded")
# If you want the raw GISTIC data (instead of thresholded)
# gistic_data_raw <- getGistic(disease, type = "data")
# Check the first few rows directly
gistic_data[1:3, 1:10]
#gistic_data contains:
#Gene Symbol: Gene name
#Locus ID: Gene locus identifier
#Cytoband: Cytogenetic band
#sample columns with GISTIC values:
#-2 → Homozygous deletion
#-1 → Hemizygous deletion
#0 → Neutral copy number
#1 → Low-level amplification
#2 → High-level amplification
# Define genes of interest
genes_of_interest_names <- c("CDKN2A", "CCND1", "PIK3CA")
# Filter rows for these genes
# Filter rows for the genes of interest
gistic_genes <- gistic_data[gistic_data$`Gene Symbol` %in% genes_of_interest_names, ]
# Check the first few rows to ensure the filter worked
gistic_genes[1:3, 1:5]
# Transpose the matrix so genes are columns and samples are rows
gistic_genes_transposed <- t(gistic_genes[, 4:ncol(gistic_genes)])
# Create a data frame with genes as rows and samples as columns
gistic_genes_final <- as.data.frame(gistic_genes_transposed)
# Shorten sample barcodes if needed
colnames(gistic_genes_final) <- substr(colnames(gistic_genes_final), 1, 15)
# View the transposed data frame
head(gistic_genes_final)
# Assign correct gene names to columns
colnames(gistic_genes_final) <- genes_of_interest_names
# Check the result
head(gistic_genes_final)
# Shorten rownames to 15 characters
rownames(gistic_genes_final) <- substr(rownames(gistic_genes_final), 1, 15)
# Check the final matrix
head(gistic_genes_final)
# Define function to categorize copy number
categorize_copy_number <- function(cnv_value) {
  if (cnv_value < 0) {
    return("Deletion")  # Deletion
  } else if (cnv_value > 0) {
    return("Amplification")  # Amplification
  } else if (cnv_value == 0) {
    return("Neutral")  # Neutral
  } else {
    return(NA)  # NA
  }
}
# Apply the function to all values in the matrix
gistic_genes_final_categorized <- as.data.frame(apply(gistic_genes_final, 2, function(x) sapply(x, categorize_copy_number)))
# Check the results
head(gistic_genes_final_categorized)
# Add shortened barcodes as a column
gistic_genes_final_categorized$barcode_short <- rownames(gistic_genes_final_categorized)
gistic_filtered <- gistic_genes_final_categorized %>%
  filter(barcode_short %in% coldata_filtered_heatmap$barcode_short)
# Join cnv_filtered with coldata_filtered_heatmap
coldata_cnv <- coldata_filtered_heatmap %>%
  left_join(gistic_filtered, by = "barcode_short")
# Verify the joined data
head(coldata_cnv$PIK3CA)
# Verify the changes
head(coldata_cnv[, c("barcode_short", "HPV_status", "TP53.m", "PIK3CA.m", "CDKN2A", "PIK3CA", "CCND1")])
# Select the desired columns
annotations_hm <- coldata_cnv[, c("barcode_short", "HPV_status","TP53.m", "PIK3CA.m", "CDKN2A", "PIK3CA", "CCND1")]
head(annotations_hm)
# Check if HPV_status is "normal" and if so, set ajcc_pathologic_stage to NA
annotations_hm$ajcc_pathologic_stage[annotations_hm$HPV_status == "normal"] <- NA
head(annotations_hm)
# Save annotations_hm as a tab-separated .txt file
write.table(annotations_hm, file = "tcga_heatmap/annotations_hm.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#keep
#cBioportal CNA
cna_cdkn2a <- read.delim("cbioportal/cna_cdkn2a.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cna_pik3ca <- read.delim("cbioportal/cna_pik3ca.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cna_ccnd1 <- read.delim("cbioportal/cna_ccnd1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Merge the datasets sequentially based on SAMPLE_ID
cna_merged <- merge(cna_ccnd1[, c("SAMPLE_ID", "CCND1")], 
                     cna_pik3ca[, c("SAMPLE_ID", "PIK3CA")], 
                     by = "SAMPLE_ID")
cna_merged <- merge(cna_merged, 
                     cna_cdkn2a[, c("SAMPLE_ID", "CDKN2A")], 
                     by = "SAMPLE_ID")
head(cna_merged)
# Define function to categorize copy number
categorize_cna <- function(cna) {
  if (cna < 0) {
    return("Deletion")  # Deletion
  } else if (cna > 0) {
    return("Amplification")  # Amplification
  } else if (cna == 0) {
    return("Neutral")  # Neutral
  } else {
    return(NA)  # NA
  }
}
# Apply the function to all relevant columns
categorized_cna <- cna_merged
# Apply to CCND1, PIK3CA, and CDKN2A
categorized_cna$CCND1 <- sapply(cna_merged$CCND1, categorize_cna)
categorized_cna$PIK3CA <- sapply(cna_merged$PIK3CA, categorize_cna)
categorized_cna$CDKN2A <- sapply(cna_merged$CDKN2A, categorize_cna)
head(categorized_cna)
# Merge coldata_filtered_heatmap and categorized_cna by barcode_short and SAMPLE_ID
coldata_cna <- merge(coldata_filtered_heatmap, 
                     categorized_cna, 
                     by.x = "barcode_short", 
                     by.y = "SAMPLE_ID", 
                     all.x = TRUE)

# Verify the changes
head(coldata_cna[, c("barcode_short", "HPV_status", "TP53.m", "PIK3CA.m", "CDKN2A", "PIK3CA", "CCND1")])

#keep
#cBioportal mutations
mutations_tp53 <- read.delim("cbioportal/mutations_tp53.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mutations_pik3ca <- read.delim("cbioportal/mutations_pik3ca.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge the datasets sequentially based on SAMPLE_ID
mut_merged <- merge(mutations_tp53[, c("SAMPLE_ID", "TP53")], 
                    mutations_pik3ca[, c("SAMPLE_ID", "PIK3CA")], 
                    by = "SAMPLE_ID")
head(mut_merged)
any(is.na(mut_merged))
# Define the function correctly
categorize_mut <- function(mut) {
  if (mut == "WT") {
    return("Wild Type")
  } else {
    return("Mutant")
  }
}
# Apply the function to all relevant columns
categorized_mut <- mut_merged
# Apply to CCND1, PIK3CA, and CDKN2A
categorized_mut$TP53 <- sapply(mut_merged$TP53, categorize_mut)
categorized_mut$PIK3CA <- sapply(mut_merged$PIK3CA, categorize_mut)
head(categorized_mut)
# Rename the columns in categorized_mut
colnames(categorized_mut)[colnames(categorized_mut) == "TP53"] <- "TP53.status"
colnames(categorized_mut)[colnames(categorized_mut) == "PIK3CA"] <- "PIK3CA.status"
head(categorized_mut)
# Merge coldata_filtered_heatmap and categorized_cna by barcode_short and SAMPLE_ID
coldata_hm <- merge(coldata_cna, 
                    categorized_mut, 
                     by.x = "barcode_short", 
                     by.y = "SAMPLE_ID", 
                     all.x = TRUE)

# Verify the changes
head(coldata_hm[, c("barcode_short", "HPV_status", "TP53.m", "TP53.status", "PIK3CA.m", "PIK3CA.status", "CDKN2A", "PIK3CA", "CCND1")])
# Select the desired columns
annotations_hm <- coldata_hm[, c("barcode_short", "HPV_status", "ajcc_pathologic_stage", "TP53.status", "PIK3CA.status", "CDKN2A", "PIK3CA", "CCND1")]
head(annotations_hm)
# Save annotations_hm as a tab-separated .txt file
write.table(annotations_hm, file = "tcga_heatmap/annotations_hm.txt", sep = "\t", row.names = FALSE, quote = FALSE)












#keep
str(annotations_hm)
str(hnsc_ddr_z_score)
str(hnsc_mage_z_score)

#keep
#HEATMAP MAGEs
# Are barcodes in same order?
if (identical(colnames(hnsc_mage_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}
#if not ordered in match then reorder to match
# Reorder the columns of hnsc_mage_z_score
hnsc_mage_z_score <- hnsc_mage_z_score[, match(annotations_hm$barcode_short, colnames(hnsc_mage_z_score))]
# Are barcodes in same order?

if (identical(colnames(hnsc_mage_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}
# Find which genes in mages$gene_name are present in rownames(hnsc_mage_z_score)
common_mages <- mages$gene_name[mages$gene_name %in% rownames(hnsc_mage_z_score)]
# Reorder rows of hnsc_tls_z_score based on tlsgenes$genes
hnsc_mage_z_score <- hnsc_mage_z_score[common_mages, , drop = FALSE]
# Convert HPV_status to a factor with the desired order
#annotations_hm$HPV_status <- factor(annotations_hm$HPV_status, levels = c("normal", "negative", "positive"))
# Convert HPV_status and Pathologic_Stage to factors with desired orders
annotations_hm$HPV_status <- factor(annotations_hm$HPV_status, levels = c("normal", "negative", "positive"))
annotations_hm$ajcc_pathologic_stage <- factor(annotations_hm$ajcc_pathologic_stage,
                                               levels = c("NA", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IVA", "Stage IVB", "Stage IVC"))

# Define the annotation
ha <- HeatmapAnnotation(
  HPV_status = annotations_hm$HPV_status,
  TP53_m = annotations_hm$TP53.status,
  PIK3CA_m = annotations_hm$PIK3CA.status,
  CDKN2A_CN = annotations_hm$CDKN2A,
  PIK3CA_CN = annotations_hm$PIK3CA,
  CCND1_CN = annotations_hm$CCND1,
  Pathologic_Stage = annotations_hm$ajcc_pathologic_stage,
  col = list(
    HPV_status = c("normal" = "gray", "negative" = "blue", "positive" = "red"),
    Pathologic_Stage = c(
      "NA" = "gray",
      "Stage 0" = "#a6cee3",
      "Stage I" = "#1f78b4",
      "Stage II" = "#6a3d9a",
      "Stage III" = "#33a02c",
      "Stage IVA" = "#e31a1c",
      "Stage IVB" = "#ff7f00",
      "Stage IVC" = "#fdbf6f"
    ),
    TP53_m = c("Wild Type" = "lightgreen", "Mutant" = "darkgreen"),
    PIK3CA_m = c("Wild Type" = "lightblue", "Mutant" = "darkblue"),
    CDKN2A_CN = c("Amplification" = "red", "Neutral" = "green", "Deletion" = "blue"),
    PIK3CA_CN = c("Amplification" = "red", "Neutral" = "green", "Deletion" = "blue"),
    CCND1_CN = c("Amplification" = "red", "Neutral" = "green", "Deletion" = "blue")
  )
)


#because of long range and most data point scewed to low values. I am using two SD as a range for scale
# Calculate mean and standard deviation
mean_value_mage <- mean(hnsc_mage_z_score)
sd_value_mage <- sd(hnsc_mage_z_score)
# Calculate two standard deviations from the mean
lower_bound_mage <- mean_value_mage - 2 * sd_value_mage
upper_bound_mage <- mean_value_mage + 2 * sd_value_mage
print(lower_bound_mage)
print(upper_bound_mage)
range(hnsc_mage_z_score)
#color scale for z-score -1, 0, 2
#A z-score difference of around 1 or more is often considered indicative 
#of a notable difference in gene expression between two samples.
#A z-score of ±1 indicates that the expression level of a gene in a sample is 1 standard deviation 
#away from the mean expression of that gene across the entire set of samples.
#Generally, a difference of ±1 SD (or a z-score difference of 2 or more between two samples) 
#can be considered a meaningful difference. This is because, in a normal distribution:
#About 68% of the data lies within ±1 SD from the mean.
#About 95% of the data lies within ±2 SD from the mean.
#Differences of 2 SD (i.e., z-score differences of 2 or more) 
#are typically considered statistically significant or large.
ht <- Heatmap(
  hnsc_mage_z_score,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = FALSE, # Disable row clustering
  row_names_gp = gpar(fontsize = 8),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status) # Explicitly set column order
)
# Create the heatmap with column splitting and ordering
ht_mage <- Heatmap(
  hnsc_mage_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = FALSE, # Disable row clustering
  row_names_gp = gpar(fontsize = 10),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage) # Order within each HPV status
)
# Draw the heatmap
draw(ht_mage)

#keep
#pdf
pdf("tcga_heatmap/heatmap_tcga_mages_final_path.pdf", width = 11, height = 8, useDingbats = FALSE) # Adjust width and height as needed
draw(ht_mage)
dev.off() # Close the graphics device





#keep everything below
#HEATMAP DDR
# Are barcodes in same order?

if (identical(colnames(hnsc_ddr_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}
#if not ordered in match then reorder to match
# Reorder the columns of hnsc_mage_z_score
hnsc_ddr_z_score <- hnsc_ddr_z_score[, match(annotations_hm$barcode_short, colnames(hnsc_ddr_z_score))]
# Are barcodes in same order?

if (identical(colnames(hnsc_ddr_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}

# Create the heatmap with explicit column order
ht_ddr <- Heatmap(
  hnsc_ddr_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Disable row clustering
  row_names_gp = gpar(fontsize = 6),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage) # Order within each HPV status
)
# Draw the heatmap
draw(ht_ddr)

pdf("tcga_heatmap/heatmap_tcga_ddr_final_path.pdf", width = 12, height = 15) # Adjust width and height as needed
draw(ht_ddr)
dev.off() # Close the graphics device






#keep everything below
#HEATMAP TLS
# Are barcodes in same order?

if (identical(colnames(hnsc_tls_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}
#if not ordered in match then reorder to match
# Reorder the columns of hnsc_mage_z_score
hnsc_tls_z_score <- hnsc_tls_z_score[, match(annotations_hm$barcode_short, colnames(hnsc_tls_z_score))]
# Are barcodes in same order?

if (identical(colnames(hnsc_tls_z_score), annotations_hm$barcode_short)) {
  print("The column names and barcode_short are in the same order.")
} else {
  print("The column names and barcode_short are NOT in the same order.")
}
# Reorder rows of hnsc_tls_z_score based on tlsgenes$genes
hnsc_tls_z_score <- hnsc_tls_z_score[tlsgenes$genes, , drop = FALSE]

# Create the heatmap with explicit column order
ht_tls <- Heatmap(
  hnsc_tls_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = FALSE, # Disable row clustering
  row_names_gp = gpar(fontsize = 10),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage) # Order within each HPV status
)
# Draw the heatmap
draw(ht_tls)

pdf("tcga_heatmap/heatmap_tcga_tls_final_path.pdf", width = 10, height = 5, useDingbats = FALSE) # Adjust width and height as needed
draw(ht_tls)
dev.off() # Close the graphics device




#experiment
#clean up ddr gene list
ddr_genes_tcga <- ddr_genes_acute
# Create a mapping between Gene_Name and DDR_pathway
gene_pathway_mapping <- ddr_genes_tcga %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  dplyr::select(Gene_Name, DDR_pathway)

# Get all genes in the heatmap
all_genes_heatmap <- rownames(hnsc_ddr_z_score)

# Create a data frame for all genes in the heatmap
all_genes_df <- tibble(Gene_Name = all_genes_heatmap)

# Merge pathway information with all genes in the heatmap
row_annotation_data_full <- all_genes_df %>%
  left_join(gene_pathway_mapping, by = "Gene_Name")

# Ensure the order matches the heatmap rows
row_annotation_data_full <- row_annotation_data_full %>%
  arrange(factor(Gene_Name, levels = all_genes_heatmap))

# Create the row annotation
row_ha <- rowAnnotation(
  DDR_pathway = row_annotation_data_full$DDR_pathway,
  col = list(
    DDR_pathway = set_names(scales::hue_pal()(nlevels(factor(row_annotation_data_full$DDR_pathway))),
                            levels(factor(row_annotation_data_full$DDR_pathway)))
  ),
  show_annotation_name = TRUE
)

# Assuming you have the 'annotations_hm' data frame and the HeatmapAnnotation 'ha' defined as in your original code

# Create the heatmap with row annotations
ht_ddr_pathway <- Heatmap(
  hnsc_ddr_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  left_annotation = row_ha, # Add the row annotation here
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Disable row clustering
  row_names_gp = gpar(fontsize = 6),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage) # Order within each HPV status
)

# Draw the heatmap
draw(ht_ddr_pathway)





#experiment 2
# Create a mapping between Gene_Name and DDR_pathway
gene_pathway_mapping <- ddr_genes_tcga %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  dplyr::select(Gene_Name, DDR_pathway)

# Get all genes in the heatmap
all_genes_heatmap <- rownames(hnsc_ddr_z_score)

# Create a data frame for all genes in the heatmap
all_genes_df <- tibble(Gene_Name = all_genes_heatmap)

# Merge pathway information with all genes in the heatmap
row_annotation_data_full <- all_genes_df %>%
  left_join(gene_pathway_mapping, by = "Gene_Name") %>%
  arrange(factor(Gene_Name, levels = all_genes_heatmap))

# Get unique DDR pathways
unique_pathways <- unique(ddr_genes_tcga$DDR_pathway)

# Create a list to store row annotation columns
annotation_list <- list()

# Create individual binary classifiers for each pathway
for (pathway in unique_pathways) {
  # Sanitize pathway names for use as column names (remove spaces and special characters)
  column_name <- gsub("[^[:alnum:]]", ".", pathway)
  annotation_list[[column_name]] <- ifelse(row_annotation_data_full$DDR_pathway == pathway, "In Pathway", "Not in Pathway")
  annotation_list[[column_name]] <- factor(annotation_list[[column_name]], levels = c("Not in Pathway", "In Pathway"))
}

# Define colors for the binary annotation
pathway_colors <- c("Not in Pathway" = "lightgray", "In Pathway" = "purple")

# Create the row annotation
row_ha <- rowAnnotation(
  df = as.data.frame(annotation_list),
  col = setNames(rep(list(pathway_colors), length(annotation_list)), names(annotation_list)),
  show_annotation_name = TRUE
  # Removed annotation_name_side here
)

# Assuming you have the 'annotations_hm' data frame and the HeatmapAnnotation 'ha' defined as in your original code

# Create the heatmap with row annotations
ht_ddr_pathway <- Heatmap(
  hnsc_ddr_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  left_annotation = row_ha, # Add the row annotation here
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Disable row clustering
  row_names_gp = gpar(fontsize = 6),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage) # Order within each HPV status
)

# Draw the heatmap
draw(ht_ddr_pathway)








#expwriment 3
# Create a mapping between Gene_Name and DDR_pathway
gene_pathway_mapping <- ddr_genes_tcga %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  dplyr::select(Gene_Name, DDR_pathway)

# Get all genes in the heatmap
all_genes_heatmap <- rownames(hnsc_ddr_z_score)

# Create a data frame for all genes in the heatmap
all_genes_df <- tibble(Gene_Name = all_genes_heatmap)

# Merge pathway information with all genes in the heatmap
row_annotation_data_full <- all_genes_df %>%
  left_join(gene_pathway_mapping, by = "Gene_Name") %>%
  arrange(factor(Gene_Name, levels = all_genes_heatmap))

# Get unique DDR pathways
unique_pathways <- unique(ddr_genes_tcga$DDR_pathway)

# Create a list to store row annotation columns
annotation_list <- list()

# Create individual binary classifiers for each pathway
for (pathway in unique_pathways) {
  # Sanitize pathway names for use as column names (remove spaces and special characters)
  column_name <- gsub("[^[:alnum:]]", ".", pathway)
  annotation_list[[column_name]] <- ifelse(row_annotation_data_full$DDR_pathway == pathway, "In Pathway", "Not in Pathway")
  annotation_list[[column_name]] <- factor(annotation_list[[column_name]], levels = c("Not in Pathway", "In Pathway"))
}

# Define colors for the binary annotation
pathway_colors <- c("Not in Pathway" = "lightgray", "In Pathway" = "purple")

# Create the row annotation
row_ha <- rowAnnotation(
  df = as.data.frame(annotation_list),
  col = setNames(rep(list(pathway_colors), length(annotation_list)), names(annotation_list)),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
)

# Assuming you have the 'annotations_hm' data frame and the HeatmapAnnotation 'ha' defined as in your original code

# Create the heatmap with row annotations and scaling adjustments
ht_ddr_pathway <- Heatmap(
  hnsc_ddr_z_score,
  name = "Z-score",
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  left_annotation = row_ha, # Add the row annotation here
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Disable row clustering
  row_names_gp = gpar(fontsize = 6),
  column_split = annotations_hm$HPV_status,
  column_order = order(annotations_hm$HPV_status, annotations_hm$ajcc_pathologic_stage), # Order within each HPV status
  width = unit(10, "cm"),  # Adjust width of the heatmap body
  height = unit(20, "cm"), # Adjust height of the heatmap body
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 8)
)

# Draw the heatmap
draw(ht_ddr_pathway,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     row_dend_width = unit(1, "cm")
)
