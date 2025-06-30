library(here)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(openxlsx)
library(circlize)
library(RColorBrewer)
library(org.Mm.eg.db)
library(readxl)
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")


#load the data
# Load count matrix
count_matrix_acute <- read.delim("Scott/counts_W_annot_AdultTg_210726.txt", 
                                 comment.char = "#",    # Skip the header lines starting with '#'
                                 header = TRUE,         # First non-comment line is the header
                                 sep = "\t",            # Tab-separated file
                                 stringsAsFactors = FALSE)
#load meatada
samples_acute <- read_xlsx("Scott/metadata.xlsx")

# Create a named vector for renaming
name_map_acute <- setNames(samples_acute$Sample_name, samples_acute$sample_id)

# Filter and rename columns in count_matrix_acute
filtered_cma <- count_matrix_acute %>%
  select(Geneid, one_of(samples_acute$sample_id)) %>%   # Filter columns
  rename_at(vars(-Geneid), ~ name_map_acute[.])               # Rename columns

# View the filtered and renamed data
View(filtered_cma)

# Extract the sample IDs
sample_acute_ids <- samples_acute$sample_id
# Specify the columns to retain (metadata columns + matching sample columns)
columns_to_keep_acute <- c(colnames(count_matrix_acute)[1:6], sample_acute_ids)
# Filter the count matrix
filtered_cma_match <- count_matrix_acute[, colnames(count_matrix_acute) %in% columns_to_keep_acute]
# Check the result
head(filtered_cma)



# load the ddr_genes list (already loaded in the project as ddr_genes)
ddr_genes_formatted <- read_csv("ddr_genes_formatted.csv")
ddr_genes_acute <- ddr_genes_formatted

# Convert Gene Names to uppercase and trim whitespaces for consistent matching
#filtered_cma$Geneid <- trimws(toupper(filtered_cma$Geneid))
#ddr_genes_acute$Gene_Name <- trimws(toupper(ddr_genes_acute$Gene_Name))

#ddr_genes_acute <- ddr_genes
#gene_mapping_acute <- data.frame(
  #column1 = c("XRCC7 (PRKDC)", "MRE11", "53BP1", "ABRAXAS1", "XPB (ERCC3)", "XPD (ERCC2)", "XPG (ERCC5)"),
  #column2 = c("PRKDC", "MRE11A", "TRP53BP1", "FAM175A", "ERCC3", "ERCC2", "ERCC5")
#)


# Replace names in ddr_genes_acute$Gene_Name using the mapping
#ddr_genes_acute$Gene_Name <- ifelse(
  #ddr_genes_acute$Gene_Name %in% gene_mapping_acute$column1,
  #gene_mapping_acute$column2[match(ddr_genes_acute$Gene_Name, gene_mapping_acute$column1)],
  #ddr_genes_acute$Gene_Name
#)
#write_csv(ddr_genes_acute, "ddr_genes_formatted.csv")
#ddr_genes_formatted <- read_csv("ddr_genes_formatted.csv")


# Filter for DDR genes using Geneid column
filtered_cma_ddr <- filtered_cma %>%
  filter(Geneid %in% ddr_genes_acute$Gene_Name)
head(filtered_cma_ddr)
#replace CCDC111 with PRIMPOL
filtered_cma_ddr$Geneid <- gsub("CCDC111", "PRIMPOL", filtered_cma_ddr$Geneid)
"CCDC111" %in% filtered_cma_ddr$Geneid
"PRIMPOL" %in% filtered_cma_ddr$Geneid

# Identify missing genes
#missing_genes_acute <- setdiff(ddr_genes$Gene_Name, filtered_cma_ddr$Geneid)
# View the missing genes
#print(missing_genes_acute)
# Optional: Save to a CSV
#write.csv(missing_genes_acute, "Scott/missing_ddr_genes.csv", row.names = FALSE)


# Set Geneid as rownames
rownames(filtered_cma_ddr) <- filtered_cma_ddr$Geneid
# Remove the Geneid column
filtered_cma_ddr$Geneid <- NULL
# View the updated data frame
head(filtered_cma_ddr)


# Log transformation
normalized_cma_ddr <- log2(filtered_cma_ddr + 1)
# Calculate Z-scores for each row (gene)
normalized_cma_ddr <- t(apply(normalized_cma_ddr, 1, function(x) (x - mean(x)) / sd(x)))


#HEATMAP TLS genes
tlsgenes <- read.xlsx("TLSgenes.xlsx")
tlsgenes
# Filter the data for TLS genes
tls_gene_names <- tlsgenes$genes
tls_gene_names
filtered_cma_tls <- normalized_cma_ddr[rownames(normalized_cma_ddr) %in% tls_gene_names, ]
filtered_cma_tls
# ensure consistent case
rownames(normalized_cma_ddr) <- trimws(toupper(rownames(normalized_cma_ddr)))

# Create the annotation for grouping
group_annotation_acute <- samples_acute$`4NQO`
names(group_annotation_acute) <- colnames(filtered_cma_tls)
# Ensure the group factor levels are set according to the desired order
group_annotation_acute <- factor(group_annotation_acute, levels = c("NONE", "2d", "14d"))
# Define the heatmap annotation
column_acute <- HeatmapAnnotation(Group = group_annotation_acute, 
                               col = list(Group = c("NONE" = "gray", "2d" = "blue", "14d" = "red")))

# Create logical indices for each group
none_cols <- group_annotation_acute == "NONE"
two_d_cols <- group_annotation_acute == "2d"
fourteen_d_cols <- group_annotation_acute == "14d"

# Order the columns based on the logical indices
ordered_cols_acute <- c(which(none_cols), which(two_d_cols), which(fourteen_d_cols))

missing_tls_genes <- !tls_gene_names %in% rownames(normalized_cma_ddr)
if (any(missing_tls_genes)) {
  print(paste("The following genes are missing from normalized_cma_ddr:", tls_gene_names[missing_tls_genes]))
  tls_gene_names <- tls_gene_names[!missing_tls_genes] # Remove missing genes
}

# Create the heatmap with the reordered columns
heatmap_tls_acute <- Heatmap(as.matrix(filtered_cma_tls),
                   name = "Expression",
                   top_annotation = column_acute,
                   cluster_rows = FALSE, # Disable row clustering
                   cluster_columns = FALSE, # Disable column clustering
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   column_title = "TLS Genes Heatmap",
                   row_title = "Genes",
                   column_order = ordered_cols_acute, # Reordering columns explicitly
                   row_order = tls_gene_names) # Ensure row order matches tlsgenes


# Draw the heatmap
draw(heatmap_tls_acute)
pdf("Scott/heatmap_4NQO_acute_tls.pdf", width = 10, height = 5, useDingbats = FALSE) # Adjust width and height as needed
draw(heatmap_tls_acute)
dev.off() # Close the graphics device






#HEATMAP of DDR genes

# Filter normalized_cma_ddr to include only genes in ddr_genes_acute
common_genes <- intersect(rownames(normalized_cma_ddr), ddr_genes_acute$Gene_Name)
ncma_ddr <- normalized_cma_ddr[common_genes, ]

# Order filtered_cma_ddr rows to match ddr_genes_acute
ncma_ddr <- ncma_ddr[ddr_genes_acute$Gene_Name[ddr_genes_acute$Gene_Name %in% common_genes], , drop = FALSE]

# 2. Prepare Pathway Annotation:
# Create a pathway annotation vector
pathway_annotation_ddr <- ddr_genes_acute$DDR_pathway[ddr_genes_acute$Gene_Name %in% common_genes]
names(pathway_annotation_ddr) <- common_genes

# Ensure the pathway annotation is in the same order as the rows in filtered_cma_ddr
pathway_annotation_ddr <- pathway_annotation_ddr[rownames(ncma_ddr)]

# Create a color mapping for pathways
unique_pathways <- unique(pathway_annotation_ddr)
pathway_colors <- colorRampPalette(brewer.pal(length(unique_pathways), "Set3"))(length(unique_pathways))
names(pathway_colors) <- unique_pathways

# Create row annotation
row_annotation <- rowAnnotation(Pathway = pathway_annotation_ddr,
                                col = list(Pathway = pathway_colors))


# 4. Create Heatmap:
heatmap_ddr_acute <- Heatmap(ncma_ddr,
                   name = "Z-score",
                   top_annotation = column_acute,
                   col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                   cluster_rows = TRUE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_title = "DDR Genes Heatmap (Z-score)",
                   row_title = "Genes",
                   row_names_gp = gpar(fontsize = 4),
                   column_order = ordered_cols_acute) # Reordering columns explicitly

# 5. Draw Heatmap:
draw(heatmap_ddr_acute)
pdf("Scott/heatmap_4NQO_acute_ddr_2.pdf", width = 10, height = 5, useDingbats = FALSE) # Adjust width and height as needed
draw(heatmap_ddr_acute)
dev.off() # Close the graphics device




#experiment
# 1. Create a mapping between Gene_Name and DDR_pathway
gene_pathway_mapping_acute <- ddr_genes_acute %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  dplyr::select(Gene_Name, DDR_pathway)
# 2. Get all genes in the heatmap
unique_gene_indices <- !duplicated(all_genes_heatmap_acute)
all_genes_heatmap_acute <- all_genes_heatmap_acute[unique_gene_indices]
ncma_ddr <- ncma_ddr[unique_gene_indices, ] # Also update the data matrix
all_genes_heatmap_acute <- rownames(ncma_ddr)
# 3. Create a data frame for all genes in the heatmap
all_genes_df_acute <- tibble(Gene_Name = all_genes_heatmap_acute)
# 4. Merge pathway information with all genes in the heatmap
row_annotation_data_full_acute <- all_genes_df_acute %>%
  left_join(gene_pathway_mapping_acute, by = "Gene_Name") %>%
  arrange(factor(Gene_Name, levels = all_genes_heatmap_acute))
# 5. Get unique DDR pathways
unique_pathways_acute <- unique(ddr_genes_acute$DDR_pathway)
# 6. Create a list to store row annotation columns
annotation_list_acute <- list()
# 7. Create individual binary classifiers for each pathway
for (pathway in unique_pathways_acute) {
  # Sanitize pathway names for use as column names (remove spaces and special characters)
  column_name <- gsub("[^[:alnum:]]", ".", pathway)
  annotation_list_acute[[column_name]] <- ifelse(row_annotation_data_full_acute$DDR_pathway == pathway, "In Pathway", "Not in Pathway")
  annotation_list_acute[[column_name]] <- factor(annotation_list_acute[[column_name]], levels = c("Not in Pathway", "In Pathway"))
}
# 8. Define colors for the binary annotation
pathway_colors_acute <- c("Not in Pathway" = "lightgray", "In Pathway" = "purple")
# 9. Create the row annotation
row_ha <- rowAnnotation(
  df = as.data.frame(annotation_list_acute),
  col = setNames(rep(list(pathway_colors_acute), length(annotation_list_acute)), names(annotation_list_acute)),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
)

# Now you can use this 'row_ha' object in your Heatmap function:
heatmap_ddr_acute <- Heatmap(
  ncma_ddr,
  name = "Z-score",
  left_annotation = row_ha, # Add the row annotation here
  top_annotation = column_acute, # Add the top annotation here
  col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "DDR Genes Heatmap (Z-score)",
  row_title = "Genes",
  row_names_gp = gpar(fontsize = 4),
  column_order = ordered_cols_acute # Reordering columns explicitly
)

# Draw the heatmap
draw(heatmap_ddr_acute,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     row_dend_width = unit(1, "cm")
)

pdf("Scott/heatmap_4NQO_acute_ddr_2.pdf", width = 12, height = 15, useDingbats = FALSE) # Adjust width and height as needed
draw(heatmap_ddr_acute)
dev.off() # Close the graphics device











#LEE et al.
#load the data
# Load count matrix
count_matrix_lee <- read.table("Lee_et_al/GSE229289_readcount_table.txt", header = TRUE, sep = "\t", row.names = 1)
# Remove everything after and including the dot in the row names and assign back to matrix 
rownames(count_matrix_lee) <- sub("\\..*", "", rownames(count_matrix_lee))
count_matrix_lee[1:10,1:10]
# convert ENSEMBL IDs to gene symbols and add as a new column to results
gene_symbols_lees <- mapIds(org.Mm.eg.db,
                             keys = rownames(count_matrix_lee),
                             keytype = 'ENSEMBL',
                             column = 'SYMBOL')
count_matrix_lee$gene_symbols <- gene_symbols_lees

#load meatada
samples_lee <- read_xlsx("Lee_et_al/lee_metadata.xlsx")
# load the ddr_genes list (already loaded in the project as ddr_genes)


# Convert Gene Names to uppercase and trim whitespaces for consistent matching
count_matrix_lee$gene_symbols <- trimws(toupper(count_matrix_lee$gene_symbols))

#ddr_genes_lee <- ddr_genes
#gene_mapping_lee <- data.frame(
#  column1 = c("XRCC7 (PRKDC)", "MRE11", "53BP1", "ABRAXAS1", "XPB (ERCC3)", "XPD (ERCC2)", "XPG (ERCC5)"),
#  column2 = c("PRKDC", "MRE11A", "TRP53BP1", "FAM175A", "ERCC3", "ERCC2", "ERCC5")
#)

# Replace names in ddr_genes_acute$Gene_Name using the mapping
#ddr_genes_lee$Gene_Name <- ifelse(
#  ddr_genes_lee$Gene_Name %in% gene_mapping_acute$column1,
#  gene_mapping_acute$column2[match(ddr_genes_acute$Gene_Name, gene_mapping_acute$column1)],
#  ddr_genes_acute$Gene_Name
#)


# Filter for DDR genes using Geneid column
#ENSMUSG00000038225 mouse ensembl id for RPIMPOL aka CCDC111
"CCDC111" %in% count_matrix_lee$gene_symbols
"PRIMPOL" %in% count_matrix_lee$gene_symbols
ddr_genes_lee <- ddr_genes_acute
"CCDC111" %in% ddr_genes_lee$Gene_Name
"PRIMPOL" %in% ddr_genes_lee$Gene_Name
ddr_genes_lee$Gene_Name <- gsub("CCDC111", "PRIMPOL", ddr_genes_lee$Gene_Name)
"CCDC111" %in% ddr_genes_lee$Gene_Name
"PRIMPOL" %in% ddr_genes_lee$Gene_Name
filtered_cml_ddr <- count_matrix_lee %>%
  filter(gene_symbols %in% ddr_genes_lee$Gene_Name)
head(filtered_cml_ddr)

# Identify missing genes
#missing_genes_acute <- setdiff(ddr_genes$Gene_Name, filtered_cma_ddr$Geneid)
# View the missing genes
#print(missing_genes_acute)
# Optional: Save to a CSV
#write.csv(missing_genes_acute, "Scott/missing_ddr_genes.csv", row.names = FALSE)

# Remove the row with ENSMUSG00000075569
filtered_cml_ddr <- filtered_cml_ddr %>%
  filter(rownames(.) != "ENSMUSG00000075569")
# Set gene_symbols as rownames
rownames(filtered_cml_ddr) <- filtered_cml_ddr$gene_symbols
head(filtered_cml_ddr)
# Remove the Geneid column
filtered_cml_ddr$gene_symbols <- NULL
# View the updated data frame
head(filtered_cml_ddr)


# Log transformation
normalized_cml_ddr <- log2(filtered_cml_ddr + 1)
# Calculate Z-scores for each row (gene)
normalized_cml_ddr <- t(apply(normalized_cml_ddr, 1, function(x) (x - mean(x)) / sd(x)))


#HEATMAP TLS genes
#tlsgenes <- read_xlsx("TLSgenes.xlsx")
# Filter the data for TLS genes
#tls_gene_names <- tlsgenes$genes
tls_gene_names
filtered_cml_tls <- normalized_cml_ddr[rownames(normalized_cml_ddr) %in% tls_gene_names, ]
# ensure consistent case
rownames(normalized_cml_ddr) <- trimws(toupper(rownames(normalized_cml_ddr)))
# Reorder columns to match metadata order
filtered_cml_tls <- filtered_cml_tls[, samples_lee$sample_id]

# Create the annotation for grouping
group_annotation_lee <- samples_lee$pathological_stage
names(group_annotation_lee) <- colnames(filtered_cml_tls)
# Ensure the group factor levels are set according to the desired order
# Ensure the group factor levels are set according to the desired order
group_annotation_lee <- factor(group_annotation_lee, levels = c("Normal", "Hyperplasia", "Dysplasia", "severe_Dysplasia", "Carcinoma"))
# Define the heatmap annotation with specified colors
column_lee <- HeatmapAnnotation(Group = group_annotation_lee,
                                col = list(Group = c("Normal" = "gray",
                                                     "Hyperplasia" = "lightgreen",
                                                     "Dysplasia" = "lightblue",
                                                     "severe_Dysplasia" = "yellow",
                                                     "Carcinoma" = "orange")))

# Create logical indices for each group
normal_cols <- group_annotation_lee == "Normal"
hyperplasia_cols <- group_annotation_lee == "Hyperplasia"
dysplasia_cols <- group_annotation_lee == "Dysplasia"
sev_dysplasia_cols <- group_annotation_lee == "severe_Dysplasia"
carcinoma_cols <- group_annotation_lee == "Carcinoma"


# Order the columns based on the logical indices
ordered_cols_lee <- c(which(normal_cols), which(hyperplasia_cols), which(dysplasia_cols), which(sev_dysplasia_cols), which(carcinoma_cols))

# Create the heatmap with the reordered columns
heatmap_lee_tls <- Heatmap(as.matrix(filtered_cml_tls),
                   name = "Z-score",
                   top_annotation = column_lee,
                   col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                   cluster_rows = FALSE, # Disable row clustering
                   cluster_columns = FALSE, # Disable column clustering
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   column_title = "TLS Genes Heatmap",
                   row_title = "Genes",
                   column_order = ordered_cols_lee, # Reordering columns explicitly
                   row_order = tls_gene_names) # Ensure row order matches tlsgenes

# Draw the heatmap
draw(heatmap_lee_tls)
pdf("Lee_et_al/heatmap_4NQO_lee_tls_2.pdf", width = 10, height = 5, useDingbats = FALSE) # Adjust width and height as needed
draw(heatmap_lee_tls)
dev.off() # Close the graphics device





#HEATMAP of DDR genes

# Filter normalized_cma_ddr to include only genes in ddr_genes_acute
common_genes_lee <- intersect(rownames(normalized_cml_ddr), ddr_genes_lee$Gene_Name)
ncml_ddr <- normalized_cml_ddr[common_genes_lee, ]

# Order filtered_cma_ddr rows to match ddr_genes_acute
ncml_ddr <- ncml_ddr[ddr_genes_lee$Gene_Name[ddr_genes_lee$Gene_Name %in% common_genes_lee], , drop = FALSE]
# Reorder columns to match metadata order
ncml_ddr <- ncml_ddr[, samples_lee$sample_id]


# 2. Prepare Pathway Annotation:
# Create a pathway annotation vector
pathway_annotation_ddr_lee <- ddr_genes_lee$DDR_pathway[ddr_genes_lee$Gene_Name %in% common_genes_lee]
names(pathway_annotation_ddr_lee) <- common_genes_lee

# Ensure the pathway annotation is in the same order as the rows in filtered_cma_ddr
pathway_annotation_ddr_lee <- pathway_annotation_ddr_lee[rownames(ncml_ddr)]

# Create a color mapping for pathways
unique_pathways_lee <- unique(pathway_annotation_ddr_lee)
pathway_colors_lee <- colorRampPalette(brewer.pal(length(unique_pathways_lee), "Set3"))(length(unique_pathways_lee))
names(pathway_colors_lee) <- unique_pathways_lee

# Create row annotation
row_annotation_lee <- rowAnnotation(Pathway = pathway_annotation_ddr_lee,
                                col = list(Pathway = pathway_colors_lee))


# 4. Create Heatmap:
heatmap_ddr_lee <- Heatmap(ncml_ddr,
                   name = "Z-score",
                   right_annotation = row_annotation_lee,
                   top_annotation = column_lee,
                   cluster_rows = TRUE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_title = "DDR Genes Heatmap (Z-score)",
                   row_title = "Genes",
                   row_names_gp = gpar(fontsize = 4),
                   column_order = ordered_cols_lee) # Reordering columns explicitly

# 5. Draw Heatmap:
draw(heatmap_ddr_lee)


#experiment
# Assuming you have ddr_genes_lee and ncml_ddr defined as in your provided data
# Assuming you have samples_lee defined elsewhere and it contains 'pathological_stage' and 'sample_id'
# 1. Create a mapping between Gene_Name and DDR_pathway
gene_pathway_mapping_lee <- ddr_genes_lee %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  dplyr::select(Gene_Name, DDR_pathway)
# 2. Get all genes in the heatmap
all_genes_heatmap_lee <- rownames(ncml_ddr)
unique_gene_indices <- !duplicated(all_genes_heatmap_lee)
all_genes_heatmap_lee <- all_genes_heatmap_lee[unique_gene_indices]
ncml_ddr <- ncml_ddr[unique_gene_indices, ] # Also update the data matrix
# 3. Create a data frame for all genes in the heatmap
all_genes_df_lee <- tibble(Gene_Name = all_genes_heatmap_lee)
# 4. Merge pathway information with all genes in the heatmap
row_annotation_data_full_lee <- all_genes_df_lee %>%
  left_join(gene_pathway_mapping_lee, by = "Gene_Name") %>%
  arrange(factor(Gene_Name, levels = all_genes_heatmap_lee))
# 5. Get unique DDR pathways
unique_pathways_lee <- unique(ddr_genes_lee$DDR_pathway)
# 6. Create a list to store row annotation columns
annotation_list_lee <- list()
# 7. Create individual binary classifiers for each pathway
for (pathway in unique_pathways_lee) {
  # Sanitize pathway names for use as column names (remove spaces and special characters)
  column_name <- gsub("[^[:alnum:]]", ".", pathway)
  annotation_list_lee[[column_name]] <- ifelse(row_annotation_data_full_lee$DDR_pathway == pathway, "In Pathway", "Not in Pathway")
  annotation_list_lee[[column_name]] <- factor(annotation_list_lee[[column_name]], levels = c("Not in Pathway", "In Pathway"))
}
# 8. Define colors for the binary annotation
pathway_colors_lee <- c("Not in Pathway" = "lightgray", "In Pathway" = "purple")
# 9. Create the row annotation
row_annotation_lee_final <- rowAnnotation(
  df = as.data.frame(annotation_list_lee),
  col = setNames(rep(list(pathway_colors_lee), length(annotation_list_lee)), names(annotation_list_lee)),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
)
# Define the heatmap annotation with specified colors
group_annotation_lee <- samples_lee$pathological_stage
names(group_annotation_lee) <- colnames(ncml_ddr) # Use ncml_ddr column names
# Ensure the group factor levels are set according to the desired order
group_annotation_lee <- factor(group_annotation_lee, levels = c("Normal", "Hyperplasia", "Dysplasia", "severe_Dysplasia", "Carcinoma"))
column_lee <- HeatmapAnnotation(Group = group_annotation_lee,
                                col = list(Group = c("Normal" = "gray",
                                                     "Hyperplasia" = "lightgreen",
                                                     "Dysplasia" = "lightblue",
                                                     "severe_Dysplasia" = "yellow",
                                                     "Carcinoma" = "orange")))
# Create logical indices for each group
normal_cols <- group_annotation_lee == "Normal"
hyperplasia_cols <- group_annotation_lee == "Hyperplasia"
dysplasia_cols <- group_annotation_lee == "Dysplasia"
sev_dysplasia_cols <- group_annotation_lee == "severe_Dysplasia"
carcinoma_cols <- group_annotation_lee == "Carcinoma"
# Order the columns based on the logical indices
ordered_cols_lee <- c(which(normal_cols), which(hyperplasia_cols), which(dysplasia_cols), which(sev_dysplasia_cols), which(carcinoma_cols))
# Create the heatmap
heatmap_ddr_lee <- Heatmap(ncml_ddr,
                           name = "Z-score",
                           left_annotation = row_annotation_lee_final, # Use the final row annotation
                           top_annotation = column_lee,
                           col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                           cluster_rows = TRUE,
                           cluster_columns = FALSE,
                           show_row_names = TRUE,
                           show_column_names = TRUE,
                           column_title = "DDR Genes Heatmap (Z-score)",
                           row_title = "Genes",
                           row_names_gp = gpar(fontsize = 4),
                           column_order = ordered_cols_lee) # Reordering columns explicitly

# Draw the heatmap
draw(heatmap_ddr_lee,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     row_dend_width = unit(1, "cm")
)
pdf("Lee_et_al/heatmap_4NQO_lee_ddr_2.pdf", width = 12, height = 15, useDingbats = FALSE) # Adjust width and height as needed
draw(heatmap_ddr_lee)
dev.off() # Close the graphics device











#pathway hetampa


#experiment 3
# Assuming you have gene_pathway_mapping_lee defined
pathway_data_lee <- ncml_ddr %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_Name") %>%
  left_join(gene_pathway_mapping_lee, by = "Gene_Name") %>%
  group_by(DDR_pathway) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  filter(!is.na(DDR_pathway)) %>% # Ensure no NA pathways are used as row names.
  column_to_rownames(var = "DDR_pathway")
# pathway_data_lee now has pathways as rows and samples as columns
# Define colors for the heatmap values (adjust as needed)
#heatmap_colors_pathway <- colorRamp2(c(min(pathway_data_lee), 0, max(pathway_data_lee)), c("blue", "white", "red"))
heatmap_ddr_pathway_lee <- Heatmap(as.matrix(pathway_data_lee), # Convert to matrix
                                   name = "Average Z-score",
                                   top_annotation = column_lee,
                                   col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   column_title = "DDR Pathways Heatmap",
                                   row_title = "DDR Pathway",
                                   row_names_gp = gpar(fontsize = 8),
                                   column_order = ordered_cols_lee)
# Draw the heatmap
draw(heatmap_ddr_pathway_lee,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     row_dend_width = unit(1, "cm"))

pdf("Lee_et_al/heatmap_4NQO__lee_ddr_pathway_2new.pdf", width = 11, height = 8)
draw(heatmap_ddr_pathway_lee)
dev.off()





#experiment 4
# Assuming you have gene_pathway_mapping_lee defined
pathway_data_acute <- ncma_ddr %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_Name") %>%
  left_join(gene_pathway_mapping_lee, by = "Gene_Name") %>%
  group_by(DDR_pathway) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames(var = "DDR_pathway")
pathway_data_acute <- ncma_ddr %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_Name") %>%
  left_join(gene_pathway_mapping_lee, by = "Gene_Name") %>%
  group_by(DDR_pathway) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  filter(!is.na(DDR_pathway)) %>% # Ensure no NA pathways are used as row names.
  column_to_rownames(var = "DDR_pathway")
# pathway_data_lee now has pathways as rows and samples as columns
# Define colors for the heatmap values (adjust as needed)
#heatmap_colors_pathway <- colorRamp2(c(min(pathway_data_lee), 0, max(pathway_data_lee)), c("blue", "white", "red"))

heatmap_ddr_pathway_acute <- Heatmap(pathway_data_acute,
                                   name = "Average Z-score", # Adjust name based on your summary metric
                                   top_annotation = column_acute,
                                   col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   column_title = "DDR Pathways Heatmap",
                                   row_title = "DDR Pathway",
                                   row_names_gp = gpar(fontsize = 8),
                                   column_order = ordered_cols_acute)
heatmap_ddr_pathway_acute <- Heatmap(as.matrix(pathway_data_acute), # Convert to matrix
                                   name = "Average Z-score",
                                   top_annotation = column_acute,
                                   col = colorRamp2(c((-2), 0, 2), c("blue", "white", "red")),
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   column_title = "DDR Pathways Heatmap",
                                   row_title = "DDR Pathway",
                                   row_names_gp = gpar(fontsize = 8),
                                   column_order = ordered_cols_acute)
# Draw the heatmap
draw(heatmap_ddr_pathway_acute,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     row_dend_width = unit(1, "cm"))

pdf("Scott/heatmap_4NQO__acute_ddr_pathway_2new.pdf", width = 11, height = 8, useDingbats = FALSE)
draw(heatmap_ddr_pathway_acute)
dev.off()
