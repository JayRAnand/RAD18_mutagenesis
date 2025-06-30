#library
library(MutationalPatterns)
library(BSgenome)
library(here)
library(readxl)
library(ggplot2)
library(openxlsx)
library(dplyr)


#load reference genome
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)




# SNVs


#load data i.e. VCF files
vcf_path_snv <- here("MAGEA4", "SNVs")
vcf_files_snv <- list.files(vcf_path_snv, pattern = '*.vcf$')
##If vcf_files is a vector of file names, file.path() will combine each name with vcf_path, producing a vector of full file paths.
vcf_dir_snv <- file.path(vcf_path_snv, vcf_files_snv)
samples <- read_excel(here("MAGEA4", "SNVs", "samples.xlsx"))
sample_names <- c(samples$name)
group <- c(samples$group)

#load the VCF files into GRangesList
snv_grl <- read_vcfs_as_granges(vcf_dir_snv, sample_names, ref_genome, type = "snv")

# base substituition types
# count mutation type occurences
type_occurrences <- mut_type_occurrences(snv_grl, ref_genome)
type_occurrences
#save file for future reference
# Create a workbook object
snv_count <- createWorkbook()
# Add a worksheet
addWorksheet(snv_count, sheetName = "Sheet1")
# Write the data to the worksheet, including row names
writeData(snv_count, sheet = "Sheet1", x = type_occurrences, rowNames = TRUE)
# Save the workbook to a file
saveWorkbook(snv_count, "output_data/mage_snv_count.xlsx")

# plot_spectrum
# Specify the path to the folder
output_folder <- "plots"
# Ensure the directory exists (create it if it doesn't)
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}
# Create the full path to the output PDF file
pdf_file_path <- file.path(output_folder, "MAGEA4_mut_spectrum.pdf")
# Open a PDF file to save the plot
pdf(pdf_file_path)
# Create the plot
MAGEA4_mut_spectrum <- plot_spectrum(type_occurrences, by = group, CT = TRUE, indv_points = TRUE, legend = TRUE)
# Print the plot (it will be saved to the PDF)
print(MAGEA4_mut_spectrum)
# Close the PDF device to save the plot
dev.off()


# 96 mutational profile
mut_mat <- mut_matrix(vcf_list = snv_grl, ref_genome = ref_genome)
head(mut_mat)
# Rearrange columns in ascending order based on column names
mut_mat_ordered <- mut_mat[, order(colnames(mut_mat))]
head(mut_mat_ordered)
#plot the 96 profile of samples
pdf_file_path_96 <- file.path(output_folder, "MAGEA4_plot_96_profile.pdf")
pdf(pdf_file_path_96, height = 14)
MAGEA4_plot_96_profile <- plot_96_profile(mut_mat_ordered)
print(MAGEA4_plot_96_profile)
dev.off()


#COSMIC signature
signatures = get_known_signatures()
# Inspect the structure of the 'signatures' object
str(signatures)
# Look for attributes or metadata that might indicate the COSMIC version
attributes(signatures)


fit_res <- fit_to_signatures(mut_mat_ordered, signatures)
#Plot
pdf_file_path_cosmic <- file.path(here("plots"), "MAGEA4_cosmic_abs.pdf")
pdf(pdf_file_path_cosmic)
MAGEA4_cosmic_abs <- plot_contribution(fit_res$contribution,
                                  coord_flip = FALSE,
                                  mode = "absolute")
# Rotate x-axis labels (sample names)
MAGEA4_cosmic_abs + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate labels at 45 degrees
)
print(MAGEA4_cosmic_abs)
dev.off()


# cosine similarity
plot_original_vs_reconstructed(mut_mat_ordered, fit_res$reconstructed, 
                              y_intercept = 0.95)

#strict
strict_refit <- fit_to_signatures_strict(mut_mat_ordered,
                                         signatures,
                                         max_delta = 0.004,
                                         method = c("backwards", "best_subset"))
fit_res_strict <- strict_refit$fit_res
#Plot
pdf_file_path_cosmic_strict <- file.path(here("plots"), "MAGEA4_cosmic_strict_abs.pdf")
pdf(pdf_file_path_cosmic_strict)
MAGEA4_cosmic_strict_abs <- plot_contribution(fit_res_strict$contribution,
                                         coord_flip = FALSE,
                                         mode = "absolute")
MAGEA4_cosmic_strict_abs + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate labels at 45 degrees
)
# Print the plot with the modified y-axis
print(MAGEA4_cosmic_strict_abs)
# Close the PDF device to save the plot
dev.off()

# Plot relative
pdf_file_path_cosmic_strict <- file.path(here("plots"), "MAGEA4_cosmic_strict_rel.pdf")
pdf(pdf_file_path_cosmic_strict)
MAGEA4_cosmic_strict_rel <- plot_contribution(fit_res_strict$contribution,
                                              coord_flip = FALSE,
                                              mode = "relative")
MAGEA4_cosmic_strict_rel + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate labels at 45 degrees
)
# Print the plot with the modified y-axis
print(MAGEA4_cosmic_strict_rel)
# Close the PDF device to save the plot
dev.off()


# pooling samples
pooled_mut_mat <- pool_mut_mat(mut_mat, grouping = group)
head(pooled_mut_mat)

#plot the 96 profile of samples
pdf_file_path_96_pooled <- file.path(here("plots"), "MAGEA4_pooled_96_profile.pdf")
pdf(pdf_file_path_96_pooled)
MAGEA4_pooled_96_profile <- plot_96_profile(pooled_mut_mat)
print(MAGEA4_pooled_96_profile)
dev.off()

#COSMIC signature
fit_res_pooled <- fit_to_signatures(pooled_mut_mat, signatures)
#Plot
pdf_file_path_pooled_cosmic <- file.path(here("plots"), "MAGEA4_pooled_cosmic_abs.pdf")
pdf(pdf_file_path_pooled_cosmic)
MAGEA4_pooled_cosmic_abs <- plot_contribution(fit_res_pooled$contribution,
                                       coord_flip = FALSE,
                                       mode = "absolute")
print(MAGEA4_pooled_cosmic_abs)
dev.off()

# cosine similarity
plot_original_vs_reconstructed(pooled_mut_mat, fit_res_pooled$reconstructed, 
                               y_intercept = 0.95)


#strict
strict_refit_pool <- fit_to_signatures_strict(pooled_mut_mat,
                                         signatures,
                                         max_delta = 0.004,
                                         method = c("backwards", "best_subset"))
fit_res_strict_pool <- strict_refit_pool$fit_res
#Plot
pdf_file_path_pooled_cosmic_strict <- file.path(here("plots"), "MAGEA4_pooled_cosmic_strict_abs.pdf")
pdf(pdf_file_path_pooled_cosmic_strict)
MAGEA4_pooled_cosmic_strict_abs <- plot_contribution(fit_res_strict_pool$contribution,
                                              coord_flip = FALSE,
                                              mode = "absolute")
print(MAGEA4_pooled_cosmic_strict_abs)
# Close the PDF device to save the plot
dev.off()

# Plot relative
pdf_file_path_cosmic_strict <- file.path(here("plots"), "MAGEA4_pooled_cosmic_strict_rel.pdf")
pdf(pdf_file_path_cosmic_strict)
MAGEA4_pooled_cosmic_strict_rel <- plot_contribution(fit_res_strict_pool$contribution,
                                              coord_flip = FALSE,
                                              mode = "relative")
print(MAGEA4_pooled_cosmic_strict_rel)
dev.off()


#by group
# Assuming 'mut_mat' is a matrix with rows = SNV types and columns = sample names
# And 'samples' is a dataframe with columns 'name' (sample names) and 'group' (sample groups)

# Step 1: Match sample names in 'mut_mat' to sample names in 'samples'
sample_names <- colnames(mut_mat)  # Sample names in the 'mut_mat' matrix
group_info <- samples$group[match(sample_names, samples$name)]  # Group info corresponding to each sample

# Step 2: Create a new matrix for average counts by group
# We will create a matrix where rows are SNVS\s types and columns are groups

# Step 3: Initialize the new matrix with zeros, with indel types as rows and group names as columns
group_names <- unique(group_info)  # Unique group names
mut_avg_by_group <- matrix(0, nrow = nrow(mut_mat), ncol = length(group_names),
                           dimnames = list(rownames(mut_mat), group_names))

# Step 4: Loop through each group and compute the average count of SNV per type
for (group_name in group_names) {
  # Find which samples belong to the current group
  group_samples <- sample_names[group_info == group_name]
  
  # Subset the 'mut_mat' matrix to include only the current group's samples
  group_mut <- mut_mat[, group_samples, drop = FALSE]
  
  # Calculate the average count for each dbs type in the current group
  mut_avg_by_group[, group_name] <- rowMeans(group_mut, na.rm = TRUE)
}

# The 'mut_avg_by_group' matrix now contains the average count per indel type per group
head(mut_avg_by_group)  # Display the first few rows


#plot
pdf_file_path_96_grouped <- file.path(here("plots"), "MAGEA4_grouped_96_profile.pdf")
pdf(pdf_file_path_96_grouped, height = 4)
MAGEA4_grouped_96_profile <- plot_96_profile(mut_avg_by_group)
print(MAGEA4_grouped_96_profile)
dev.off()





#DSB

dbs_grl <- read_vcfs_as_granges(vcf_dir_snv, sample_names, ref_genome, type = "dbs")
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)
#Plot
pdf_file_path_dsb <- file.path(here("plots"), "MAGEA4_dsb.pdf")
pdf(pdf_file_path_dsb)
MAGEA4_dsb <- plot_dbs_contexts(dbs_counts,  same_y = TRUE, condensed = TRUE)
print(MAGEA4_dsb)
dev.off()



# INDELS

#load data i.e. VCF files
vcf_path_indel <- here("MAGEA4", "INDELs")
vcf_files_indel <- list.files(vcf_path_indel, pattern = '*.vcf$')
##If vcf_files is a vector of file names, file.path() will combine each name with vcf_path, producing a vector of full file paths.
vcf_dir_indel <- file.path(vcf_path_indel, vcf_files_indel)
samples_indel <- read_excel(here("MAGEA4", "INDELS", "samples.xlsx"))
sample_names_indel <- c(samples_indel$name)
group_indel <- c(samples_indel$group)

indel_grl <- read_vcfs_as_granges(vcf_dir_indel, sample_names_indel, ref_genome, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
head(indel_grl)

indel_counts <- count_indel_contexts(indel_grl)
head(indel_counts)


#Plot
pdf_file_path_indel <- file.path(here("plots"), "MAGEA4_indel.pdf")
pdf(pdf_file_path_indel, width = 12, height = 20)
MAGEA4_indel <- plot_indel_contexts(indel_counts,  same_y = TRUE, condensed = TRUE, extra_labels = FALSE)
print(MAGEA4_indel)
dev.off()

# pooling samples
pooled_indel_counts <- pool_mut_mat(indel_counts, grouping = group_indel)
head(pooled_indel_counts)

#Plot
pdf_file_path_pooled_indel <- file.path(here("plots"), "MAGEA4_pooled_indel.pdf")
pdf(pdf_file_path_pooled_indel, width = 12, height = 4.5)
MAGEA4_pooled_indel <- plot_indel_contexts(pooled_indel_counts,  same_y = TRUE, condensed = TRUE, extra_labels = FALSE)
print(MAGEA4_pooled_indel)
dev.off()

#by group
sample_names <- colnames(indel_counts)  # Sample names in the 'indel_counts' matrix
group_info <- samples$group[match(sample_names, samples$name)]  # Group info corresponding to each sample
group_names <- unique(group_info)  # Unique group names
indel_avg_by_group <- matrix(0, nrow = nrow(indel_counts), ncol = length(group_names),
                             dimnames = list(rownames(indel_counts), group_names))
for (group_name in group_names) {
  group_samples <- sample_names[group_info == group_name]
  group_indels <- indel_counts[, group_samples, drop = FALSE]
  indel_avg_by_group[, group_name] <- rowMeans(group_indels, na.rm = TRUE)
}
head(indel_avg_by_group)  # Display the first few rows


#plot
pdf_file_path_grouped_indel <- file.path(here("plots"), "MAGEA4_grouped_indel.pdf")
pdf(pdf_file_path_grouped_indel, width = 12, height = 4.5)
MAGEA4_grouped_indel <- plot_indel_contexts(indel_avg_by_group,  same_y = TRUE, condensed = TRUE)
print(MAGEA4_grouped_indel)
dev.off()






signatures_indel = get_known_signatures(muttype = "indel")
signatures_indel

signatures_head = get_known_signatures(source = "SIGNAL", sig_type = "tissue", tissue_type = "Head")
signatures_head [1:5, 1:5]
fit_res_head <- fit_to_signatures(mut_mat_ordered, signatures_head)
fit_res_head$contribution[1:5, 1:5]

#plot absolute
pdf_file_path_head <- file.path(here("plots"), "MAGEA4_cosmic_abs_head.pdf")
pdf(pdf_file_path_head)
MAGEA4_cosmic_abs_head <- plot_contribution(fit_res_head$contribution,
                                       coord_flip = FALSE,
                                       mode = "absolute")
MAGEA4_cosmic_abs_head + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(MAGEA4_cosmic_abs_head)
dev.off()


#plot relative
pdf_file_path_head <- file.path(here("plots"), "MAGEA4_cosmic_rel_head.pdf")
pdf(pdf_file_path_head)
MAGEA4_cosmic_rel_head <- plot_contribution(fit_res_head$contribution,
                                            coord_flip = FALSE,
                                            mode = "relative")
MAGEA4_cosmic_rel_head + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(MAGEA4_cosmic_rel_head)
dev.off()


#convert to reference signature
fit_res_headref <- convert_sigs_to_ref(fit_res_head)
fit_res_headref$contribution[1:5, 1:5]

#plot absolute
pdf_file_path_headref <- file.path(here("plots"), "MAGEA4_cosmic_abs_headref.pdf")
pdf(pdf_file_path_headref)
MAGEA4_cosmic_abs_headref <- plot_contribution(fit_res_headref$contribution,
                                            coord_flip = FALSE,
                                            mode = "absolute")
MAGEA4_cosmic_abs_headref + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(MAGEA4_cosmic_abs_headref)
dev.off()


#plot relative
pdf_file_path_headref <- file.path(here("plots"), "MAGEA4_cosmic_rel_headref.pdf")
pdf(pdf_file_path_headref)
MAGEA4_cosmic_rel_headref <- plot_contribution(fit_res_headref$contribution,
                                            coord_flip = FALSE,
                                            mode = "relative")
MAGEA4_cosmic_rel_headref + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(MAGEA4_cosmic_rel_headref)
dev.off()








#Convert Absolute to Relative Contributions
#first normalize the contributions per sample
# Transpose, normalize columns, transpose back
rel_contrib <- t(apply(fit_res_strict$contribution, 2, function(x) x / sum(x)))
rel_contrib[, 1:5]
# Define sample groups (must match sample order in your matrix)
group_info <- data.frame(
  sample = colnames(fit_res_strict$contribution),
  group = c("MAGE-A4", "MAGE-A4", "MAGE-A4", "MAGE-A4", "MAGE-A4", "MAGE-A4", "MAGE-A4", "MAGE-A4", 
            "WT", "WT", "WT", "WT", "WT")
)
# Match group info to normalized matrix
rel_contrib_df <- as.data.frame(rel_contrib)
rel_contrib_df$sample <- rownames(rel_contrib_df)
rel_contrib_df <- left_join(rel_contrib_df, group_info, by = "sample")


#Perform t-test + Wilcoxon Per Signature
#Use two statistical tests:
#A Student's t-test: checks if the means of the two groups are significantly different.
#A Wilcoxon rank-sum test (also called Mann–Whitney U): a non-parametric test that compares medians/ranks, useful if data isn't normally distributed.
#Take the smallest p-value between the two tests as the most conservative/significant result for that signature.
#Repeat for each signature (i.e., one test per row in your contribution matrix).
#Apply FDR correction to adjust for multiple testing (since you test many signatures).

# Function to test each signature
test_signature <- function(sig, df) {
  x <- df %>% filter(group == "WT") %>% pull(sig)
  y <- df %>% filter(group == "MAGE-A4") %>% pull(sig)
  
  t_p <- tryCatch(t.test(x, y)$p.value, error = function(e) NA)
  w_p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA)
  min_p <- min(t_p, w_p, na.rm = TRUE)
  
  return(c(t_p = t_p, w_p = w_p, min_p = min_p))
}


# Run for each signature
#pval_results <- sapply(colnames(rel_contrib)[1:ncol(rel_contrib)], function(sig) {
#  test_signature(sig, rel_contrib_df)
#})

#We just need to skip signatures that have no variation (e.g., all 0s or same value across both groups).
test_signature <- function(sig, df) {
  x <- df %>% filter(group == "MAGE-A4") %>% pull(sig)
  y <- df %>% filter(group == "WT") %>% pull(sig)
  
  # Skip if all values are the same (e.g., all zero)
  if (length(unique(c(x, y))) <= 1) {
    return(c(t_p = NA, w_p = NA, min_p = NA))
  }
  
# Try tests, use NA if they fail
  t_p <- tryCatch(t.test(x, y)$p.value, error = function(e) NA)
  w_p <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA)  # use approx p-value
  min_p <- min(t_p, w_p, na.rm = TRUE)
  
  return(c(t_p = t_p, w_p = w_p, min_p = min_p))
}
pval_results <- sapply(colnames(rel_contrib)[1:ncol(rel_contrib)], function(sig) {
  test_signature(sig, rel_contrib_df)
})
# Convert and process
pval_df <- as.data.frame(t(pval_results))
pval_df$adj_p <- p.adjust(pval_df$min_p, method = "fdr")
pval_df$signature <- rownames(pval_df)

# Add significance labels
pval_df$signif <- cut(pval_df$adj_p,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", "ns"))

# View sorted results
pval_df <- pval_df %>% arrange(adj_p)
print(pval_df)

write.csv(pval_df, file = "pval_results.csv", row.names = TRUE)
