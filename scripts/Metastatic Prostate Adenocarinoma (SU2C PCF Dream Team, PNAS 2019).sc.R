################################################################################
# Install all packages and libraries
install.packages("BiocManager")
install.packages("dplyr")

# Install Bioconductor packages
BiocManager::install("preprocessCore")
BiocManager::install("limma")

library(preprocessCore)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(tidyverse)

################################################################################
#Process mRNA expression data 

# Assign the FPKM count mRNA expression data
fpkm <- read.delim("data_mrna_seq_fpkm_capture.tsv", check.names = FALSE)

#check if there a duplicate  gene names  in column 1, then change to a unique identifier 

any(duplicated(fpkm[, 1]))
duplicates <- fpkm[duplicated(fpkm[, 1]), 1]
unique_duplicates <- unique(duplicates)
print(unique_duplicates)

# Make any of the duplicate names unique using in to make compatible downstream with limma. 
rownames(fpkm) <- make.unique(as.character(fpkm[, 1]))
fpkm<- fpkm[, -1]  

#convert FPKM to TPM
fpkm_to_tpm <- function(fpkm_mat) {
  apply(fpkm_mat, 2, function(x) {
    (x / sum(x, na.rm = TRUE)) * 1e6
  })
}
tpm_matrix <- fpkm_to_tpm(fpkm)

################################################################################
# Normalise data using quantile based normalisation  

log2_tpm <- log2(tpm_matrix + 0.01)
qn_matrix <- normalize.quantiles(log2_tpm)
colnames(qn_matrix) <- colnames(log2_tpm)
rownames(qn_matrix) <- rownames(log2_tpm)



################################################################################
#Z-score and file download 
#Z-score normalization (Per Gene)
z_score<- t(scale(t(qn_matrix)))


#Write all the stages of normalisation QC 
write.csv(tpm_matrix, "tpm_normalized.csv")
write.csv(log2_tpm, "log2_tpm.csv")
write.csv(qn_matrix, "quantile_normalized.csv")
write.csv(z_score, "zscore_normalized.csv")

################################################################################
#Load the altered and unaltered data from cBioPortal website (pre-filtered) and create groups

pre_alterations <- read.delim("alterations_across_samples.tsv", stringsAsFactors = FALSE)

# ERG alteration column, then convert the vector into a factor 
alterations <- pre_alterations %>%
  mutate(TMPRSS2ERG_Fusion = if_else(grepl("TMPRSS2-ERG", `ERG..FUSION`, ignore.case = TRUE), 1, 0))

group <- ifelse(alterations$TMPRSS2ERG_Fusion == 1, "Fusion_Positive", "Fusion_Negative")
group <- factor(group)

# load the normalised expression file into R 
expr_data <- read.csv("quantile_normalized.csv", check.names = FALSE)

#name a new column name and make that the rownames
colnames(expr_data)[1] <- "GeneSymbol"
rownames(expr_data) <- expr_data[[1]]


# Extract sample IDs from the alterations data and filter expression data  
sample_ids <- alterations$Sample.ID 
expr_filtered <- expr_data[, colnames(expr_data) %in% sample_ids] 

# Match order between expression and group 

matched_indices <- match(colnames(expr_filtered), alterations$Sample.ID) 
group_matched <- group[matched_indices] 

# Check that the samples match the alterations file to the mRNA expression data
all(colnames(expr_filtered) == names(group_matched))
length(group_matched) == ncol(expr_filtered)
################################################################################
#Check if the data is normalised before using limma and other QC checks 

expr_t <- as.data.frame(t(expr_filtered))


summary_stats <- data.frame(
  Sample = rownames(expr_t),
  Mean = rowMeans(expr_t, na.rm = TRUE),
  Median = apply(expr_t, 1, median, na.rm = TRUE),
  SD = apply(expr_t, 1, sd, na.rm = TRUE),
  Max = apply(expr_t, 1, max, na.rm = TRUE),
  Min = apply(expr_t, 1, min, na.rm = TRUE)
)

print(head(summary_stats))

summary(as.numeric(as.matrix(expr_filtered)))
# Densitiy plot of normalised expression 
plotDensities(z_score, 
              group=group_matched, 
              main="Density Plot of Z-Score Normalised Expression 
   by TMPRSS2-ERG Fusion Status",
              col=c("#2C2CFD", "firebrick"))

################################################################################
# design matrix for limma and fit expression data to a linear model  
design <- model.matrix(~0+ group_matched)
colnames(design) <- levels(group_matched)
colnames(design) 
table(group_matched)

fit <- lmFit(expr_filtered, design)
fit <- eBayes(fit)

#Contrast the results between TMPRSS2-ERG fusion Status 
contrast.matrix <- makeContrasts(ERG_Fusion_vs_NEG = Fusion_Positive - Fusion_Negative, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = 'ERG_Fusion_vs_NEG', number = Inf) 
results

################################################################################
#Format the data and filter data 
gene_info <- do.call(rbind, strsplit(rownames(results), "\\|"))
colnames(gene_info) <- c("GeneSymbol")
results <- cbind(gene_info, results)

#Filter results to have the adjusted P.Value cut off of 0.05 
filtered_results <- results[results$adj.P.Val <0.05, ]
#sort results from largest to smallest LogFC
sorted_filtered_results <- filtered_results[order(-filtered_results$logFC), ]
print(sorted_filtered_results)

################################################################################

#write a gene expression table and a pre-ranked gene list
write.table( sorted_filtered_results, file = "Metastatic Prostate Adenocarinoma (SU2C PCF Dream Team, PNAS 2019) TMPRSS2-ERG DEGs 
             ranked_gene_list.tsv",
             sep = "\t",
             row.names = FALSE,
             quote = FALSE)


# write a pre-ranked gene list using .rnk and in a txt file type 
write.table(sorted_filtered_results[, 1:2],
            file = "Metastatic Prostate Adenocarinoma (SU2C PCF Dream Team, PNAS 2019) TMPRSS2-ERG DEGs 
             Pre-ranked_gene_list.rnk",
            sep = "\t",         
            row.names = FALSE,  
            col.names = FALSE,
            quote = FALSE)      
################################################################################
#Volcano Plots 

# Assign significance categories
df <- results %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -0.5~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# ERG + EMT genes from GSEA - hallmark 
highlight_genes <- c(
  "ERG", "ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COMP","COPA","CRLF1","CCN2","CTHRC1","CXCL1","CXCL12","CXCL6","CCN1","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","COLGALT1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","CXCL8","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","P3H1","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A"
)

# Filter genes to label 
label_genes <- df %>%
  filter(
    GeneSymbol %in% highlight_genes,
    abs(logFC) > 0.5,
    adj.P.Val < 0.05
  )

ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = label_genes,
    aes(label = GeneSymbol),
    color = "black",
    size = 2.5,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.1,
    point.padding = 0.0,
    segment.color = "black",
    segment.size = 0.3,
    nudge_x = -0.7,
    direction = "y",
    force = 0.09
  ) +
  scale_color_manual(values = c("Upregulated" = "lightcoral", 
                                "Downregulated" = "skyblue", 
                                "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Differential gene expression of EMT related genes 
       in TMPRSS2-ERG fusion Metastatic Prostate Adenocarcinoma ",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
