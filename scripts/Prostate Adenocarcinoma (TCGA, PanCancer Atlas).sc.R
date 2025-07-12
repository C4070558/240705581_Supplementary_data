################################################################################
# Install all packages and libraries
install.packages("BiocManager")
install.packages("dplyr")

# Install Bioconductor packages
BiocManager::install("preprocessCore")
BiocManager::install("limma")

library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Load files ################################################################################
# load the expression file and altered and unaltered data 
expr_data <- read.delim("TCGA.PRAD.sampleMap_HiSeqV2.tsv", row.names = 1, check.names = FALSE)
pre_alterations <- read.delim("alterations_across_samples.tsv", stringsAsFactors = FALSE)

# ERG alteration column then convert the vector into a factor 
alterations <- pre_alterations %>%
  mutate(TMPRSS2ERG_Fusion = if_else(grepl("TMPRSS2-ERG fusion", `ERG..FUSION`, ignore.case = TRUE), 1, 0))

group <- ifelse(alterations$TMPRSS2ERG_Fusion == 1, "Fusion_Positive", "Fusion_Negative")
group <- factor(group)

# Extract sample IDs from the alterations data and format
sample_ids <- alterations$Sample.ID 
short_sample_ids <- substr(colnames(expr_data), 1, 15) 
colnames(expr_data) <- short_sample_ids 
# Filter expression data ################################################################################
# Filter expression data with alterations and order samples to match expression in group
expr_filtered <- expr_data[, colnames(expr_data) %in% sample_ids] 
matched_indices <- match(colnames(expr_filtered), alterations$Sample.ID) 
group_matched <- group[matched_indices] 

# design matrix for limma 
design <- model.matrix(~0+ group_matched)
colnames(design) <- levels(group_matched)
colnames(design) 
table(group_matched)
length(group_matched) == ncol(expr_filtered)

# Check data quality ################################################################################
# Check data quality and normalisation / distribution 

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

# Calculate Z-Score and plot Density plot 
z_score<- t(scale(t(expr_data), center = TRUE, scale = TRUE))
plotDensities(z_score, 
              group=group_matched, 
              main="Density Plot of Z-Score Normalised Expression 
   by TMPRSS2-ERG Fusion Status in PCa samples",
              col=c("#2C2CFD", "firebrick"))

# Fit expression data to a linear model ################################################################################
# Fit expression data to a linear model 
fit <- lmFit(expr_filtered, design)
fit <- eBayes(fit)

# Contrast the results between TMPRSS2-ERG fusion Status 
contrast.matrix <- makeContrasts(ERG_Fusion_vs_NEG = Fusion_Positive - Fusion_Negative, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef = 'ERG_Fusion_vs_NEG', number = Inf)  # coef=2 

results
results <- rownames_to_column(results, var = "GeneSymbol")  
rownames(results) <- NULL
#Filter result data  ################################################################################
#Filter result data to be significany with the adj.P.Value <0.05 and sort with Log Fold Chage
filtered_results <- results[results$adj.P.Val < 0.05, ]
sorted_filtered_results <- filtered_results[order(-filtered_results$logFC), ]


# Make results file ################################################################################
# Make results file 
write.table(sorted_filtered_results,file = "Prostate Adenocarcinoma (TCGA, PanCancer Atlas) TMPRSS2-ERG 
            ranked_gene_list.tsv",
            sep = "\t", 
            row.names = FALSE,
            quote = FALSE)

write.table(sorted_filtered_results[, 1:2],
            file = "Prostate Adenocarcinoma (TCGA, PanCancer Atlas) TMPRSS2-ERG Pre-ranked_gene_list.rnk",
            sep = "\t",         
            row.names = FALSE,  
            col.names = FALSE,
            quote = FALSE)      
#Volcano Plots ################################################################################# 
#Volcano Plots 

# Assign significance categories
df <- results %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -0.5~ "Downregulated",
    TRUE ~ "Not Significant"
  ))


# Filter genes to label 
label_genes <- df %>%
  filter(
    GeneSymbol %in% highlight_genes,
    abs(logFC) > 0.5,
    adj.P.Val < -log10(0.04)
  )


highlight_genes <- c(
    "ERG", "ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COMP","COPA","CRLF1","CCN2","CTHRC1","CXCL1","CXCL12","CXCL6","CCN1","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","COLGALT1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","CXCL8","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","P3H1","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A"
  )

  # Filter genes to label 
label_genes <- df %>%
  filter(
    GeneSymbol %in% highlight_genes,
    abs(logFC) > 0.5,
    adj.P.Val < 0.05)


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
  labs(title = "Differential gene expression of EMT related genes in 
       TMPRSS2-ERG fusion Prostate Adenocarcinoma ",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

