#TGCA PanCancer Atlas PCa ######################################################
#Load TGCA PanCancer Atlas GSEA enriched pathways list 
library(ggplot2)
library(readr)
# Read txt file and sort and re-order 
gsea_data <- read_delim("Pan Cancer Atlas GSEA.txt", delim = "\t")
colnames(gsea_data) <- c("pathway", "NES")
gsea_data$pathway <- factor(gsea_data$pathway, levels = gsea_data$pathway[order(gsea_data$NES)])

# GSEA plot
ggplot(gsea_data, aes(x = pathway, y = NES, fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_gradient2(low = "darkblue", mid = "dodgerblue", high = "red", midpoint = 0) +
  labs(
    title = "Hallmark Pathways Enrichment in  
TMPSS2-ERG fusion Prostate Adenocarcinoma ",
    x = " ",
    y = "Normalized Enrichment Score (NES)",
    fill = "NES"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # Centered title
    axis.title.y = element_text(size = 14, face = "bold"),             # Y-axis title
    axis.title.x = element_text(size = 20, face = "bold"),             # X-axis title (optional)
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size = 18, face ="bold"), # Y-axis tick labels
    legend.title = element_text(size = 16, face = "bold"),             # Legend title
    legend.text = element_text(size = 14),                             # Legend labels
    legend.key.size = unit(1.2, "cm"),                                 # Bigger legend boxes
    panel.grid.major.y = element_blank()
  )




#SU2C mPCa######################################################################
# Read txt file and sort and re-order 
gsea_data2 <- read_delim("SU2C_GSEA.txt", delim = "\t")
colnames(gsea_data2) <- c("pathway", "NES")
gsea_data2$pathway <- factor(gsea_data2$pathway, levels = gsea_data2$pathway[order(gsea_data2$NES)])

# GSEA plot

ggplot(gsea_data2, aes(x = pathway, y = NES, fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_gradient2(low = "darkblue", mid = "blue", high = "dodgerblue", midpoint = -1.6) +
  labs(
    title = "Hallmark Pathways Enriched in TMPRSS2-ERG fusion in
    Metastatic Prostate Adenocarcinoma
    ",
    x = " ",
    y = "Normalized Enrichment Score (NES)",
    fill = "NES"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # Centered title
    axis.title.y = element_text(size = 14, face = "bold"),             # Y-axis title
    axis.title.x = element_text(size = 20, face = "bold"),             # X-axis title (optional)
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size = 18, face ="bold"), # Y-axis tick labels
    legend.title = element_text(size = 16, face = "bold"),             # Legend title
    legend.text = element_text(size = 14),                             # Legend labels
    legend.key.size = unit(1.2, "cm"),                                 # Bigger legend boxes
    panel.grid.major.y = element_blank()
  )








