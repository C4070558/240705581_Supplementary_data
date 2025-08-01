# Unravelling the role of ERG in EMT related genes in a prostate cancer _in vitro_ model 240705581_Supplementary_data
## Author 
Student No. 240705581
## Process
Scripts associated to each public data set can be found in Scripts folder. Before running script download data package ZIP file and place in R directory. 
## Data 
All the required data for the differentally expressed gene analysis can be sourced from using the reference links below from the public data sets.
For Gene Set Enrichment analysis files previously created can be found in he data folder.

### Prostate Adenocarcinoma (TCGA, PanCancer Atlas):
#### RSEM normalised count data  
https://xenabrowser.net/datapages/?dataset=TCGA.PRAD.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#### Sample matrix (Alterations file)
https://www.cbioportal.org/results/download?cancer_study_list=prad_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants%2Cgistic&case_set_id=prad_tcga_pan_can_atlas_2018_cnaseq&gene_list=ERG&geneset_list=%20&tab_index=tab_visualize&Action=Submit
### Metastatic Prostate Adenocarinoma (SU2C/PCF Dream Team, PNAS 2019): 
#### mRNA expression (FPKM capture) and Sample Matrix (Alterations file)
https://www.cbioportal.org/results/download?cancer_study_list=prad_su2c_2019&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants%2Cgistic&case_set_id=prad_su2c_2019_cnaseq&gene_list=ERG&geneset_list=%20&tab_index=tab_visualize&Action=Submit

## Results 
Result output for GSEA, pre-ranked gene list, differentally expressed gene list  has been attached in results folder 


## Software version 
* R version 4.5.1
* GSEA version 4.4.0 
* tidyverse version 2.0.0
* BiocManager version 1.30.26
* Bioconductor version 3.20 
* limma version 3.62.2
* ggplot2 version 3.5.2
* ggrepel version 0.9.6
* preprocessCore 1.66.0.
