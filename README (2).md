# Transcriptomic Analysis of Breast Cancer Subtypes Using TCGA RNA-seq Data

This project investigates differential gene expression between Basal and Luminal A subtypes in breast cancer using RNA-seq data from TCGA-BRCA. We performed preprocessing, filtering, DESeq2-based differential expression analysis, gene annotation using Ensembl and HGNC symbols, and visualization using volcano plots. A special focus was given to immune-related genes.

## Objective

To identify differentially expressed genes (DEGs) between Basal and Luminal A subtypes of breast cancer, and assess the immune-related transcriptomic differences.

## Data Source

- Project: TCGA-BRCA
- Data Type: RNA-seq gene expression quantification (STAR - Counts)
- Samples: Primary Tumor
- Subtypes: Basal (197 samples), Luminal A (571 samples)

## Tools and Packages Used

- **R packages**: TCGAbiolinks, DESeq2, SummarizedExperiment, biomaRt, EnhancedVolcano, tidyverse
- **Visualization**: ggplot2, EnhancedVolcano
- **Data annotation**: biomaRt

## Workflow Overview

1. **Data Retrieval**  
   Retrieved STAR-counts RNA-seq data from TCGA-BRCA project using `TCGAbiolinks`.

2. **Preprocessing & Filtering**  
   - Filtered samples by subtype (Basal vs Luminal A)
   - Cleaned metadata and matched with gene count matrix

3. **Differential Expression Analysis**  
   - Used DESeq2 for DE analysis between subtypes
   - Applied filters (adjusted p-value < 0.05)

4. **Gene Annotation**  
   - Converted Ensembl IDs to gene symbols using biomaRt

5. **Visualization**  
   - Created volcano plot of all DEGs
   - Highlighted immune-related genes (e.g., CD8A, IFNG, PDCD1, CD274, CTLA4, LAG3, etc.)

6. **Key Results**  
   - Identified 32,060 DEGs (padj < 0.05)
   - Multiple immune genes were significantly upregulated in Basal subtype

## Key Output Files

- `Basal_vs_Luminal_Annotated_DEGs.csv`: Annotated DEGs with gene symbols
- `Immune_DEGs_Table.png`: Highlighted immune DEGs
- `VolcanoPlot_AllGenes.png`: All DEGs in volcano plot

## Example Code Snippet

```r
EnhancedVolcano(res_annotated,
                lab = res_annotated$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Basal vs Luminal - TCGA BRCA',
                pCutoff = 0.05,
                FCcutoff = 1.0)
```

## Author

Renuka Reddy Namala  
Master’s in Bioinformatics and Computational Biology  
University of South Florida

## License

This repository is for educational and academic project demonstration purposes only.