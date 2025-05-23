# -------------------------------
# TCGA BRCA Basal vs Luminal DEG Analysis
# Author: Renuka Reddy Namala
# -------------------------------

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "DESeq2", "SummarizedExperiment", "EnhancedVolcano", "biomaRt"))

library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(biomaRt)

# -------------------------------
# STEP 1: Download and Prepare TCGA Data
# -------------------------------
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

GDCdownload(query)
data <- GDCprepare(query)

# -------------------------------
# STEP 2: Filter Basal and Luminal A Samples
# -------------------------------
metadata <- colData(data)

# Grouping subtypes
metadata$subtype_group <- NA
metadata$subtype_group[metadata$paper_BRCA_Subtype_PAM50 == "Basal"] <- "Basal"
metadata$subtype_group[metadata$paper_BRCA_Subtype_PAM50 == "LumA"] <- "Luminal"

# Clean barcodes for matching
metadata$barcode <- substr(metadata$barcode, 1, 12)
rownames(metadata) <- metadata$barcode

# Extract and align count matrix
count_matrix <- assay(data)
colnames(count_matrix) <- substr(colnames(count_matrix), 1, 12)

# Subset metadata and count matrix
filtered_meta <- metadata[metadata$subtype_group %in% c("Basal", "Luminal"), ]
filtered_counts <- count_matrix[, colnames(count_matrix) %in% filtered_meta$barcode]
filtered_meta <- filtered_meta[match(colnames(filtered_counts), filtered_meta$barcode), ]

# -------------------------------
# STEP 3: Run DESeq2 Analysis
# -------------------------------
filtered_meta$subtype_group <- factor(filtered_meta$subtype_group)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = filtered_meta,
  design = ~ subtype_group
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("subtype_group", "Basal", "Luminal"))

# Save significant DEGs
sig_res <- res[which(res$padj < 0.05), ]
sig_res <- sig_res[order(sig_res$padj), ]
write.csv(as.data.frame(sig_res), "DEGs_Basal_vs_Luminal.csv")

# -------------------------------
# STEP 4: Annotate with Gene Names
# -------------------------------
res_cleaned <- as.data.frame(res)
res_cleaned$ensembl_id <- gsub("\\..*", "", rownames(res_cleaned))

# Use biomaRt to map Ensembl to HGNC symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = res_cleaned$ensembl_id,
  mart = mart
)

# Merge annotations
res_annotated <- merge(res_cleaned, gene_map, by.x = "ensembl_id", by.y = "ensembl_gene_id")
write.csv(res_annotated, "Basal_vs_Luminal_Annotated_DEGs.csv", row.names = FALSE)

# -------------------------------
# STEP 5: Volcano Plot (All Genes)
# -------------------------------
EnhancedVolcano(res_annotated,
                lab = res_annotated$hgnc_symbol,
                x = "log2FoldChange",
                y = "padj",
                title = "Basal vs Luminal - TCGA BRCA",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.5
)

# -------------------------------
# STEP 6: Highlight Immune Markers
# -------------------------------
immune_genes <- c("CD8A", "GZMB", "CXCL10", "PDCD1", "CD274", "IFNG", "CTLA4", "LAG3")
res_immune <- res_annotated[res_annotated$hgnc_symbol %in% immune_genes, ]
print(res_immune)

EnhancedVolcano(res_annotated,
                lab = res_annotated$hgnc_symbol,
                selectLab = immune_genes,
                x = "log2FoldChange",
                y = "padj",
                title = "Immune Gene Upregulation in Basal vs Luminal - TCGA BRCA",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.5
)
 