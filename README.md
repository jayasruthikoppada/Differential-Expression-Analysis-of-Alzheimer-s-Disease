# **Project Title: Differential Expression Analysis of Alzheimer’s Disease**

# **Overview**
This project performs a differential expression analysis (DEA) on RNA-Seq data to identify differentially expressed genes (DEGs) between:
1. **Healthy vs. Mild Cognitive Impairment (MCI)**
2. **Healthy vs. Alzheimer’s Disease (AD)**

The analysis uses the DESeq2 package in R, statistical models, and pathway analysis to uncover critical genes and biological mechanisms associated with Alzheimer’s Disease and its precursor, MCI.

### **Motivation**
The motivation of this project is to identify highly expressed genes in Alzheimer’s and Mild Cognitive Impairment conditions compared to healthy states to improve our understanding of disease progression and potential therapeutic targets.

# **Goals**
- Identify DEGs in **Healthy vs. MCI** and **Healthy vs. AD** datasets.
- Perform functional pathway analysis to highlight affected biological processes.
- Visualize results using **volcano plots**, **heatmaps**, and **pathway analysis plots**.


### **Dependencies**
This project uses R packages:
- `DESeq2`
- `dplyr`
- `reshape`
- `ggplot2`
- `ggrepel`
- `RColorBrewer`
- `pheatmap`

Install the required packages using:
```R
install.packages(c("dplyr", "reshape", "ggplot2", "ggrepel", "RColorBrewer", "pheatmap"))
BiocManager::install("DESeq2")
```
# **Data**
# Input Files:
1. **coldata.HvsM.txt**: Metadata for Healthy vs. Mild Cognitive Impairment.
2. **rawcounts.HvsM.txt**: Count data for Healthy vs. Mild Cognitive Impairment.
3. **coldata.HvsAD.txt**: Metadata for Healthy vs. Alzheimer's Disease.
4. **rawcounts.HvsAD.txt**: Count data for Healthy vs. Alzheimer's Disease.

# Output Files:
1. **DEGs_HvsM_results.csv**: Differential expression results for Healthy vs. MCI.
2. **DEGs_HvsAD_results.csv**: Differential expression results for Healthy vs. AD.
3. Heatmaps and Volcano plots:
   - `H_MCI_volcano.png`
   - `H_AD_volcano.png`
   - `H_MCI_AD_logFC.png`

# **Pipeline**

# 1. **Data Preprocessing**
- Load count and metadata files.
- Filter out low-count genes (genes with fewer than 5 reads across all samples).

# 2. **Differential Expression Analysis**
Using DESeq2:
- Model gene expression levels for **Healthy vs. MCI** and **Healthy vs. AD**.
- Apply thresholds:
  - Adjusted p-value (`padj`) < 0.05.
  - |log2 Fold Change| > 1 for MCI and > 2 for AD.
- Save DEGs to `.csv` files.

#### 3. **Visualization**
- Volcano plots to show significant DEGs.
- Heatmaps for top 20 DEGs.
- Pathway analysis using EnrichR with KEGG.

---

### **Key Scripts**

#### R Code for Healthy vs. MCI
```R
# Load data
coldata.1 <- read.table("coldata.HvsM.txt", header = TRUE, row.names = 1)
cts.1 <- read.table("rawcounts.HvsM.txt", header = TRUE, row.names = 1)

# Create DESeq2 object and filter low-count genes
dds.1 <- DESeqDataSetFromMatrix(countData = cts.1, colData = coldata.1, design = ~ condition)
dds.1 <- dds.1[rowSums(counts(dds.1)) >= 5,]

# Perform DEA
dds.1$condition <- factor(dds.1$condition, levels = c("Healthy", "Mild_cog"))
dds.1 <- DESeq(dds.1)
res.1 <- results(dds.1)

# Export results
degs.HvsM <- as.data.frame(subset(res.1, abs(log2FoldChange) >= 1 & padj <= 0.05))
write.csv(degs.HvsM, "DEGs_HvsM_results.csv")
```

#### R Code for Visualization: Volcano Plot
```R
library(ggplot2)
library(ggrepel)

# Volcano plot
ggplot(res.1, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)) +
  geom_point() +
  geom_text_repel(data = subset(res.1, padj < 0.05 & abs(log2FoldChange) > 1),
                  aes(label = rownames(subset(res.1, padj < 0.05 & abs(log2FoldChange) > 1)))) +
  labs(title = "Healthy vs Mild Cognitive Impairment",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_classic()
```

#### R Code for Heatmap
```R
library(pheatmap)

# Top 20 DEGs for Healthy vs MCI
top20_HvsM <- head(degs.HvsM[order(degs.HvsM$padj),], 20)
normalized_counts <- counts(dds.1, normalized = TRUE)
pheatmap(normalized_counts[rownames(top20_HvsM), ], scale = "row")
```

---

### **Results**
- **Healthy vs MCI**: Identified **944 DEGs**.
- **Healthy vs AD**: Identified **739 DEGs**.
- Significant pathways:
  - Cytokine-cytokine receptor interaction (Healthy vs. MCI).
  - Rheumatoid arthritis pathway (Healthy vs. AD).

---

### **Figures**
1. **Volcano Plots**:
   - `H_MCI_volcano.png`: Healthy vs. MCI.
   - `H_AD_volcano.png`: Healthy vs. AD.
2. **Heatmap**:
   - `H_MCI_AD_logFC.png`: Top 20 DEGs combined.

---

### **Contact**
For questions or collaboration opportunities, contact:
- **Jaya Sruthi Koppada**: sruthikoppada@gmail.com
