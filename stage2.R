###########################################
#### Stage 2 task - HackBio Internship #### 
###########################################

# set your working directory 
setwd("C:/Users/agdad/Desktop/HackBio/stage2")

# load the glioblastoma dataset
gbm <- read.csv(file = 'glioblastoma.csv', header = TRUE, row.names = 1)

# 1.generate a heatmap for the entire dataset ####
library(gplots)
library(RColorBrewer) # for color palettes 

# diverging color palette -> for normalized data
div_color_palette <- rev(brewer.pal(11, "RdBu")) # red for up-regulation and blue for down-regulation

# sequential color palette -> for raw counts
seq_color_palette <- brewer.pal(9, "Blues") # blue gradient where light blue represents low counts


# heatmap with sequential palette 
heatmap.2(as.matrix(gbm),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap with diverging palette 
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by genes
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = T, Colv = F, dendrogram = "row",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by samples
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = F, Colv = T, dendrogram = "column",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by genes and samples
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = T, Colv = T, dendrogram = "both",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))


# 2.subset genes that are significantly up- and down-regulated ####
# first 5 samples (02) are recurrent tumors, last 5 samples (01) are primary tumors
# PCA plot for confirmation 
# change names of the samples in a duplicated df for simplicity
gbm_copy <- gbm
colnames(gbm_copy) <- c("02A_1","02A_2","02A_3","02A_4","02A_5","01A_1","01A_2","01B_3","01A_4","01A_5")

pca <- prcomp(t(gbm_copy), center = TRUE, scale. = TRUE)
plot(pca) # how much variation is explained by each component
PC1_and_PC2 <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], type = rownames(pca$x))
library(ggplot2)
ggplot(PC1_and_PC2,
       aes(x = PC1, y = PC2, col = type)) +
  geom_point() +
  geom_text(aes(label = type),
            hjust = 0, vjust = 0)

# create the df for "metadata"
sample_code <- c("TCGA.19.4065.02A.11R.2005.01", "TCGA.19.0957.02A.11R.2005.01", "TCGA.06.0152.02A.01R.2005.01", "TCGA.14.1402.02A.01R.2005.01", "TCGA.14.0736.02A.01R.2005.01", "TCGA.06.5410.01A.01R.1849.01", "TCGA.19.5960.01A.11R.1850.01", "TCGA.14.0781.01B.01R.1849.01", "TCGA.02.2483.01A.01R.1849.01", "TCGA.06.2570.01A.01R.1849.01")
tumor_stage <- c("recurrent", "recurrent", "recurrent", "recurrent", "recurrent", "primary", "primary", "primary", "primary", "primary")
meta <- data.frame(sample_code, tumor_stage)
meta
meta$tumor_stage <- as.factor(meta$tumor_stage)

# differential expression analysis
library(DESeq2)

dds <- DESeqDataSetFromMatrix( # dds is the DESeq2 dataset object
  countData = gbm,
  colData = meta,
  design = ~ tumor_stage)

dds <- DESeq(dds) # differential expression analysis

res <- results(dds) # results of the differential expression analysis 
res

# filtering
up_genes <- res[which(res$padj < 0.05 & res$log2FoldChange > 1), ] # 2 times change in expression
up_genes@rownames # list of upregulated genes
down_genes <- res[which(res$padj < 0.05 & res$log2FoldChange < -1), ]
down_genes@rownames # list of downregulated genes 

background <- as.factor(res@rownames) 
background # list of background genes
