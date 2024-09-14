###########################################
#### Stage 2 task - HackBio Internship #### 
###########################################

# load the glioblastoma dataset
gbm <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv", row.names = 1)

# 1.generate a heatmap for the entire dataset ####
library(gplots) # for heatmap generation
library(RColorBrewer) # for color palettes 

# diverging color palette 
div_color_palette <- rev(brewer.pal(11, "RdBu")) # red for up-regulation and blue for down-regulation

# sequential color palette 
seq_color_palette <- brewer.pal(9, "Blues") # blue gradient where light blue represents low counts


# heatmap with sequential palette 
heatmap.2(as.matrix(gbm),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top\n 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap with diverging palette 
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top\n 500+ differentially expressed genes in Glioblastoma",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by genes
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = T, Colv = F, dendrogram = "row",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top\n 500+ differentially expressed genes in Glioblastoma\n clustering by genes",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by samples
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = F, Colv = T, dendrogram = "column",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top\n 500+ differentially expressed genes in Glioblastoma\n clustering by samples",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))

# heatmap clustering by genes and samples
heatmap.2(as.matrix(gbm),
          col = div_color_palette,
          Rowv = T, Colv = T, dendrogram = "both",
          trace = "none",
          key = T, scale = "row",
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top\n 500+ differentially expressed genes in Glioblastoma\n clustering by genes and samples",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))


# 2.subset genes that are significantly up- and down-regulated ####
# getting column names
colnames(gbm) 
# selecting the two groups using index position 
group1 <- c(1, 2, 3, 4, 5)
group2 <- c(6, 7, 8, 9, 10)

# df with the groups data
group1_gbm <- gbm[, group1]
group2_gbm <- gbm[, group2]

# means of the two groups
group1_mean <- rowMeans(group1_gbm)
group2_mean <- rowMeans(group2_gbm)

# fold change 
gbm$log2fold_change <- log2(group2_mean) - log2(group1_mean)

# p-values
gbm$pvalues <- apply(gbm, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
})

# visualize the fold change and -log10pvalue
plot(gbm$log2fold_change, -log10(gbm$pvalues))
abline(h = 0.6, col = "red", lty = 2)
abline(v = 1.5, col = "blue", lty = 2)
abline(v = -1.5, col = "blue", lty = 2)

# subsetting of up- and down-regulated genes 
up_genes <- gbm[gbm$log2fold_change > 1.5 & gbm$pvalues < 0.25, ]
rownames(up_genes) # list of up-regulated genes
down_genes <- gbm[gbm$log2fold_change < -1.5 & gbm$pvalues < 0.25, ]
rownames(down_genes) # list of down-regulated genes 

# read enrichment.csv downloaded from ShinyGO after functional enrichment analysis using GO biological processes (FDR cutoff 0.05)
# the enrichment analysis was done on all deregulated genes (57 up-regulated and 100 down-regulated) and sorted by fold enrichment
enrichment <- read.csv("enrichment.csv")

# select the top 5 pathways by fold enrichment 
top5_index <- c(1, 2, 3, 4, 5)
top5 <- enrichment[top5_index, c("Enrichment.FDR", "nGenes", "Fold.Enrichment", "Pathway")]

# calculate -log10pvalue
top5$Neglog10FDR <- -log10(top5$Enrichment.FDR)

# visualizing the top5 down-regulated pathways using a dot plot
library(ggplot2)
# factorizing Pathway column to ensure that ggplot2 uses the order that I want for the plot
top5$Pathway <- factor(top5$Pathway, levels = rev(top5$Pathway)) 

# dot plot code
ggplot(top5, aes(x = nGenes, y = Pathway)) + 
  geom_point(aes(size = Neglog10FDR), color = "darkgreen") +
  geom_segment(aes(x = 0, xend = nGenes, yend = Pathway), color = "lightgrey", linetype = "dotted") +
  labs(x = "Number of Genes", y = "Pathways", size = "-log10(p-value)") + 
  theme_minimal() + ggtitle("Top 5 Enriched Pathways \n in Glioblastoma") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.3))
