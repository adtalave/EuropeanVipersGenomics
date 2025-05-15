library(DESeq2)
library(ggplot2)
library(tximport)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


cts <- as.matrix(read.table("all_counts_f_alliso_strict.txt",
                            sep="\t",
                            header = TRUE,
                            row.names=1,
                            stringsAsFactors = F))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = df,
                              design = ~ venom)
cts <- cts[,2:(ncol(cts))]

df$venom <- factor(df$venom, levels=c("venom", "non_venom")) #to classify libraries belonging to venomous or non-venomous organs
df[,2]=c("venom","non_venom","non_venom","non_venom","non_venom","non_venom","non_venom","venom","non_venom","non_venom","non_venom","non_venom","venom","non_venom","non_venom","venom","non_venom","non_venom","non_venom")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = df,
                              design = ~ venom)


dds@colData
saveRDS(dds, file = "DESeq2.rds")

dds <- estimateSizeFactors(dds)
# extract size factors
sizeFactors(dds)
# plot histogram
hist(sizeFactors(dds),
     breaks=6, col = "cornflowerblue",
     xlab="Size factors", ylab="No. of samples",
     main= "Size factor distribution over samples")

# calculate normalized counts
counts_norm <- counts(dds, normalized=TRUE)

###DESeq2 normalized read counts are not normalized for gene length, so cannot be used to compare expression levels between genes within the same sample.

# drop genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# check the new dimensions
dim(dds)
rld <- rlog(dds, blind = FALSE)


########################

# set plotting window to 1 row vs 2 columns
par(mfrow=c(1,1))

# plot standard log counts
cts <- counts(dds, normalized=FALSE)
plot(log2(cts[,1]+1), log2(cts[,2]+1), col = "cornflowerblue", xlab = "Sample 1", ylab = "Sample 2", main = "Log2 + 1")

# plot rlog counts
plot(assay(rld)[,2], assay(rld)[,2], col = "indianred", xlab = "Sample 2", ylab = "Sample 2", main = "rlog")


####PCA


# calculate gene expression level variance between samples
var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])

# reset plotting window to 1 row vs 1 columns
par(mfrow=c(1,1))

# plot variance for genes accross samples
plot(var,
     las = 1,
     main="Sample gene expression variance",
     xlab = "Gene", ylab = "Variance")

# add vertical lines at specific gene number indexes
abline(v=1000, col="red")
abline(v=600, col="green")
abline(v=300, col="blue")

# modify variable feature number to be used in PCA and hierarchical clustering based on no. of most variable features
var_feature_n <- 1000

# calculate the row variance
rv <- rowVars(assay(rld))

# order variance by size and select top 1000 with most variance
select <- order(rv, decreasing = TRUE)[1:1000]

# subset rlog values for genes with top variance ranks
rld_sub <- assay(rld)[select, ]

# transpose the matrix (rows to columns and columns to rows)
rld_sub <- t(rld_sub)

# run principal components analysis
pca <- prcomp(rld_sub)

# extract the variance explained by each PC
percentVar <- pca$sdev^2/sum(pca$sdev^2)

# subset for first 5 elements
percentVar <- percentVar[1:5]

# give the string names
names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# plot variance for top 10 PCs
barplot(percentVar, col = "indianred", las = 1, ylab = "% Variance", cex.lab = 1.2)

# construct data frame w/ PC loadings and add sample labels
pca_df <- as.data.frame(pca$x)

# add a column containing venom
pca_df$venom <- dds@colData$venom

# add column containing sample IDs
pca_df$sample_ids <- colnames(dds)

# add colors for plotting to df
pca_df$col <- NA
for(i in 1:length(levels(pca_df$venom))){
  ind1 <- which(pca_df$venom == levels(pca_df$venom)[i])
  pca_df$col[ind1] <- i
}

# plot PC1 vs PC2
plot(pca_df[, 1], pca_df[, 2],
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1,
     panel.first = grid(),
     col=pca_df$col)

# add sample names to data points
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$sample, cex=0.6, font=2, pos=1)




#####heatmap

var_feature_n <- 1000

# calculate the row variance
rv <- rowVars(assay(rld))

# order variance by size and select top 1000 with most variance
select <- order(rv, decreasing = TRUE)[1:1000]

# select top X no. of variable genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), var_feature_n)

# set up gene expression matrix
mat1 <- assay(rld)[topVarGenes,]

# set up colors for heatmap
col = colorRamp2(c(0,7), c("white","orange"))
cols1 <- brewer.pal(11, "Paired")

# set up annotation bar for samples
ha1 = HeatmapAnnotation(Group = colData(dds)$venom,
                        col = list(Group = c("venom" = "grey40",
                                             "non_venom" = "grey90")),
                        show_legend = TRUE)

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData(dds)$sample,
                                    which="column", rot = 45,
                                    gp = gpar(fontsize = 10)))
# generate heatmap object
ht1 = Heatmap(mat1,
              name = "Expression",
              col = col,
              top_annotation = c(ha1),
              bottom_annotation = c(ha),
              show_row_names = FALSE,
              show_column_names = FALSE)

# plot the heatmap
draw(ht1, row_title = "Genes", column_title = "Top 1000 most variable genes")

#####scaling with Zscore

# scale matrix by each row
mat_scaled = t(apply(mat1, 1, scale))

# set up colors for heatmap
col = colorRamp2(c(-1,-0.5,0, 1,2,4), c("#281579","#2670D1", "white", "#FFBE5C","#E57B24","#E13535"))

# generate heatmap object
ht1 = Heatmap(mat_scaled, name = "Expression", col = col,
              top_annotation = c(ha1),
              bottom_annotation = c(ha),
              show_row_names = FALSE,
              show_column_names = FALSE)

draw(ht1, row_title = "Genes", column_title = "Top 1000 most variable genes")
