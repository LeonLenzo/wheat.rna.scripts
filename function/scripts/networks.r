library(WGCNA)
library(DESeq2)
library(tximport)

# Load metadata
metadata <- read.csv("quant/cs_cds/metadata.csv")
metadata$condition <- factor(metadata$condition)

# Create named vector of file paths
files <- metadata$path
names(files) <- metadata$sample

# Import Kallisto quantifications
txi <- tximport(files, type = "kallisto", txOut = TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)

# Check total genes before filtering
nrow(dds)
# Apply filtering
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
# Check total genes after filtering
nrow(dds)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Load DESeq2 object (dds)
vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
datExpr <- assay(vsd)  # Extract transformed expression matrix

# Transpose to ensure genes are columns, samples are rows
datExpr <- t(datExpr)





# Ensure genes and samples have no missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
    # Remove problematic genes/samples
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Check sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub="", xlab="", cex.lab=0.8)

powers <- c(1:20)  # Test a range of powers
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-Free Topology Model Fit",
     type="o", main = "Scale-Free Topology Fit")
abline(h=0.9, col="red")




softPower <- 17  # Chosen power

datExpr <- datExpr[, geneVar > 0.01]  # Remove low-variance genes
dim(datExpr)  # Check new dimensions

# Calculate variance for each gene
geneVar <- apply(datExpr, 2, var)

# Select the top 5000 most variable genes
numGenesToKeep <- 10000  # Adjust based on memory limits
topGenes <- order(geneVar, decreasing = TRUE)[1:numGenesToKeep]

# Subset the dataset
datExpr <- datExpr[, topGenes]
dim(datExpr)  # Should now be (samples, ~5000 genes)




# Construct adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Convert adjacency to Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Perform hierarchical clustering of genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot dendrogram to visualize gene clustering
plot(geneTree, main = "Gene Clustering Dendrogram", sub="", xlab="")




# Identify modules using dynamic tree cutting
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

# Convert module numbers to colors for visualization
moduleColors <- labels2colors(dynamicMods)

# Plot the dendrogram with module colors
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




# Calculate eigengenes (summary expression patterns of modules)
MEList <- moduleEigengenes(datExpr, colors = moduleColors)
MEs <- MEList$eigengenes

# Calculate correlation between modules
moduleDissimilarity <- 1 - cor(MEs)

# Perform hierarchical clustering on module eigengenes
mergeTree <- hclust(as.dist(moduleDissimilarity), method = "average")

# Plot the module clustering tree
plot(mergeTree, main = "Clustering of Module Eigengenes", xlab="", sub="")

# Merge modules that have high correlation (>0.75)
mergeThreshold <- 0.25  # Adjust if needed
mergedColors <- mergeCloseModules(datExpr, moduleColors, cutHeight = mergeThreshold, verbose = 3)$colors

# Update module colors after merging
moduleColors <- mergedColors



# Load your sample trait data (ensure row names match sample names)
traitData <- read.csv("networks/treatments.csv", row.names = 1)

# Calculate module eigengenes again
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# Compute correlation between module eigengenes and traits
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

png("Module_Trait_Heatmap.png", width = 1200, height = 800)  # Adjust size
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               main = "Module-Trait Relationships")
dev.off()  # Closes the plotting device