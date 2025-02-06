#Differential expression analysis for WT and PAR1-LratCre mice treated with Vehicle and CCl4
# read in gene counts data
counts <- read.csv('gene_count.csv', row.names = 1)
head(counts)

# read in sample data
colData <- read.csv('RNAmetadata.csv', row.names = 1)
#View(colData)

# make sure sample names and orders match in counts and colData
all(colnames(counts) %in% rownames(colData))
all(colnames(counts) == rownames(colData))

# Create a DESeqDataSet object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~Genotype+Treatment+Genotype:Treatment)

# Relevel genotype and treatment factors
colData(dds)$Genotype <- relevel(colData(dds)$Genotype, ref = "WT")
colData(dds)$Treatment <- relevel(colData(dds)$Treatment, ref = "Vehicle")

# Check levels of Genotype factor - first level should be "WT"
levels(dds$Genotype)

# Check levels of Treatment factor - first level should be "Vehicle"
levels(dds$Treatment)

# Update the design formula with the interaction term
design(dds) <- ~ Genotype + Treatment + Genotype:Treatment

# pre-filtering to remove rows with read counts below 10
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

# Run DESeq analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)
plotMA(res)
write.csv(as.data.frame(res), file = "results.csv")
