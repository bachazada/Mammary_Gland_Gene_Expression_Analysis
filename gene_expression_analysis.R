#This study examines the expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Each group contains two biological replicates. We will first use the counts file as a starting point for our analysis. This data has already been aligned to the mouse genome. The command line tool featureCounts (Liao, Smyth, and Shi 2014) was used to count reads mapped to mouse genes from Refseq annotation (see the paper for details).

#### exercise 1 ####
rm(list=ls())
library(edgeR)
library(limma)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

library(Glimma)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Glimma")

library(gplots)
#install.packages("gplots")

library(org.Mm.eg.db)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")

library(RColorBrewer)
#install.packages("RColorBrewer")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("GO.db")
library(GO.db)

#### exercise 2 ####
# Read the data into R
seqdata <- read.delim("./tutorial_4_data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
dim(seqdata)
head(seqdata)
seqdata[1:5,1:5]
# Read the sample information into R
sampleinfo <- read.delim("./tutorial_4_data/SampleInfo.txt", stringsAsFactors = FALSE)
head(sampleinfo)
sampleinfo$CellType <- as.factor(sampleinfo$CellType)
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)

rownames(countdata) <- paste("ID_",seqdata[,1],sep="")
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
table(colnames(countdata)==sampleinfo$SampleName)

#### exercise 3 ####
# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)

dim(counts.keep)

#### exercise 4 ####
y <- DGEList(counts.keep)     #used by edge R
# have a look at y
y

# See what slots are stored in y
names(y)

# Library size information is stored in the samples slot
y$samples

#### Quality controls ####
y$samples$lib.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
pdf(file="bar_plot_library_sizes.pdf", width = 7, height=5)
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
#log2(cpm(y)+.5)[1:5,1:5]

# Check distributions of samples using boxplots
pdf(file = "distribution_log2_CPMs.pdf", width = 7, height = 5)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

# I'm going to write over the sampleinfo object with the corrected sample info
sampleinfo <- read.delim("data 2/SampleInfo_Corrected.txt")
sampleinfo
# Redo the MDSplot with corrected information
sampleinfo$CellType <- as.factor(sampleinfo$CellType)
sampleinfo$Status <- as.factor(sampleinfo$Status)
pdf(file = "MDS.pdf", width=10, height=5)
par(mfrow=c(1,2), mai=c(1,1,.3,.3))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")
dev.off()


#### exercise 5 ####
#### Hierarchical clustering with heatmaps ####

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
pdf("heatmap_most_variable.pdf", width=10, height=10)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
# Save the heatmap
#png(file="./plots/High_var_genes.heatmap.png")
#heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()
#dir.create("./plots")
#### Normalisation for composition bias ####

heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell,scale="row",distfun =  function(x){dist(x, method="manhattan")})
?dist


