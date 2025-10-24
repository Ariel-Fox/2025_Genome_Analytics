

#load required libraries
library(data.table)
library(DESeq2)

# load counts table from GEO2R url
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

#define path for access to RNA-seq count data for GEO accession GSE128749 [file imports the raw 
#count data in a TSV format]. While this command concatenates the separate parameters with "&". 
path <- paste(urld, "acc=GSE128749", "file=GSE128749_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")

#reads in large data file from path variable defined upstream, accounting for column name headers.
#All columns are read in as integers [these are counts of reads mapped to target genes]
tbl <- as.matrix(fread(path, header=T, colClasses="integer"), rownames="GeneID")

#getting a sense of what the data looks like
head(tbl)

#create concatenated variable name for for GEO2R to access GRCh38.p13 version of 
#human genome annotations in a tsv format
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")

#Read in file from GRCh38.p13 file path, accounting for header, treating quotes 
#as text characters, not convert strings into factors
annot <- fread(apath, header=T, quote="", stringsAsFactors=F)
#convert format from default data table to a data frame.
annot <- as.data.frame(annot)

#define row names as GeneID column, matching row names from tbl samples
head(annot)
rownames(annot) <- annot$GeneID

# select both R1881 [synthetic androgen methyltrienolone] and control [96% Ethanol] in LNCaP cells
gsms <- "XXXXX000111"
#split selection in a vector c(,)
sml <- strsplit(gsms, split="")[[1]]


# exclude first 5 samples LAPC4 cells from split vector and select samples that are LNCaP cells
sel <- which(sml != "X")
sml <- sml[sel]
#select raw data from table for LNCaP cells
tbl <- tbl[ ,sel]

# create factors for sml (000 and 111 )
gs <- factor(sml)
#converts vector into names ( not necessary in this case; no additional formatting needed)
groups <- make.names(c("control","test"))
#apply (replace) control and test to 000 and 111 samples
levels(gs) <- groups
#create data frame labeling samples with either test or control group
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

#creates a DEseqdataset object fro prep in differential analysis of the filtered LNCaP cells
#falling under test or control sample 
ds <- DESeqDataSetFromMatrix(countData=tbl, 
                             colData=sample_info, 
                             design= ~Group)

#using the DeSeqDataSeq object containing filtered read counts with assigned 
#test or control grouping. Using the default statistical testing algorithm doing 
#a simple comparison of the two levels for the treated and untreated LNCap using
#positive count normalization
ds <- DESeq(ds, test="Wald", sfType="poscount")


#provide statistical summary of Differential Sequencing and compares the control and test groups
#the significance threshold is 0.05 and the p-adjusted value uses FDR for multiple testing corrections
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")


#order the results for the top 10 genes based on the adjusted p-value
tT <- r[order(r$padj)[1:10],] 
#merge the ordered results with annotation data for top 10 genes
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)
#select annotation columns of interest in addition to statistical test values
tT <- subset(tT, select=c("Symbol","GeneID","padj","pvalue","stat","log2FoldChange","Description"))

print(tT)


all_tT <- r[order(r$padj),] 


# Assuming r is the DESeq2 results object with log2 fold changes
positive_fc_genes <- all_tT[which(all_tT$padj<=0.05), ]
positive_fc_genes <- positive_fc_genes[which(positive_fc_genes$log2FoldChange>=1), ]

negative_fc_genes <- all_tT[which(all_tT$padj<=0.05), ]
negative_fc_genes <- negative_fc_genes[which(negative_fc_genes$log2FoldChange<= -1), ]


# Count the number of genes with positive log2 fold change
num_positive_fc_genes <- nrow(positive_fc_genes)
num_negative_fc_genes <- nrow(negative_fc_genes)

# Print the result
print(num_positive_fc_genes)
print(num_negative_fc_genes)


head(all_tT)









write.table(tT, file=stdout(), row.names=F, sep="\t")

plotDispEsts(ds, main="GSE128749 Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE128749 Frequencies of padj-values")

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(groups[1], "vs", groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

# MD plot
par(mar=c(4,4,2,1), cex.main=1.5)
plot(log10(r$baseMean), r$log2FoldChange, main=paste(groups[1], "vs", groups[2]),
     xlab="log10(mean of normalized counts)", ylab="log2FoldChange", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log10(baseMean), log2FoldChange, pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)
abline(h=0)
palette(old.pal) # restore palette

################################################################
#   General expression data visualization
dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts

# box-and-whisker plot
lbl <- "log10(raw counts + 1)"
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
boxplot(dat[,ord], boxwex=0.6, notch=T, main="GSE128749", ylab="lg(norm.counts)", outline=F, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# UMAP plot (multi-dimensional scaling)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
par(mar=c(3,3,2,6), xpd=TRUE, cex.main=1.5)
ump <- umap(t(dat), n_neighbors = 3, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=groups, pch=20,
       col=1:length(groups), title="Group", pt.cex=1.5)