# Log file open
log_file <- snakemake@log[[1]]
log_stream <- file(log_file, open="wt")
sink(log_stream, append=T)
sink(log_stream, append=T, type="message")

#### DEPENDENCIES ####
# Package for differential analysis
library(DESeq2)
library(regionReport)

# Package annotation
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# Package vizualisation
library(regionReport)
library(EnhancedVolcano)

#### FROM CALLISTO H5 TO COUNT MATRIX with tximport ####
# Sample data to df with rownames set to sample
samples <- read.table(file = snakemake@config[["samples"]],
                      header = TRUE, 
                      sep = "\t",
                      )

# Set condition as factor
samples$condition <- factor(samples$condition)

# Load counts
load(file=snakemake@input[["counts"]])

#### DESeq2 Diffrential gene expression ####
# Set samples meta data rownames as same count matrix
rownames(samples) <- colnames(txi$counts)

# Build DESeq2 data structure from tximport count matrix.
dds <- DESeqDataSetFromTximport(txi,
                                colData=samples,
                                design=~condition)

# Reference level to BSA (control)
dds$condition <- relevel(dds$condition, ref="BSA")

# Differential expression analysis following parameters recomended for single cell
dds <- DESeq(object=dds,
             test="LRT",
             reduced=~1,
             useT = T,
             minmu=1e-6,
             minReplicatesForReplace=Inf)

# Dispersion estimae for model fitting assessment
dir.create(path = snakemake@params[["report_dir"]], recursive = T)
png(file=paste0(snakemake@params[["report_dir"]],"/disp.png"))
plotDispEsts(dds)
dev.off()

# Result table from dds object
res <- results(dds, name="condition_INFy_vs_BSA")
res <- na.omit(res)
res$GENENAME <- select(EnsDb.Hsapiens.v86,
                     keys=rownames(res),
                     columns="GENENAME",
                     keytype = "GENEID")$GENENAME

# Plot volcano with symbols
png(file=paste0(snakemake@params[["report_dir"]],"/volcano.png"))
EnhancedVolcano(res, x="log2FoldChange",y="padj", lab=res$GENENAME)
dev.off()

# Order by smallest pvalues
ordered_res <- res[order(res$pvalue),]
# Log FC shrinkage for ranking-sum
LFC_shrink <- lfcShrink(dds, coef="condition_INFy_vs_BSA", type="apeglm")
LFC_shrink <- na.omit(LFC_shrink)
LFC_shrink$GENENAME <- select(EnsDb.Hsapiens.v86,
                            keys=rownames(LFC_shrink),
                            columns = "GENENAME",
                            keytype = "GENEID")$GENENAME

# Differnetially expressed GENES set
TopGenes <- LFC_shrink[LFC_shrink$padj < 0.1,]




# Save files
save(ordered_res,
     file=snakemake@output[["raw"]])
save(LFC_shrink,
     file=snakemake@output[["shrink"]])
write.table(TopGenes,
            file=snakemake@output[["TopGenes"]],
            sep = "\t",
            row.names = T,
            col.names = T)

# Report of Differential Gene Epression analysis
DESeq2Report(dds=dds,
             project="DESeq2 Differential Gene Expression Analysis",
             intgroup="condition",
             outdir=snakemake@params[["report_dir"]])

# Log file close
sink()
sink(type = "message")
close(log_stream)