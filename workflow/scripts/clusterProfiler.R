# Log file open
log_file <- snakemake@log[[1]]
log_stream <- file(log_file, open="wt")
sink(log_stream, append=T)
sink(log_stream, append=T, type="message")

#### DEPENDENCIES ####
# Package GO enrichment analysis
library(clusterProfiler)

# Package for annotation
library(AnnotationDbi)
library(org.Hs.eg.db)

#### GO Enrichment Analysis ####$
# Read Top Differentially Expressed Genes
TopGenes <- read.table(file=snakemake@input[["TopGenes"]], header=T, sep="\t")

# GO enrichement over Ensembl Gene ID for Biological Process with strict parameters
TopGO <- enrichGO(gene=rownames(TopGenes),
                OrgDb=org.Hs.eg.db,
                keyType="ENSEMBL",
                ont="BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = T)

# Plot results
dir.create(snakemake@params[["report_dir"]], recursive = TRUE, showWarnings = FALSE)
png(file=paste0(snakemake@params[["report_dir"]],"/dotplot.png"))
dotplot(TopGO)
dev.off()

png(file=paste0(snakemake@params[["report_dir"]],"/goplot.png"))
goplot(TopGO)
dev.off()

png(file=paste0(snakemake@params[["report_dir"]],"/cnetplot.png"))
cnetplot(TopGO)
dev.off()

# Save result GO enrichment
write.table(as.data.frame(TopGO),
            file=snakemake@output[["TopGO"]],
            row.names=T,
            col.names=T,
            sep="\t")


# Log file close
sink()
sink(type = "message")
close(log_stream)