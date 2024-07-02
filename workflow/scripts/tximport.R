# Log file open
log_file <- snakemake@log[[1]]
log_stream <- file(log_file, open="wt")
sink(log_stream, append=T)
sink(log_stream, append=T, type="message")

#### DEPENDENCIES ####
# Packages for h5 into count matrix
library(tximport)
library(rhdf5)
library(readr)

# Package for annotation
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#### FROM CALLISTO H5 TO COUNT MATRIX with tximport ####
# Sample data to df with rownames set to sample
samples <- read.table(file = snakemake@config[["samples"]],
                      header = TRUE,
                      sep = "\t")

# Set condition as factor
samples$condition <- factor(samples$condition)

# H% files path
fastq_files <- snakemake@input[["all_h5"]]
# Name files by sample name
names(fastq_files) <- samples$sample

# Transcript ID from Ensembl transcriptome to Gene ID with org.Hs
tx_keys <- keys(EnsDb.Hsapiens.v86, keytype = "TXID")
tx2gene <- select(EnsDb.Hsapiens.v86, keys=tx_keys, columns="GENEID", keytype="TXID")

# tximport count table for DESeq2 analysis with Ensembl Gene ID
txi <- tximport(fastq_files,
                type="kallisto",
                txOut=FALSE,
                tx2gene=tx2gene,
                ignoreTxVersion = T)

save(txi, file=snakemake@output[["counts"]])

# Log file close
sink()
sink(type = "message")
close(log_stream)