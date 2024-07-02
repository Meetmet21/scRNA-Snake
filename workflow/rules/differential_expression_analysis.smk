rule tximport_counts:
    """
    Run R script to generate counts from kallisto h5 files and 
    analysis count matrix for gene expression.
    """
    input:
        all_h5=collect("results/{sample}/kallisto/abundance.h5",
            sample=samples["sample"])
    output:
        counts="results/tximport/txi.RData"
    params:
        samples_tsv=config["samples"]
    threads: 8
    conda:
        "../envs/r.yaml"
    log:
        "workflow/logs/tximport_counts.log"
    script:
        "../scripts/tximport.R"

rule deseq2_differential_expression:
    """
    Differential expression analysis, retrieving topGenes, normalized DE and raw DE.
    """
    input:
        counts=rules.tximport_counts.output.counts
    output:
        raw="results/DESeq2/DESeq2.RData",
        shrink="results/DESeq2/DESeq2_shrink.RData",
        TopGenes="results/DESeq2/TopGenes.tsv",
    conda:
        "../envs/r.yaml"
    params:
        report_dir="workflow/report/DESeq2"
    threads: 8
    log:
        "workflow/logs/deseq2_differential_expression.log"
    script:
        "../scripts/DESeq2.R"

rule clusterProfiler_GOenrich_KGBB:
    """
    Running Go Enrichment Analysis and Pathway Enrichment Analysis using
    clusterProfiler::enrichGO. Retrieve reports from analysis and processed data.
    """
    input:
        TopGenes=rules.deseq2_differential_expression.output.TopGenes
    output:
        TopGO="results/clusterProfiler/TopGO.tsv"
    params:
        report_dir="workflow/report/clusterProfiler"
    conda:
        "../envs/r.yaml"
    log:
        "workflow/logs/clusterProfiler_GOenrich_KGBB"
    script:
        "../scripts/clusterProfiler.R"
