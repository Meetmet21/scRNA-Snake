rule Ensembl_get_ref_transcriptoome:
    """
    Download the reference transcriptome GRCh38 from Ensembl FTP.
    """
    output:
        fasta=config["ref"]["fasta_path"]
    params:
        fasta_url=config["ref"]["fasta_url"]
    log:
        "workflow/logs/Ensembl_ref_transcriptoome.log"
    shell:
        "wget --quiet --timestamping -O {output.fasta} {params.fasta_url} 2>&1 | "
        "tee -a {log}"

rule Ensembl_get_ref_gene_annotation:
    """
    Download the reference gene annotation for GRCh38 from Ensembl FTP.
    """
    output:
        gtf=config["ref"]["gtf_path"]
    params:
        gtf_url=config["ref"]["gtf_url"],
    log:
        "workflow/logs/Ensembl_ref_gene_annotation.log"
    shell:
        "wget --quiet --timestamping -O {output.gtf} {params.gtf_url} 2>&1 | "
        "tee -a {log}"

rule ENA_get_sample_scRNAseq:
    """
    Download samples fastq.gz files from ENA. 
    """
    output:
        fastq="resources/data/raw_fastq/{sample}.fastq.gz"
    log:
        "workflow/logs/{sample}/ENA_get_sample_scRNAseq.log"
    params:
        fastq_url=lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, "fastq_url"].values[0]
    shell:
        "wget --quiet --timestamping -O {output.fastq} {params.fastq_url} 2>&1 | "
        "tee -a {log}"

