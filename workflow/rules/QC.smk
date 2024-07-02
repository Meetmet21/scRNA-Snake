rule fastqc:
    """
    Quality control report for raw sequencing data
    """
    input:
        fastq="resources/data/raw_fastq/{sample}.fastq.gz"
    output:
        html_report="workflow/report/fastqc/{sample}_fastqc.html"
    log:
        "workflow/logs/{sample}/fastqc.log"
    params:
        outdir=lambda wildcards, output: dirname(output.html_report),
    threads: 4
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --outdir {params.outdir} {input.fastq} 2>&1 | "
        "tee {log}"

rule multiqc:
    """
    Aggregate multi sample fastqc and alignment information into one html
    """
    input:
        # wait for files
        fastqc=collect("workflow/report/fastqc/{sample}_fastqc.html", sample=samples["sample"]),
        align=collect("workflow/logs/{sample}/kallisto_pseudo_alignement.log", sample=samples["sample"])
    output:
        multiqc="workflow/report/multiqc/multiqc_report.html"
    params:
        qc_dir=lambda wildcards, input: dirname(input.fastqc[0]),
        align_dir=lambda wildcards, input: [dirname(path) for path in input.align],
        outdir=lambda wildcards, output: dirname(output.multiqc)
    conda:
        "../envs/multiqc.yaml"
    log:
        "workflow/logs/multiqc/multiqc.log"
    shell:
        "multiqc --outdir {params.outdir} {params.qc_dir} {params.align_dir}"