
rule kallisto_index_ref:
    input:
        ref=config["ref"]["fasta_path"]
    output:
        index="resources/data/kallisto/"+re.sub(
            pattern="cdna.all.fa.gz",
            repl="index",
            string=basename(config["ref"]["fasta_path"]))
    threads: 8
    log:
        "workflow/logs/kallisto_index_ref.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index --index {output.index} --threads {threads} {input.ref} 2>&1 | "
        "tee {log}"

rule kallisto_pseudo_alignement:
    """
    Pseudo-allign and quantify reads from single read data.
    """
    input:
        fasta="resources/data/raw_fastq/{sample}.fastq.gz",
        index=rules.kallisto_index_ref.output.index
    output:
        tsv="results/{sample}/kallisto/abundance.h5"
    params:
        outdir=lambda wildcards, output: dirname(output.tsv)
    log:
        "workflow/logs/{sample}/kallisto_pseudo_alignement.log"
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--index {input.index} "
        "--output-dir {params.outdir} "
        "--single "
        "--fragment-length 100 "
        "--sd 1 "
        "--threads {threads} "
        "--verbose 2>&1 "
        "{input.fasta} | "
        "tee {log} "
