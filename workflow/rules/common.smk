def get_fastqc_reports(wildcards):
    """
    Get FASTQC output path for each sample.
    :param wildcards:
    :return: list of paths
    """
    return ["workflow/report/fastqc/"+sample+"_fastqc.html"
            for sample in samples["sample"]]
