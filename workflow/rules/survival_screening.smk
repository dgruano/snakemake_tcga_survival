rule all_survival:
    input:
        survival_output

rule survival_screening:
    input:
        signatures_file = os.path.abspath(config["signatures_file"]),
        deseq2_file = "<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds"
    params:
        THRESHOLD = config["THRESHOLD"],
        DPI = config["DPI"],
        PERCENTILES = ",".join(map(str, config["PERCENTILES"]))
    output:
        "<results>/screening/survival/{project}/survival_pval_filtered.tsv"
    threads:
        config["resources"]["survival_screening"]["threads"]
    resources:
        mem = config["resources"]["survival_screening"]["mem"],
        time = config["resources"]["survival_screening"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/survival_screening_{project}.log"
    shell:
        """
        Rscript scripts/survival_screening.R {params.DPI} {params.THRESHOLD} {input.signatures_file} {wildcards.project} {output} {params.PERCENTILES} > {log} 2>&1
        """
