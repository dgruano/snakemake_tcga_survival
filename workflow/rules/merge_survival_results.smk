rule all_merge:
    input:
        "<results>/screening/survival/merged_per_signature.xlsx"

rule merge_survival_results:
    input:
        tsv_files = lambda wildcards: survival_output(wildcards)
    output:
        "<results>/screening/survival/merged_per_signature.xlsx"
    threads:
        config["resources"]["merge_survival_results"]["threads"]
    resources:
        mem = config["resources"]["merge_survival_results"]["mem"],
        time = config["resources"]["merge_survival_results"]["time"]
    conda:
        "../../envs/merge_smk.yml"
    log:
        "<logs>/merge_survival_results.log"
    params:
        cohorts = lambda wildcards: ",".join(TCGA_PROJECTS)
    shell:
        """
        Rscript scripts/merge_survival_results.R {params.cohorts} {output} > {log} 2>&1
        """
