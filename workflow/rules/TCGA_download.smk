rule all_tcga:
    input:
        tcga_download_outputs,


rule TCGA_download:
    output:
        "<results>/rds/{cohort}_STAR_Counts.rds",
    log:
        "<logs>/TCGA_download_{cohort}.log",
    retries: 3
    conda:
        "../../envs/tcga_smk.yml"
    threads: config["resources"]["TCGA_download"]["threads"]
    resources:
        mem=config["resources"]["TCGA_download"]["mem"],
        time=config["resources"]["TCGA_download"]["time"],
    shell:
        """
        Rscript scripts/TCGA_download.R {wildcards.cohort} {output} > {log} 2>&1
        """


rule rule_separate_cohorts:
    input:
        "<results>/rds/TCGA-SKCM_STAR_Counts.rds",
    output:
        outfile_prim="<results>/rds/TCGA-SKCM_prim_STAR_Counts.rds",
        outfile_met="<results>/rds/TCGA-SKCM_met_STAR_Counts.rds",
    log:
        "<logs>/TCGA_split_SKCM.log",
    conda:
        "../../envs/tcga_smk.yml"
    threads: config["resources"]["TCGA_download"]["threads"]
    resources:
        mem=config["resources"]["TCGA_download"]["mem"],
        time=config["resources"]["TCGA_download"]["time"],
    shell:
        """
        Rscript scripts/TCGA_split_SKCM.R {input} {output.outfile_prim} {output.outfile_met} > {log} 2>&1
        """
