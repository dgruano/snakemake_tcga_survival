rule all_survival:
    input:
        survival_output,


rule survival_screening:
    input:
        signatures_file=os.path.abspath(config["signatures_file"]),
        deseq2_file="<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
    output:
        "<results>/screening/survival/{project}/survival_pval_filtered.tsv",
    log:
        "<logs>/survival_screening_{project}.log",
    benchmark:
        "<logs>/benchmarks/survival_screening_{project}.tsv",
    conda:
        "../../envs/tcga_smk.yml"
    params:
        THRESHOLD=config["THRESHOLD"],
        DPI=config["DPI"],
        PERCENTILES=",".join(map(str, config["PERCENTILES"])),
    script:
        "../scripts/survival_screening.R"
