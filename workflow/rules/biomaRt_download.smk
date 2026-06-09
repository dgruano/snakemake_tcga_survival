rule biomaRt_download:
    output:
        "<results>/biomart/biomart_protein_coding_genes.csv",
    log:
        "<logs>/biomaRt_download.log",
    retries: 3
    conda:
        "../../envs/tcga_smk.yml"
    threads: config["resources"]["biomaRt_download"]["threads"]
    resources:
        mem=config["resources"]["biomaRt_download"]["mem"],
        time=config["resources"]["biomaRt_download"]["time"],
    script:
        "../scripts/biomaRt_download.R"
