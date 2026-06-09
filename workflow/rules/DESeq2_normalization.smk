rule all_deseq2:
    input:
        deseq2_output,


rule DESeq2_normalization:
    input:
        rds_file="<results>/rds/{project}_STAR_Counts.rds",
        biomart_file="<results>/biomart/biomart_protein_coding_genes.csv",
    output:
        "<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
    log:
        "<logs>/DESeq2_normalization_{project}.log",
    conda:
        "../../envs/tcga_smk.yml"
    script:
        "../scripts/DESeq2_normalization.R"
