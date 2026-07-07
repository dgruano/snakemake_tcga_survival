rule all_survival:
    input:
        survival_output,


rule compute_scores:
    input:
        signatures_file=os.path.abspath(config["signatures_file"]),
        deseq2_file="<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
    output:
        scores="<results>/screening/survival/{project}/patient_scores.tsv",
        categorical="<results>/screening/survival/{project}/patient_scores_categorical.tsv",
        clinical="<results>/screening/survival/{project}/clinical_data.tsv",
    log:
        "<logs>/compute_scores_{project}.log",
    benchmark:
        "<logs>/benchmarks/compute_scores_{project}.tsv",
    conda:
        "../../envs/tcga_smk.yml"
    threads: 4
    params:
        PERCENTILES=",".join(map(str, config["PERCENTILES"])),
    script:
        "../scripts/compute_scores.R"


rule survival_per_signature:
    input:
        scores="<results>/screening/survival/{project}/patient_scores.tsv",
        categorical="<results>/screening/survival/{project}/patient_scores_categorical.tsv",
        clinical="<results>/screening/survival/{project}/clinical_data.tsv",
        deseq2_file="<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
    output:
        "<results>/screening/survival/{project}/{signature}_pval.tsv",
    log:
        "<logs>/survival_{project}_{signature}.log",
    benchmark:
        "<logs>/benchmarks/survival_{project}_{signature}.tsv",
    conda:
        "../../envs/tcga_smk.yml"
    wildcard_constraints:
        signature="|".join(re.escape(s) for s in SIGNATURES),
    params:
        THRESHOLD=config["THRESHOLD"],
        DPI=config["DPI"],
        PERCENTILES=",".join(map(str, config["PERCENTILES"])),
    script:
        "../scripts/survival_per_signature.R"


rule gather_survival:
    input:
        lambda wc: expand(
            "<results>/screening/survival/{project}/{signature}_pval.tsv",
            project=wc.project,
            signature=SIGNATURES,
        ),
    output:
        "<results>/screening/survival/{project}/survival_pval_filtered.tsv",
    run:
        import pandas as pd

        frames = [
            pd.read_csv(f, sep="\t", index_col=0)
            for f in input
            if os.path.getsize(f) > 0
        ]
        result = pd.concat(frames) if frames else pd.DataFrame()
        result.to_csv(output[0], sep="\t")
