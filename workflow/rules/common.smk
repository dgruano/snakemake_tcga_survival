configfile: "config.yaml"


def get_tcga_projects():
    """Load TCGA projects from config and expand SKCM variants."""
    with open(config["TCGA_cohorts"]) as f:
        projects = [line.strip() for line in f if line.strip()]
    # Automatically add SKCM_prim and SKCM_met
    if "TCGA-SKCM" in projects:
        projects.extend(["TCGA-SKCM_prim", "TCGA-SKCM_met"])
    return projects


TCGA_PROJECTS = get_tcga_projects()


def tcga_download_outputs(wildcards):
    return expand("<results>/rds/{cohort}_STAR_Counts.rds", cohort=TCGA_PROJECTS)


def deseq2_output(wildcards):
    return expand(
        "<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
        project=TCGA_PROJECTS,
    )


def survival_output(wildcards):
    return expand(
        "<results>/screening/survival/{project}/survival_pval_filtered.tsv",
        project=TCGA_PROJECTS,
    )
