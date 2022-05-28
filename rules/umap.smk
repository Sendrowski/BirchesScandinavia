include: "header.smk"
include: "vcf_to_plink.smk"

# target output files
rule run_umap:
    input:
        expand("results/{sample_class}/graphs/png/umap/{sample_set}/biallelic/umap.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/umap/{sample_set}/{flag}/umap_marked.png",sample_set=['pendula_pubescens', 'birch'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/umap/{sample_set}/{flag}/umap_tight_layout.png",sample_set=['pendula_pubescens', 'birch', 'pendula', 'pubescens'],flag=['biallelic', 'synonymous'],sample_class=['default'])

# resource VCF files
rule resource_vcf_umap:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# resource for list of sample names
rule resource_sample_sets_umap:
    output:
        protected("results/{sample_class}/sample_sets/{sample_set}.args")

# generate a UMAP graph
rule plot_umap_tight_layout:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/umap/{sample_set}/{flag}/umap_tight_layout.svg"
    conda:
        "../envs/umap.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        tight_layout=True
    script:
        "../scripts/plot_umap.py"

# generate a UMAP graph
rule plot_umap:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/umap/{sample_set}/{flag}/umap.svg"
    conda:
        "../envs/umap.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        min_dist=0
    script:
        "../scripts/plot_umap.py"

# generate a UMAP graph with labeled subpopulations
rule plot_umap_marked:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_marked_subpopulations(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/umap/{sample_set}/{flag}/umap_marked.svg"
    conda:
        "../envs/umap.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        subsample_sets=lambda w: get_marked_subpopulations(w.sample_set),
        tight_layout=True
    script:
        "../scripts/plot_umap.py"
