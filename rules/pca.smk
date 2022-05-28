include: "header.smk"
include: "vcf_to_plink.smk"

# target output files
rule run_pca:
    input:
        expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_tight_layout.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'synonymous', 'nonsynonymous'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_labeled.png",sample_set=['left_cluster', 'right_cluster', 'pendula_luis', 'birch', 'pendula_pubescens'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_marked.png",sample_set=['pendula_pubescens', 'birch'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        expand("results/{sample_class}/graphs/png/pca/{sample_set}/clustered_locations.2.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens'],sample_class=['default'])

# resource VCF files
rule resource_vcf_pca:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# resource for list of sample names
rule resource_sample_sets_pca:
    output:
        protected("results/{sample_class}/sample_sets/{sample_set}.args")

# generate PCA plot
rule plot_pca:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/pca/{sample_set}/{flag}/pca.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        add_names=False,
        subsample_sets=[]
    script:
        "../scripts/plot_pca.py"

# generate PCA plot
rule plot_pca_tight_layout:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/pca/{sample_set}/{flag}/pca_tight_layout.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        add_names=False,
        subsample_sets=[],
        tight_layout=True
    script:
        "../scripts/plot_pca.py"

# generate PCA plot with labeled subpopulations
rule plot_pca_marked:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_marked_subpopulations(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/pca/{sample_set}/{flag}/pca_marked.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        add_names=False,
        subsample_sets=lambda w: get_marked_subpopulations(w.sample_set),
        tight_layout=True
    script:
        "../scripts/plot_pca.py"

# generate labeled PCA plot
rule plot_pca_labeled:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_marked_subpopulations(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/pca/{sample_set}/{flag}/pca_labeled.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}",
        add_names=True,
        subsample_sets=[]
    script:
        "../scripts/plot_pca.py"

# plot sample locations with labeled subpopulations
rule plot_sample_locations_subsets:
    input:
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_subpopulations(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/pca/{sample_set}/clustered_locations.2.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        subsample_sets=lambda w: get_subpopulations(w.sample_set)
    script:
        "../scripts/plot_clusters_location.py"