include: "header.smk"
include: "vcf_to_plink.smk"


# get the ADMIXTURE subpopulations for the given sample set
def get_subpopulations_admixture(sample_set, n=2):
    return get_subpopulations(sample_set + '_admixture',n)


# sort from south to north in ascending order
# the location is inferred from the sample set name
def sort_from_south_to_north(sample_sets):
    return sorted(sample_sets,key=lambda x: 2 * int('pubescens' in x) + int('north' in x))


# R2 thresholds depending on sample set
# we need to make sure enough sites remain
R2_thresholds = {
    # we need almost all sites to discover
    # population structure in pendula
    'pendula': 0.99,
    'default': 0.2
}


# the R2 threshold per sample set
def get_R2(sample_set):
    if sample_set in R2_thresholds:
        return R2_thresholds[sample_set]

    return R2_thresholds['default']


rule run_admixture_all:
    input:
        expand("results/default/admixture/pendula_pubescens/{flag}/snps.admixture.{K}.P",K=[1, 2, 3, 4, 5, 6, 7],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/admixture/pendula/{flag}/snps.admixture.{K}.P",K=[1, 2, 3, 4],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/admixture/pubescens/{flag}/snps.admixture.{K}.P",K=[1, 2, 3, 4],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/graphs/png/admixture/{sample_set}/{flag}/cv_error.4.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/graphs/png/admixture/{sample_set}/{flag}/cv_error.7.png",sample_set=['pendula_pubescens'],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/graphs/png/admixture/{sample_set}/{flag}/synonymous/pca.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/graphs/png/admixture/{sample_set}/{flag}/locations.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'synonymous', 'nonsynonymous']),
        expand("results/default/graphs/png/admixture/{sample_set}/biallelic/barplot.{K}.2.png",sample_set=['pendula', 'pubescens'],K=[2, 3, 4]),
        expand("results/default/graphs/png/admixture/{sample_set}/biallelic/barplot.{K}.4.png",sample_set=['pendula_pubescens'],K=[2, 3, 4, 5, 6, 7])

# unfiltered vcf resource file
rule resource_vcf_admixture:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# recode the ids of the variants
rule recode_variant_ids:
    input:
        bed="results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        bim="results/{sample_class}/snps/{sample_set}/{flag}/snps.bim",
        fam="results/{sample_class}/snps/{sample_set}/{flag}/snps.fam"
    output:
        bed=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bed"),
        bim=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bim"),
        fam=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.fam"),
        log=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.log")
    conda:
        "../envs/plink2.yaml"
    script:
        "../scripts/recode_variant_ids_plink.py"

# create a pair of files listing sites that are below
# or above the biggest R2 statistics using PLINK
rule create_pairwise_ld_lists:
    input:
        bed="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bed",
        bim="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bim",
        fam="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.fam"
    output:
        temp("results/{sample_class}/snps/{sample_set}/{flag}/plink.prune.in"),
        temp("results/{sample_class}/snps/{sample_set}/{flag}/plink.prune.out")
    params:
        window_size=500,
        step_size=50,
        R2=lambda w: get_R2(w.sample_set)
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/create_pairwise_ld_lists.py"

# extract sites with the ids of given list
rule extract_sites_plink:
    input:
        bed="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bed",
        bim="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.bim",
        fam="results/{sample_class}/snps/{sample_set}/{flag}/snps.recoded_ids.fam",
        ids="results/{sample_class}/snps/{sample_set}/{flag}/plink.prune.in"
    output:
        bed=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.bed"),
        bim=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.bim"),
        fam=temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.fam")
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/extract_sites_plink.py"

# run ADMIXTURE
rule run_admixture:
    input:
        bed="results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.bed",
        bim="results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.bim",
        fam="results/{sample_class}/snps/{sample_set}/{flag}/snps.admixture.fam"
    output:
        P="results/{sample_class}/admixture/{sample_set}/{flag}/snps.admixture.{K}.P",
        Q="results/{sample_class}/admixture/{sample_set}/{flag}/snps.admixture.{K}.Q",
        log="results/{sample_class}/admixture/{sample_set}/{flag}/snps.admixture.{K}.log"
    params:
        K="{K}"
    conda:
        "../envs/admixture.yaml"
    script:
        "../scripts/run_admixture.py"

# plot the CV error for different values of K
rule plot_cv_error:
    input:
        lambda w: [f"results/{{sample_class}}/admixture/{{sample_set}}/{{flag}}/snps.admixture.{n}.log" for n in range(1,int(w.K_max) + 1)]
    output:
        "results/{sample_class}/graphs/svg/admixture/{sample_set}/{flag}/cv_error.{K_max}.svg"
    params:
        K=lambda w: range(1,int(w.K_max) + 1)
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_cv_error_admixture.py"

# create the two cluster sample sets based on ADMIXTURE
rule create_sample_set_ADMIXTURE:
    input:
        lambda w: f"results/default/admixture/{get_superpopulation(w.sample_set)}/biallelic/snps.admixture.2.Q"
    output:
        "results/{sample_class}/sample_sets/{sample_set}_admixture.args"
    params:
        sample_set="{sample_set}_admixture",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# generate PCA plot with labeled subpopulations according
# to the two cluster subpopulation division from ADMIXTURE
rule plot_pca_marked_ADMIXTURE_subsets:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_subpopulations_admixture(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/admixture/{sample_set}/{flag}/{flag_pca}/pca.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag_pca}",
        add_names=False,
        subsample_sets=lambda w: get_subpopulations_admixture(w.sample_set),
        demarcation_type='color',
        tight_layout=True,
        scaling=2 / 3,
        marker_size=10
    script:
        "../scripts/plot_pca.py"

# plot sample locations with labeled subpopulations according
# to the two cluster subpopulation division from ADMIXTURE
rule plot_sample_locations_ADMIXTURE_subsets:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_subpopulations_admixture(w.sample_set)],
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/admixture/{sample_set}/{flag}/locations.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        subsample_sets=lambda w: get_subpopulations_admixture(w.sample_set)
    script:
        "../scripts/plot_clusters_location.py"

# plot sample locations with labeled subpopulations according
# to the two cluster subpopulation division from ADMIXTURE
rule create_barplot_ADMIXTURE:
    input:
        "results/{sample_class}/sample_sets/{sample_set}.args",
        lambda w: ["results/{{sample_class}}/sample_sets/{}.args".format(s) for s in get_subpopulations_admixture(w.sample_set,w.N)],
        Q="results/{sample_class}/admixture/{sample_set}/{flag}/snps.admixture.{K}.Q"
    output:
        "results/{sample_class}/graphs/svg/admixture/{sample_set}/{flag}/barplot.{K}.{N}.svg"
    conda:
        "../envs/base.yaml"
    params:
        sample_set="{sample_set}",
        flag="{flag}",
        sample_class="{sample_class}",
        subsample_sets=lambda w: sort_from_south_to_north(get_subpopulations_admixture(w.sample_set,w.N))
    script:
        "../scripts/create_barplot_admixture.py"
