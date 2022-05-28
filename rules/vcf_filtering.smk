include: "header.smk"
include: "sample_sets.smk"
include: "vcf_to_plink.smk"

ruleorder:
    generate_vcf_sample_set_biallelic >
    generate_vcf_sample_set_no_missing >
    generate_vcf_sample_set >
    generate_vcf_sample_set_nonsynonymous >
    generate_vcf_sample_set_synonymous >
    generate_vcf_sample_set_no_low_freqs >
    generate_vcf_sample_set_no_very_low_freqs >
    create_tbi_index_vcf

# target output files
rule run_all_vcf_filtering:
    input:
        expand("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],flag=flags),
        expand("results/{sample_class}/graphs/png/synonymy/{sample_set}/synonymy.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"]),
        expand("results/{sample_class}/graphs/png/degeneracy/{sample_set}/degeneracy.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"]),
        expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.1D.{n_proj}.{unfolded}.{scale}.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],flag=flags,n_proj=[20],unfolded=['unfolded'],scale=['linear']),
        expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.2D.{n_proj}.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],flag=flags,n_proj=[20],scale=['linear']),
        expand("results/{sample_class}/graphs/png/sfs/{sample_set}/synonymous_nonsynonymous/sfs.1D.{n_proj}.{unfolded}.{scale}.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],n_proj=[20],unfolded=['unfolded'],scale=['linear']),
        expand("results/{sample_class}/graphs/png/sfs/{sample_set}/0fold_4fold/sfs.1D.{n_proj}.{unfolded}.{scale}.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],n_proj=[20],unfolded=['unfolded'],scale=['linear']),
        expand("results/{sample_class}/graphs/png/sfs/{sample_set}/0fold_4fold/sfs.1D.{n_proj}.{unfolded}.{scale}_tight_layout.png",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],n_proj=[20],unfolded=['unfolded'],scale=['linear'])


rule resource_raw_vcf:
    output:
        protected("results/{sample_class}/snps/raw/snps.vcf.gz")

# create vcf files containing the samples of the specified sample set
rule generate_vcf_sample_set:
    input:
        vcf="results/{sample_class}/snps/raw/snps.vcf.gz",
        tbi="results/{sample_class}/snps/raw/snps.vcf.gz.tbi",
        samples="results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz.tbi"
    params:
        max_nocall_fraction=0.5
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/generate_vcf_sample_set.py"

# create vcf files containing the samples of the specified sample set
# allow no missing values
rule generate_vcf_sample_set_no_missing:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz.tbi",
        samples="results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/snps/{sample_set}/no_missing/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/no_missing/snps.vcf.gz.tbi"
    params:
        max_nocall_fraction=0
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/generate_vcf_sample_set.py"

# create vcf files containing the samples of the specified sample set
# restrict to biallelic alleles
rule generate_vcf_sample_set_biallelic:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz.tbi",
        samples="results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz.tbi"
    params:
        max_nocall_fraction=0.5,
        biallelic=True
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/generate_vcf_sample_set.py"

# remove all sites that have a minor or major allele frequency lower than equal 0.05
rule generate_vcf_sample_set_no_low_freqs:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz.tbi",
        samples="results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/snps/{sample_set}/no_low_freqs/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/no_low_freqs/snps.vcf.gz.tbi"
    params:
        filter='0.05 <= AF && AF <= 0.95'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# remove all sites that have a minor or major allele frequency lower than equal 0.01
rule generate_vcf_sample_set_no_very_low_freqs:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz.tbi",
        samples="results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/snps/{sample_set}/no_very_low_freqs/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/no_very_low_freqs/snps.vcf.gz.tbi"
    params:
        filter='0.01 <= AF && AF <= 0.99'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# create a vcf file containing only the synonymous sites
rule generate_vcf_sample_set_synonymous:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz.tbi"
    params:
        filter='Synonymy != "." && Synonymy == 1'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# create a vcf file containing only the non-synonymous sites
rule generate_vcf_sample_set_nonsynonymous:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/biallelic/snps.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz.tbi"
    params:
        filter='Synonymy != "." && Synonymy == 0'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# restrict sample set to 4-fold degenerate sites
rule filter_nfold_sites:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/{n_fold_deg}fold/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/{n_fold_deg}fold/snps.vcf.gz.tbi"
    params:
        filter='Degeneracy != "." && Degeneracy == {n_fold_deg}'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# restrict the set of synonymous sites to 4-fold degenerate sites
rule filter_4fold_degenerate_sites:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/synonymous_4fold/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/synonymous_4fold/snps.vcf.gz.tbi"
    params:
        filter='Degeneracy != "." && Degeneracy == 4'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# restrict the set of non-synonymous sites to 0-fold degenerate sites
rule filter_0fold_degenerate_sites:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/nonsynonymous_0fold/snps.vcf.gz",
        "results/{sample_class}/snps/{sample_set}/nonsynonymous_0fold/snps.vcf.gz.tbi"
    params:
        filter='Degeneracy != "." && Degeneracy == 0'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"

# plot a 1D SFS
rule plot_sfs_1D:
    input:
        vcf=["results/{sample_cass}/snps/{sample_set}/{flag}/snps.vcf.gz"],
        pops="results/{sample_cass}/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    output:
        "results/{sample_cass}/graphs/svg/sfs/{sample_set}/{flag}/sfs.1D.{n_proj}.{unfolded}.{scale}.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        unfolded=lambda w: w.unfolded == "unfolded",
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_1D.py"

# plot a 2D SFS
rule plot_sfs_2D:
    input:
        vcf="results/{sample_cass}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops="results/{sample_cass}/sample_sets/subpopulations/{sample_set}/2_pops.txt"
    output:
        "results/{sample_cass}/graphs/svg/sfs/{sample_set}/{flag}/sfs.2D.{n_proj}.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        sample_set=lambda w: w.sample_set
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_2D.py"

# plot the SFS
rule plot_sfs_synonymous_nonsynonymous:
    input:
        vcf=["results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
             "results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz"],
        pops="results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    output:
        "results/{sample_class}/graphs/svg/sfs/{sample_set}/synonymous_nonsynonymous/sfs.1D.{n_proj}.{unfolded}.{scale}.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        unfolded=lambda w: w.unfolded == "unfolded",
        log_scale=lambda w: w.scale == "log",
        labels=['synonymous', 'nonsynonymous']
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_1D.py"

# plot the SFS using only 0 and 4-fold degenerate sites
rule plot_sfs_0fold_4fold:
    input:
        vcf=["results/{sample_class}/snps/{sample_set}/synonymous_4fold/snps.vcf.gz",
             "results/{sample_class}/snps/{sample_set}/nonsynonymous_0fold/snps.vcf.gz"],
        pops="results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    output:
        "results/{sample_class}/graphs/svg/sfs/{sample_set}/0fold_4fold/sfs.1D.{n_proj}.{unfolded}.{scale}.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        unfolded=lambda w: w.unfolded == "unfolded",
        log_scale=lambda w: w.scale == "log",
        labels=['4fold', '0fold']
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_1D.py"

# plot the SFS using only 0 and 4-fold degenerate sites
rule plot_sfs_0fold_4fold_tight_layout:
    input:
        vcf=["results/{sample_class}/snps/{sample_set}/synonymous_4fold/snps.vcf.gz",
             "results/{sample_class}/snps/{sample_set}/nonsynonymous_0fold/snps.vcf.gz"],
        pops="results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    output:
        "results/{sample_class}/graphs/svg/sfs/{sample_set}/0fold_4fold/sfs.1D.{n_proj}.{unfolded}.{scale}_tight_layout.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        unfolded=lambda w: w.unfolded == "unfolded",
        log_scale=lambda w: w.scale == "log",
        labels=['4fold', '0fold'],
        tight_layout=True
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_1D.py"

# plot pie chart of degeneracy classes
rule plot_degeneracy:
    input:
        sites_all="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        synonymous="results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
        nonsynonymous="results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz"
    output:
        "results/{sample_class}/graphs/svg/degeneracy/{sample_set}/degeneracy.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_degeneracy.py"

# plot pie chart of synonymy classes
rule plot_synonymy:
    input:
        synonymous="results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
        nonsynonymous="results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz"
    output:
        "results/{sample_class}/graphs/svg/synonymy/{sample_set}/synonymy.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_synonymy.py"

# add the R2 info field to the info field using BCFtools denoting the largest R2 values
rule annotate_ld_R2:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz"
    output:
        temp("results/{sample_class}/snps/{sample_set}/{flag}/snps.ld_annot.vcf.gz")
    params:
        window_size_kb=50
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/annotate_ld.py"

# filter sites that are likely in LD
rule filter_ld_R2:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.ld_annot.vcf.gz",
        tbi="results/{sample_class}/snps/{sample_set}/{flag}/snps.ld_annot.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.ld.vcf.gz"
    params:
        filter='R2 <= 0.2'
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/filter_by_info.py"
