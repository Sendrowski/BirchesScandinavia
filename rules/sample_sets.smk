include: "header.smk"
include: "vcf_filtering.smk"

# target output files
rule run_all_sample_sets:
    input:
        expand("results/{sample_class}/sample_sets/{sample_set}.args",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"])

# create a file containing all birch sample names
rule create_sample_set_all:
    input:
        config['birch_samples'],
        config['outgroup_samples']
    output:
        "results/{sample_class}/sample_sets/all.args"
    params:
        sample_set="all",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create a file containing all birch sample names
rule create_sample_set_birch:
    input:
        config['birch_samples']
    output:
        "results/{sample_class}/sample_sets/birch.args"
    params:
        sample_set="birch",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create file containing the samples of the left initial PCA clusters
rule create_sample_set_left_cluster:
    input:
        "results/{sample_class}/snps/birch/biallelic/snps.bed"
    output:
        "results/{sample_class}/sample_sets/left_cluster.args"
    params:
        sample_set="left_cluster",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create file containing the samples of the right initial PCA clusters
rule create_sample_set_right_cluster:
    input:
        "results/{sample_class}/snps/birch/biallelic/snps.bed"
    output:
        "results/{sample_class}/sample_sets/right_cluster.args"
    params:
        sample_set="right_cluster",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set only containing putative B. pendula or B. pubescens species
rule create_sample_set_pendula_pubescens:
    input:
        birch=config['birch_samples']
    output:
        "results/{sample_class}/sample_sets/pendula_pubescens.args"
    params:
        sample_set="pendula_pubescens",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set with all putative
# B. pendula species based on PCA
rule create_sample_set_pendula:
    input:
        "results/{sample_class}/sample_sets/right_cluster.args"
    output:
        "results/{sample_class}/sample_sets/pendula.args"
    params:
        sample_set="pendula",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set including all southern
# B. pendula species based on PCA
rule create_sample_set_pendula_south:
    input:
        "results/{sample_class}/sample_sets/pendula.args"
    output:
        "results/{sample_class}/sample_sets/pendula_south.args"
    params:
        sample_set="pendula_south",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set including all northern
# B. pendula species based on PCA
rule create_sample_set_pendula_north:
    input:
        "results/{sample_class}/sample_sets/pendula.args"
    output:
        "results/{sample_class}/sample_sets/pendula_north.args"
    params:
        sample_set="pendula_north",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set including all southern
# B. pubescens species based on PCA
rule create_sample_set_pubescens_south:
    input:
        "results/{sample_class}/sample_sets/pubescens.args"
    output:
        "results/{sample_class}/sample_sets/pubescens_south.args"
    params:
        sample_set="pubescens_south",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set including all northern
# B. pubescens species based on PCA
rule create_sample_set_pubescens_north:
    input:
        "results/{sample_class}/sample_sets/pubescens.args"
    output:
        "results/{sample_class}/sample_sets/pubescens_north.args"
    params:
        sample_set="pubescens_north",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set with all putative
# B. pendula species based the list obtained from Luis
rule create_sample_set_pendula_luis:
    input:
        config['pendula_samples_luis']
    output:
        "results/{sample_class}/sample_sets/pendula_luis.args"
    conda:
        "../envs/base.yaml"
    shell:
        "cp {input} {output}"

# create the subsample set with all putative
# B. pubescens species based on PCA
checkpoint create_sample_set_pubescens:
    input:
        "results/{sample_class}/sample_sets/left_cluster.args"
    output:
        "results/{sample_class}/sample_sets/pubescens.args"
    params:
        sample_set="pubescens",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create files containing subpopulation labels used for dadi
rule create_subpopulation_files_1_pop:
    input:
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        n_pops=1
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_subpopulation_files.py"

# create files containing subpopulation labels used for dadi
rule create_subpopulation_files_2_pops:
    input:
        "results/{sample_class}/sample_sets/{sample_set}.args",
        lambda w: expand("results/{{sample_class}}/sample_sets/{sample_set}.args",sample_set=get_subpopulations(w.sample_set))
    output:
        "results/{sample_class}/sample_sets/subpopulations/{sample_set}/2_pops.txt"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        n_pops=2
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_subpopulation_files.py"
