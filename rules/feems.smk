import numpy as np

include: "header.smk"
include: "vcf_to_plink.smk"

# lambda values to use for cross-validation procedure
lambdas_cv = np.geomspace(1e-8,1e2,11)

# target output files
rule run_feems:
    input:
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/cv_error_warm_start.png",sample_set=['pendula', 'pubescens'],flag=['no_low_freqs', 'no_very_low_freqs']),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/cv_error_warm_start_buffer_1.png",sample_set=['pendula', 'pubescens'],flag=['no_low_freqs', 'no_very_low_freqs']),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_cv.png",sample_set=['pendula', 'pubescens'],flag=['no_low_freqs', 'no_very_low_freqs']),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_cv_buffer_1.png",sample_set=['pendula', 'pubescens'],flag=['no_low_freqs', 'no_very_low_freqs']),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/cv_error_cold_start.png",sample_set=['pendula', 'pubescens'],flag=['no_low_freqs', 'no_very_low_freqs']),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_lambda_{n}.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs'],n=[10, 1, 0.1, 0.01]),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_lambda_{n}_buffer_1.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs'],n=[10, 1, 0.1, 0.01]),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_lambda_{n}_buffer_0.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs'],n=[10, 1, 0.1, 0.01]),
        expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_lambda_{n}.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs'],n=lambdas_cv),
        expand("results/default/graphs/png/feems/{sample_set}/allocation.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens']),
        expand("results/default/graphs/png/feems/{sample_set}/locations.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'])

# resource VCF files
rule resource_vcf_feems:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# resource for list of sample names
rule resource_sample_sets_feems:
    output:
        protected("results/{sample_class}/sample_sets/{sample_set}.args")

# calculate the migration surfaces and determine the best lambda using cross-validation
rule plot_cv_error_feems:
    input:
        expand("results/{{sample_class}}/feems/{{sample_set}}/{{flag}}/cv_error_{lamb}.txt",lamb=lambdas_cv)
    output:
        "results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/cv_error_cold_start.svg"
    params:
        lambdas=lambdas_cv
    conda:
        "../envs/feems.yaml"
    script:
        "../scripts/plot_cv_error_feems.py"

# Calculate the migration surfaces for lambda
# providing the lowest CV error.
# Not parallelizing the calculations, we can use warm start
# which perform considerably better according to the paper.
rule calculate_migration_surfaces_cv:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        surfaces="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/migration_surfaces_cv.svg",
        cv="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/cv_error_warm_start.svg"
    params:
        lambdas=lambdas_cv,
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}"
    conda:
        "../envs/feems.yaml"
    script:
        "../scripts/calculate_migration_surfaces_cv.py"

rule calculate_migration_surfaces_cv_buffer:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        surfaces="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/migration_surfaces_cv_buffer_{buffer}.svg",
        cv="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/cv_error_warm_start_buffer_{buffer}.svg"
    params:
        lambdas=lambdas_cv,
        buffer= lambda w: float(w.buffer),
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}"
    conda:
        "../envs/feems.yaml"
    script:
        "../scripts/calculate_migration_surfaces_cv.py"

# calculate the migration surfaces for the given value of lambda
rule calculate_migration_surfaces_lambda:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        surfaces="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/migration_surfaces_lambda_{lamb}.svg",
        cv_error="results/{sample_class}/feems/{sample_set}/{flag}/cv_error_{lamb}.txt"
    params:
        lamb=lambda w: float(w.lamb),
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}"
    conda:
        "../envs/feems.yaml"
    script:
        "../scripts/calculate_migration_surfaces_lambda.py"

# calculate the migration surfaces for the given value of lambda
rule calculate_migration_surfaces_lambda_buffer:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        surfaces="results/{sample_class}/graphs/svg/feems/{sample_set}/{flag}/migration_surfaces_lambda_{lamb}_buffer_{buffer}.svg",
        cv_error="results/{sample_class}/feems/{sample_set}/{flag}/cv_error_{lamb}_buffer_{buffer}.txt"
    params:
        lamb=lambda w: float(w.lamb),
        buffer=lambda w: float(w.buffer),
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="{flag}"
    conda:
        "../envs/feems.yaml"
    script:
        "../scripts/calculate_migration_surfaces_lambda.py"

# plot the sample allocation for the migration surfaces
rule plot_migration_surfaces_sample_allocation:
    input:
        "results/{sample_class}/snps/{sample_set}/4fold/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/feems/{sample_set}/allocation.svg"
    conda:
        "../envs/feems.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}",
        flag="4fold"
    script:
        "../scripts/plot_sample_allocation_feems.py"

# plot the sample location for the migration surfaces
rule plot_migration_surfaces_sample_locations:
    input:
        "results/{sample_class}/snps/{sample_set}/4fold/snps.bed",
        "results/{sample_class}/sample_sets/{sample_set}.args"
    output:
        "results/{sample_class}/graphs/svg/feems/{sample_set}/locations.svg"
    conda:
        "../envs/feems.yaml"
    params:
        sample_set="{sample_set}",
        sample_class="{sample_class}"
    script:
        "../scripts/plot_sample_locations.py"
