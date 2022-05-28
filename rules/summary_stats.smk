include: "header.smk"
include: "vcf_to_plink.smk"

# target output files
rule run_all_summary_stats:
    input:
        expand("results/default/graphs/png/pi/{sample_set}/{flag}/pi.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic', '0fold', '4fold'],scale=['log', 'linear']),
        expand("results/default/graphs/png/tajimas_d/{sample_set}/{flag}/tajimas_d.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_admixture', 'pubescens_admixture', 'pendula_pubescens_admixture'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.windowed.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        expand("results/default/graphs/png/missingness/{sample_set}/{flag}/missingness.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear'],type=['l', 'i']),
        expand("results/{sample_class}/graphs/png/hwe/{sample_set}/{flag}/{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['p', 'midp']),
        expand("results/{sample_class}/graphs/png/heterozygosity/{sample_set}/{flag}/het.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['expected', 'observed']),
        expand("results/{sample_class}/graphs/png/heterozygosity/{sample_set}/{flag}/het.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'])

# resource VCF files
rule resource_vcf_summary_stats:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# resource for list of sample names
rule resource_sample_sets_summary_stats:
    output:
        protected("results/{sample_class}/sample_sets/{sample_set}.args")

# calculate the Fst using vcftools
rule calculate_Fst:
    input:
        vcf="results/default/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: expand("results/default/sample_sets/{sample_set}.args",sample_set=get_subpopulations(w.sample_set))
    output:
        "results/default/stats/{sample_set}/{flag}/fst.weir.fst"
    params:
        prefix="results/default/stats/{sample_set}/{flag}/fst"
    conda:
        "../envs/vcftools.yaml"
    script:
        "../scripts/calculate_fst.py"

# calculate the windowed Fst using vcftools
rule calculate_Fst_windowed:
    input:
        vcf="results/default/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: expand("results/default/sample_sets/{sample_set}.args",sample_set=get_subpopulations(w.sample_set))
    output:
        "results/default/stats/{sample_set}/{flag}/fst.windowed.weir.fst"
    params:
        prefix="results/default/stats/{sample_set}/{flag}/fst"
    conda:
        "../envs/vcftools.yaml"
    script:
        "../scripts/calculate_fst_windowed.py"

# calculate Pi using vcftools
rule calculate_pi:
    input:
        vcf="results/default/snps/{sample_set}/{flag}/snps.vcf.gz"
    output:
        "results/default/stats/{sample_set}/{flag}/pi.sites.pi"
    params:
        prefix="results/default/stats/{sample_set}/{flag}/pi"
    conda:
        "../envs/vcftools.yaml"
    script:
        "../scripts/calculate_pi.py"

# calculate Tajima's D using vcftools
rule calculate_Tajimas_D:
    input:
        vcf="results/default/snps/{sample_set}/{flag}/snps.vcf.gz"
    output:
        "results/default/stats/{sample_set}/{flag}/D.Tajima.D"
    params:
        prefix="results/default/stats/{sample_set}/{flag}/D"
    conda:
        "../envs/vcftools.yaml"
    script:
        "../scripts/calculate_tajimas_d.py"

# plot frequency distribution of Fst values
rule plot_Fst:
    input:
        "results/default/stats/{sample_set}/{flag}/fst.weir.fst"
    output:
        "results/default/graphs/svg/fst/{sample_set}/{flag}/fst.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_fst.py"

# plot windowed frequency distribution of Fst values
rule plot_Fst_windowed:
    input:
        "results/default/stats/{sample_set}/{flag}/fst.windowed.weir.fst"
    output:
        "results/default/graphs/svg/fst/{sample_set}/{flag}/fst.windowed.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_fst_windowed.py"

# plot frequency distribution of Tajima's D values
rule plot_Tajimas_D:
    input:
        "results/default/stats/{sample_set}/{flag}/D.Tajima.D"
    output:
        "results/default/graphs/svg/tajimas_d/{sample_set}/{flag}/tajimas_d.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_tajimas_d.py"

# plot frequency distribution of Pi values
rule plot_pi:
    input:
        "results/default/stats/{sample_set}/{flag}/pi.sites.pi"
    output:
        "results/default/graphs/svg/pi/{sample_set}/{flag}/pi.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_pi.py"

# calculate number of sites as well as
# number of samples for all sample sets
rule collect_sample_set_stats:
    input:
        vcfs=expand("results/{{sample_class}}/snps/{sample_set}/all/snps.vcf.gz",sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"])
    output:
        "results/{sample_class}/stats/sample_sets.csv"
    params:
        sample_sets=sample_sets,
        sample_class="{sample_class}"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/collect_sample_set_stats.py"

# rule to determine expected and observed heterozygosity
# as well as a mid-p value for HWE
rule calculate_hwe_midp:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed"
    output:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.midp.hwe"
    params:
        input_prefix="results/{sample_class}/snps/{sample_set}/{flag}/snps",
        output_prefix="results/{sample_class}/stats/{sample_set}/{flag}/hwe.midp",
        midp=True
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/calculate_hwe.py"

# rule to determine expected and observed heterozygosity
# as well as a p value for HWE
rule calculate_hwe_p:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed"
    output:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.hwe"
    params:
        input_prefix="results/{sample_class}/snps/{sample_set}/{flag}/snps",
        output_prefix="results/{sample_class}/stats/{sample_set}/{flag}/hwe",
        midp=False
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/calculate_hwe.py"

# plot frequency distribution of HWE p-values
rule plot_hwe_midp:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.midp.hwe"
    output:
        "results/{sample_class}/graphs/svg/hwe/{sample_set}/{flag}/midp.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_hwe.py"

# plot frequency distribution of HWE p-values
rule plot_hwe_p:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.hwe"
    output:
        "results/{sample_class}/graphs/svg/hwe/{sample_set}/{flag}/p.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_hwe.py"

# plot frequency distribution of observed or expected heterozygosity
rule plot_heterozygosity:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.hwe"
    output:
        "results/{sample_class}/graphs/svg/heterozygosity/{sample_set}/{flag}/het.{type}.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log",
        observed=lambda w: w.type == "observed"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_heterozygosity.py"

# plot overlaid frequency distribution of observed or expected heterozygosity
rule plot_heterozygosity_comparison:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/hwe.hwe"
    output:
        "results/{sample_class}/graphs/svg/heterozygosity/{sample_set}/{flag}/het.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_heterozygosity_comparison.py"

# calculate missingness over sites and samples
rule calculate_missingness:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed"
    output:
        "results/{sample_class}/stats/{sample_set}/{flag}/missingness.imiss",
        "results/{sample_class}/stats/{sample_set}/{flag}/missingness.lmiss"
    params:
        input_prefix="results/{sample_class}/snps/{sample_set}/{flag}/snps",
        output_prefix="results/{sample_class}/stats/{sample_set}/{flag}/missingness"
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/calculate_missingness.py"

# plot distribution of missingness over sites
rule plot_missingness_sites:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/missingness.lmiss"
    output:
        "results/{sample_class}/graphs/svg/missingness/{sample_set}/{flag}/missingness.l.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_missingness_sites.py"

# plot distribution of missingness over sites
rule plot_missingness_individuals:
    input:
        "results/{sample_class}/stats/{sample_set}/{flag}/missingness.imiss"
    output:
        "results/{sample_class}/graphs/svg/missingness/{sample_set}/{flag}/missingness.i.{scale}.svg"
    params:
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_missingness_individuals.py"

# calculate linkage over sites
rule calculate_linkage:
    input:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed"
    output:
        "results/{sample_class}/stats/{sample_set}/{flag}/linkage.ld"
    params:
        input_prefix="results/{sample_class}/snps/{sample_set}/{flag}/snps",
        output_prefix="results/{sample_class}/stats/{sample_set}/{flag}/linkage"
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/calculate_linkage.py"
