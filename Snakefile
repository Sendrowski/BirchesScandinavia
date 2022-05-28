import numpy as np

include: "rules/header.smk"
include: "rules/snp_calling.smk"
include: "rules/snp_calling_tetraploid.smk"
include: "rules/sample_sets.smk"
include: "rules/vcf_filtering.smk"
include: "rules/summary_stats.smk"
include: "rules/feems.smk"
include: "rules/pca.smk"
include: "rules/umap.smk"
include: "rules/polydfe.smk"
include: "rules/dadi.smk"
include: "rules/admixture.smk"

schematics = ['pipeline', 'sample_sets', 'dadi']

components = [
    'snp_calling',
    'snp_calling_tetraploid',
    'sample_sets',
    'vcf_filtering',
    'summary_stats',
    'feems',
    'pca',
    'umap',
    'polydfe',
    'dadi'
]

ruleorder:
    create_sample_set_birch >
    create_sample_set_right_cluster >
    create_sample_set_left_cluster >
    create_sample_set_pendula >
    create_sample_set_pubescens >
    create_sample_set_ingroups >
    create_sample_set_outgroups >
    create_sample_set_pendula_luis >
    create_sample_set_pendula_pubescens >
    create_sample_set_all >
    create_sample_set_pendula_north >
    create_sample_set_pendula_south >
    create_sample_set_pubescens_north >
    create_sample_set_pubescens_south >
    create_sample_set_ADMIXTURE >
    create_subpopulation_files_1_pop >
    create_subpopulation_files_2_pops >
    generate_vcf_sample_set >
    generate_vcf_sample_set_biallelic >
    generate_vcf_sample_set_synonymous >
    generate_vcf_sample_set_nonsynonymous >
    generate_vcf_sample_set_no_missing >
    generate_vcf_sample_set_no_low_freqs >
    generate_vcf_sample_set_no_very_low_freqs >
    filter_nfold_sites >
    filter_0fold_degenerate_sites >
    filter_4fold_degenerate_sites >
    filter_snps >
    create_tbi_index_vcf >
    resource_sample_sets_pca >
    resource_sample_sets_feems >
    resource_sample_sets_summary_stats >
    resource_sample_sets_umap >
    resource_vcf_polydfe >
    resource_vcf_umap >
    resource_vcf_polydfe >
    resource_vcf_dadi >
    resource_vcf_pca >
    resource_vcf_summary_stats >
    resource_vcf_feems >
    resource_vcf_admixture >
    resource_ingroup_samples >
    resource_outgroup_samples >
    resource_raw_vcf >
    resource_reads_birch >
    resource_targeted_genes >
    resource_reads_outgroup >
    resource_annotation_file >
    resource_reference_genome >
    resource_subpopulation_files_dadi >
    resource_subpopulation_files_polydfe


def load_target_files(wildcards):
    targets_remote = [
        expand("results/default/trimmed_reads/{name}.{side}.fastq.gz",name=sample_names,side=[1, 2]),

        # expand("results/default/quality_checks/{name}.1_fastqc.html",name=sample_names),
        # expand("results/default/quality_checks/{name}.2_fastqc.html",name=sample_names),
        # expand("results/default/quality_checks/multiqc_report.html"),

        # expand("results/default/mapped_reads/{name}.bam",name=sample_names),
        # expand("results/default/mapped_reads/{name}.marked_dups.bam",name=sample_names),
        # expand("results/default/mapped_reads/{name}.marked_dups.bai",name=sample_names),

        # expand("results/default/stats/coverage/{name}.txt",name=sample_names),
        # expand("results/default/stats/coverage/all.txt"),
        # expand("results/default/stats/{sample_set}/stats.txt",sample_set=['pendula', 'pubescens', 'pendula_pubescens']),
        # expand("results/{sample_class}/stats/stats.txt",sample_class=sample_classes),

        # expand("results/default/variants/haplotypes/{name}/all.sorted.g.vcf.gz",name=sample_names),
        # expand("results/default/variants/intervals/intervals{n}.bed",n=n_genomic_intervals),
        # expand("results/default/variants/genomics_dbs/db{n}",n=n_genomic_intervals),
        # expand("results/default/variants/vcfs/variants{n}.vcf.gz",n=n_genomic_intervals),
        # expand("results/default/variants/vcfs/variants{n}.quality.tagged.vcf.gz",n=n_genomic_intervals),

        # expand("results/default/snps/raw/intervals/snps{n}.vcf.gz",n=n_genomic_intervals),
        # expand("results/default/snps/raw/intervals/snps{n}.vcf.gz",n=n_genomic_intervals),
        # expand("results/default/snps/raw/intervals/snps{n}.ancestral.vcf.gz",n=n_genomic_intervals),
        # expand("results/{sample_class}/snps/raw/snps.vcf.gz",sample_class=sample_classes),
        # expand("results/default/snps/all/biallelic/snps.vcf.gz"),

        # expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_cv.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs', 'no_missing']),
        # expand("results/default/graphs/png/feems/{sample_set}/{flag}/migration_surfaces_lambda_{n}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'no_low_freqs', 'no_very_low_freqs', 'no_missing'],n=[10, 1, 0.1, 0.01]),
        # expand("results/default/graphs/png/feems/{sample_set}/allocation.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens']),
        # expand("results/default/graphs/png/feems/{sample_set}/locations.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens']),

        # expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens', 'pendula_luis'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        # expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_tight_layout.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens', 'pendula_luis'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        # expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_labeled.png",sample_set=['left_cluster', 'right_cluster', 'pendula_luis', 'birch', 'pendula_pubescens'],flag=['biallelic', 'synonymous'],sample_class=['default']),
        # expand("results/{sample_class}/graphs/png/pca/{sample_set}/{flag}/pca_marked.png",sample_set=['pendula_pubescens', 'birch'],flag=['biallelic', 'synonymous'],sample_class=['default']),

        # expand("results/{sample_class}/graphs/png/umap/{sample_set}/biallelic/umap.png",sample_set=['birch', 'left_cluster', 'right_cluster', 'pendula', 'pubescens', 'pendula_pubescens'],sample_class=['default']),
        # expand("results/{sample_class}/graphs/png/umap/{sample_set}/{flag}/umap_marked.png",sample_set=['pendula_pubescens', 'birch'],flag=['biallelic', 'synonymous'],sample_class=['default']),

        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/synonymous_nonsynonymous/sfs.1D.100.unfolded.linear.png",sample_set=['pendula', 'pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/synonymous_nonsynonymous/sfs.1D.20.unfolded.linear.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/synonymous_nonsynonymous/sfs.1D.20.folded.linear.png",sample_set=['pendula', 'pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/0fold_4fold/sfs.1D.20.unfolded.linear.png",sample_set=['pendula', 'pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/0fold_4fold/sfs.1D.20.unfolded.linear_tight_layout.png",sample_set=['pendula', 'pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.1D.20.unfolded.linear.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs', '0fold', '2fold', '4fold'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.1D.100.unfolded.{scale}.png",sample_set=['pendula', 'pubescens'],flag=['biallelic', 'no_low_freqs'],scale=['linear', 'log'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.1D.198.unfolded.log.png",sample_set=['birch'],flag=['biallelic', 'all'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/{sample_set}/{flag}/sfs.2D.20.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens'],flag=['all', 'biallelic', '0fold', '4fold'],sample_class=sample_classes),

        # expand("results/{sample_class}/graphs/png/sfs/est-sfs.{n_proj}.unfolded.{scale}.png",n_proj=[20, config['est_sfs_n_samples']],scale=['log', 'linear'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/sfs/est-sfs.{n_proj}.folded.{scale}.png",n_proj=[20, int(config['est_sfs_n_samples'] / 2)],scale=['log', 'linear'],sample_class=sample_classes),

        # expand("results/default/graphs/png/pi/{sample_set}/{flag}/pi.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic', '0fold', '4fold'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/tajimas_d/{sample_set}/{flag}/tajimas_d.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.windowed.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/missingness/{sample_set}/{flag}/missingness.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear'],type=['l', 'i']),
        # expand("results/{sample_class}/graphs/png/hwe/{sample_set}/{flag}/{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['p', 'midp']),
        # expand("results/{sample_class}/graphs/png/heterozygosity/{sample_set}/{flag}/het.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['expected', 'observed']),

        # expand("results/default/graphs/png/rulegraphs/schematics/{schematics}.png",schematics=schematics),

        # expand("results/{sample_class}/graphs/png/degeneracy/{sample_set}/degeneracy.png",sample_set=['pendula', 'pubescens', 'birch'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/synonymy/{sample_set}/synonymy.png",sample_set=['pendula', 'pubescens', 'birch'],sample_class=sample_classes),

        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.A.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.B.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens', 'pubescens_north', 'pubescens_south'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/dfe_subpopulations_bs_{bs_type}.C.full_anc.png",bs_type=['percentile', 'bca'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.C.png",sample_set=['pubescens', 'pendula_pubescens', 'pendula', 'pendula_south', 'pendula_north', 'pubescens_south', 'pubescens_north'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_bs_{bs_type}.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],bs_type=['percentile', 'bca'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_bs_{bs_type}.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],bs_type=['percentile', 'bca'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_grad_error.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],bs_type=['percentile', 'bca'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_param_dist.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],sample_class=sample_classes),
        # expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_probs_nested.C.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],sample_class=sample_classes),

        # expand("results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.out",sample_set=['pendula'],flag=['synonymous'],pop_scenario=pop_scenarios,sample_class=['tetraploid']),
        # expand("results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.out",sample_set=['pendula_pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_since_ice_age + ['split_asymmetric_migration'],sample_class=['tetraploid']),

        # expand("results/{sample_class}/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula'],flag=['synonymous'],pop_scenario=pop_scenarios_1d,sample_class=sample_classes,type=['trajectory', 'sfs', 'residuals']),
        # expand("results/{sample_class}/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula_pubescens', 'pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_since_ice_age_1d,sample_class=sample_classes,type=['trajectory', 'sfs', 'residuals']),

        # expand("results/{sample_class}/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula'],flag=['synonymous'],pop_scenario=pop_scenarios_2d,sample_class=['tetraploid'],type=['data', 'model', 'residuals']),
        # expand("results/{sample_class}/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula_pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_since_ice_age_2d + ['split_asymmetric_migration'],sample_class=['tetraploid'],type=['data', 'model', 'residuals']),

        # "results/default/graphs/png/rulegraphs/until_raw_vcf.png",
        # "results/default/graphs/png/rulegraphs/sfs.png",
        # "results/default/graphs/png/rulegraphs/all.png"

        # expand("results/default/graphs/png/pi/{sample_set}/{flag}/pi.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic', '0fold', '4fold'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/tajimas_d/{sample_set}/{flag}/tajimas_d.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/fst/{sample_set}/{flag}/fst.windowed.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear']),
        # expand("results/default/graphs/png/missingness/{sample_set}/{flag}/missingness.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['all', 'biallelic'],scale=['log', 'linear'],type=['l', 'i']),
        # expand("results/{sample_class}/graphs/png/hwe/{sample_set}/{flag}/{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['p', 'midp']),
        # expand("results/{sample_class}/graphs/png/heterozygosity/{sample_set}/{flag}/het.{type}.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'],type=['expected', 'observed']),
        # expand("results/{sample_class}/graphs/png/heterozygosity/{sample_set}/{flag}/het.{scale}.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],flag=['biallelic', 'all'],sample_class=['default'],scale=['log', 'linear'])
    ]

    targets_local = [
        "results/default/graphs/png/rulegraphs/schematics/dadi.png"
    ]

    targets = targets_local if is_local else targets_remote

    return np.hstack(targets)


# target output files
rule all:
    input: load_target_files

# generate an image depicting the rule graph up to the SNP calling
rule visualize_rulegraph_all:
    output:
        "results/default/graphs/svg/rulegraphs/all.svg"
    params:
        path=lambda w: ' '.join(load_target_files(w))
    conda:
        "envs/snakemake.yaml"
    script:
        "scripts/visualize_complete_rulegraph.py"

# visualize the pipeline schematic which was manually created
rule visualize_pipeline_schematic:
    input:
        "resources/schematics/{schematics}.dot"
    output:
        "results/default/graphs/svg/rulegraphs/schematics/{schematics}.svg"
    script:
        "scripts/visualize_dotfile.py"

# visualize dag
rule visualize_component_dags:
    output:
        "results/default/graphs/svg/rulegraphs/{component}.svg"
    params:
        snakefile="{component}.smk"
    conda:
        "envs/snakemake.yaml"
    script:
        "scripts/visualize_complete_rulegraph.py"
