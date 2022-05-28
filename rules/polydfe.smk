include: "header.smk"

polydfe_models = ['A', 'B', 'C', 'D']

# The type of DFEs. The order of the list is important here
# as we use it to pass the right id to polyDFE
polydfe_types = ['full', 'deleterious', 'full_anc', 'deleterious_anc']

polydfe_n_bootstraps = range(1,config['polydfe_n_bootstraps'] + 1)


# conform with polyDFE's output file numbering format
def pad_number(n):
    return str(n).zfill(len(str(config['polydfe_n_bootstraps'])))


wildcard_constraints:
    polydfe_model=make_wildcard(polydfe_models),
    polydfe_type=make_wildcard(polydfe_types),

# rules to execute locally
localrules:
    fetch_postprocessing_script_polyDFE

# target output files
rule run_all_polydfe:
    input:
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.A.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens'],sample_class=sample_classes),
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.B.png",sample_set=['pendula_pubescens', 'pendula', 'pubescens', 'pubescens_north', 'pubescens_south'],sample_class=sample_classes),
        expand("results/default/graphs/png/polydfe/dfe_subpopulations_bs_{bs_type}.C.{dfe_type}.png",bs_type=['percentile', 'bca'],dfe_type=['full_anc']),
        expand("results/tetraploid/graphs/png/polydfe/dfe_subpopulations_bs_{bs_type}.C.{dfe_type}.png",bs_type=['percentile', 'bca'],dfe_type=['full_anc', 'deleterious_anc']),
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_type.C.png",sample_set=['pubescens', 'pendula_pubescens', 'pendula', 'pendula_south', 'pendula_north', 'pubescens_south', 'pubescens_north'],sample_class=sample_classes),
        expand("results/tetraploid/graphs/png/polydfe/{sample_set}/dfe_type_bs_bca.C.png",sample_set=['pubescens', 'pendula']),
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_bs_{bs_type}.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],bs_type=['percentile', 'bca'],sample_class=sample_classes),
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_grad_error.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],bs_type=['percentile', 'bca'],sample_class=sample_classes),
        expand("results/default/graphs/png/polydfe/{sample_set}/dfe_param_dist.C.full_anc.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south']),
        expand("results/tetraploid/graphs/png/polydfe/{sample_set}/dfe_param_dist.C.deleterious_anc.png",sample_set=['pendula', 'pubescens']),
        expand("results/{sample_class}/graphs/png/polydfe/{sample_set}/dfe_probs_nested.C.png",sample_set=['pendula', 'pubescens', 'pendula_pubescens', 'pendula_north', 'pendula_south', 'pubescens_north', 'pubescens_south'],sample_class=sample_classes)

# annotated resource VCF that can be obtain from the SNP calling
# part of the pipeline
rule resource_vcf_polydfe:
    output:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz"

# subpopulation file
rule resource_subpopulation_files_polydfe:
    output:
        "results/{sample_class}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt"

# fetch the postprocessing script from the GitHub repository
rule fetch_postprocessing_script_polyDFE:
    output:
        config['polydfe_postprocessing_source']
    params:
        url=config['polydfe_postprocessing_source_url']
    shell:
        "wget {params.url} -O {output}"

# create polyDFE SFS file
rule prepare_input_polyDFE:
    input:
        synonymous="results/{sample_class}/snps/{sample_set}/synonymous/snps.vcf.gz",
        nonsynonymous="results/{sample_class}/snps/{sample_set}/nonsynonymous/snps.vcf.gz",
        pops="results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt",
        vcf="results/{sample_class}/snps/{sample_set}/all/snps.vcf.gz",
        targets=config['targeted_genes']
    output:
        "results/{sample_class}/polydfe/{sample_set}/spectra/sfs.txt"
    params:
        n_pops=1,
        unfolded=True
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/prepare_input_polydfe.py"


# run polyDFE on one population
rule run_polyDFE:
    input:
        sfs="results/{sample_class}/polydfe/{sample_set}/spectra/sfs.txt",
        init="resources/polydfe/init_model_{polydfe_model}.txt"
    output:
        "results/{sample_class}/polydfe/{sample_set}/out/{polydfe_model}.{polydfe_type}.out"
    params:
        m="{polydfe_model}",
        id=lambda w: polydfe_types.index(w.polydfe_type) + 1
    conda:
        "../envs/polydfe.yaml"
    script:
        "../scripts/run_polydfe.py"

# run polyDFE on the bootstraps
rule run_polydfe_n_bootstraps:
    input:
        sfs="results/{sample_class}/polydfe/{sample_set}/spectra/bootstraps/bs_{n}",
        init="results/{sample_class}/polydfe/{sample_set}/init/{polydfe_model}.{polydfe_type}_init"
    output:
        "results/{sample_class}/polydfe/{sample_set}/out/bootstraps/{polydfe_model}.{polydfe_type}.{n}.out"
    params:
        m="{polydfe_model}",
        id=1
    conda:
        "../envs/polydfe.yaml"
    script:
        "../scripts/run_polydfe.py"

# extract the discretized DFE from the polyDFE output
rule discretize_dfe:
    input:
        "results/{sample_class}/polydfe/{sample_set}/out/{polydfe_model}.{polydfe_type}.out",
        config['polydfe_postprocessing_source']
    output:
        "results/{sample_class}/polydfe/{sample_set}/out/dfe_{polydfe_model}.{polydfe_type}.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/parse_dfe_output.R"

# plot discretized DFE values
rule plot_polyDFE_types:
    input:
        unpack(lambda w: {t: "results/{{sample_class}}/polydfe/{{sample_set}}/out/dfe_{{polydfe_model}}.{}.txt".format(t)
                          for t in ['full', 'full_anc', 'deleterious', 'deleterious_anc']})
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_type.{polydfe_model}.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_types.py"

# plot discretized DFE values
rule plot_polyDFE_types_bootstrapped:
    input:
        full=["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_C.full.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        full_anc=["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_C.full_anc.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        deleterious=["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_C.deleterious.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        deleterious_anc=["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_C.deleterious_anc.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps]
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_type_bs_{bs_type}.C.svg"
    params:
        bs_type="{bs_type}",
        label_type="polydfe_type"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_bootstrapped.py"

# create bootstrapped spectra
rule create_bootstrapped_spectra_polyDFE:
    input:
        "results/{sample_class}/polydfe/{sample_set}/spectra/sfs.txt",
        config['polydfe_postprocessing_source']
    output:
        ["results/{{sample_class}}/polydfe/{{sample_set}}/spectra/bootstraps/bs_{}".format(pad_number(n)) for n in polydfe_n_bootstraps]
    params:
        n=config['polydfe_n_bootstraps'],
        prefix="results/{sample_class}/polydfe/{sample_set}/spectra/bootstraps/bs"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/create_bootstraps_dfe.R"

# create init file for bootstraps
rule create_init_bootstraps_polyDFE:
    input:
        "results/{sample_class}/polydfe/{sample_set}/out/C.{polydfe_type}.out",
        config['polydfe_postprocessing_source']
    output:
        "results/{sample_class}/polydfe/{sample_set}/init/C.{polydfe_type}_init"
    params:
        prefix="results/{sample_class}/polydfe/{sample_set}/init/C.{polydfe_type}",
        type="{polydfe_type}"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/create_init_bootstraps_dfe.R"

# extract the discretized DFE from the bootstrapped polyDFE output
rule discretize_dfe_bootstrapped:
    input:
        "results/{sample_class}/polydfe/{sample_set}/out/bootstraps/{polydfe_model}.{polydfe_type}.{n}.out",
        config['polydfe_postprocessing_source']
    output:
        "results/{sample_class}/polydfe/{sample_set}/out/bootstraps/dfe_{polydfe_model}.{polydfe_type}.{n}.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/parse_dfe_output.R"

# plot discretized DFE with confidence intervals
rule plot_polyDFE_bootstrapped:
    input:
        unpack(lambda w: {w.sample_set: ["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps]})
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_bs_{bs_type}.{polydfe_model}.{polydfe_type}.svg"
    params:
        bs_type="{bs_type}",
        label_type="sample_set"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_bootstrapped.py"

# plot distribution of bootstrapped interval values
rule plot_polyDFE_param_dist:
    input:
        ["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps]
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_param_dist.{polydfe_model}.{polydfe_type}.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_param_dist.py"

# plot DFE values for different subpopulations including their confidence intervals
rule plot_polyDFE_subpopulations:
    input:
        pendula_north=["results/{{sample_class}}/polydfe/pendula_north/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        pendula_south=["results/{{sample_class}}/polydfe/pendula_south/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        pubescens_north=["results/{{sample_class}}/polydfe/pubescens_north/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps],
        pubescens_south=["results/{{sample_class}}/polydfe/pubescens_south/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps]
    output:
        "results/{sample_class}/graphs/svg/polydfe/dfe_subpopulations_bs_{by_type}.{polydfe_model}.{polydfe_type}.svg"
    params:
        bs_type="{by_type}",
        label_type="sample_set"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_bootstrapped.py"

# plot discretized DFE vs gradient
rule plot_polyDFE_bootstrapped_gradient:
    input:
        bootstrapped=["results/{{sample_class}}/polydfe/{{sample_set}}/out/bootstraps/dfe_{{polydfe_model}}.{{polydfe_type}}.{}.txt".format(
            pad_number(n)) for n in polydfe_n_bootstraps]
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_grad_error.{polydfe_model}.{polydfe_type}.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_dfe_gradient_error.py"

# compare the nested DFE types using LRTs
rule compare_nested_DFE_types:
    input:
        unpack(lambda w: {t: "results/{{sample_class}}/polydfe/{{sample_set}}/out/{{polydfe_model}}.{}.out".format(t)
                          for t in polydfe_types}),
        postprocessing=config['polydfe_postprocessing_source']
    output:
        "results/{sample_class}/polydfe/{sample_set}/type_comparison.{polydfe_model}.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/compare_dfe_types.R"

# plot tiles indicating the p-values of the nested models
rule plot_props_nested_DFE_types:
    input:
        "results/{sample_class}/polydfe/{sample_set}/type_comparison.{polydfe_model}.txt"
    output:
        "results/{sample_class}/graphs/svg/polydfe/{sample_set}/dfe_probs_nested.{polydfe_model}.svg"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_props_nested_dfe_types.py"
