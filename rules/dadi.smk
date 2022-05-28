include: "header.smk"

# 1D population scenarios
pop_scenarios_1d = [
    'constant_pop_size',
    'lin_pop_growth',
    'exp_pop_growth',
    'one_pop_size_change',
    'two_pop_size_change',
    'constant_pop_size_since_ice_age',
    'lin_pop_growth_since_ice_age',
    'exp_pop_growth_since_ice_age',
    'one_pop_size_change_since_ice_age',
    'two_pop_size_change_since_ice_age'
]

# 2D population scenarios
pop_scenarios_2d = [
    'split_no_migration',
    'split_no_migration_since_ice_age',
    'split_no_migration_one_pop_size_change',
    'split_no_migration_one_pop_size_change_since_ice_age',
    'split_symmetric_migration',
    'split_asymmetric_migration',
    'split_symmetric_migration_since_ice_age',
    'split_asymmetric_migration_since_ice_age',
    'split_symmetric_migration_one_pop_size_change',
    'split_symmetric_migration_one_pop_size_change_since_ice_age',
    'split_asymmetric_migration_one_pop_size_change',
    'split_asymmetric_migration_one_pop_size_change_since_ice_age',
    'split_unidirectional_migration_from_two_to_one_since_ice_age',
    'split_unidirectional_migration_from_two_to_one',
    'split_unidirectional_migration_from_one_to_two_since_ice_age',
    'split_unidirectional_migration_from_one_to_two'
]

pop_scenarios = pop_scenarios_1d + pop_scenarios_2d

pop_scenarios_since_ice_age_1d = [s for s in pop_scenarios_1d if 'since_ice_age' in s]
pop_scenarios_since_ice_age_2d = [s for s in pop_scenarios_2d if 'since_ice_age' in s]
pop_scenarios_since_ice_age = pop_scenarios_since_ice_age_1d + pop_scenarios_since_ice_age_2d

optimization_modes_dadi = ['random_hopping', 'local_optimization']

dadi_n_batches_global_optimization = range(1,config['dadi_n_batches_global_optimization'] + 1)
dadi_n_bootstraps = range(1,config['dadi_n_bootstraps'] + 1)


# get the number of populations for the given population scenario
def get_n_pops(scenario):
    if scenario in pop_scenarios_1d:
        return 1

    return 2


wildcard_constraints:
    pop_scenario_1d=make_wildcard(pop_scenarios_1d),
    pop_scenario_2d=make_wildcard(pop_scenarios_2d),
    pop_scenario=make_wildcard(pop_scenarios),
    optimization_mode_dadi=make_wildcard(optimization_modes_dadi),
    _since_ice_age='(\\_since\\_ice\\_age)*',
    _ci='(\\_ci)*'


# target output files
rule run_all_dadi:
    input:
        expand("results/tetraploid/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.out",sample_set=['pendula', 'pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_1d),
        expand("results/tetraploid/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.out",sample_set=['pendula_pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_2d),

        expand("results/{sample_class}/stats/{sample_set}/stats.txt",sample_set=['pendula', 'pubescens', 'pendula_pubescens'],sample_class=sample_classes),
        expand("results/{sample_class}/stats/stats.txt",sample_class=sample_classes),

        expand("results/tetraploid/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula', 'pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_1d,type=['trajectory', 'sfs', 'residuals']),

        expand("results/tetraploid/graphs/png/dadi/{pop_scenario}/{sample_set}/{flag}/{type}.png",sample_set=['pendula_pubescens'],flag=['synonymous'],pop_scenario=pop_scenarios_2d,type=['data', 'model', 'residuals']),

        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_1d.png",sample_set=['pendula', 'pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_1d_ci.png",sample_set=['pendula', 'pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_1d_since_ice_age.png",sample_set=['pendula', 'pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_1d_since_ice_age_ci.png",sample_set=['pendula', 'pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_lgm_1d.png",sample_set=['pendula', 'pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_lgm_1d_ci.png",sample_set=['pendula', 'pubescens']),

        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_2d.png",sample_set=['pendula_pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_2d_ci.png",sample_set=['pendula_pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_lgm_2d.png",sample_set=['pendula_pubescens']),
        expand("results/tetraploid/graphs/png/dadi/{sample_set}/nested_models_lgm_2d_ci.png",sample_set=['pendula_pubescens']),

        expand("results/tetraploid/dadi/{sample_set}/results_1d.{ext}",sample_set=['pendula', 'pubescens'],ext=['out', 'tex']),
        expand("results/tetraploid/dadi/{sample_set}/results_2d.{ext}",sample_set=['pendula_pubescens'],ext=['out', 'tex'])

# annotated resource VCF that can be obtain from the SNP calling
# part of the pipeline
rule resource_vcf_dadi:
    output:
        protected("results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz")

# subpopulation file
rule resource_subpopulation_files_dadi:
    output:
        protected("results/{sample_class}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt")

# fit various population scenarios to the given SFS using dadi
rule fit_pop_scenarios_dadi_batches:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt".format(
            sample_set=w.sample_set,n=get_n_pops(w.pop_scenario))
    output:
        params="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/batches/{n}/dadi.csv",
        params_pretty="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/batches/{n}/dadi.out"
    params:
        n_pops=lambda w: get_n_pops(w.pop_scenario),
        scenario="{pop_scenario}",
        sample_set="{sample_set}",
        mode="random_hopping"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/fit_pop_scenarios_dadi.py"

# determine the batch with the highest likelihood
rule determine_best_batch_dadi:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt".format(
            sample_set=w.sample_set,n=get_n_pops(w.pop_scenario)),
        results=expand("results/{{sample_class}}/dadi/{{pop_scenario}}/{{sample_set}}/{{flag}}/batches/{n}/dadi.csv",
            n=dadi_n_batches_global_optimization)
    output:
        params="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi.csv",
        params_pretty="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi.out"
    params:
        n_pops=lambda w: get_n_pops(w.pop_scenario),
        scenario="{pop_scenario}",
        sample_set="{sample_set}"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/determine_best_batch_dadi.py"

# plot the trajectory of 1D population scenarios
rule plot_trajectory_1D_dadi:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops="results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt",
        result="results/{sample_class}/dadi/{pop_scenario_1d}/{sample_set}/{flag}/dadi.csv"
    output:
        "results/{sample_class}/graphs/svg/dadi/{pop_scenario_1d}/{sample_set}/{flag}/trajectory.svg"
    params:
        n_pops=1,
        scenario="{pop_scenario_1d}",
        sample_set="{sample_set}"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_trajectory_dadi.py"

# plot data SFS vs model SFS
rule plot_sfs_comparison_dadi_2D:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/2_pops.txt".format(
            sample_set=w.sample_set),
        result="results/{sample_class}/dadi/{pop_scenario_2d}/{sample_set}/{flag}/dadi.csv"
    output:
        data="results/{sample_class}/graphs/svg/dadi/{pop_scenario_2d}/{sample_set}/{flag}/data.svg",
        model="results/{sample_class}/graphs/svg/dadi/{pop_scenario_2d}/{sample_set}/{flag}/model.svg",
        residuals="results/{sample_class}/graphs/svg/dadi/{pop_scenario_2d}/{sample_set}/{flag}/residuals.svg"
    params:
        n_pops=2,
        scenario="{pop_scenario_2d}",
        sample_set="{sample_set}"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_comparison_2D.py"

# plot data SFS vs model SFS
rule plot_sfs_comparison_dadi_1D:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/1_pops.txt".format(
            sample_set=w.sample_set),
        result="results/{sample_class}/dadi/{pop_scenario_1d}/{sample_set}/{flag}/dadi.csv"
    output:
        sfs="results/{sample_class}/graphs/svg/dadi/{pop_scenario_1d}/{sample_set}/{flag}/sfs.svg",
        residuals="results/{sample_class}/graphs/svg/dadi/{pop_scenario_1d}/{sample_set}/{flag}/residuals.svg"
    params:
        n_pops=1,
        scenario="{pop_scenario_1d}",
        sample_set="{sample_set}"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_comparison_1D.py"

# create bootstraps for the dadi simulations
rule bootstrap_pop_scenarios_dadi:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt".format(
            sample_set=w.sample_set,n=get_n_pops(w.pop_scenario)),
        result="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi.csv"
    output:
        params="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/bootstraps/{n}/dadi.csv",
        params_pretty="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/bootstraps/{n}/dadi.out"
    params:
        n_pops=lambda w: get_n_pops(w.pop_scenario),
        scenario="{pop_scenario}",
        sample_set="{sample_set}",
        mode="local_optimization"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/bootstrap_pop_scenarios_dadi.py"

# determine confidence intervals using the bootstraps
rule determine_cis_dadi_1D:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz",
        pops=lambda w: "results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/{n}_pops.txt".format(
            sample_set=w.sample_set,n=get_n_pops(w.pop_scenario)),
        results=expand("results/{{sample_class}}/dadi/{{pop_scenario}}/{{sample_set}}/{{flag}}/bootstraps/{n}/dadi.csv",
            n=dadi_n_bootstraps),
        result_original="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi.csv"
    output:
        params="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.csv",
        params_pretty="results/{sample_class}/dadi/{pop_scenario}/{sample_set}/{flag}/dadi_ci.out"
    params:
        n_pops=lambda w: get_n_pops(w.pop_scenario),
        scenario="{pop_scenario}",
        sample_set="{sample_set}",
        mode="local_optimization"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/determine_cis_dadi.py"

# plot a graph visualizing the p-values for various
# nested fixed-time vs. variable-time models
rule plot_models_LRT_dadi_LGM_1d:
    input:
        constant_pop_size="results/tetraploid/dadi/constant_pop_size/{sample_set}/synonymous/dadi{_ci}.csv",
        lin_pop_growth="results/tetraploid/dadi/lin_pop_growth/{sample_set}/synonymous/dadi{_ci}.csv",
        exp_pop_growth="results/tetraploid/dadi/exp_pop_growth/{sample_set}/synonymous/dadi{_ci}.csv",
        one_pop_size_change="results/tetraploid/dadi/one_pop_size_change/{sample_set}/synonymous/dadi{_ci}.csv",
        two_pop_size_change="results/tetraploid/dadi/two_pop_size_change/{sample_set}/synonymous/dadi{_ci}.csv",
        constant_pop_size_since_ice_age="results/tetraploid/dadi/constant_pop_size_since_ice_age/{sample_set}/synonymous/dadi{_ci}.csv",
        lin_pop_growth_since_ice_age="results/tetraploid/dadi/lin_pop_growth_since_ice_age/{sample_set}/synonymous/dadi{_ci}.csv",
        exp_pop_growth_since_ice_age="results/tetraploid/dadi/exp_pop_growth_since_ice_age/{sample_set}/synonymous/dadi{_ci}.csv",
        one_pop_size_change_since_ice_age="results/tetraploid/dadi/one_pop_size_change_since_ice_age/{sample_set}/synonymous/dadi{_ci}.csv",
        two_pop_size_change_since_ice_age="results/tetraploid/dadi/two_pop_size_change_since_ice_age/{sample_set}/synonymous/dadi{_ci}.csv"
    output:
        "results/tetraploid/graphs/svg/dadi/{sample_set}/nested_models_lgm_1d{_ci}.svg"
    params:
        comparisons=[
            ['constant_pop_size_since_ice_age', 'constant_pop_size', 1],
            ['lin_pop_growth_since_ice_age', 'lin_pop_growth', 1],
            ['exp_pop_growth_since_ice_age', 'exp_pop_growth', 1],
            ['one_pop_size_change_since_ice_age', 'one_pop_size_change', 1],
            ['two_pop_size_change_since_ice_age', 'two_pop_size_change', 1]
        ],
        labels=[
            'constant',
            'lin growth',
            'exp growth',
            '1-size change',
            '2-size change',
            'constant LGM',
            'lin growth LGM',
            'exp growth LGM',
            '1-size change LGM',
            '2-size change LGM',
        ],
        scaling=1
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_models_lrt_dadi.py"

# plot a graph visualizing the p-values for various
# nested fixed-time vs. variable-time models
rule plot_models_LRT_dadi_LGM_2d:
    input:
        split_no_migration="results/tetraploid/dadi/split_no_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_two_to_one="results/tetraploid/dadi/split_unidirectional_migration_from_two_to_one/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_one_to_two="results/tetraploid/dadi/split_unidirectional_migration_from_one_to_two/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_symmetric_migration="results/tetraploid/dadi/split_symmetric_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration="results/tetraploid/dadi/split_asymmetric_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration_one_pop_size_change="results/tetraploid/dadi/split_asymmetric_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_no_migration_since_ice_age="results/tetraploid/dadi/split_no_migration_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_two_to_one_since_ice_age="results/tetraploid/dadi/split_unidirectional_migration_from_two_to_one_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_one_to_two_since_ice_age="results/tetraploid/dadi/split_unidirectional_migration_from_one_to_two_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_symmetric_migration_since_ice_age="results/tetraploid/dadi/split_symmetric_migration_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration_since_ice_age="results/tetraploid/dadi/split_asymmetric_migration_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration_one_pop_size_change_since_ice_age="results/tetraploid/dadi/split_asymmetric_migration_one_pop_size_change_since_ice_age/pendula_pubescens/synonymous/dadi{_ci}.csv"
    output:
        "results/tetraploid/graphs/svg/dadi/pendula_pubescens/nested_models_lgm_2d{_ci}.svg"
    params:
        comparisons=[
            ['split_no_migration_since_ice_age', 'split_no_migration', 1],
            ['split_unidirectional_migration_from_two_to_one_since_ice_age', 'split_unidirectional_migration_from_two_to_one', 1],
            ['split_unidirectional_migration_from_one_to_two_since_ice_age', 'split_unidirectional_migration_from_one_to_two', 1],
            ['split_symmetric_migration_since_ice_age', 'split_symmetric_migration', 1],
            ['split_asymmetric_migration_since_ice_age', 'split_asymmetric_migration', 1],
            ['split_asymmetric_migration_one_pop_size_change_since_ice_age', 'split_asymmetric_migration_one_pop_size_change', 1]
        ],
        labels=[
            'no migration',
            'unilateral migr\n$pen$ to $pub$',
            'unilateral migr\n$pub$ to $pen$',
            'sym migration',
            'asym migration',
            'asym migr + \n1-size change',
            'no migration LGM',
            'unilateral migr\n$pen$ to $pub$ LGM',
            'unilateral migr\n$pub$ to $pen$ LGM',
            'sym migr LGM',
            'asym migr LGM',
            'asym migr + \n1-size change LGM'
        ],
        scaling=1.2
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_models_lrt_dadi.py"

# plot a graph visualizing the p-values for various nested models
rule plot_models_LRT_dadi_1d:
    input:
        constant_pop_size="results/tetraploid/dadi/constant_pop_size{_since_ice_age}/{sample_set}/synonymous/dadi{_ci}.csv",
        lin_pop_growth="results/tetraploid/dadi/lin_pop_growth{_since_ice_age}/{sample_set}/synonymous/dadi{_ci}.csv",
        exp_pop_growth="results/tetraploid/dadi/exp_pop_growth{_since_ice_age}/{sample_set}/synonymous/dadi{_ci}.csv",
        one_pop_size_change="results/tetraploid/dadi/one_pop_size_change{_since_ice_age}/{sample_set}/synonymous/dadi{_ci}.csv",
        two_pop_size_change="results/tetraploid/dadi/two_pop_size_change{_since_ice_age}/{sample_set}/synonymous/dadi{_ci}.csv"
    output:
        "results/tetraploid/graphs/svg/dadi/{sample_set}/nested_models_1d{_since_ice_age}{_ci}.svg"
    params:
        comparisons=[
            # the likelihood of the by_slope scenarios is indistinguishable from
            # the other parameterization so that we can regard them as nested
            ['constant_pop_size', 'lin_pop_growth', 1],
            ['constant_pop_size', 'exp_pop_growth', 1],
            ['constant_pop_size', 'one_pop_size_change', 2],
            ['one_pop_size_change', 'two_pop_size_change', 2],
            ['constant_pop_size', 'two_pop_size_change', 4],
        ],
        labels=[
            'constant',
            'lin growth',
            'exp growth',
            '1-size change',
            '2-size change',
        ],
        scaling=1
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_models_lrt_dadi.py"

# plot a graph visualizing the p-values for various nested models
rule plot_models_LRT_dadi_2d:
    input:
        split_no_migration="results/tetraploid/dadi/split_no_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_no_migration_one_pop_size_change="results/tetraploid/dadi/split_no_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_two_to_one="results/tetraploid/dadi/split_unidirectional_migration_from_two_to_one/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_unidirectional_migration_from_one_to_two="results/tetraploid/dadi/split_unidirectional_migration_from_one_to_two/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_symmetric_migration="results/tetraploid/dadi/split_symmetric_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration="results/tetraploid/dadi/split_asymmetric_migration/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_symmetric_migration_one_pop_size_change="results/tetraploid/dadi/split_symmetric_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi{_ci}.csv",
        split_asymmetric_migration_one_pop_size_change="results/tetraploid/dadi/split_asymmetric_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi{_ci}.csv"
    output:
        "results/tetraploid/graphs/svg/dadi/pendula_pubescens/nested_models_2d{_ci}.svg"
    params:
        comparisons=[
            # the likelihood of the by_slope scenarios is indistinguishable from
            # the other parameterization so that we can regard them as nested
            ['split_no_migration', 'split_symmetric_migration', 1],
            ['split_no_migration', 'split_unidirectional_migration_from_two_to_one', 1],
            ['split_no_migration', 'split_unidirectional_migration_from_one_to_two', 1],
            ['split_no_migration', 'split_asymmetric_migration', 2],
            ['split_unidirectional_migration_from_two_to_one', 'split_asymmetric_migration', 1],
            ['split_unidirectional_migration_from_one_to_two', 'split_asymmetric_migration', 1],
            ['split_no_migration_one_pop_size_change', 'split_symmetric_migration_one_pop_size_change', 1],
            ['split_no_migration_one_pop_size_change', 'split_asymmetric_migration_one_pop_size_change', 2],
            ['split_symmetric_migration', 'split_symmetric_migration_one_pop_size_change', 2],
            ['split_asymmetric_migration', 'split_asymmetric_migration_one_pop_size_change', 2],
        ],
        labels=[
            'no migration',
            'no mig + \n1-size change',
            'unilateral migr\n$pen$ to $pub$',
            'unilateral migr\n$pub$ to $pen$',
            'sym migration',
            'asym migration',
            'sym migr + \n1-size change',
            'asym migr + \n1- size change'
        ],
        scaling=1.2,
        transpose=True
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/plot_models_lrt_dadi.py"

# summarize the 1D results
rule tabulate_results_dadi_1d:
    input:
        constant_pop_size="results/tetraploid/dadi/constant_pop_size/{sample_set}/synonymous/dadi_ci.csv",
        lin_pop_growth="results/tetraploid/dadi/lin_pop_growth/{sample_set}/synonymous/dadi_ci.csv",
        exp_pop_growth="results/tetraploid/dadi/exp_pop_growth/{sample_set}/synonymous/dadi_ci.csv",
        one_pop_size_change="results/tetraploid/dadi/one_pop_size_change/{sample_set}/synonymous/dadi_ci.csv",
        constant_pop_size_since_ice_age="results/tetraploid/dadi/constant_pop_size_since_ice_age/{sample_set}/synonymous/dadi_ci.csv",
        lin_pop_growth_since_ice_age="results/tetraploid/dadi/lin_pop_growth_since_ice_age/{sample_set}/synonymous/dadi_ci.csv",
        exp_pop_growth_since_ice_age="results/tetraploid/dadi/exp_pop_growth_since_ice_age/{sample_set}/synonymous/dadi_ci.csv",
        one_pop_size_change_since_ice_age="results/tetraploid/dadi/one_pop_size_change_since_ice_age/{sample_set}/synonymous/dadi_ci.csv"
    output:
        pretty="results/tetraploid/dadi/{sample_set}/results_1d.out",
        tex="results/tetraploid/dadi/{sample_set}/results_1d.tex"
    params:
        labels=[
            'constant size',
            'lin. growth',
            'exp. growth',
            '1-size change',
            'constant size LGM',
            'lin. growth LGM',
            'exp. growth LGM',
            '1-size change LGM'
        ]
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/tabulate_results_dadi.py"

# summarize the 2D results
rule tabulate_results_dadi_2d:
    input:
        split_no_migration="results/tetraploid/dadi/split_no_migration/pendula_pubescens/synonymous/dadi_ci.csv",
        split_unidirectional_migration_from_two_to_one="results/tetraploid/dadi/split_unidirectional_migration_from_two_to_one/pendula_pubescens/synonymous/dadi_ci.csv",
        split_unidirectional_migration_from_one_to_two="results/tetraploid/dadi/split_unidirectional_migration_from_one_to_two/pendula_pubescens/synonymous/dadi_ci.csv",
        split_symmetric_migration="results/tetraploid/dadi/split_symmetric_migration/pendula_pubescens/synonymous/dadi_ci.csv",
        split_asymmetric_migration="results/tetraploid/dadi/split_asymmetric_migration/pendula_pubescens/synonymous/dadi_ci.csv",
        split_symmetric_migration_one_pop_size_change="results/tetraploid/dadi/split_symmetric_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi_ci.csv",
        split_asymmetric_migration_one_pop_size_change="results/tetraploid/dadi/split_asymmetric_migration_one_pop_size_change/pendula_pubescens/synonymous/dadi_ci.csv",
    output:
        pretty="results/tetraploid/dadi/{sample_set}/results_2d.out",
        tex="results/tetraploid/dadi/{sample_set}/results_2d.tex"
    params:
        labels=[
            'no migration',
            'unilateral migr pen to pub',
            'unilateral migr pub to pen',
            'sym migration',
            'asym migration',
            'sym migr + 1-size change',
            'asym migr + 1-size change'
        ]
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/tabulate_results_dadi.py"

# calculate various population genetics
# statistics over different sample set flags
rule calculate_basic_stats_types:
    input:
        vcfs=expand("results/{{sample_class}}/snps/{{sample_set}}/{flag}/snps.vcf.gz",flag=flags),
        pops_1=lambda w: ["results/{sample_class}/sample_sets/subpopulations/{sample_set}/1_pops.txt"] * len(flags),
        pops_2=lambda w: ["results/{sample_class}/sample_sets/subpopulations/{sample_set}/2_pops.txt"] * len(flags)
    output:
        txt="results/{sample_class}/stats/{sample_set}/stats.txt",
        csv="results/{sample_class}/stats/{sample_set}/stats.csv"
    params:
        names=flags
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/calculate_basic_stats.py"

# calculate various population genetics
# statistics over different sample sets
rule calculate_basic_stats_sample_sets:
    input:
        vcfs=expand("results/{{sample_class}}/snps/{sample_set}/all/snps.vcf.gz",
            sample_set=['pendula', 'pubescens', 'pendula_pubescens']),
        pops_1=expand("results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/1_pops.txt",
            sample_set=['pendula', 'pubescens', 'pendula_pubescens']),
        pops_2=expand("results/{{sample_class}}/sample_sets/subpopulations/{sample_set}/2_pops.txt",
            sample_set=['pendula', 'pubescens', 'pendula_pubescens'])
    output:
        txt="results/{sample_class}/stats/stats.txt",
        csv="results/{sample_class}/stats/stats.csv"
    params:
        names=['pendula', 'pubescens', 'pendula_pubescens']
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/calculate_basic_stats.py"
