"""
Plot 2D SFS comparison.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import pandas as pd

import demographic_scenarios
import plot_utils

try:
    test_mode = False
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    result = snakemake.input.result
    scenario_name = snakemake.params.scenario
    sample_set = snakemake.params.sample_set
    n_pops = snakemake.params.n_pops

    out_data = snakemake.output.data
    out_model = snakemake.output.model
    out_residuals = snakemake.output.residuals

    generation_time = snakemake.config['generation_time']
    mu = snakemake.config['mu']
    N_e = snakemake.config['N_e'][sample_set]
except NameError:
    # testing
    test_mode = True
    sample_set = "pendula"
    n_pops = 2
    vcf = "output/tetraploid/snps/pendula/synonymous/snps.vcf.gz"
    pops = f"output/tetraploid/sample_sets/subpopulations/pendula/{n_pops}_pops.txt"
    scenario_name = 'split_asymmetric_migration'
    result = f"output/tetraploid/dadi_old/{scenario_name}/pendula/synonymous/dadi.csv"

    out_data = "scratch/data.svg"
    out_model = "scratch/model.svg"
    out_residuals = "scratch/residuals.svg"

    generation_time = 20
    mu = 1e-9
    N_e = 675000

# load params from file
params = pd.read_csv(result, index_col=None, header=0).iloc[0]

# determine initial value
param_names = demographic_scenarios.get_class(scenario_name).params
p0 = [params[name] for name in param_names]

scenario = demographic_scenarios.restore(scenario_name, p0, vcf, pops,
                                         generation_time, n_pops, N_e, mu,
                                         sample_set=sample_set)

scenario.plot_data_2d()
plot_utils.save_fig(out_data, tight_layout=True, show=test_mode)

scenario.plot_model_2d()
plot_utils.save_fig(out_model, tight_layout=True, show=test_mode)

scenario.plot_residuals_2d()
plot_utils.save_fig(out_residuals, tight_layout=True, show=test_mode)
