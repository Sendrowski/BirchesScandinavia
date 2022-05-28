import numpy as np
import pandas as pd

import demographic_scenarios
import utils

try:
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    results = snakemake.input.results
    result = snakemake.input.result_original
    scenario_name = snakemake.params.scenario
    sample_set = snakemake.params.sample_set
    n_pops = snakemake.params.n_pops
    params_out = snakemake.output.params
    params_pretty_out = snakemake.output.params_pretty

    generation_time = snakemake.config['generation_time']
    mu = snakemake.config['mu']
    N_e = snakemake.config['N_e'][sample_set]
except NameError:
    # testing
    vcf = "output/tetraploid/snps/pendula/synonymous/snps.vcf.gz"
    """pops = "output/tetraploid/sample_sets/subpopulations/pendula/2_pops.txt"
    scenario_name = 'split_symmetric_migration'
    n_pops = 2"""
    scenario_name = 'exp_pop_growth'
    sample_set = "pendula"
    pops = f"output/tetraploid/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    results = [f"output/tetraploid/dadi/{scenario_name}/{sample_set}/synonymous/bootstraps/{n}/dadi.csv" for n in range(1, 11)]
    result = f"output/tetraploid/dadi/{scenario_name}/pendula/synonymous/dadi.csv"
    n_pops = 1
    params_out = "scratch/dadi.csv"
    params_pretty_out = "scratch/dadi.out"

    generation_time = 20
    mu = 1e-9
    N_e = 675000

# load boostrap results
bootstraps = []
for file in results:
    bootstraps.append(pd.read_csv(file, index_col=None, header=0).head(1))

bootstraps = pd.concat(bootstraps, axis=0, ignore_index=True)

# load params from file
params = pd.read_csv(result, index_col=None, header=0).iloc[0]

# determine initial value
param_names = demographic_scenarios.get_class(scenario_name).params
p0 = [params[name] for name in param_names]

scenario = demographic_scenarios.restore(scenario_name, p0, vcf, pops,
                                         generation_time, n_pops, N_e, mu)

# calculate standard deviation from GIM
scenario.perform_uncertainty_analysis()

# get dataframe of params
data = scenario.to_dataframe()

# confidence interval threshold
a = 0.05

# calculate additional statistics based on the bootstraps
cols = [
    'cis_bs_percentile',  # percentile bootstrap confidence interval
    'cis_bs_bca',  # bca bootstrap confidence interval
    'std_bs',  # bca bootstrap standard deviation
    'mean_bs',  # bca bootstrap mean
]

col_data = {key: [] for key in cols}
for name in ['lnl', 'theta'] + scenario.params:
    values = bootstraps[name].values.tolist()

    bs_percentile = utils.get_ci_percentile_bootstrap(values, a)
    bs_bca = utils.get_ci_bca(values, params[name], a)

    col_data['cis_bs_percentile'].append(bs_percentile)
    col_data['cis_bs_bca'].append(bs_bca)
    col_data['std_bs'].append(np.std(values, ddof=1))
    col_data['mean_bs'].append(np.mean(values))

# add to data frame
for col in cols:
    data.loc[col] = col_data[col]

# save to file
open(params_out, 'w').write(data.to_csv())
open(params_pretty_out, 'w').write(data.to_markdown())
