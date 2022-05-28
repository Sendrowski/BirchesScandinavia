import pandas as pd

import demographic_scenarios

try:
    test_mode = False
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    results = snakemake.input.results
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
    test_mode = True
    vcf = "output/tetraploid/snps/pendula/synonymous/snps.vcf.gz"
    sample_set = 'pendula_pubescens'
    scenario_name = 'split_symmetric_migration'
    n_pops = 2
    """scenario_name = 'one_pop_size_change_since_ice_age'
    n_pops = 1"""
    results = [f"output/tetraploid/dadi/{scenario_name}/{sample_set}/synonymous/batches/{n}/dadi.csv"
               for n in range(1, 101)]
    pops = f"output/tetraploid/sample_sets/subpopulations/pendula/{n_pops}_pops.txt"
    params_out = "scratch/dadi.csv"
    params_pretty_out = "scratch/dadi.out"
    generation_time = 20
    mu = 1e-9
    N_e = 675000

res = []
for file in results:
    df = pd.read_csv(file, index_col=None, header=0)
    res.append(df)

res = pd.concat(res, axis=0, ignore_index=True)

# determine result with highest likelihood
best_params = res[res.lnl == res.lnl.max()].iloc[0]

# determine initial value
param_names = demographic_scenarios.get_class(scenario_name).params
p0 = [best_params[name] for name in param_names]

scenario = demographic_scenarios.restore(scenario_name, p0, vcf, pops,
                                         generation_time, n_pops, N_e, mu)

# calculate standard deviation from GIM
scenario.perform_uncertainty_analysis()

# save to file
data = scenario.to_dataframe()
open(params_pretty_out, 'w').write(data.to_markdown())
open(params_out, 'w').write(data.to_csv())
