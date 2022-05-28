import numpy as np
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
    out_sfs = snakemake.output.sfs
    out_residuals = snakemake.output.residuals
    scaling = [1, 0.5]

    generation_time = snakemake.config['generation_time']
    mu = snakemake.config['mu']
    N_e = snakemake.config['N_e'][sample_set]
except NameError:
    # testing
    test_mode = True
    n_pops = 1
    vcf = "output/tetraploid/snps/pendula/synonymous/snps.vcf.gz"
    pops = f"output/tetraploid/sample_sets/subpopulations/pendula/{n_pops}_pops.txt"
    scenario_name = 'one_pop_size_change_since_ice_age'
    result = f"output/tetraploid/dadi/{scenario_name}/pendula/synonymous/dadi.csv"
    out_sfs = "scratch/sfs.png"
    out_residuals = "scratch/residuals.png"
    scaling = np.array([1, 1.15]) * 0.6

    generation_time = 20
    mu = 1e-9
    N_e = 945883

# load params from file
params = pd.read_csv(result, index_col=None, header=0).iloc[0]

# determine initial value
param_names = demographic_scenarios.get_class(scenario_name).params
p0 = [params[name] for name in param_names]

scenario = demographic_scenarios.restore(scenario_name, p0, vcf, pops,
                                         generation_time, n_pops, N_e, mu)

# generate SFS plot
scenario.plot_sfs(scaling_1d=scaling)
plot_utils.save_fig(out_sfs, tight_layout=True, show=test_mode, clear=True)

# generate plot of residuals
scenario.plot_residuals_1d()
plot_utils.save_fig(out_residuals, tight_layout=True, show=test_mode)
