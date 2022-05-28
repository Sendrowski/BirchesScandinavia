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
    out = snakemake.output[0]

    plot_opts = {}
    generation_time = snakemake.config['generation_time']
    mu = snakemake.config['mu']
    N_e = snakemake.config['N_e'][sample_set]
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    vcf = f"output/tetraploid/snps/{sample_set}/synonymous/snps.vcf.gz"
    pops = f"output/tetraploid/sample_sets/subpopulations/{sample_set}/1_pops.txt"
    scenario_name = 'one_pop_size_change'
    result = f"output/tetraploid/dadi/{scenario_name}/{sample_set}/synonymous/dadi.csv"
    n_pops = 1
    out = "scratch/trajectory.png"

    plot_opts = {'scaling': np.array([1, 1 / 2]) * 0.8}
    generation_time = 20
    mu = 1e-9
    N_e = 675000

# load params from file
params = pd.read_csv(result, index_col=None, header=0).iloc[0]

# determine initial value
param_names = demographic_scenarios.get_class(scenario_name).params
p0 = [params[name] for name in param_names]

scenario = demographic_scenarios.restore(scenario_name, p0, vcf, pops,
                                         generation_time, n_pops, N_e, mu)

# plot trajectory
scenario.plot_trajectory(opts=plot_opts)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
