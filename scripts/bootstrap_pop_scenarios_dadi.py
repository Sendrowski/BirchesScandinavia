"""
Create bootstraps for a certain demographic scenario in dadi.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import dadi
import pandas as pd

import demographic_scenarios
import sfs_utils

try:
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    result = snakemake.input.result
    scenario_name = snakemake.params.scenario
    sample_set = snakemake.params.sample_set
    n_pops = snakemake.params.n_pops
    params_out = snakemake.output.params
    params_pretty_out = snakemake.output.params_pretty
    mode = snakemake.params.mode

    n_iter = snakemake.config['dadi_n_iter_' + mode]

    generation_time = snakemake.config['generation_time']
    mu = snakemake.config['mu']
    N_e = snakemake.config['N_e'][sample_set]
except NameError:
    # testing
    vcf = "output/default/snps/pendula/synonymous/snps.vcf.gz"
    """pops = "output/default/sample_sets/subpopulations/pendula/2_pops.txt"
    scenario_name = 'split_symmetric_migration'
    n_pops = 2"""
    pops = "output/default/sample_sets/subpopulations/pendula/1_pops.txt"
    result = "output/default/dadi/one_pop_size_change_since_ice_age/pendula/synonymous/batches/1/dadi.csv"
    scenario_name = 'one_pop_size_change_since_ice_age'
    n_pops = 1
    sample_set = "pendula"
    params_out = "scratch/dadi.csv"
    params_pretty_out = "scratch/dadi.out"
    mode = 'local_optimization'

    n_iter = 10

    generation_time = 20
    mu = 1e-9
    N_e = 675000

n_proj = 20
dd = dadi.Misc.make_data_dict_vcf(vcf, pops)

# divide data into chunks of size chunk_size_bs
# we then subsample from these chunks with replacement
# to create the bootstraps
# https://dadi.readthedocs.io/en/latest/user-guide/bootstrapping/
chunks = dadi.Misc.fragment_data_dict(dd, demographic_scenarios.DemographicScenario.chunk_size_bs)

# obtain one bootstrap sample
fs_bs = dadi.Misc.bootstraps_from_dd_chunks(chunks, 1, sfs_utils.get_pop_ids(n_pops), [n_proj] * n_pops)[0]

# we pass the original data dict here as it is only
# used to create parametric bootstraps
params = dict(fs=fs_bs, dd=dd, n_pops=n_pops, n_proj=n_proj, n_iter=n_iter,
              N_e=N_e, mu=mu, generation_time=generation_time, mode=mode)

scenario = demographic_scenarios.create(scenario_name, params)

# load params from simulation result
params = pd.read_csv(result, index_col=None, header=0).iloc[0]

# determine initial value
scenario.p0 = [params[name] for name in scenario.params]

# repeat the simulation with the appropriate initial value for one
# iteration for convenience reasons
scenario.simulate()

# save to file
data = scenario.to_dataframe()
open(params_pretty_out, 'w').write(data.to_markdown())
open(params_out, 'w').write(data.to_csv())
