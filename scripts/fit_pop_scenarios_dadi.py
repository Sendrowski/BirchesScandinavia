"""
Run MLE for specified demographic scenario using dadi.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import dadi

import demographic_scenarios
import sfs_utils

try:
    test_mode = False
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
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
    test_mode = True
    sample_set = "pendula"
    vcf = f"output/default/snps/{sample_set}/synonymous/snps.vcf.gz"
    n_pops = 2
    pops = f"output/default/sample_sets/subpopulations/{sample_set}/{n_pops}_pops.txt"
    scenario_name = 'split_unidirectional_migration_from_two_to_one_since_ice_age'
    params_out = "scratch/dadi.csv"
    params_pretty_out = "scratch/dadi.out"
    # mode = 'local_optimization'
    mode = 'random_hopping'

    n_iter = 1

    generation_time = 20
    mu = 1e-9
    N_e = 675000

n_proj = 20
dd = dadi.Misc.make_data_dict_vcf(vcf, pops)
fs = sfs_utils.get_spectrum_from_dd(dd, n_proj=n_proj, n_pops=n_pops)

params = dict(fs=fs, dd=dd, n_pops=n_pops, n_proj=n_proj, n_iter=n_iter,
              N_e=N_e, mu=mu, generation_time=generation_time, mode=mode)

scenario = demographic_scenarios.create(scenario_name, params)

# perform simulation
scenario.simulate()

# save to file
data = scenario.to_dataframe()
open(params_pretty_out, 'w').write(data.to_markdown())
open(params_out, 'w').write(data.to_csv())
