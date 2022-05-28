import os
import utils

sample_set = snakemake.params.sample_set
sample_class = snakemake.params.sample_class
n_pops = snakemake.params.n_pops
out = snakemake.output[0]

samples = utils.get_samples(sample_set, sample_class)
samples['pop'] = 'pop0'


def assign_pop1(subset_name, samples):
    subset = utils.get_samples(subset_name, sample_class)

    samples.loc[samples.name.isin(subset.name), 'pop'] = 'pop1'

# subpopulations:
# pubescens: pop0, pendula: pop1
# northern: pop0, southern: pop1
if n_pops == 2:
    if sample_set == 'pendula_pubescens':
        assign_pop1('pendula', samples)
    elif sample_set == 'pendula':
        assign_pop1('pendula_south', samples)
    elif sample_set == 'pubescens':
        assign_pop1('pubescens_south', samples)

# write to file
with open(out, 'w') as f:
    f.write(os.linesep.join(samples[['name', 'pop']].agg('\t'.join, axis=1)) + os.linesep)
