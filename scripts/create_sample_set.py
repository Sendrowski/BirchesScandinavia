import pandas as pd

import pca_utils
import utils

# The derivation of some sample sets depends on the previous deviation of
# other, more general sample sets. This is because the utils module will
# look for previously generated files matching the sample set name.
# The order of derivation is thus important, which has to be from less
# specific to more specific.

try:
    input = snakemake.input[0]
    sample_set = snakemake.params.sample_set
    sample_class = snakemake.params['sample_class']

    pendula_latitudinal_barrier = snakemake.config['latitudinal_barriers']['pendula']
    pubescens_latitudinal_barrier = snakemake.config['latitudinal_barriers']['pubescens']

    out = snakemake.output[0]
except NameError:
    # testing
    input = 'output/default/admixture/pubescens/biallelic/snps.admixture.2.Q'
    sample_set = 'pubescens_north_admixture'
    sample_class = 'default'

    pendula_latitudinal_barrier = 64
    pubescens_latitudinal_barrier = 65

    out = f'scratch/{sample_set}.args'


# get the samples to be included depending on ADMIXTURE's
# Q matrix denoting the provenance of an individual
def get_samples_admixture(condition, supersample_set):
    Q = pd.read_csv(input, index_col=None, header=None, sep=' ')

    mask = condition(Q.iloc[:, 0]).values

    return utils.get_samples(supersample_set, sample_class)[mask]


# based on PCA
pendula_pubescens_barrier = 20

birch_outliers = ['GNA02', 'GNA04', 'BUR27a', 'BUR27b']

birch_samples = pd.read_csv(utils.config['birch_samples'], index_col=0)

if sample_set == 'all':

    outgroup_samples = pd.read_csv(utils.config['outgroup_samples'], index_col=0)

    # concatenate outgroup and birch samples
    samples = pd.concat([birch_samples, outgroup_samples])

elif sample_set == 'birch':

    samples = birch_samples

elif sample_set == 'left_cluster':

    pca, pc, _ = pca_utils.get('birch', 'biallelic', sample_class)

    samples = birch_samples[pc[:, 0] < pendula_pubescens_barrier]

elif sample_set == 'right_cluster':

    pca, pc, _ = pca_utils.get('birch', 'biallelic', sample_class)

    samples = birch_samples[pc[:, 0] >= pendula_pubescens_barrier]

elif sample_set == 'pendula_pubescens':

    # remove outliers
    samples = birch_samples[~birch_samples.name.isin(birch_outliers)]

elif sample_set == 'pendula':

    samples = utils.get_samples('right_cluster', sample_class)

    samples = samples[~samples.name.isin(birch_outliers)]

elif sample_set == 'pendula_south':

    samples = utils.get_samples('pendula', sample_class)

    samples = samples[samples.latitude <= pendula_latitudinal_barrier]

elif sample_set == 'pendula_north':

    samples = utils.get_samples('pendula', sample_class)

    samples = samples[samples.latitude > pendula_latitudinal_barrier]

elif sample_set == 'pubescens':

    samples = utils.get_samples('left_cluster', sample_class)

    # remove outliers
    samples = samples[~samples.name.isin(birch_outliers)]

elif sample_set == 'pubescens_south':

    samples = utils.get_samples('pubescens', sample_class)

    samples = samples[samples.latitude <= pubescens_latitudinal_barrier]

elif sample_set == 'pubescens_north':

    samples = utils.get_samples('pubescens', sample_class)

    samples = samples[samples.latitude > pubescens_latitudinal_barrier]

elif sample_set == 'outgroups':

    samples = pd.read_csv(utils.config['outgroup_samples'], index_col=0)

elif sample_set == 'pendula_north_admixture':

    samples = get_samples_admixture(lambda Q: Q > 0.5, 'pendula')

elif sample_set == 'pendula_south_admixture':

    samples = get_samples_admixture(lambda Q: Q <= 0.5, 'pendula')

elif sample_set == 'pubescens_north_admixture':

    samples = get_samples_admixture(lambda Q: Q > 0.5, 'pubescens')

elif sample_set == 'pubescens_south_admixture':

    samples = get_samples_admixture(lambda Q: Q <= 0.5, 'pubescens')

elif sample_set == 'pubescens_admixture':

    samples = get_samples_admixture(lambda Q: Q > 0.5, 'pendula_pubescens')

elif sample_set == 'pendula_admixture':

    samples = get_samples_admixture(lambda Q: Q <= 0.5, 'pendula_pubescens')

# write samples to file
samples.name.to_csv(out, index=False, header=False)
