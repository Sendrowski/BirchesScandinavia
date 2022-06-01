"""
Plot sample locations colored by specified subpopulations.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import matplotlib.pyplot as plt
import numpy as np

import utils

try:
    test_mode = False
    sample_set = snakemake.params.sample_set
    subsample_sets = snakemake.params.subsample_sets
    sample_class = snakemake.params['sample_class']
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    subsample_sets = ['pendula_south_admixture', 'pendula_north_admixture']
    sample_class = 'default'
    out = 'scratch/admixture_locations.png'

# fetch all samples
samples = utils.get_samples(sample_set=sample_set, sample_class=sample_class)

# fetch subsamples and choose colors according to subsample set
# only two subsample sets are currently supported
subsamples = utils.get_samples(subsample_sets[0], sample_class=sample_class)
mask = samples.name.isin(subsamples.name)
colors = [int(m) for m in mask]

# perturb the sample locations to make individual samples visible
np.random.seed(seed=0)
samples.latitude += np.random.normal(size=samples.shape[0], scale=0.3)
samples.longitude += np.random.normal(size=samples.shape[0], scale=0.3)

# plot samples
scatter = plt.scatter(samples.longitude, samples.latitude, c=colors, s=3.5)
plt.gca().axis('square')
plt.xlabel('longitude')
plt.ylabel('latitude')

# add legend
labels = [utils.get_proper_name(s) for s in subsample_sets]
plt.legend(handles=scatter.legend_elements()[0], labels=labels)

utils.scale_cf(2 / 3)

utils.save_fig(out, show=test_mode, tight_layout=True)
