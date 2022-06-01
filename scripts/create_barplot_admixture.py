"""
Create a barplot for ADMIXTURE output.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utils

try:
    test_mode = False
    Q = snakemake.input.Q
    sample_set = snakemake.params.sample_set
    flag = snakemake.params.flag
    sample_class = snakemake.params.sample_class
    subsample_sets = snakemake.params.subsample_sets
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula_pubescens'
    flag = 'biallelic'
    sample_class = 'default'
    subsample_sets = [
        'pendula_north_admixture',
        'pendula_south_admixture',
        'pubescens_north_admixture',
        'pubescens_south_admixture'
    ]
    Q = f"output/{sample_class}/admixture/{sample_set}/{flag}/snps.admixture.4.Q"
    out = 'scratch/admixture_locations.png'

# load data
data = pd.read_csv(Q, sep=' ', header=None)
K = data.columns.shape[0]

data[['name', 'latitude']] = utils.get_samples(sample_set, sample_class=sample_class)[['name', 'latitude']]
data['sample_set'] = ''

# assign subsample sets to samples
for subsample_set in subsample_sets:
    subsamples = utils.get_samples(subsample_set, sample_class=sample_class)
    data.sample_set[data.name.isin(subsamples.name)] = subsample_set

# sort individuals
# sort sample sets in the specified order
data['sample_set_order'] = [subsample_sets.index(s) for s in data.sample_set]
data = data.sort_values(['sample_set_order', 'latitude']).reset_index()

# plot data
ax = data[range(K)].plot(kind='bar', stacked=True, xticks=[], yticks=[0, 1], width=1, legend=None)

# define ticks
labels = [utils.get_proper_name(s) for s in subsample_sets]
tick_locs = [np.median(data[data.sample_set == s].index) for s in subsample_sets]
labels2 = [data[data.sample_set == s].sample_set.unique() for s in subsample_sets]
plt.xticks(ticks=tick_locs, labels=labels)

# plot vertical lines separating the subsample sets
plt.vlines([np.max(data[data.sample_set == s].index) for s in subsample_sets[:-1]], 0, 1, linestyles='dashed', color='black')

ax.set_aspect(data.shape[0] / 5)
utils.scale_cf(4 / 3)
ax.margins(0)

utils.save_fig(out, tight_layout=True, show=test_mode)
