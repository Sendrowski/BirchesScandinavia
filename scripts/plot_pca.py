import matplotlib.cm as cm
import matplotlib.pyplot as plt

import pca_utils
import utils

try:
    test_mode = False
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']
    flag = snakemake.params['flag']
    out = snakemake.output[0]
    add_names = snakemake.params.get('add_names')
    subsample_sets = snakemake.params.get('subsample_sets', [])
    demarcation_type = snakemake.params.get('demarcation_type', 'ellipse')
    scaling = snakemake.params.get('scaling', 1)
    marker_size = snakemake.params.get('marker_size', 20)
    tight_layout = snakemake.params.get('tight_layout', False)
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    sample_class = 'default'
    flag = 'no_missing'
    out = "scratch/all_sites_no_missing.png"
    add_names = False
    subsample_sets = ['pendula_north', 'pendula_south']
    demarcation_type = 'gradient'
    scaling = 1
    marker_size = 20
    tight_layout = False

pca, pc, samples = pca_utils.get(sample_set, flag, sample_class)

# determine color of dots
if demarcation_type == 'color':
    # we only support two colors here
    subsamples = utils.get_samples(subsample_sets[0], sample_class=sample_class)
    mask = samples.name.isin(subsamples.name)
    colors = [int(m) for m in mask]
else:
    colors = samples.latitude

# plot the 2 components
scatter = plt.scatter(x=pc[:, 0], y=pc[:, 1], c=colors, s=marker_size)

# draw ellipses around subsample sets if specified
if demarcation_type == 'ellipse':
    for i, subsample_set in enumerate(subsample_sets):
        pca_utils.encircle_subset(sample_set, subsample_set, sample_class, pc, i)

# add sample names if specified
if add_names:
    eps = 1
    min_lat = samples.latitude.min()
    max_lat = samples.latitude.max()

    for i in range(len(samples)):
        c = cm.viridis((samples.latitude[i] - min_lat) / (max_lat - min_lat))

        plt.text(x=pc[i, 0] + eps, y=pc[i, 1] + eps, s=samples.name[i], size=6, c=c)

plt.gca().axis('square')

# add axis labels
v1 = round(pca.explained_variance_ratio_[0] * 100, 2)
v2 = round(pca.explained_variance_ratio_[1] * 100, 2)
plt.xlabel(f'{v1}%')
plt.ylabel(f'{v2}%')

# only show legend and colorbar if we don't
# color the samples according to their specified subsets
if demarcation_type == 'color':
    labels = [utils.get_proper_name(s) for s in subsample_sets]
    plt.legend(handles=scatter.legend_elements()[0], labels=labels)
else:
    plt.colorbar().set_label('latitude')

    if len(subsample_sets):
        plt.legend()

# increase scaling
utils.scale_cf(scaling)

utils.save_fig(out, tight_layout=tight_layout, show=test_mode)
