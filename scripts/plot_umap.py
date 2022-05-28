import matplotlib.pyplot as plt
import umap

import pca_utils
import utils

try:
    sample_set = snakemake.params['sample_set']
    sample_class = snakemake.params['sample_class']
    flag = snakemake.params['flag']
    out = snakemake.output[0]
    subsample_sets = snakemake.params.get('subsample_sets', [])
    min_dist = snakemake.params.get('min_dist', 1.5)
    spread = snakemake.params.get('spread', 2)
    tight_layout = snakemake.params.get('tight_layout', False)
except NameError:
    # testing
    sample_set = 'birch'
    sample_class = 'default'
    flag = 'biallelic'
    out = "scratch/pca.svg"
    subsample_sets = ['pubescens', 'pendula']
    min_dist = 1.5
    spread = 2
    tight_layout = False

# fetch genotypes
genotypes, samples, bim, fam = utils.get_genotypes(sample_set, flag, sample_class)

# create umap object
reducer = umap.UMAP(n_neighbors=samples.shape[0] - 1, metric='euclidean',
                    min_dist=min_dist, random_state=4, spread=spread)

# fit embedding using umap
embedding = reducer.fit_transform(genotypes)

# plot and save embedding
plt.scatter(x=embedding[:, 0], y=embedding[:, 1], c=samples.latitude, s=20)
plt.gca().set_aspect('equal', 'datalim')

# draw ellipses around subsample sets if specified
for i, subsample_set in enumerate(subsample_sets):
    pca_utils.encircle_subset(sample_set, subsample_set, sample_class, embedding, i, buf=10)

plt.xlabel(f'x')
plt.ylabel(f'y')

if len(subsample_sets):
    plt.legend()

plt.colorbar().set_label('latitude')

if tight_layout:
    plt.savefig(out, bbox_inches='tight', pad_inches=0)
else:
    plt.savefig(out)
