"""
Create a pie chart for the sites' degeneracy classes.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import matplotlib.pyplot as plt

import numpy as np
import utils
import warnings

try:
    sites_all = snakemake.input.sites_all
    vcf_synonymous = snakemake.input.synonymous
    vcf_nonsynonymous = snakemake.input.nonsynonymous
    out = snakemake.output[0]
except NameError:
    # testing
    sites_all = "output/default/snps/pendula/snps.vcf.gz"
    vcf_synonymous = "output/default/snps/pendula_synonymous/snps.vcf.gz"
    vcf_nonsynonymous = "output/default/snps/pendula_nonsynonymous/snps.vcf.gz"
    out = 'scratch/degeneracy.svg'

# the degeneracy values
ns = [0, 2, 4]

labels_degeneracy = ['non-coding'] + [f"{n}-fold" for n in ns]
labels_synonymy = ['_', 'non-synonymous', '_', 'non-synonymous', 'synonymous', '_', 'synonymous', '_']

n_lines_total = utils.count_lines_vcf(sites_all)
degeneracy_all = {n: utils.get_n_degenerate(sites_all, n) for n in ns}
n_non_coding = n_lines_total - sum(degeneracy_all.values())

# the degeneracy classes with respect to synonymy
synonymous = {n: utils.get_n_degenerate(vcf_synonymous, n) for n in [2, 4]}
nonsynonymous = {n: utils.get_n_degenerate(vcf_nonsynonymous, n) for n in [0, 2]}

# the degeneracy classes with synonymy subclasses
degeneracy_synonymy = {
    0: [nonsynonymous[0], degeneracy_all[0] - nonsynonymous[0]],
    2: [nonsynonymous[2], synonymous[2], degeneracy_all[2] - nonsynonymous[2] - synonymous[2]],
    4: [synonymous[4], degeneracy_all[4] - synonymous[4]]
}

# all values in a nested array
degeneracy_synonymy = [[n_non_coding]] + list(degeneracy_synonymy.values())
degeneracy = [sum(s) for s in degeneracy_synonymy]

# prepare colors of disks
cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(4) * 4)
inner_colors = cmap([1, 5, 6, 9, 10, 11, 13, 14])

fig, ax = plt.subplots()


# get wedge labels of outer disk
def autopict(percentage):
    n_sites = int(np.round(percentage * n_lines_total / 100))

    return str(n_sites) + '\n' + str(np.round(percentage, 1)) + '%'


widths = [0.5, 0.1]
# plot outer disk showing the degeneracy classes
ax.pie(degeneracy, labels=labels_degeneracy, autopct=autopict, startangle=90,
       wedgeprops=dict(edgecolor='w', linewidth=1, antialiased=True, width=widths[0]),
       pctdistance=0.75, colors=outer_colors, radius=1)

# plot inner disk showing the synonymy classes
wedges_synonymy, _ = ax.pie(utils.flatten(degeneracy_synonymy), startangle=90,
                            wedgeprops=dict(edgecolor='w', linewidth=1, antialiased=True, width=widths[1]),
                            colors=inner_colors, radius=1 - widths[0])

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.legend(wedges_synonymy, labels_synonymy, loc="upper left")

ax.axis('equal')
utils.save_fig(out, tight_layout=True)
