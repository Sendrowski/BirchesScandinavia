import matplotlib.pyplot as plt
import numpy as np
import utils

try:
    vcf_synonymous = snakemake.input.synonymous
    vcf_nonsynonymous = snakemake.input.nonsynonymous
    out = snakemake.output[0]
except NameError:
    # testing
    vcf = 'testing/snps1.ancestral.vep.annotated.vcf.gz'
    out = 'scratch/synonymy.svg'

# the labels
labels_degeneracy = [f"{n}-fold" for n in [2, 4] + [0, 2]]
labels_synonymy = ['synonymous', 'non-synonymous']

# the degeneracy values for each synonymy class
synonymous = [utils.get_n_degenerate(vcf_synonymous, n) for n in [2, 4]]
nonsynonymous = [utils.get_n_degenerate(vcf_nonsynonymous, n) for n in [0, 2]]

# the nested array containing the values
values = [synonymous, nonsynonymous]

# prepare colors for inner and outer disk
cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(2) * 5)
inner_colors = cmap([1, 2, 6, 7])

fig, ax = plt.subplots()


# get wedge labels of outer disk
def autopict(percentage):
    n_sites = int(np.round(percentage * sum(utils.flatten(values)) / 100))

    return str(n_sites) + '\n' + str(np.round(percentage, 1)) + '%'


# plot outer disk showing the synonymy classes
ax.pie([sum(s) for s in values], startangle=90, labels=labels_synonymy, autopct=autopict,
       wedgeprops=dict(edgecolor='w', linewidth=1, antialiased=True, width=0.4),
       colors=outer_colors, radius=1, pctdistance=0.8)

# plot inner disk showing the degeneracy classes
ax.pie(utils.flatten(values), startangle=90, labels=labels_degeneracy,
       wedgeprops=dict(edgecolor='w', linewidth=1, antialiased=True, width=0.2),
       colors=inner_colors, radius=0.6, labeldistance=0.25)

ax.axis('equal')
utils.save_fig(out, tight_layout=True)
