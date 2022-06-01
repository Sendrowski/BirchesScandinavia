"""
Visualize LRT comparisons of different DFE models.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cycler

import plot_utils
import polydfe_utils

plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('Set2').colors)

try:
    testing = False
    input = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    input = "output/tetraploid/polydfe/pendula/type_comparison.C.txt"
    out = "scratch/probs_nested.svg"

# read data containing p-values
# columns: 'type1': complex model name, 'type2': simple model name, 'p.value': LRT p-value
data = pd.read_csv(input, sep=' ')

types = ['full', 'deleterious', 'full_anc', 'deleterious_anc']
labels = [polydfe_utils.get_short_label(t) for t in types]

plot_utils.plot_lrt_p_values(data.type1, data.type2, data['p.value'], types, labels)

plt.tight_layout()
plot_utils.scale_cf(0.75)
plot_utils.save_fig(out, show=testing, tight_layout=True)
