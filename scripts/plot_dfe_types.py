"""
Compare the DFE of different types using histogram plots.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import matplotlib.pyplot as plt
import numpy as np
import polydfe_utils
import plot_utils
from matplotlib import cycler

plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('Set2').colors)

try:
    test_mode = False
    files = snakemake.input
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    files = {'full_anc': "output/default/polydfe/pendula/output/dfe_C.full_anc.txt",
             'full': "output/default/polydfe/pendula/output/dfe_C.full.txt",
             'deleterious': "output/default/polydfe/pendula/output/dfe_C.deleterious.txt"}
    out = "scratch/dfe.svg"

width_total = 0.9
width = width_total / len(files)

# iterate over files and a bars
for i, (name, file) in enumerate(files.items()):
    values = polydfe_utils.parse_output(file)['discretized']
    indices = np.arange(len(values)) + i * width

    label = polydfe_utils.get_label(name)

    plt.bar(indices, values, width=width, label=label)

# adjust x-axis
plt.xlabel('$4N_e s$')
plt.xticks([r + (width_total - width) / 2 for r in range(len(values))], polydfe_utils.get_interval_names(file))

plt.autoscale(tight=True)
plt.legend()

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
