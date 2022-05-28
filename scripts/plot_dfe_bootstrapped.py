import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cycler

import plot_utils
import polydfe_utils
import utils

plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('Set2').colors)

try:
    test_mode = False
    input = snakemake.input
    bs_type = snakemake.params.bs_type
    label_type = snakemake.params.label_type
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    input = {}
    for key in ['pendula_south', 'pendula_north', 'pubescens_south', 'pubescens_north']:
        files = [f"output/default/polydfe/{key}/output/bootstraps/dfe_C.full_anc.{str(n).zfill(3)}.txt"
                 for n in range(1, 101)]
        input[key] = files

    bs_type = 'percentile'
    label_type = 'sample_set'
    out = "scratch/dfe.svg"

width_total = 0.9
width = width_total / len(input.keys())

for i, (key, files) in enumerate(input.items()):
    # parse bootstrapped discretized DFE values
    dfe = []
    for bs in files:
        dfe.append(polydfe_utils.parse_output(bs)['discretized'])

    n = len(dfe[0])

    # determine mean and 95% CIs
    a = 0.05
    values = np.mean(dfe, axis=0)
    dfe = np.array(dfe)

    if bs_type == 'percentile':
        ci = [utils.get_ci_percentile_bootstrap(dfe[:, i], a) for i in range(n)]
    elif bs_type == 'bca':
        ci = [utils.get_ci_bca(dfe[:, i], values[i], a) for i in range(n)]

    yerr = np.array([[values[i] - ci[i][0], ci[i][1] - values[i]] for i in range(n)]).T

    indices = np.array(range(n)) + i * width

    if label_type == 'sample_set':
        label = utils.get_proper_name(key)
    elif label_type == 'polydfe_type':
        label = polydfe_utils.get_label(key)

    # plot bar
    plt.bar(indices, values, width=width, label=label, yerr=yerr, error_kw=dict(capsize=3))

# adjust x-axis
plt.xlabel('$4N_e s$')
plt.xticks([r + (width_total - width) / 2 for r in range(n)], polydfe_utils.get_interval_names(files[0]))

plt.autoscale(tight=True)

eps = 0.01
plt.ylim(min(values - yerr[0]), max(values + yerr[1]) + eps)

plt.legend()

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
