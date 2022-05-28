import matplotlib.pyplot as plt
import numpy as np

import polydfe_utils
import utils

try:
    input = snakemake.input
    out = snakemake.output[0]
except NameError:
    # testing
    input = [f"output/default/polydfe/pendula_south/output/bootstraps/dfe_C.full_anc.{str(n).zfill(3)}.txt"
             for n in range(1, 101)]
    out = "scratch/param_dist.svg"

# number of intervals
n = len(polydfe_utils.parse_output(input[0])['discretized'])

fig, axes = plt.subplots(n)

# parse bootstrapped discretized DFE values
values = []
for bs in input:
    values.append(polydfe_utils.parse_output(bs)['discretized'])

# determine mean and 95% CIs
a = 0.05
mean = np.mean(values, axis=0)
values = np.array(values)

for i, ax in enumerate(axes):
    ax.hist(values[:, i], bins=50)

    ci_bca = utils.get_ci_bca(values[:, i], mean[i], a)
    ci_percentile = utils.get_ci_percentile_bootstrap(values[:, i], a)

    for x in ci_bca:
        ax.axvline(x=x, c='red', linewidth=1)

    for x in ci_percentile:
        ax.axvline(x=x, c='black', linewidth=1)

plt.tight_layout(pad=0)
fig.savefig(out)
