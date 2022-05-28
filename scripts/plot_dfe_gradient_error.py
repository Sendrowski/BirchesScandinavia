import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import polydfe_utils

try:
    bootstrapped = snakemake.input.bootstrapped
    out = snakemake.output[0]
except NameError:
    # testing
    bootstrapped = [f"output/default/polydfe/pendula/output/bootstraps/dfe_C.full_anc.{str(n).zfill(3)}.txt"
                    for n in range(1, 101)]
    out = "scratch/dfe.svg"

data = polydfe_utils.parse_output_files(bootstrapped)
data = data.sort_values(by=['gradient'])

names = polydfe_utils.get_interval_names(bootstrapped[0])
x = data.gradient.values

# prepare axes and colors
axs = plt.subplots(2, 3)[1].flatten()
colors = plt.get_cmap('Set2').colors

# iterate over intervals
for i, name in enumerate(names):
    y = data[i].values

    # plot points
    axs[i].scatter(x, y, c=[colors[i]] * len(x), s=10, alpha=0.5)

    # plot linear regression line
    m, b = np.polyfit(np.log(x), y, 1)
    axs[i].plot(x, m * np.log(x) + b, c='black')

    # plot mean per bin
    bins = stats.binned_statistic(np.log(x), y, 'mean', bins=10)
    x_bins = np.exp(bins.bin_edges[1:] + (bins.bin_edges[:-1] - bins.bin_edges[1:]) / 2)
    y_bins = bins.statistic
    axs[i].plot(x_bins, y_bins, c='red', alpha=0.5)

    axs[i].set_xlabel("log($\delta$)")
    axs[i].set_ylabel(name)
    axs[i].set_xscale('log')

plt.tight_layout(pad=1)
plt.savefig(out)
