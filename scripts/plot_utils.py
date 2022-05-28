import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

import utils

save_fig = utils.save_fig
scale_cf = utils.scale_cf


def get_color(cm, index):
    return plt.cm.get_cmap(cm).colors[index]


def plot_lrt_p_values(complex, simple, p_values, types, labels, transpose=False):
    n = len(types)
    complex, simple = np.array(complex), np.array(simple)
    p_values = np.array(p_values)

    labels_x, labels_y = np.array(labels), np.array(labels)

    # create matrix
    m = np.full((n, n), -1, np.float)
    for i in range(n):
        for j in range(n):
            # check if p-value is present for combination and determine position
            if sum(mask := (complex == types[i]) & (simple == types[j])):
                m[i, j] = p_values[mask][0]

    # remove empty columns
    nonempty_cols = np.sum(m != -1, axis=0) != 0
    m = m[:, nonempty_cols]
    labels_x = labels_x[nonempty_cols]

    # remove empty rows
    nonempty_rows = np.sum(m != -1, axis=1) != 0
    m = m[nonempty_rows]
    labels_y = labels_y[nonempty_rows]

    def format_number(x):
        if x == 0:
            return 0

        if x < 0.0001:
            return "{:.1e}".format(x)

        return np.round(x, 4)

    # determine values to display
    annot = np.vectorize(lambda x: str(format_number(x)))(m)
    annot[m == -1] = '-'

    # change to 1 to get a nicer color
    m[m == -1] = 1

    if transpose:
        m = m.T
        annot = annot.T
        labels_x, labels_y = labels_y, labels_x

    fig, ax = plt.subplots()

    # make the cbar have the same height as the heatmap
    cbar_ax = make_axes_locatable(ax).new_horizontal(size="4%", pad=0.15)
    fig.add_axes(cbar_ax)

    # plot heatmap
    sns.heatmap(m, ax=ax, cbar_ax=cbar_ax, cmap="inferno", annot=annot, fmt="",
                square=True, cbar_kws={'label': 'p-value'}, vmin=0, vmax=1,
                linewidths=0.5, linecolor='#cccccc')

    ax.set_xticklabels(labels_x, rotation=45)
    ax.set_yticklabels(labels_y, rotation=0)


# plot a histogram plot of a basic statistic calculated with VCFtools
def plot_hist_basic_stat(file, colname, xlabel, ylabel, color, log_scale=False, weights=None,
                         sep='\t', delim_whitespace=False, bins=50, opts={}, data_filter=None):
    # load data from csv files and remove empty rows
    if delim_whitespace:
        data = pd.read_csv(file, delim_whitespace=True)
    else:
        data = pd.read_csv(file, sep=sep)

    data = data[data[colname].notnull()]

    # apply filter to data if specified
    if data_filter:
        data = data_filter(data)

    # whether to weight the dat~ using specified column name
    if weights:
        weights = data[weights]

    # plot histogram
    opts = opts | dict(grid=False, bins=bins, color=color, weights=weights,
                       edgecolor=sns.set_hls_values(color, s=0.2))

    data[colname].hist(**opts)

    # add vertical line to indicate mean
    mean = np.average(data[colname], weights=weights)
    plt.axvline(mean, color=color, linestyle='dashed', linewidth=1, label='$\\bar{x}=' + str(np.round(mean, 4)) + '$')

    # plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if log_scale:
        plt.gca().set_yscale('log')
        plt.ylim(bottom=1)

    plt.legend()
    plt.tight_layout()

    return plt
