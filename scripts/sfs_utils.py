import dadi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utils


# load the given population file in dadi format
def load_samples_pop_file(file):
    return pd.read_csv(file, names=['name', 'pop'], header=None, sep='\t')


# obtain the unprojected spectrum from the given vcf file
def get_full_spectrum(vcf_file, pops):
    # create data dict from vcf file
    dd = dadi.Misc.make_data_dict_vcf(vcf_file, pops)

    # read sample file
    samples = load_samples_pop_file(pops)

    # determine pop ids
    pop_ids = sorted(samples['pop'].unique())

    # do not down-project
    # we assume all samples to be diploid here
    projections = [2 * sum(samples['pop'].values == pop_id) for pop_id in pop_ids]

    return dadi.Spectrum.from_data_dict(dd, pop_ids, projections=projections, polarized=True, mask_corners=False)


# load SFS from vcf file
def get_spectrum(vcf_file, pops, n_proj=20, polarized=True, mask_corners=True, n_pops=1):
    # dadi can handle polyploidy
    dd = dadi.Misc.make_data_dict_vcf(vcf_file, pops)

    return get_spectrum_from_dd(dd, n_proj=n_proj, polarized=polarized,
                                mask_corners=mask_corners, n_pops=n_pops)


# load the spectrum object from the given data dict
def get_spectrum_from_dd(dd, n_proj=20, polarized=True, mask_corners=True, n_pops=1):
    pop_ids = get_pop_ids(n_pops)

    return dadi.Spectrum.from_data_dict(dd, pop_ids, projections=[n_proj] * n_pops,
                                        polarized=polarized, mask_corners=mask_corners)


def get_pop_ids(n_pops):
    return ['pop' + str(i) for i in range(n_pops)]


# plot the sfs of the given vcf files
def plot_sfs_from_file(files, pops, n_proj, unfolded=True, log_scale=False, labels=[], sample_set=None):
    # double projection value if SFS is not specified as folded
    # like this we can better compare folded and unfolded spectra
    n_proj = n_proj if unfolded else n_proj * 2 - 1

    n_pops = len(load_samples_pop_file(pops)['pop'].unique())

    # load spectra from files
    spectra = []
    for file in files:
        spectra.append(get_spectrum(file, pops, n_proj, unfolded, n_pops=n_pops))

    plot_sfs(spectra, log_scale, labels, sample_set=sample_set)


# plot the sfs of the given spectra
def plot_sfs(spectra, log_scale=False, labels=[], sample_set=None):
    if spectra[0].ndim == 1:
        plot_sfs_1D(spectra, log_scale, labels)
    elif spectra[0].ndim == 2:
        plot_sfs_2D(spectra, sample_set)


# get proper sample set names for the pops ids used
def get_proper_names(sample_set):
    names = utils.get_name_for_pop_ids(sample_set)
    return [utils.get_proper_name(name) for name in names]


def plot_sfs_2D(spectra, sample_set=None):
    fig, axs = plt.subplots(ncols=len(spectra))

    # pyplot does not return an array if only one plot is requested
    if len(spectra) == 1:
        axs = [axs]

    pop_ids = get_proper_names(sample_set) if sample_set else None

    for i, ax in enumerate(axs):
        dadi.Plotting.plot_single_2d_sfs(spectra[i], vmin=1, ax=ax, pop_ids=pop_ids)


def plot_sfs_1D(spectra, log_scale=False, labels=[]):
    width_total = 0.9
    width = width_total / len(spectra)

    for i, fs in enumerate(spectra):
        heights = fs.compressed()
        n = len(heights)
        indices = np.arange(1, n + 1) + i * width

        plt.bar(indices, heights, width=width, label=labels[i] if len(labels) else None)

    # adjust ticks
    indices_ticks = np.arange(1, n)
    indices_ticks = indices_ticks[indices_ticks % int(n / 7) == 1]
    plt.xticks([i + (width_total - width) / 2 for i in indices_ticks], indices_ticks)

    plt.xlabel('Allele count')
    plt.autoscale(tight=True)

    if log_scale:
        plt.gca().set_yscale('log')

    if len(labels):
        plt.legend()
