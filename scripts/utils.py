"""
Main utilities.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import gzip
import os.path
import re
import warnings
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from pandas_plink import read_plink
from scipy.stats import chi2
from scipy.stats import norm
from sklearn.impute import SimpleImputer

pd.options.mode.chained_assignment = None

# load config file
if os.path.exists('config.yaml'):
    config = yaml.load(open('config.yaml'), Loader=yaml.FullLoader)

sample_sets_proper_names = {
    'pendula': 'B. pendula',
    'pubescens': 'B. pubescens',
    'pendula_pubescens': 'B. pendula & pubescens',
    'birch': 'Betula spp.',
    'pendula_south': 'B. pendula south',
    'pendula_north': 'B. pendula north',
    'pubescens_south': 'B. pubescens south',
    'pubescens_north': 'B. pubescens north',
    'pendula_admixture': 'B. pendula',
    'pubescens_admixture': 'B. pubescens',
    'pendula_south_admixture': 'B. pendula south',
    'pendula_north_admixture': 'B. pendula north',
    'pubescens_south_admixture': 'B. pubescens south',
    'pubescens_north_admixture': 'B. pubescens north'
}

pop_id_names = {
    'pendula_pubescens': ['pubescens', 'pendula'],
    'pendula': ['pendula_north', 'pendula_south'],
    'pubescens': ['pubescens_north', 'pubescens_south']
}


# get proper name for the given sample set
def get_proper_name(sample_set):
    if sample_set in sample_sets_proper_names:
        return sample_sets_proper_names[sample_set]

    return sample_set.replace('_', ' ')


# We used generic names for the population ids in dadi and DILS.
# This function returns the proper names based on the
# assignment in the script 'create_subpopulation_files'.
# The returned names are in the order ['pop0', 'pop1']
def get_name_for_pop_ids(sample_set):
    if sample_set in pop_id_names:
        return pop_id_names[sample_set]


# get sample list
def get_samples(sample_set='all', sample_class='default'):
    birch_samples = pd.read_csv(config['birch_samples'], index_col=0)
    outgroup_samples = pd.read_csv(config['outgroup_samples'], index_col=0)

    # concatenate outgroups with all birch samples
    samples = pd.concat([birch_samples, outgroup_samples])

    path = config['sample_sets'].format(sample_set=sample_set, sample_class=sample_class)
    subsamples = pd.read_csv(path, header=None, index_col=False)[0]

    # reduce sample set
    samples = samples[samples.name.isin(subsamples)]
    samples = samples.reset_index(drop=True)

    return samples


# get sample list from file
def get_samples_from_file(path):
    return pd.read_csv(path, header=None, index_col=False)[0].tolist()


# get the genotypes of specified sample set in PLINK format
def get_genotypes(sample_set, flag, sample_class='default'):
    path = config['plink_data'].format(sample_set=sample_set, flag=flag, sample_class=sample_class)
    (bim, fam, genotypes) = read_plink(path)

    samples = get_samples(sample_set)

    # make sure the samples have the same order
    assert np.all(samples.name.values == fam.iid.values), "Samples are in different order!"

    # Silence FutureWarning: The `numpy.may_share_memory` function is not implemented by Dask array.
    # You may want to use the da.map_blocks function or something similar to silence this warning.
    # Your code may stop working in a future release
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # impute missing values
        genotypes = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(genotypes)

    # transpose
    genotypes = np.array(genotypes).T

    return genotypes, samples, bim, fam


# return the samples of 'ref' that are not in 'comp'
def not_in_sample_set(ref, comp):
    sample_set = get_samples(ref)

    return sample_set[~sample_set.name.isin(get_samples(comp).name)]


# remove monomorphic sites from plink genotypes
def remove_monomorphic(genotypes):
    sums = sum(genotypes)
    return genotypes[:, np.logical_and(sums != len(genotypes[:, 0]) * 2, sums != 0)]

# remove missing sites from plink genotypes
def remove_missing(genotypes):
    return genotypes[:, ~np.isnan(genotypes).any(axis=0)]


# load given GFF file into a Dataframe
def load_gff(file):
    cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attribute']
    return pd.read_csv(file, sep='\t', header=None, comment='#', names=cols)


# open potentially gzipped file for reading
def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, "rt")

    return open(file, 'r')


# count number of lines the given pattern was matched in file
def count_matches(file, patterns: List):
    patterns = [re.compile(p) for p in patterns]

    i = 0
    with open_file(file) as f:
        for line in f:
            matches = [bool(pattern.search(line)) for pattern in patterns]
            # add 1 if all patterns were matched
            if sum(matches) == len(patterns):
                i += 1

    return i


# get the number of sites in given vcf file
def count_lines_vcf(file):
    i = 0
    with open_file(file) as f:
        for line in f:
            if not line.startswith('#'):
                i += 1

    return i


# get the number of n-fold degenerate lines
def get_n_degenerate(vcf, degeneracy):
    return count_matches(vcf, [f'Degeneracy={degeneracy}'])


# get the number of n-fold synonymous lines
def get_n_synonymous(vcf, synonymy):
    return count_matches(vcf, [f'Synonymy={synonymy}'])


# get the number of sites that are (non)-synonymous with given degeneracy
def get_n_degenerate_synonymous(vcf, degeneracy, synonymy):
    return count_matches(vcf, [f'Degeneracy={degeneracy}', f'Synonymy={synonymy}'])


# flatten the given 2D array
def flatten(values):
    return [item for sub in values for item in sub]


def get_bounds(data, a1, a2, n):
    if np.isnan(a1) or np.isnan(a2):
        return [None, None]

    return [data[max(int(a1 * n), 0)], data[min(int(a2 * n), n - 1)]]


# get the (1 - a)% confidence intervals using the percentile bootstrap
def get_ci_percentile_bootstrap(bootstraps: List, a: float):
    n = len(bootstraps)
    data = np.sort(bootstraps)

    return get_bounds(data, a, 1 - a, n)


# get the (1 - a)% confidence intervals using the BCa method
# cf. An Introduction to the Bootstrap, Bradley Efron, Robert J. Tibshirani, section 14.2
def get_ci_bca(bootstraps: List, original: float, a: float):
    if sum(bootstraps) == 0:
        return [0, 0]

    n = len(bootstraps)
    data = np.sort(np.array(bootstraps))

    theta_hat_i = np.array([np.var(np.delete(data, i)) for i in range(n)])
    theta_hat = sum(theta_hat_i) / n

    a_hat = sum((theta_hat - theta_hat_i) ** 3) / (6 * (sum((theta_hat - theta_hat_i) ** 2)) ** (3 / 2))

    # we add epsilon here to avoid getting -inf when sum(data < original) / n is 0
    z0_hat = norm.ppf(sum(data < original) / n + np.finfo(float).eps)
    z_a = norm.ppf(a)

    a1 = norm.cdf(z0_hat + (z0_hat + z_a) / (1 - a_hat * (z0_hat + z_a)))
    a2 = norm.cdf(z0_hat + (z0_hat - z_a) / (1 - a_hat * (z0_hat - z_a)))

    return get_bounds(data, a1, a2, n)


# save current plot
def save_fig(out, tight_layout=False, show=False, pad=0, clear=False):
    if tight_layout:
        plt.savefig(out, bbox_inches='tight', pad_inches=pad)
    else:
        plt.savefig(out)

    if show:
        plt.show()

    if clear:
        plt.clf()


# scale the current figure by the given multiplier
def scale_cf(m):
    plt.gcf().set_size_inches(np.array(m) * np.array(plt.rcParams["figure.figsize"]))


# perform an LRT
def lrt(lnl_simple, lnl_complex, df=1):
    lr = -2 * (lnl_simple - lnl_complex)

    return chi2.sf(lr, df)
