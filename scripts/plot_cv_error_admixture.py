"""
Plot the cross-validation error for ADMIXTURE.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import plot_utils

try:
    test_mode = False
    log_files = snakemake.input
    K = snakemake.params.K
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    log_files = [f"output/default/admixture/pendula/biallelic/snps.admixture.{n}.log" for n in range(1, 5)]
    K = range(1, 5)
    out = 'scratch/cv_error_admixture.svg'

cv_errors = {}
for K, file in zip(K, log_files):
    # get file contents
    contents = open(file, 'r').read()

    # match cv error
    match = re.search("CV error \(K=\d\): (\d+\.\d+)", contents)

    cv_errors[K] = float(match.group(1))

# plot cv errors
plt.plot(cv_errors.keys(), cv_errors.values(), marker='o')
plt.xlabel('K')
plt.ylabel('CV error')

# mark lambda with lowest cv error
best_K = list(cv_errors.keys())[list(cv_errors.values()).index(min(cv_errors.values()))]
plt.axvline(best_K, color="orange")

# only show integers on x-axis
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

# increase scaling
plot_utils.scale_cf(2 / 3)

# save plot
plot_utils.save_fig(out, show=test_mode, tight_layout=True)
