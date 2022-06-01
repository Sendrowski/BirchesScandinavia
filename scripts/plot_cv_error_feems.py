"""
Plot the cross-validation error for FEEMS.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import feems_utils
import utils

try:
    test_mode = False
    files = snakemake.input
    lambdas = snakemake.params.lambdas
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    lambdas = [0.01, 1]
    files = [f"output/default/feems/pendula/biallelic/cv_error_{n}.txt" for n in lambdas]

    out = 'scratch/cv_error_admixture.svg'

cv_errors = []
for file in files:
    # get file contents
    cv_errors.append(float(open(file, 'r').read()))

feems_utils.plot_cv_errors(lambdas, cv_errors)

# save plot
utils.save_fig(out, show=test_mode, tight_layout=True)
