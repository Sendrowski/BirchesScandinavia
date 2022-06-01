"""
Plot SFS obtained from EST-SFS.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import numpy as np
from matplotlib import pyplot as plt
import sfs_utils
import dadi

try:
    input = snakemake.input[0]
    n_proj = snakemake.params.n_proj
    unfolded = snakemake.params.unfolded
    log_scale = snakemake.params.log_scale
    out = snakemake.output[0]
except NameError:
    # testing
    input = "output_old/est-sfs/sfs/1.txt"
    n_proj = 100
    unfolded = False
    log_scale = True
    out = "scratch/pendula.pendula_one_pop.svg"

fs = dadi.Spectrum([int(np.round(float(n))) for n in open(input, 'r').readline().split(',')])

if not unfolded:
    fs = fs.fold()
    n_proj *= 2

fs = fs.project([n_proj])

sfs_utils.plot_sfs([fs], log_scale=log_scale)
plt.savefig(out)
