"""
Plot 1D SFS.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import sfs_utils
import plot_utils

try:
    test_mode = False
    files = snakemake.input.vcf
    pops = snakemake.input.pops
    n_proj = snakemake.params.n_proj
    unfolded = snakemake.params.unfolded
    log_scale = snakemake.params.log_scale
    labels = snakemake.params.get('labels', [])
    tight_layout = snakemake.params.get('tight_layout', False)
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    files = ["output/tetraploid/snps/pendula/biallelic/snps.vcf.gz"]
    pops = "output/tetraploid/sample_sets/subpopulations/pendula/1_pops.txt"
    n_proj = 20
    unfolded = True
    log_scale = False
    labels = ['biallelic']
    tight_layout = False
    out = "scratch/pendula.pendula_one_pop.svg"

sfs_utils.plot_sfs_from_file(files, pops, n_proj, unfolded, log_scale, labels)

plot_utils.save_fig(out, tight_layout=tight_layout, show=test_mode)
