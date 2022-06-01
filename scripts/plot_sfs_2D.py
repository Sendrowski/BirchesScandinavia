"""
Plot 2D SFS.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import sfs_utils
import plot_utils

try:
    test_mode = False
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    n_proj = snakemake.params.n_proj
    sample_set = snakemake.params.get('sample_set', None)
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    sample_set = 'pendula'
    vcf = f"output/default/snps/{sample_set}/4fold/snps.vcf.gz"
    pops = f"output/default/sample_sets/subpopulations/{sample_set}/2_pops.txt"
    n_proj = 20
    out = "scratch/sfs_{sample_set}.svg"

sfs_utils.plot_sfs_from_file([vcf], pops, n_proj, unfolded=True, sample_set=sample_set)

plot_utils.scale_cf(0.65)
plot_utils.save_fig(out, tight_layout=True, show=test_mode)
