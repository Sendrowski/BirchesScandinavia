"""
Plot a frequency distribution site-wise heterozygosity.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    observed = snakemake.params.observed
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/biallelic/hwe.hwe'
    out = 'scratch/hwe.svg'
    observed = False
    log_scale = True

colname = 'O(HET)' if observed else 'E(HET)'

plot_utils.plot_hist_basic_stat(input, colname=colname, xlabel="heterozygosity", delim_whitespace=True,
                                ylabel="n SNPs", color=plot_utils.get_color('tab20b', 18), log_scale=log_scale)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
