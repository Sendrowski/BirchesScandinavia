"""
Plot a frequency distribution for both expected and observed site-wise heterozygosity.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/all/hwe.hwe'
    out = 'scratch/hwe.svg'
    log_scale = True

opts = dict(xlabel="heterozygosity", delim_whitespace=True,
            ylabel="n SNPs", log_scale=log_scale)

opts_observed = dict(colname='O(HET)', color=plot_utils.get_color('tab20b', 18),
                     bins=50, opts=dict(alpha=0.3, label='observed'))

plot_utils.plot_hist_basic_stat(input, **(opts_observed | opts))

opts_expected = dict(colname='E(HET)', color='cornflowerblue',
                     bins=25, opts=dict(alpha=0.3, label='expected'))

plot_utils.plot_hist_basic_stat(input, **(opts_expected | opts))

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
