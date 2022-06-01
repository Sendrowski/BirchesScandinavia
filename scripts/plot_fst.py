"""
Plot a frequency distribution for site-wise Fst values.
The input file was generated with VCFtools.
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
    input = 'output/default/stats/pendula/biallelic/fst.weir.fst'
    out = 'scratch/fst.svg'
    log_scale = True

plot_utils.plot_hist_basic_stat(input, colname='WEIR_AND_COCKERHAM_FST', xlabel="$F_{st}$",
                                ylabel="n SNPs", color=plot_utils.get_color('tab20c', 1),
                                log_scale=log_scale)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
