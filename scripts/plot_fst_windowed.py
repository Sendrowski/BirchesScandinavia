import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/biallelic/fst.windowed.weir.fst'
    out = 'scratch/fst_windowed.svg'
    log_scale = True

plot_utils.plot_hist_basic_stat(input, colname='WEIGHTED_FST', xlabel="$F_{st}$", ylabel="n windows",
                                color="turquoise", log_scale=log_scale, weights='N_VARIANTS')

plot_utils.save_fig(out, tight_layout=True, show=test_mode)