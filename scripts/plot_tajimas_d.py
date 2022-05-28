import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/all/D.Tajima.D'
    out = 'scratch/tajimas_d.svg'
    log_scale = False

plot_utils.plot_hist_basic_stat(input, colname='TajimaD', xlabel="Tajima's D", weights='N_SNPS',
                                ylabel="n windows", color=plot_utils.get_color('tab20b', 10),
                                log_scale=log_scale)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
