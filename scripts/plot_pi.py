import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/biallelic/pi.sites.pi'
    out = 'scratch/pi.svg'
    log_scale = True

plot_utils.plot_hist_basic_stat(input, colname='PI', xlabel="$\pi$", ylabel="n SNPs",
                                color=plot_utils.get_color('tab20b', 14), log_scale=log_scale)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
