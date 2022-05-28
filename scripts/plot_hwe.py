import plot_utils

try:
    test_mode = False
    input = snakemake.input[0]
    out = snakemake.output[0]
    log_scale = snakemake.params.log_scale
except NameError:
    # testing
    test_mode = True
    input = 'output/default/stats/pendula/biallelic/hwe.hwe'
    out = 'scratch/hwe.svg'
    log_scale = False

# filtering dependent on allele frequency to check for bias
'''
def filter_ac(data):
    data['AC'] = data.GENO.str.split('/').apply(lambda x: float(x[1])) / \
                 data.GENO.str.split('/').apply(lambda x: float(x[2]))

    return data[~data.AC.between(0.001, 0.999)]
'''

plot_utils.plot_hist_basic_stat(input, colname='P', xlabel="p-value", delim_whitespace=True,
                                ylabel="n SNPs", color=plot_utils.get_color('Paired', 0),
                                log_scale=log_scale)

plot_utils.save_fig(out, tight_layout=True, show=test_mode)
