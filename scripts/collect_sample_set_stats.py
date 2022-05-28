import pandas as pd

import utils

try:
    vcfs = snakemake.input.vcfs
    sample_sets = snakemake.params.sample_sets
    sample_class = snakemake.params['sample_class']
    out = snakemake.output[0]
except NameError:
    sample_sets = ['pendula', 'pubescens', 'pendula_pubescens']
    sample_class = 'default'
    vcfs = ['output/default/snps/pendula/all/snps.vcf.gz',
            'output/default/snps/pubescens/all/snps.vcf.gz',
            'output/default/snps/pendula_pubescens/all/snps.vcf.gz'
            ]
    out = 'scratch/sample_set_stats.csv'

n_samples, n_sites = [], []
for i, sample_set in enumerate(sample_sets):
    n_samples.append(len(utils.get_samples(sample_set, sample_class)))
    n_sites.append(utils.count_lines_vcf(vcfs[i]))

data = pd.DataFrame(columns=[s.replace('_', ' ') for s in sample_sets])
data.loc['n samples'] = n_samples
data.loc['n sites'] = n_sites

data.T.to_csv(out)
