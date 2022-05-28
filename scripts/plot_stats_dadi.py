import matplotlib.pyplot as plt

import pandas as pd

try:
    stats = snakemake.input
    out = snakemake.output[0]
except NameError:
    # testing
    stats = {
        'pendula': 'output/default/stats/pendula/stats.csv',
        'pubescens': 'output/default/stats/pubescens/stats.csv',
        'pendula_pubescens': 'output/default/stats/pendula_pubescens/stats.csv'
    }

    out = 'scratch/stats.svg'

df = {}
for name, file in stats.items():
    data = pd.read_csv(file, header='infer', index_col=0).T
    data = data.filter(items=['all'], axis=0)
    data['name'] = name

    df[name] = data

df = pd.concat(df.values()).reset_index(drop=True)



plt.bar(df.index.values, df.Tajima_D.values, label='Tajima_D')
plot.bar(df.index.values, df.pi.values, label='$\pi$')

ax0.get_xaxis().set_visible(False)
plt.gca().set_yscale('log')

plt.autoscale(tight=True)
plt.savefig(out)
