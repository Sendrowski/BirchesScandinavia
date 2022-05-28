import pandas as pd

try:
    test_mode = False
    results = snakemake.input
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    results = {
        "constant_pop_size_since_ice_age": "output/tetraploid/dadi/constant_pop_size_since_ice_age/pendula/synonymous/dadi_ci.csv",
        "exp_pop_growth_since_ice_age": "output/tetraploid/dadi/exp_pop_growth_since_ice_age/pendula/synonymous/dadi_ci.csv",
        "lin_pop_growth_since_ice_age": "output/tetraploid/dadi/lin_pop_growth_since_ice_age/pendula/synonymous/dadi_ci.csv",
        "lin_pop_growth_by_slope_since_ice_age": "output/tetraploid/dadi/lin_pop_growth_by_slope_since_ice_age/pendula/synonymous/dadi_ci.csv"
    }
    out = 'scratch/likelihoods_dadi.csv'

# parse files and fetch likelihood
data = pd.DataFrame(columns=['scenario', 'lnl'])
for i, (scenario, file) in enumerate(results.items()):
    lnl = pd.read_csv(file, index_col=None, header=0).iloc[0]['lnl']

    data.loc[i] = [scenario, lnl]

# write to file
open(out, 'w').write(data.to_csv())
