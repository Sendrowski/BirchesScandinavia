import matplotlib.pyplot as plt
import pandas as pd

import plot_utils
import utils

try:
    test_mode = False
    results = snakemake.input
    df = snakemake.params.get('df', 1)
    comparisons = snakemake.params['comparisons']
    labels = snakemake.params.get('labels', list(results.keys()))
    scaling = snakemake.params.get('scaling', 1)
    transpose = snakemake.params.get('transpose', False)
    out = snakemake.output[0]
except NameError:
    # testing
    test_mode = True
    scenarios = [
        "split_no_migration",
        "split_unidirectional_migration_from_two_to_one",
        "split_unidirectional_migration_from_one_to_two",
        "split_symmetric_migration",
        "split_asymmetric_migration",
        "split_asymmetric_migration_one_pop_size_change",
        "split_no_migration_since_ice_age",
        "split_unidirectional_migration_from_two_to_one_since_ice_age",
        "split_unidirectional_migration_from_one_to_two_since_ice_age",
        "split_symmetric_migration_since_ice_age",
        "split_asymmetric_migration_since_ice_age",
        "split_asymmetric_migration_one_pop_size_change_since_ice_age"
    ]
    sample_set = "pendula_pubescens"
    results = {scenario: f"output/tetraploid/dadi/{scenario}/{sample_set}/synonymous/dadi.csv" for scenario in scenarios}
    # each element is a pair of models to compare
    # where the first element is the more simple model
    comparisons = [
        ['split_no_migration_since_ice_age', 'split_no_migration', 1],
        ['split_unidirectional_migration_from_two_to_one_since_ice_age', 'split_unidirectional_migration_from_two_to_one', 1],
        ['split_unidirectional_migration_from_one_to_two_since_ice_age', 'split_unidirectional_migration_from_one_to_two', 1],
        ['split_symmetric_migration_since_ice_age', 'split_symmetric_migration', 1],
        ['split_asymmetric_migration_since_ice_age', 'split_asymmetric_migration', 1],
        ['split_asymmetric_migration_one_pop_size_change_since_ice_age', 'split_asymmetric_migration_one_pop_size_change', 1]
    ]
    labels = [key.replace('_', ' ') for key in results.keys()]
    scaling = 1
    df = 1
    transpose = False
    out = 'scratch/dadi_nested.svg'


def parse_lnl(file):
    # take the last row which if the file contains bootstrap ci estimates
    # and take the first row otherwise
    n_row = -1 if '_ci' in file else 0

    return float(pd.read_csv(file, index_col=None, header=0).iloc[n_row]['lnl'])


simple, complex, p_values = [], [], []
for comp in comparisons:
    simple.append(comp[0])
    complex.append(comp[1])
    df = comp[2] if len(comp) == 3 else 1

    lnl_simple = parse_lnl(results[comp[0]])
    lnl_complex = parse_lnl(results[comp[1]])
    p = utils.lrt(lnl_simple, lnl_complex, df)

    p_values.append(p)

types = list(results.keys())
plot_utils.plot_lrt_p_values(simple, complex, p_values, types, labels, transpose=transpose)

plt.tight_layout()
plot_utils.scale_cf(scaling)
plot_utils.save_fig(out, show=test_mode, tight_layout=True)
