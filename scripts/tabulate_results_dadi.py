import math

import numpy as np
import pandas
import pandas as pd

try:
    test_mode = False
    results = snakemake.input
    labels = snakemake.params.labels
    out_pretty = snakemake.output.pretty
    out_tex = snakemake.output.tex
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
    results = {scenario: f"output/tetraploid/dadi/{scenario}/{sample_set}/synonymous/dadi_ci.csv" for scenario in scenarios}
    labels = [key.replace('_', ' ') for key in results.keys()]
    out_pretty = 'scratch/dadi_tabulated.out'
    out_tex = 'scratch/dadi_tabulated.tex'


# format given float string
def format_number(x):
    x = float(x)

    if x == 0:
        return str(0)

    magnitude = math.floor(math.log(abs(float(x)), 10))

    if magnitude > 0:
        return str(int(x))
    else:
        return str(np.round(x, abs(magnitude) + 1))


def parse_row(file, delimiter):
    data = pd.read_csv(file, index_col=0, header=0)

    mean = data.loc['mean_bs']
    std = data.loc['std_bs']

    for key in mean.index:
        if mean[key]:
            mean[key] = format_number(mean[key])

            if float(std[key]) != 0:
                mean[key] += delimiter + format_number(std[key])

    return mean


def prepare_table(delimiter):
    data = pandas.DataFrame()
    for i, (key, file) in enumerate(results.items()):
        row = parse_row(file, delimiter)
        row.name = labels[i]
        data = data.append(row)

    # remove unnamed columns
    data = data.loc[:, ~data.columns.str.contains('unnamed')]
    data = data.loc[:, ~(data.columns == 'theta')]
    return data.fillna('')


open(out_pretty, 'w').write(prepare_table(" ± ").to_markdown())
open(out_tex, 'w').write(prepare_table("\,$\pm$\,").to_latex(escape=False))
