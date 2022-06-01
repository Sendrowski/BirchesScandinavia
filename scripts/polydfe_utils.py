"""
Utilities for polyDFE.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import pandas as pd


# parse the discretized DFE value from the given file
def parse_output(file):
    with open(file) as f:
        lines = f.readlines()

        return {
            'discretized': list(map(float, lines[1].split(' '))),
            'gradient': float(lines[2])
        }


# load several polydfe files into a data frame
def parse_output_files(files):
    data = pd.DataFrame(columns=list(range(6)) + ['gradient'])

    for f in files:
        output = parse_output(f)
        col = {i: v for i, v in enumerate(output['discretized'])} | {'gradient': output['gradient']}
        data = data.append(col, ignore_index=True)

    return data


# get a list of interval names used for the x-axis
def get_interval_names(file):
    with open(file) as f:
        intervals = f.read().split('\n')[0].split(' ')

    intervals = ['-$\infty$'] + intervals + ['$\infty$']

    names = [f"({intervals[i]}, {intervals[i + 1]}]" for i in range(len(intervals) - 1)]
    names[-1] = names[-1].replace(']', ')')

    return names


labels = {
    'full_anc': 'full DFE + ancestral misidentification',
    'full': 'full DFE',
    'deleterious_anc': 'deleterious DFE + ancestral misidentification',
    'deleterious': 'deleterious DFE'
}

short_labels = {
    'full_anc': 'full anc',
    'full': 'full',
    'deleterious_anc': 'del anc',
    'deleterious': 'del'
}


# extract plot label from file name
def get_label(name):
    return labels[name]


def get_short_label(name):
    return short_labels[name]
