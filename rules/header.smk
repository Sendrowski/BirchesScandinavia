import re

import pandas as pd
from snakemake.utils import min_version

min_version("6.0")

configfile: "config.yaml"

# update with testing environment if specified
if config.get('testing',False):
    print("Loading config for test environment.")

    configfile: "config.testing.yaml"

# whether we run snakemake locally
is_local = config.get('local',False)

# currently we are using macOS locally
is_macos = is_local

# whether we are using macOS
is_macos = config.get('macos',is_macos)

birch_samples = pd.read_csv(config['birch_samples'],index_col=0)
outgroup_samples = pd.read_csv(config['outgroup_samples'],index_col=0)

# merge birch and outgroup samples
samples = pd.concat([birch_samples, outgroup_samples])

sample_names = samples.name.values

# we replicate part of the pipeline for different sample classes
# this is useful for fundamental differences close to the resource files
sample_classes = [
    # here we treat all samples as diploid
    'default',

    # here we treat all pubescens sample as tetraploid 
    # and integrate them with the remaining diploid samples
    'tetraploid'
]

sample_sets = [
    'all',
    'birch',
    'left_cluster',
    'right_cluster',
    'pendula',
    'pubescens',
    'pendula_pubescens',
    'pendula_south',
    'pendula_north',
    'pubescens_south',
    'pubescens_north',
    'pendula_luis',
    'pendula_pubescens_admixture',
    'pendula_admixture',
    'pubescens_admixture',
    'pendula_south_admixture',
    'pendula_north_admixture',
    'pubescens_south_admixture',
    'pubescens_north_admixture'
]

# flags for pseudo sample sets
flags = [
    'all',
    'biallelic',
    'no_missing',
    'synonymous',
    'nonsynonymous',
    'no_low_freqs',
    'no_very_low_freqs',
    '0fold',
    '2fold',
    '4fold',
    'nonsynonymous_0fold',
    'synonymous_4fold'
]

n_folds_degeneracy = [0, 2, 4]

# immediate subpopulations for 2n
subpopulations_2n = {
    'birch': ['right_cluster', 'left_cluster'],
    'pendula_pubescens': ['pendula', 'pubescens'],
    'pendula': ['pendula_north', 'pendula_south'],
    'pubescens': ['pubescens_north', 'pubescens_south'],
    'pendula_pubescens_admixture': ['pendula_admixture', 'pubescens_admixture'],
    'pendula_admixture': ['pendula_north_admixture', 'pendula_south_admixture'],
    'pubescens_admixture': ['pubescens_north_admixture', 'pubescens_south_admixture']
}

# immediate superpopulations
superpopulations = {
    'pendula': 'pendula_pubescens',
    'pubescens': 'pendula_pubescens',
    'pendula_south': 'pendula',
    'pendula_north': 'pendula',
    'pubescens_south': 'pubescens',
    'pubescens_north': 'pubescens',
    'pendula_admixture': 'pendula_pubescens_admixture',
    'pubescens_admixture': 'pendula_pubescens_admixture',
    'pendula_south_admixture': 'pendula_admixture',
    'pendula_north_admixture': 'pendula_admixture',
    'pubescens_south_admixture': 'pubescens_admixture',
    'pubescens_north_admixture': 'pubescens_admixture',
}


# Get the subpopulations of the given sample set
# where the number of subpopulations is specified by n.
# Simply return the sample set if n equals 1.
def get_subpopulations(sample_set, n=2):
    n = int(n)

    if n == 1:
        return [sample_set]

    if n == 2 and sample_set in subpopulations_2n:
        return subpopulations_2n[sample_set]

    # if n = 4 we try to nest two calls with n = 2
    if n == 4:
        sets = []
        for s in get_subpopulations(sample_set,2):
            sets += get_subpopulations(s,2)

        return sets

    return []


# Get the immediate superpopulation of given sample set
def get_superpopulation(sample_set):
    if sample_set in superpopulations:
        return superpopulations[sample_set]


# get subpopulations for sample sets for which we mark
# the subpopulations in PCA and UMAP plots by ellipses
def get_marked_subpopulations(sample_set):
    if sample_set in ['pendula_pubescens', 'birch']:
        return get_subpopulations('pendula_pubescens')

    return []


# remove metacharacter used for optional wildcards
def strip_wildcard(name):
    return re.sub(r"..","",name)


# compile regex to only match with the given values
def make_wildcard(values):
    return "|".join(map(str,values))


# compile regex to optionally match with the given values
def make_optional_wildcard(values):
    return "|\.".join([".{0}"] + list(map(str,values)))


wildcard_constraints:
    name="\w+",
    path=".+",
    component="\w+",
    n="\d+",
    optional_flags='.*',
    sample_class=make_wildcard(sample_classes),
    sample_set=make_wildcard(sample_sets),
    flag=make_wildcard(flags),
    n_fold_deg=make_wildcard(n_folds_degeneracy),
    scale=make_wildcard(['log', 'linear']),
    bs_type=make_wildcard(['percentile', 'bca']),
    lamb="\d*\.?\d+([eE][-+]?\d+)?"

# rules to execute locally
localrules:
    convert_svg_to_png

# convert SVG files to PNG format
rule convert_svg_to_png:
    input:
        "results/{sample_class}/graphs/svg/{path}.svg"
    output:
        "results/{sample_class}/graphs/png/{path}.png"
    conda:
        "../envs/cairosvg.yaml"
    shell:
        "cairosvg {input} -d 200 -o {output}"

# create a tbi index for a vcf file using the GATK's IndexFeatureFile
rule create_tbi_index_vcf:
    input:
        "{path}.vcf.gz"
    output:
        "{path}.vcf.gz.tbi"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk IndexFeatureFile -I {input}"

# file with ingroup sample information
rule resource_ingroup_samples:
    output:
        protected(config['birch_samples'])

# file with outgroup sample information
rule resource_outgroup_samples:
    output:
        protected(config['outgroup_samples'])
