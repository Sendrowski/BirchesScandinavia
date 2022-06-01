"""
Prepare the input for est-sfs.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import os

import utils
from pyvcf_utils import *

try:
    vcf_file = snakemake.input.vcf
    ingroups_file = snakemake.input.ingroups
    outgroups_file = snakemake.input.outgroups
    out_data = snakemake.output.data
    n_max_subsamples = snakemake.config['est_sfs_n_samples']
except NameError:
    # testing
    vcf_file = "testing/snps/pendula/synonymous/snps.vcf.gz"
    ingroups_file = "output/default/sample_sets/ingroups.args"
    outgroups_file = "output/default/sample_sets/outgroups.args"
    out_data = "scratch/pendula.synonymous.data.txt"
    n_max_subsamples = 50

# load ingroup and outgroup samples
ingroups = utils.get_samples_from_file(ingroups_file)
outgroups = utils.get_samples_from_file(outgroups_file)

# seed rng
np.random.seed(seed=0)

# number of subsamples to sample from haplotypes
n_subsamples = min(n_max_subsamples, len(ingroups))

# write to data file
with open(out_data, 'w') as f:
    i = 0
    for record in vcf.Reader(filename=vcf_file):

        # only do inference for polymorphic sites
        if not record.is_monomorphic:

            ingroup = count(subsample(haplotypes(restrict(record.samples, outgroups, True)), n_subsamples))
            outgroup1 = count(subsample(haplotypes(restrict(record.samples, [outgroups[0]])), 1))
            outgroup2 = count(subsample(haplotypes(restrict(record.samples, [outgroups[1]])), 1))

            f.write(' '.join([base_dict_to_string(r) for r in [ingroup, outgroup1, outgroup2]]) + os.linesep)

        i += 1
        if i % 1000 == 0: print(f"Processed sites: {i}", flush=True)
