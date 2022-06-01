"""
Annotate the VCF with regards to the ancestral alleles.
An AA info tag will added in addition to further information.
This is same tag dadi looks for when parsing a unfolded SFS.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from collections import Counter

import vcf
from vcf.parser import _Info as Info

import pyvcf_utils
import utils

try:
    vcf_file = snakemake.input.vcf
    probs = snakemake.input.probs
    data = snakemake.input.data
    samples_file = snakemake.input.samples
    out = snakemake.output.vcf
except NameError:
    # testing
    vcf_file = "output/default/snps/raw/intervals/snps1.vcf.gz"
    probs = "output/default/est-sfs/probs/1.txt"
    data = "output/default/est-sfs/data/1.txt"
    samples_file = "output/default/sample_sets/ingroups.args"
    out = "scratch/1.ancestral.vcf"

samples = utils.get_samples_from_file(samples_file)

vcf_reader = vcf.Reader(filename=vcf_file)

# add AA info field to header
vcf_reader.infos['AA'] = Info('AA', 1, 'String', 'Ancestral Allele', None, None)
vcf_reader.infos['EST_SFS_probs'] = Info('EST_SFS_probs', 1, 'String', 'EST-SFS probabilities', None, None)
vcf_reader.infos['EST_SFS_input'] = Info('EST_SFS_input', 1, 'String', 'EST-SFS input', None, None)

probs_reader = open(probs, 'r')
data_reader = open(data, 'r')
writer = vcf.Writer(open(out, 'w'), vcf_reader)

# write to data file
i = 0
for record in vcf_reader:

    # simply assign the ancestral allele to be the reference allele
    # if the record is monomorphic
    if record.is_monomorphic:
        record.INFO['AA'] = record.REF
        record.INFO['EST_SFS_probs'] = None
        record.INFO['EST_SFS_input'] = None
    else:

        line = probs_reader.readline()

        # read est-sfs input data
        est_input = data_reader.readline()

        # get probability of major allele
        prob_major_allele = float(line.split(' ')[2])

        # restrict haplotypes to the non-outgroup birch samples
        hps = pyvcf_utils.haplotypes(pyvcf_utils.restrict(record.samples, samples))

        major_allele = '.'

        if hps:
            # determine major allele
            major_allele = Counter(hps).most_common()[0][0]

            # take the major allele to be the ancestral allele
            # if its probability is greater than equal 0.5
            if prob_major_allele >= 0.5:
                ancestral_allele = major_allele
            else:
                # there are exactly two alleles
                for allele in record.alleles:
                    if str(allele) != major_allele:
                        ancestral_allele = str(allele)
        else:
            ancestral_allele = '.'

        # add ancestral allele annotation to record
        record.INFO['AA'] = ancestral_allele

        # add additional information from est-sfs
        record.INFO['EST_SFS_probs'] = line.replace('\n', '').replace(' ', '|')
        record.INFO['EST_SFS_input'] = est_input.replace('\n', '').replace(' ', '|').replace(',', ':')

        if not (major_allele == ancestral_allele == record.REF):
            print(f"site: {record.CHROM}:{record.POS}, ancestral allele: {ancestral_allele}, "
                  f"major allele: {major_allele}, reference: {record.REF}, "
                  f"prob major allele: {prob_major_allele}", flush=True)

    writer.write_record(record)

    i += 1
    if i % 1000 == 0: print(f"Processed sites: {i}", flush=True)

# raise error if there are lines left from the EST-SFS output
if next(probs_reader, None) is not None:
    raise AssertionError("Number of sites don't match.")

probs_reader.close()
