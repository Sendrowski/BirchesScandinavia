from collections import Counter

import numpy as np

import scripts.pyvcf_utils as pyvcf_utils
import scripts.utils as utils

sample_set = 'birch'
outgroups = utils.get_samples('outgroups')
samples = utils.get_samples(sample_set)
vcf = "testing/snps/snps1.ancestral.test.vcf.gz"
n_sites = utils.count_lines_vcf(vcf)


# check how many to derived alleles the ancestral state has been assigned
def test_count_minor_ancestral_alleles():
    def count(record, _):
        return record.INFO['AA'] != record.REF

    n = pyvcf_utils.count_sites(vcf, count)

    # output/snps/all/biallelic/snps.vcf.gz
    np.testing.assert_equal(n / n_sites, 0.0470540186733888)


# check how many major alleles don't equal the reference alleles
def test_count_minor_reference_alleles():
    def count(record, _):
        # restrict haplotypes to the non-outgroup birch samples
        hps = pyvcf_utils.haplotypes(pyvcf_utils.restrict(record.samples, samples.name.values))

        if hps:
            # determine major allele
            major_allele = Counter(hps).most_common()[0][0]

            if str(major_allele) != record.INFO['AA']:
                ingroup = pyvcf_utils.count(pyvcf_utils.haplotypes(
                    pyvcf_utils.restrict(record.samples, outgroups.name.values, True)))

                outgroup = pyvcf_utils.count(pyvcf_utils.haplotypes(
                    pyvcf_utils.restrict(record.samples, outgroups.name.values)))

                b = max(outgroup, key=outgroup.get)
                if b != record.INFO['AA'] or outgroup[b] == 0:
                    print(f"ingroup: {ingroup}, outgroup: {outgroup}, AA: {record.INFO['AA']}, "
                          f"major allele: {major_allele}, REF: {record.REF}")
                    return True

        return False

    n = pyvcf_utils.count_sites(vcf, count)

    # output/snps/all/biallelic/snps.vcf.gz
    np.testing.assert_equal(n / n_sites, 0.005577133212883123)
