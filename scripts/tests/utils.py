import scripts.utils as utils

import numpy as np


def test_get_n_degenerate():
    file = 'testing/snps1.ancestral.vep.annotated.vcf.gz'
    np.testing.assert_equal(utils.get_n_degenerate(file, 4), 2312)
    np.testing.assert_equal(utils.get_n_degenerate(file, 2), 3604)
    np.testing.assert_equal(utils.get_n_degenerate(file, 0), 10705)


def test_get_n_synonymous():
    file = 'testing/snps1.ancestral.vep.annotated.vcf.gz'
    np.testing.assert_equal(utils.get_n_synonymous(file, 0), 1343)
    np.testing.assert_equal(utils.get_n_synonymous(file, 1), 919)

def test_get_ci_percentile_bootstrap():
    a = 0.05
    data = np.arange(100)
    np.random.shuffle(data)
    ci = utils.get_ci_percentile_bootstrap(data.T, a)
    np.testing.assert_equal(ci, [5, 95])

def test_get_ci_bca():
    a = 0.05
    data = np.arange(100)
    np.random.shuffle(data)
    mean = np.mean(data, axis=0)
    ci = utils.get_ci_bca(data.T, mean, a)
    np.testing.assert_equal(ci, [5, 95])
