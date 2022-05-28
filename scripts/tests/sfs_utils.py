import matplotlib.pyplot as plt
import scipy.stats

import scripts.sfs_utils as sfs_utils


def test_plot_log_dist():
    y = scipy.stats.uniform.rvs(0.01, 100, size=10000)
    plt.hist(y, bins=100)
    plt.show()


def test_plot_spectrum_2D():
    vcf = "../../testing/snps/pendula/synonymous/snps.vcf.gz"
    pops = "../../testing/sample_sets/subpopulations/pendula/2_pops.txt"

    fs = sfs_utils.get_spectrum(vcf, n_proj=20, pops=pops, n_pops=2)

    sfs_utils.plot_sfs([fs], log_scale=False, labels=[])
    plt.show()
