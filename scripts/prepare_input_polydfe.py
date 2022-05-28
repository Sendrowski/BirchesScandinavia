import os

import sfs_utils
import utils

try:
    synonymous = snakemake.input.synonymous
    nonsynonymous = snakemake.input.nonsynonymous
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    n_pops = snakemake.params.n_pops
    out = snakemake.output[0]
except NameError:
    # testing
    synonymous = "output/default/snps/pendula/synonymous/snps.vcf.gz"
    nonsynonymous = "output/default/snps/pendula/nonsynonymous/snps.vcf.gz"
    vcf = "output/default/snps/pendula/all/snps.vcf.gz"
    pops = "output/default/sample_sets/subpopulations/pendula/1_pops.txt"
    n_pops = 1
    out = "scratch/pendula.1_pops.txt"

# number of two fold degenerate sites
n_2fold = utils.get_n_degenerate(vcf, 2)

# polyDFE fails to converge for large sample sizes
# n = 20 was used in their own simulations
n_proj = 20
freqs_neut = sfs_utils.get_spectrum(synonymous, pops, n_proj)
freqs_neut = freqs_neut[~freqs_neut.mask].tolist()
n_sites_neut = utils.get_n_degenerate(vcf, 4) + 1 / 3 * n_2fold

freqs_sel = sfs_utils.get_spectrum(nonsynonymous, pops, n_proj)
freqs_sel = freqs_sel[~freqs_sel.mask].tolist()
n_sites_sel = utils.get_n_degenerate(vcf, 0) + 2 / 3 * n_2fold


def write_line(file, arr, sep):
    file.write(sep.join(map(str, arr)) + os.linesep)


m_neut = 1
m_sel = 1
with open(out, 'w') as f:
    write_line(f, [m_neut, m_sel, n_proj], ' ')

    write_line(f, freqs_neut + [n_sites_neut], '\t')
    write_line(f, freqs_sel + [n_sites_sel], '\t')
