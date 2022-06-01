"""
Calculate some basic summary statistics of the SFS computed from the specified VCF file.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import pandas as pd

import sfs_utils
import utils

try:
    vcfs = snakemake.input.vcfs
    pops_1 = snakemake.input.pops_1
    pops_2 = snakemake.input.pops_2
    names = snakemake.params.names
    out_txt = snakemake.output.txt
    out_csv = snakemake.output.csv
except NameError:
    # testing
    sample_sets = ["pubescens"]
    flag = "0fold"
    vcfs = [f"output/tetraploid/snps/{sample_set}/{flag}/snps.vcf.gz" for sample_set in sample_sets]
    # vcf_all = [f"output/tetraploid/snps/{sample_set}/all/snps.vcf.gz" for sample_set in sample_sets]
    pops_1 = [f"output/tetraploid/sample_sets/subpopulations/{sample_set}/1_pops.txt" for sample_set in sample_sets]
    pops_2 = [f"output/tetraploid/sample_sets/subpopulations/{sample_set}/2_pops.txt" for sample_set in sample_sets]
    names = sample_sets
    out_txt = "scratch/stats.txt"
    out_csv = "scratch/stats.csv"

stats = pd.DataFrame(index=['n sites', "Tajima's D", 'pi', 'pi per site', 'S', 'FST', 'FST_0'])

for i, name in enumerate(names):
    print(f'Calculating for subset: {name}', flush=True)

    # load spectrum for subset
    # down-projecting does not seem to affect any of
    # the statistics calculated however
    fs_1D = sfs_utils.get_full_spectrum(vcfs[i], pops_1[i])
    fs_2D = sfs_utils.get_full_spectrum(vcfs[i], pops_2[i])

    n_sites = utils.count_lines_vcf(vcfs[i])
    # n_sites_all = utils.count_lines_vcf(vcf_all[i])

    stats[name] = [
        n_sites,
        fs_1D.Tajima_D(),
        fs_1D.pi(),
        fs_1D.pi() / n_sites,
        fs_1D.S(),
        fs_2D.Fst(),
        fs_2D.scramble_pop_ids().Fst()
    ]

open(out_txt, 'w').write(stats.to_markdown())
stats.to_csv(out_csv)
