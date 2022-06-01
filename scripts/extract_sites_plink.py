"""
Extract a specified set of sites using PLINK.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

try:
    ids = snakemake.input.ids
    input = snakemake.input.bed
    out = snakemake.output.bed
except NameError:
    # testing
    ids = "output/default/snps/pendula/biallelic/plink.prune.in"
    input = "output/default/snps/pendula/biallelic/snps.bed"
    out = "output/default/snps/pendula/biallelic/snps.admixture.bed"

input_prefix = input.replace('.bed', '')
out_prefix = out.replace('.bed', '')

shell(f"plink --bfile {input_prefix} --extract {ids} --make-bed --out {out_prefix} --allow-extra-chr '0'")
