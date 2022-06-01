"""
Test for Hardy-Weinberg equilibrium of sites in given PLINK files.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

file = snakemake.params.input_prefix
out = snakemake.params.output_prefix
midp = snakemake.params.midp

flag = 'midp' if midp else ''

# execute command
shell(f"plink --bfile {file} --hardy {flag} --hwe 0 --allow-extra-chr --out {out}")
