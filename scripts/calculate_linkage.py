"""
Calculate linkage of sites in given PLINK files.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

file = snakemake.params.input_prefix
out = snakemake.params.output_prefix

# execute command
# LD window options obtained from Luis Leal
shell(f"plink --bfile {file} --r2 --ld-window-kb 50 --ld-window 500 "
      f"--ld-window-r2 0 --allow-extra-chr --out {out}")
