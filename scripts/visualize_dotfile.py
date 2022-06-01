"""
Visualize a dotfile.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

input = snakemake.input[0]
out = snakemake.output[0]

shell(f"dot {input} -Tsvg > {out}")
