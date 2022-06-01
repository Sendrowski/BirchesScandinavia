"""
Visualize the complete rulegraph using snakemake.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

snakefile = snakemake.params.snakefile
out = snakemake.output[0]

shell(f"snakemake --snakefile {snakefile} --rulegraph | dot -Tsvg > {out}")
