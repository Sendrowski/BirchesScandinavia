"""
Sanitize given annotation file using AGAT.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

input = snakemake.input[0]
out = snakemake.output[0]

# cleanup file
shell(f"agat_convert_sp_gxf2gxf.pl -g {input} -o {out} --gvi 3")
