"""
Run polyDFE.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

sfs = snakemake.input.sfs
init = snakemake.input.init
out = snakemake.output[0]
id = snakemake.params.id
m = snakemake.params.m

shell(f"polyDFE -d {sfs} -m {m} -i {init} {id} -v 100 > {out}")
#shell(f"/Users/janek/polydfe/bin/polydfe -d {sfs} -m {m} -i {init} {id} -v 100 > {out}")
