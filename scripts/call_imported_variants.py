"""
Perform joint genotyping on genomic database generated with HaplotypeCaller.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

ref = snakemake.input.ref
db = snakemake.input.db
intervals = snakemake.input.intervals
java_opts = get_java_opts(snakemake)
out = snakemake.output[0]

shell(f"gatk GenotypeGVCFs --java-options '{java_opts}' -V gendb://{db} \
    -R {ref} -O {out} --include-non-variant-sites -L {intervals}")






