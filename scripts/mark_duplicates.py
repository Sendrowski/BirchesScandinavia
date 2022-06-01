"""
Mark duplicates using GATK.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

in_bam = snakemake.input.bam
java_opts = get_java_opts(snakemake)
out_bam = snakemake.output.bam
out_metric = snakemake.output.metrics

shell(f"gatk MarkDuplicates --java-options '{java_opts}' -I {in_bam} -O {out_bam} -M {out_metric}")





