"""
Call the haplotypes util GATK's HaplotypeCaller in GVCF mode.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

ref = snakemake.input.ref
bam = snakemake.input.bam
intervals = snakemake.input.intervals
java_opts = get_java_opts(snakemake)
tmp_dir = snakemake.resources.tmpdir
ploidy = snakemake.params.ploidy
out = snakemake.output[0]

shell(f"gatk --java-options '{java_opts}' HaplotypeCaller "
      f"-R {ref} -I {bam} -O {out} --tmp-dir {tmp_dir} "
      f"--interval-padding 0 -L {intervals} -ERC GVCF "
      f"--sample-ploidy {ploidy}")
