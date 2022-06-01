"""
Restrict given VCF to the specified samples.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

vcf = snakemake.input.vcf
samples = snakemake.input.samples
max_nocall = snakemake.params.max_nocall_fraction
java_opts = get_java_opts(snakemake)
out = snakemake.output[0]
biallelic = snakemake.params.get('biallelic')

flag_biallelic = '--select-type-to-include SNP --restrict-alleles-to BIALLELIC'
optional_flags = flag_biallelic if biallelic else ''

shell(f"gatk --java-options '{java_opts}' SelectVariants -V {vcf} "
      f"-O {out} --sample-name {samples} --max-nocall-fraction {max_nocall} "
      f"--remove-unused-alternates {optional_flags}")
