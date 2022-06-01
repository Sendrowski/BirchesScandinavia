"""
Apply quality filters to SNPs using GATK.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

java_opts = get_java_opts(snakemake)
ref = snakemake.input.ref
vcf = snakemake.input.vcf
out = snakemake.output[0]

# the filters are based on the GATK recommendations for hard filtering
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
shell(f"gatk --java-options '{java_opts}' VariantFiltration -R {ref} -V {vcf} -O {out} \
    --filter-name 'filter1' \
    --filter-expression 'QD < 2.0' \
    --filter-name 'filter2' \
    --filter-expression 'SOR > 3.0' \
    --filter-name 'filter3' \
    --filter-expression 'MQ < 40.0' \
    --filter-name 'filter4' \
    --filter-expression 'MQRankSum < -12.5' \
    --filter-name 'filter5' \
    --filter-expression 'ReadPosRankSum < -8.0'")
