from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

java_opts = get_java_opts(snakemake)
vcf = snakemake.input.vcf
filter = snakemake.params.filter
out = snakemake.output[0]

shell(f"gatk --java-options '{java_opts}' SelectVariants -V {vcf} -O {out} "
      f"--selectExpressions '{filter}'")
