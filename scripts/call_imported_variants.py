from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

ref = snakemake.input.ref
db = snakemake.input.db
intervals = snakemake.input.intervals
java_opts = get_java_opts(snakemake)
out = snakemake.output[0]

shell(f"gatk GenotypeGVCFs --java-options '{java_opts}' -V gendb://{db} \
    -R {ref} -O {out} --include-non-variant-sites -L {intervals}")






