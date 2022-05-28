from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

in_bam = snakemake.input.bam
java_opts = get_java_opts(snakemake)
out_bam = snakemake.output.bam
out_metric = snakemake.output.metrics

shell(f"gatk MarkDuplicates --java-options '{java_opts}' -I {in_bam} -O {out_bam} -M {out_metric}")





