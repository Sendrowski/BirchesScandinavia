from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

intervals = snakemake.input.intervals
variants = ' '.join([f"-V {vcf}" for vcf in snakemake.input.vcfs])
java_opts = get_java_opts(snakemake)
db_dir = snakemake.output[0]
tmp_dir = snakemake.resources.tmpdir
batch_size = snakemake.params.batch_size

shell(f"gatk --java-options '{java_opts}' GenomicsDBImport {variants} \
    --genomicsdb-workspace-path {db_dir} --tmp-dir {tmp_dir} \
    --interval-padding 100 -L {intervals} --batch-size {batch_size}")
