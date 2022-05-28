from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

vcfs = ' '.join([f"-I {vcf}" for vcf in snakemake.input.vcfs])
java_opts = get_java_opts(snakemake)
tmp_dir=snakemake.resources.tmpdir
out = snakemake.output[0]

shell(f"gatk --java-options '{java_opts}' GatherVcfs {vcfs} -O {out} \
    --TMP_DIR {tmp_dir} --REORDER_INPUT_BY_FIRST_VARIANT")
