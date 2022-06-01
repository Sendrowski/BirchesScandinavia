"""
Gather multiple VCF files using Picard's GatherVcfs.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

vcfs = ' '.join([f"-I {vcf}" for vcf in snakemake.input.vcfs])
java_opts = get_java_opts(snakemake)
tmp_dir=snakemake.resources.tmpdir
out = snakemake.output[0]

shell(f"gatk --java-options '{java_opts}' GatherVcfs {vcfs} -O {out} \
    --TMP_DIR {tmp_dir} --REORDER_INPUT_BY_FIRST_VARIANT")
