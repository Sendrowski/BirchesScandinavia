import os

from snakemake.shell import shell

try:
    bed = snakemake.input.bed
    out_log = snakemake.output.log
    K = snakemake.params.K
except NameError:
    # testing
    bed = "output/default/snps/pendula/biallelic/snps.admixture.bed"
    out_log = "scratch/admixture"
    K = 4

# We can't specify the path of the output files
# so we change the working directory to the output
# path of the given log file.
bed_abs_path = os.path.abspath(bed)
out_path = os.path.dirname(out_log)
log = os.path.basename(out_log)

shell(f"""
    mkdir -p {out_path}
    cd {out_path}
    admixture --cv {bed_abs_path} {K} | tee {log}
""")
