from snakemake.shell import shell

ref = snakemake.input.ref
left = snakemake.input.left
right = snakemake.input.right
name = snakemake.wildcards.name
out = snakemake.output.bam
tmp_dir = snakemake.resources.tmpdir
threads = snakemake.threads

# only map paired end reads which account for most of the data
shell(f"samtools sort -T {tmp_dir} <(bwa mem {ref} {left} {right} -R '@RG\\tID:{name}\\tSM:{name}' -t {threads}) > {out}")
