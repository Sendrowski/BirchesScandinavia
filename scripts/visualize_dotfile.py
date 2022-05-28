from snakemake.shell import shell

input = snakemake.input[0]
out = snakemake.output[0]

shell(f"dot {input} -Tsvg > {out}")
