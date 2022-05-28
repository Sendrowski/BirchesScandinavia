from snakemake.shell import shell

snakefile = snakemake.params.snakefile
out = snakemake.output[0]

shell(f"snakemake --snakefile {snakefile} --rulegraph | dot -Tsvg > {out}")
