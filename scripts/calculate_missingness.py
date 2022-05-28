from snakemake.shell import shell

file = snakemake.params.input_prefix
out = snakemake.params.output_prefix

# execute command
shell(f"plink --bfile {file} --missing --allow-extra-chr --out {out}")
