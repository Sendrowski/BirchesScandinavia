from snakemake.shell import shell

file = snakemake.params.input_prefix
out = snakemake.params.output_prefix
midp = snakemake.params.midp

flag = 'midp' if midp else ''

# execute command
shell(f"plink --bfile {file} --hardy {flag} --hwe 0 --allow-extra-chr --out {out}")
