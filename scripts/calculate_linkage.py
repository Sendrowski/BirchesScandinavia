from snakemake.shell import shell

file = snakemake.params.input_prefix
out = snakemake.params.output_prefix

# execute command
# LD window options obtained from Luis Leal
shell(f"plink --bfile {file} --r2 --ld-window-kb 50 --ld-window 500 "
      f"--ld-window-r2 0 --allow-extra-chr --out {out}")
