from snakemake.shell import shell

input = snakemake.input[0]
out = snakemake.output[0]

# cleanup file
shell(f"agat_convert_sp_gxf2gxf.pl -g {input} -o {out} --gvi 3")
