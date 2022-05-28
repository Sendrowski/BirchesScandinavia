from snakemake.shell import shell

vcf = snakemake.input.vcf
prefix = snakemake.params.prefix

# create bed, bim and fam files
shell(f"plink --vcf {vcf} --out {prefix} --allow-extra-chr")
