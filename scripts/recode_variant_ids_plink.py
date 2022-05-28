from snakemake.shell import shell

try:
    input = snakemake.input.bed
    out = snakemake.output.bed
except NameError:
    # testing
    input = "output/default/snps/pendula/biallelic/snps.bed"
    out = "output/default/snps/pendula/biallelic/snps.recoded_ids.bed"

input_prefix = input.replace('.bed', '')
out_prefix = out.replace('.bed', '')

shell(f"plink2 --bfile {input_prefix} --set-all-var-ids @_# --make-bed --out {out_prefix} --allow-extra-chr")
