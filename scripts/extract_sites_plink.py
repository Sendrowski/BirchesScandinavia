from snakemake.shell import shell

try:
    ids = snakemake.input.ids
    input = snakemake.input.bed
    out = snakemake.output.bed
except NameError:
    # testing
    ids = "output/default/snps/pendula/biallelic/plink.prune.in"
    input = "output/default/snps/pendula/biallelic/snps.bed"
    out = "output/default/snps/pendula/biallelic/snps.admixture.bed"

input_prefix = input.replace('.bed', '')
out_prefix = out.replace('.bed', '')

shell(f"plink --bfile {input_prefix} --extract {ids} --make-bed --out {out_prefix} --allow-extra-chr '0'")
