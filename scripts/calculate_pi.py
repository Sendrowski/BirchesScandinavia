from snakemake.shell import shell

vcf = snakemake.input.vcf
out = snakemake.params.prefix

shell(f"vcftools --site-pi --gzvcf {vcf} --out {out}")
