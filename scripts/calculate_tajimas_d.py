from snakemake.shell import shell

vcf = snakemake.input.vcf
out = snakemake.params.prefix

shell(f"vcftools --gzvcf {vcf} --TajimaD 1000 --out {out}")
