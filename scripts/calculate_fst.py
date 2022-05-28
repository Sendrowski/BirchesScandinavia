from snakemake.shell import shell

vcf = snakemake.input.vcf
pops = snakemake.input.pops
out = snakemake.params.prefix

shell(f"vcftools --gzvcf {vcf} --weir-fst-pop {pops[0]} --weir-fst-pop {pops[1]} --out {out}")
