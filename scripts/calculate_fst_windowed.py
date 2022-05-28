from snakemake.shell import shell

vcf = snakemake.input.vcf
pops = snakemake.input.pops
out = snakemake.params.prefix

shell(f"vcftools --gzvcf {vcf} --weir-fst-pop {pops[0]} --weir-fst-pop {pops[1]} --fst-window-size 1000000000 --out {out}")
