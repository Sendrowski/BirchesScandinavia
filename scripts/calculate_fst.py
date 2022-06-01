"""
Calculate the Fst from given VCF file.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

vcf = snakemake.input.vcf
pops = snakemake.input.pops
out = snakemake.params.prefix

shell(f"vcftools --gzvcf {vcf} --weir-fst-pop {pops[0]} --weir-fst-pop {pops[1]} --out {out}")
