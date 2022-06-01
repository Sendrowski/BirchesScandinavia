"""
Calculate site-wise nucleotide diversity for given VCF file.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

vcf = snakemake.input.vcf
out = snakemake.params.prefix

shell(f"vcftools --site-pi --gzvcf {vcf} --out {out}")
