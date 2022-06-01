"""
Predict the variant effects with VEP.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

ref = snakemake.input.ref
vcf = snakemake.input.vcf
gff = snakemake.input.gff
out = snakemake.output[0]

shell(f"vep --gff {gff} --fasta {ref} -i {vcf} --vcf --species betula_pendula --allow_non_variant \
    -o {out} --fields 'Consequence,Codons' --compress_output bgzip")
