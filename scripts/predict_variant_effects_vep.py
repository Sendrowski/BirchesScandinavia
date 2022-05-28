from snakemake.shell import shell

ref = snakemake.input.ref
vcf = snakemake.input.vcf
gff = snakemake.input.gff
out = snakemake.output[0]

shell(f"vep --gff {gff} --fasta {ref} -i {vcf} --vcf --species betula_pendula --allow_non_variant \
    -o {out} --fields 'Consequence,Codons' --compress_output bgzip")
