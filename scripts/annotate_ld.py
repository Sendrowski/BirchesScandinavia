from snakemake.shell import shell

try:
    vcf = snakemake.input[0]
    window_size_kb = snakemake.params.get('window_size_kb', 50)
    out = snakemake.output[0]
except NameError:
    # testing
    vcf = 'output/default/snps/pendula/biallelic/snps.vcf.gz'
    window_size_kb = 50
    out = 'scratch/vcf.r2.gz'

shell(f"bcftools +prune --annotate r2 {vcf} --window {window_size_kb}kb -Oz -o {out}")
