from snakemake.shell import shell

# from snakemake_wrapper_utils.java import get_java_opts

vcf = snakemake.input.vcf
# samples = snakemake.input.samples
# java_opts = get_java_opts(snakemake)
out = snakemake.output[0]

# keep only mono and bi-allelic SNPs
# we retain non-variants sites here which we will annotate
# we do this as we need the total number of surveyed sites for predicting the DFE
shell(f"bcftools view {vcf} -Oz --apply-filters .,PASS --max-alleles 2 --exclude-types indels,mnps,other > {out}")

# this didn't work for some reasons. There were still indels among the filtered variants
"""shell(f"gatk --java-options '{java_opts}' SelectVariants -V {vcf} -O {out} \
    --select-type-to-include SNP --select-type-to-include NO_VARIATION \
    --select-type-to-exclude INDEL --select-type-to-exclude SYMBOLIC \
    --select-type-to-exclude MNP --select-type-to-exclude MIXED \
    --sample-name {samples}")"""
