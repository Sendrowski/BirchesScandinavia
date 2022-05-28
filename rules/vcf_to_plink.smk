include: "header.smk"

# target output files
rule run_all_vcf_to_plink:
    input:
        expand("results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",sample_class=sample_classes,sample_set=["all", "birch", "left_cluster", "right_cluster", "pendula", "pubescens", "pendula_pubescens", "pendula_south", "pendula_north", "pubescens_south", "pubescens_north", "pendula_luis"],flag=flags),

# convert vcf file to PLINK format
rule convert_plink_sample_set:
    input:
        vcf="results/{sample_class}/snps/{sample_set}/{flag}/snps.vcf.gz"
    output:
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bed",
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.bim",
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.fam",
        "results/{sample_class}/snps/{sample_set}/{flag}/snps.log"
    params:
        prefix="results/{sample_class}/snps/{sample_set}/{flag}/snps"
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/convert_plink.py"
