import numpy as np

include: "header.smk"
include: "snp_calling.smk"
include: "sample_sets.smk"

ruleorder:
    sort_variants_all >
    resource_raw_vcf


# get the names of all tetraploid samples
def get_tetraploid_names():
    # make sure the sample set file exists
    file = checkpoints.create_sample_set_pubescens.get(sample_class='default').output[0]
    return pd.read_csv(file,header=None)[0].values


# get the names of all diploid samples
def get_diploid_names():
    return sample_names[~np.isin(sample_names,get_tetraploid_names())]


# target output files
rule run_all_snp_calling_tetraploid:
    input:
        "results/tetraploid/snps/raw/snps.vcf.gz"

# call the variants of the aligned sequences using GVCF mode
# this step is dependent on the derivation of the pubescens sample set
rule call_haplotypes_tetraploid:
    input:
        "results/default/sample_sets/pubescens.args",
        ref=config['reference_genome'],
        bam="results/default/mapped_reads/{name}.marked_dups.bam",
        dict="resources/reference/genome.dict",
        fai="resources/reference/genome.fasta.fai",
        bai="results/default/mapped_reads/{name}.marked_dups.bai",
        intervals="results/default/variants/haplotypes/intervals{n}.bed",
    output:
        "results/tetraploid/variants/haplotypes/{name}/{n}.g.vcf.gz",
        "results/tetraploid/variants/haplotypes/{name}/{n}.g.vcf.gz.tbi"
    params:
        ploidy=4
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/call_haplotypes.py"


# import diploid and tetraploid samples into the database
rule import_variants_tetraploid_diploid:
    input:
        vcfs=lambda w: expand("results/tetraploid/variants/haplotypes/{name}/all.sorted.g.vcf.gz",name=get_tetraploid_names()) + \
                       expand("results/default/variants/haplotypes/{name}/all.sorted.g.vcf.gz",name=get_diploid_names()),
        intervals="results/default/variants/intervals/intervals{n}.bed"
    output:
        directory("results/tetraploid/variants/genomics_dbs/db{n}")
    params:
        batch_size=config['genomics_import_batch_size']
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/import_variants.py"
