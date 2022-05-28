include: "header.smk"

n_genomic_intervals = range(1,config['n_genomic_intervals'] + 1)
n_files_haplotype_calling = range(1,config['n_files_haplotype_calling'] + 1)


# get input function to determine path of fastq sequences
def get_path_loader(paths):
    def get_path(wildcards):
        path = [s for s in paths if wildcards.name in s][0]
        return "resources/" + path

    return get_path


# get resource path of specified outgroup
def get_url_outgroup(w):
    sample = outgroup_samples[outgroup_samples.name == w.name].iloc[0]
    return sample.url_left_read if w.n == 1 else sample.url_right_read

ruleorder:
    quality_tag_variants >
    call_imported_variants >
    create_tbi_index_vcf

# target output files
rule run_all_snp_calling:
    input:
        expand("results/default/trimmed_reads/{name}.{side}.fastq.gz",name=sample_names,side=[1, 2]),

        expand("results/default/quality_checks/{name}.1_fastqc.html",name=sample_names),
        expand("results/default/quality_checks/{name}.2_fastqc.html",name=sample_names),
        expand("results/default/quality_checks/multiqc_report.html"),

        expand("results/default/mapped_reads/{name}.bam",name=sample_names),
        expand("results/default/mapped_reads/{name}.marked_dups.bam",name=sample_names),
        expand("results/default/mapped_reads/{name}.marked_dups.bai",name=sample_names),

        expand("results/default/stats/coverage/{name}.txt",name=sample_names),
        expand("results/default/stats/coverage/all.txt"),

        expand("results/default/variants/haplotypes/{name}/all.sorted.g.vcf.gz",name=sample_names),
        expand("results/default/variants/intervals/intervals{n}.bed",n=n_genomic_intervals),
        expand("results/default/variants/genomics_dbs/db{n}",n=n_genomic_intervals),
        expand("results/default/variants/vcfs/variants{n}.vcf.gz",n=n_genomic_intervals),
        expand("results/default/variants/vcfs/variants{n}.quality.tagged.vcf.gz",n=n_genomic_intervals),

        expand("results/default/snps/raw/intervals/snps{n}.vcf.gz",n=n_genomic_intervals),
        expand("results/default/snps/raw/intervals/snps{n}.vcf.gz",n=n_genomic_intervals),
        expand("results/default/snps/raw/intervals/snps{n}.ancestral.vcf.gz",n=n_genomic_intervals),
        expand("results/default/snps/raw/snps.vcf.gz",sample_class=sample_classes)

# rules to execute locally
localrules:
    resource_reads_outgroup

# pseudo rule for origin of short reads
rule resource_reads_birch:
    output:
        protected('resources/' + birch_samples.path_left_read),
        protected('resources/' + birch_samples.path_right_read)

# fetch the outgroup reads from its online resource
rule resource_reads_outgroup:
    output:
        protected("resources/outgroup/{name}.{n}.fastq.gz")
    params:
        url=get_url_outgroup
    shell:
        "wget {params.url} -O {output}"

# pseudo rule for origin of reference genome
rule resource_reference_genome:
    output:
        protected(config['reference_genome'])

# pseudo rule for origin of annotation file
rule resource_annotation_file:
    output:
        protected(config['annotation_file'])

# pseudo rule for origin of list of targeted genes
rule resource_targeted_genes:
    output:
        protected(config['targeted_genes'])

# create a file containing all ingroup sample names
rule create_sample_set_ingroups:
    input:
        config['birch_samples']
    output:
        "results/{sample_class}/sample_sets/ingroups.args"
    params:
        sample_set="birch",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# create the subsample set with all putative B. pubescens species
rule create_sample_set_outgroups:
    input:
        config['outgroup_samples']
    output:
        "results/{sample_class}/sample_sets/outgroups.args"
    params:
        sample_set="outgroups",
        sample_class="default"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_sample_set.py"

# trim the Illumina reads using Trimmomatic
rule trim_reads:
    input:
        get_path_loader(samples.path_left_read.tolist()),
        get_path_loader(samples.path_right_read.tolist())
    output:
        "results/default/trimmed_reads/{name}.1.fastq.gz",
        "results/default/trimmed_reads/{name}.1.unpaired.fastq.gz",
        "results/default/trimmed_reads/{name}.2.fastq.gz",
        "results/default/trimmed_reads/{name}.2.unpaired.fastq.gz"
    conda:
        "../envs/trimmomatic.yaml"
    script:
        "../scripts/trim_reads.py"

# perform a quality check on the trimmed reads using FastQC
rule quality_check_trimmed_reads:
    input:
        "results/default/trimmed_reads/{name}.1.fastq.gz",
        "results/default/trimmed_reads/{name}.2.fastq.gz"
    output:
        "results/default/quality_checks/{name}.1_fastqc.html",
        "results/default/quality_checks/{name}.2_fastqc.html",
        "results/default/quality_checks/{name}.1_fastqc.zip",
        "results/default/quality_checks/{name}.2_fastqc.zip"
    params:
        outdir="results/default/quality_checks"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} --outdir {params.outdir}"

# aggregate the quality check using MultiQC
rule aggregate_quality_checks_trimmed_reads:
    input:
        expand("results/default/quality_checks/{name}.1_fastqc.html",name=sample_names),
        expand("results/default/quality_checks/{name}.2_fastqc.html",name=sample_names)
    output:
        "results/default/quality_checks/multiqc_report.html",
        directory("results/default/quality_checks/multiqc_data")
    params:
        outdir="results/default/quality_checks"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {params.outdir} -o {params.outdir}"

# index reference genome use BWA
rule index_reference_genome_bwa:
    input:
        config['reference_genome']
    output:
        "resources/reference/genome.fasta.bwt",
        "resources/reference/genome.fasta.amb",
        "resources/reference/genome.fasta.ann",
        "resources/reference/genome.fasta.pac",
        "resources/reference/genome.fasta.sa"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input}"

# map the reads to the reference genome
rule map_reads:
    input:
        ref=config['reference_genome'],
        bwt="resources/reference/genome.fasta.bwt",
        left="results/default/trimmed_reads/{name}.1.fastq.gz",
        right="results/default/trimmed_reads/{name}.2.fastq.gz"
    output:
        bam="results/default/mapped_reads/{name}.bam"
    conda:
        "../envs/read_mapping.yaml"
    script:
        "../scripts/map_reads.py"

# calculate the coverage of the alignments
rule calculate_coverage:
    input:
        "results/default/mapped_reads/{name}.bam"
    output:
        "results/default/stats/coverage/{name}.txt"
    group:
        "calculate_coverage"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools coverage -m {input} > {output}"

# calculate the coverage of all alignments
rule calculate_aggregate_coverage:
    input:
        expand("results/default/mapped_reads/{name}.bam",name=sample_names)
    output:
        "results/default/stats/coverage/all.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools coverage -m {input} > {output}"

# index reference genome using samtools
rule index_reference_genome_samtools:
    input:
        config['reference_genome']
    output:
        "resources/reference/genome.fasta.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

# mark the duplicates in the BAM files
rule mark_duplicates:
    input:
        bam="results/default/mapped_reads/{name}.bam"
    output:
        bam="results/default/mapped_reads/{name}.marked_dups.bam",
        metrics="results/default/mapped_reads/{name}.marked_dup_metrics.txt"
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/mark_duplicates.py"

# index BAM files which were marked for duplicates
rule index_mapped_reads:
    input:
        "results/default/mapped_reads/{name}.marked_dups.bam"
    output:
        "results/default/mapped_reads/{name}.marked_dups.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} {output}"

# create dict index for reference genome
rule index_reference_genome_dict_gatk:
    input:
        config['reference_genome']
    output:
        "resources/reference/genome.dict"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

# generate the files containing the intervals to operate over
# for calling the haplotypes
rule generate_intervals_haplotype_calling:
    input:
        targets=config['targeted_genes'],
        gff="resources/reference/genome.corrected.gff.gz"
    output:
        expand("results/default/variants/haplotypes/intervals{n}.bed",n=n_files_haplotype_calling)
    conda:
        "../envs/base.yaml"
    params:
        n_files=config['n_files_haplotype_calling']
    script:
        "../scripts/generate_genomic_intervals.py"

# call the variants of the aligned sequences using GVCF mode
rule call_haplotypes_diploid:
    input:
        ref=config['reference_genome'],
        bam="results/default/mapped_reads/{name}.marked_dups.bam",
        dict="resources/reference/genome.dict",
        fai="resources/reference/genome.fasta.fai",
        bai="results/default/mapped_reads/{name}.marked_dups.bai",
        intervals="results/default/variants/haplotypes/intervals{n}.bed"
    output:
        "results/default/variants/haplotypes/{name}/{n}.g.vcf.gz",
        "results/default/variants/haplotypes/{name}/{n}.g.vcf.gz.tbi"
    params:
        ploidy=2
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/call_haplotypes.py"

# gather the variants of the different intervals
rule gather_haplotypes:
    input:
        vcfs=expand("results/{{sample_class}}/variants/haplotypes/{{name}}/{n}.g.vcf.gz",n=n_files_haplotype_calling)
    output:
        "results/{sample_class}/variants/haplotypes/{name}/all.g.vcf.gz"
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/gather_variants.py"

# sort gvcf file
rule sort_haplotypes:
    input:
        "results/{sample_class}/variants/haplotypes/{name}/all.g.vcf.gz"
    output:
        "results/{sample_class}/variants/haplotypes/{name}/all.sorted.g.vcf.gz"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SortVcf -I {input} -O {output} --MAX_RECORDS_IN_RAM 50000"

# generate the files containing the genomic intervals
rule generate_genomic_intervals:
    input:
        targets=config['targeted_genes'],
        gff="resources/reference/genome.corrected.gff.gz"
    output:
        expand("results/default/variants/intervals/intervals{n}.bed",n=n_genomic_intervals)
    conda:
        "../envs/base.yaml"
    params:
        n_files=config['n_genomic_intervals']
    script:
        "../scripts/generate_genomic_intervals.py"

# import into a database the combined gvcf files
rule import_variants_diploid:
    input:
        vcfs=expand("results/default/variants/haplotypes/{name}/all.sorted.g.vcf.gz",name=sample_names),
        intervals="results/default/variants/intervals/intervals{n}.bed"
    output:
        directory("results/default/variants/genomics_dbs/db{n}")
    params:
        batch_size=config['genomics_import_batch_size']
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/import_variants.py"

# call the variants
rule call_imported_variants:
    input:
        ref=config['reference_genome'],
        db="results/{sample_class}/variants/genomics_dbs/db{n}",
        intervals="results/default/variants/intervals/intervals{n}.bed"
    output:
        "results/{sample_class}/variants/vcfs/variants{n}.vcf.gz",
        "results/{sample_class}/variants/vcfs/variants{n}.vcf.gz.tbi"
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/call_imported_variants.py"

# add quality filters to variants
rule quality_tag_variants:
    input:
        ref=config['reference_genome'],
        vcf="results/{sample_class}/variants/vcfs/variants{n}.vcf.gz",
        tbi="results/{sample_class}/variants/vcfs/variants{n}.vcf.gz.tbi"
    output:
        "results/{sample_class}/variants/vcfs/variants{n}.quality.tagged.vcf.gz",
        "results/{sample_class}/variants/vcfs/variants{n}.quality.tagged.vcf.gz.tbi"
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/quality_filter_snps.py"

# restrict to SNPs and remove sites that failed the quality filter
rule filter_snps:
    input:
        vcf="results/{sample_class}/variants/vcfs/variants{n}.quality.tagged.vcf.gz",
        tbi="results/{sample_class}/variants/vcfs/variants{n}.quality.tagged.vcf.gz.tbi"
    output:
        "results/{sample_class}/snps/raw/intervals/snps{n}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/filter_snps.py"

# download and install EST-SFS
rule setup_est_sfs:
    input:
        config['est_sfs_config']
    output:
        config['est_sfs_bin']
    conda:
        "../envs/est_sfs.yaml"
    script:
        "../scripts/setup_est_sfs.py"

# prepare the input files for EST-SFS
rule prepare_input_est_sfs:
    input:
        ingroups="results/{sample_class}/sample_sets/ingroups.args",
        outgroups="results/{sample_class}/sample_sets/outgroups.args",
        vcf="results/{sample_class}/snps/raw/intervals/snps{n}.vcf.gz"
    output:
        data="results/{sample_class}/est-sfs/data/{n}.txt"
    params:
        sample_class="{sample_class}"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/prepare_input_est_sfs.py"

# EST-SFS needs the information of all sites to correctly
# predict the ancestral state. If used on a subset, the high variation
# of frequencies seems to inflate the number of high frequency derived alleles.
# We thus gather all subsets here and split the result again for
# further processing.
rule gather_input_est_sfs:
    input:
        data=expand("results/{{sample_class}}/est-sfs/data/{n}.txt",n=n_genomic_intervals)
    output:
        data="results/{sample_class}/est-sfs/data/all.txt",
        seed="results/{sample_class}/est-sfs/seed.txt"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/gather_input_est_sfs.py"

# determine the ancestral alleles using the outgroups
rule derive_ancestral_alleles:
    input:
        data="results/{sample_class}/est-sfs/data/all.txt",
        seed=ancient("results/{sample_class}/est-sfs/seed.txt"),
        config=config['est_sfs_config'],
        bin=config['est_sfs_bin']
    output:
        sfs="results/{sample_class}/est-sfs/sfs.txt",
        probs="results/{sample_class}/est-sfs/probs/all.txt"
    conda:
        "../envs/est_sfs.yaml"
    script:
        "../scripts/derive_ancestral_alleles.py"

# split EST-SFS output into chunks
rule scatter_output_est_sfs:
    input:
        probs="results/{sample_class}/est-sfs/probs/all.txt",
        data=expand("results/{{sample_class}}/est-sfs/data/{n}.txt",n=n_genomic_intervals)
    output:
        probs=expand("results/{{sample_class}}/est-sfs/probs/{n}.txt",n=n_genomic_intervals)
    conda:
        "../envs/est_sfs.yaml"
    script:
        "../scripts/scatter_output_est_sfs.py"

# correct the ancestral alleles
rule recode_ancestral_alleles_vcf:
    input:
        vcf="results/{sample_class}/snps/raw/intervals/snps{n}.vcf.gz",
        probs="results/{sample_class}/est-sfs/probs/{n}.txt",
        data="results/{sample_class}/est-sfs/data/{n}.txt",
        samples="results/{sample_class}/sample_sets/ingroups.args"
    output:
        vcf=temp("results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vcf")
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/recode_ancestral_alleles_vcf.py"

# compress the vcf file
rule compress_ancestral_alleles_vcf:
    input:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vcf"
    output:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vcf.gz"
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip {input} -c > {output}"

# clean up the annotation file using AGAT
rule cleanup_annotation_file:
    input:
        config['annotation_file']
    output:
        temp("resources/reference/genome.corrected.gff"),
        temp("genome.agat.log")
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/cleanup_annotation_file.py"

# sort the annotation in chromosomal order and compress it for use with VEP
rule sort_and_compress_gff:
    input:
        "resources/reference/genome.corrected.gff"
    output:
        "resources/reference/genome.corrected.gff.gz"
    conda:
        "../envs/tabix.yaml"
    shell:
        "grep -v '#' {input} | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output}"

# create a tbi index for the annotation file
rule index_gff:
    input:
        "resources/reference/genome.corrected.gff.gz"
    output:
        "resources/reference/genome.corrected.gff.gz.tbi"
    conda:
        "../envs/tabix.yaml"
    shell:
        "tabix -p gff {input}"

# predict the variants' effects with VEP
rule predict_variant_effects_vep:
    input:
        ref=config['reference_genome'],
        gff="resources/reference/genome.corrected.gff.gz",
        tbi="resources/reference/genome.corrected.gff.gz.tbi",
        vcf="results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vcf.gz"
    output:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.vcf.gz",
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.vcf.gz_summary.html"
    conda:
        "../envs/vep.yaml"
    script:
        "../scripts/predict_variant_effects_vep.py"

# determine the degeneracy for all coding sites
rule determine_degeneracy:
    input:
        vcf="results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.vcf.gz",
        gff="resources/reference/genome.corrected.gff.gz",
        ref=config['reference_genome']
    output:
        temp("results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.vcf")
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/determine_degeneracy.py"

# compress the vcf file
rule compress_degeneracy_vcf:
    input:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.vcf"
    output:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.vcf.gz"
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip {input} -c > {output}"

# determine the synonymy of the variant sites
rule determine_synonymy:
    input:
        vcf="results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.vcf.gz",
        gff="resources/reference/genome.corrected.gff.gz",
        ref=config['reference_genome']
    output:
        temp("results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.synonymy.vcf")
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/determine_synonymy.py"

# compress the vcf file
rule compress_synonymy_vcf:
    input:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.synonymy.vcf"
    output:
        "results/{sample_class}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.synonymy.vcf.gz"
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip {input} -c > {output}"

# gather the variants of the different intervals
rule gather_variants:
    input:
        vcfs=expand("results/{{sample_class}}/snps/raw/intervals/snps{n}.ancestral.vep.degeneracy.synonymy.vcf.gz",n=n_genomic_intervals)
    output:
        temp("results/{sample_class}/snps/raw/snps.unsorted.vcf.gz")
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/gather_variants.py"

# sort variants file
rule sort_variants_all:
    input:
        "results/{sample_class}/snps/raw/snps.unsorted.vcf.gz"
    output:
        "results/{sample_class}/snps/raw/snps.vcf.gz"
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SortVcf -I {input} -O {output} --MAX_RECORDS_IN_RAM 50000"

# plot the spectrum from EST-SFS
rule plot_sfs_est:
    input:
        "results/{sample_cass}/est-sfs/sfs.txt"
    output:
        "results/{sample_cass}/graphs/svg/sfs/est-sfs.{n_proj}.{unfolded}.{scale}.svg"
    params:
        n_proj=lambda w: int(w.n_proj),
        unfolded=lambda w: w.unfolded == "unfolded",
        log_scale=lambda w: w.scale == "log"
    conda:
        "../envs/dadi.yaml"
    script:
        "../scripts/plot_sfs_est.py"
