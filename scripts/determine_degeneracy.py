import vcf
from Bio import SeqIO
from vcf.parser import _Info as Info

import codon_utils

try:
    vcf_file = snakemake.input.vcf
    gff = snakemake.input.gff
    ref = snakemake.input.ref
    out = snakemake.output[0]
except NameError:
    # testing
    vcf_file = 'testing/snps1.ancestral.vep.annotated.vcf.gz'
    gff = 'resources/reference/genome.corrected.gff.gz'
    ref = 'resources/reference/genome.fasta'
    out = 'scratch/annotated.vcf'

vcf_reader = vcf.Reader(filename=vcf_file)

# add new fields to header
vcf_reader.infos['Degeneracy'] = Info('Degeneracy', 1, 'Integer', 'n-fold degeneracy', None, None)
vcf_reader.infos['Degeneracy_Info'] = Info('Degeneracy_Info', 1, 'String',
                                           'Additional information about custom degeneracy annotation:', None, None)

vcf_writer = vcf.Writer(open(out, 'w'), vcf_reader)
ref_reader = SeqIO.parse(ref, 'fasta')


# add degeneracy annotation for given record
def annotate_degeneracy(record, cd, codon, codon_pos, codon_start, pos_codon, contig):
    degeneracy = None
    if 'N' not in codon:
        degeneracy = codon_utils.get_degeneracy(codon, pos_codon)

    record.INFO['Degeneracy'] = degeneracy
    record.INFO['Degeneracy_Info'] = f"{pos_codon},{cd.strand},{codon}"

    print(f'pos codon: {pos_codon}, pos abs: {record.POS}, '
          f'codon start: {codon_start}, codon: {codon}, '
          f'strand: {cd.strand}, ref allele: {contig[record.POS - 1]}, '
          f'degeneracy: {degeneracy}, codon pos: {str(codon_pos)}, '
          f'rec allele: {record.REF}', flush=True)

    return record


info_fields = {'Degeneracy': None, 'Degeneracy_Info': None}

# perform annotation
codon_utils.map_sites(vcf_reader, vcf_writer, gff, ref_reader, annotate_degeneracy, info_fields=info_fields)
