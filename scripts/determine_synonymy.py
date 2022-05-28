import numpy as np
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
vcf_reader.infos['Synonymy'] = Info('Synonymy', 1, 'Integer', 'Whether coding site is synonymous', None, None)
vcf_reader.infos['Synonymy_Info'] = Info('Synonymy_Info', 1, 'String',
                                         'Additional information about custom synonymy annotation:', None, None)

vcf_writer = vcf.Writer(open(out, 'w'), vcf_reader)
ref_reader = SeqIO.parse(ref, 'fasta')

# Maximum number of assertion errors until raising an exception.
# We introduce some leniency here as VEP produced different
# codons in a very cases for which no explanation could be found.
n_errors_tol = 3

# the erroneous records
records_err = []


# add synonymy annotation for given record
def annotate_synonymy(record, cd, codon, codon_pos, codon_start, pos_codon, contig):
    # fetch the alternative allele if present
    alt = codon_utils.get_alt_allele(record, cd.strand)

    info = ''
    synonymy, alt_codon, codons_vep = None, None, None
    if alt:
        # alternative codon
        # 'n' might not be in uppercase
        alt_codon = codon_utils.mutate(codon, alt, pos_codon).upper()

        # whether the alternative codon is synonymous
        if 'N' not in codon and 'N' not in alt_codon:
            synonymy = int(codon_utils.is_synonymous(codon, alt_codon))

        info += f'{alt_codon}'

        if alt_codon in codon_utils.stop_codons: info += ',stop_gained'
        if alt_codon in codon_utils.start_codons: info += ',start_gained'

        # fetch the codons from the VEP annotation
        codons_vep = codon_utils.parse_codons_vep(record)

        # make sure the codons determined by VEP are the same as our codons
        # we can only do the comparison for variant sites
        # only do the comparison if we could find exactly one pair of VEP codons
        # if there are several pairs then this might be because of adjacent SNPs
        # whose alternative codons are not supported at this point
        if len(codons_vep) == 1:

            # check for equality
            if not np.array_equal(codons_vep[0], [codon, alt_codon]):
                records_err.append(record)

            if len(records_err) > n_errors_tol:
                raise AssertionError(f"VEP Codons don't match at {codon_utils.format_errors(records_err)}. "
                                     f"Number of errors ({n_errors_tol}) exceeded.")

    # add to info field
    record.INFO['Synonymy'] = synonymy
    record.INFO['Synonymy_Info'] = info

    print(f'pos codon: {pos_codon}, pos abs: {record.POS}, '
          f'codon start: {codon_start}, strand: {cd.strand}, '
          f'ref allele: {contig[record.POS - 1]}, rec allele: {record.REF}, '
          f'codon pos: {str(codon_pos)}, codons: {str([codon, alt_codon])}, '
          f'VEP: {str(codons_vep)}', flush=True)

    return record


info_fields = {'Synonymy': None, 'Synonymy_Info': None}

# perform annotation
codon_utils.map_sites(vcf_reader, vcf_writer, gff, ref_reader, annotate_synonymy, info_fields=info_fields)

print(f"In total, {len(records_err)} codon assertion errors "
      f"occurred at {str([str(r) for r in records_err])}.")
