"""
Codon utils for annotating sites with regards to synonymy and degeneracy.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import Bio.Data.CodonTable
import numpy as np
import pandas as pd
import re
import utils
from Bio.Seq import Seq

bases = np.array(['G', 'A', 'T', 'C'])
codon_table = Bio.Data.CodonTable.standard_dna_table.forward_table
stop_codons = Bio.Data.CodonTable.standard_dna_table.stop_codons
start_codons = ['ATG']

# include stop codons
for codon in stop_codons:
    codon_table[codon] = 'Σ'

# The degeneracy of the site according to how many unique amino acids
# are code for when change the site within the codon.
# We count the third position of the isoleucine codon as 2-fold degenerate.
# This is the only site that would normally have 3-fold degeneracy
# https://en.wikipedia.org/wiki/Codon_degeneracy
unique_to_degeneracy = {0: 0, 1: 2, 2: 2, 3: 4}


def load_cds(gff):
    # load the coding sequences and genes into a DataFrame
    annotation = utils.load_gff(gff).sort_values('start')

    # prepare list of coding sequences
    cds = annotation[annotation.feature == 'CDS']
    cds['phase'] = pd.to_numeric(cds['phase'], downcast='integer')
    cds = cds.drop_duplicates(subset=['seqname', 'start'])
    return cds.reset_index(drop=True)


def load_genes(gff):
    # load the coding sequences and genes into a DataFrame
    annotation = utils.load_gff(gff).sort_values('start')

    # prepare list of genes
    genes = annotation[annotation.feature == 'gene']
    return genes.reset_index(drop=True)


# get the gene enclosing given contig and position
# returns None if no gene could be found
def get_gene(record, genes):
    rows = genes[genes.seqname == record.CHROM]
    rows = rows[(rows.start <= record.POS) & (record.POS <= rows.end)]

    if len(rows) == 0:
        raise LookupError('No gene found')

    return rows.iloc[0]


# get the coding sequence among the given coding sequences
# that encloses the given position
# returns None if no coding sequence could be found
def get_cd(record, cds):
    rows = cds[(cds.seqname == record.CHROM) & (cds.start <= record.POS) & (record.POS <= cds.end)]

    if len(rows) == 0:
        raise LookupError('No coding sequence found')

    return rows.iloc[0]


# get the coding sequence on same gene preceding the given coding sequence
def get_previous_cd(cd, cds):
    rows = cds[(cds.seqname == cd.seqname) & (cds.end < cd.start)]

    return rows.tail(1).iloc[0] if len(rows) else None


# get the coding sequence on same gene subsequent to the given coding sequence
def get_next_cd(cd, cds):
    rows = cds[(cds.seqname == cd.seqname) & (cds.start > cd.end)]

    return rows.head(1).iloc[0] if len(rows) else None


# get current, previous and subsequence coding sequence
def get_cd_context(record, cds):
    cd = get_cd(record, cds)
    return get_previous_cd(cd, cds), cd, get_next_cd(cd, cds)


# fetch the contig the record is on
# we assume that the records are passed in ascending order
# as we use a stream for accessing the contig sequences
def get_contig(record, ref_reader):
    contig = next(ref_reader)
    seq = str(contig.seq)

    while contig.id != record.CHROM:
        contig = next(ref_reader)
        seq = str(contig.seq)

    return contig, seq


# initialize new record from vcf file
def init_record(record, info_fields):
    for name, value in info_fields.items():
        record.INFO[name] = value

    return record


# write given record to file and initialize new record
def write_and_init(record, writer, reader, info_fields):
    writer.write_record(record)
    return init_record(next(reader), info_fields)


# parse forward coding sequence
def parse_codon_forward(record, seq, cd_current, cd_previous, cd_next):
    # position relative to start of coding sequence
    pos_rel = record.POS - (cd_current.start + cd_current.phase)

    # position relative to codon
    pos_codon = pos_rel % 3

    # inclusive codon start, 1-based
    codon_start = record.POS - pos_codon

    # the codon positions
    codon_pos = [codon_start, codon_start + 1, codon_start + 2]

    if cd_previous is None and codon_pos[0] < cd_current.start:
        raise IndexError(f'Codon at site {record.CHROM}:{record.POS} '
                         f'starts before current CDS while no previous CDS was given.')

    # we assume here that cd_previous and cd_next have the same orientation
    # use final positions from previous coding sequence if current codon
    # starts before start position of current coding sequence
    if codon_pos[1] == cd_current.start:
        codon_pos[0] = cd_previous.end
    elif codon_pos[2] == cd_current.start:
        codon_pos[1] = cd_previous.end
        codon_pos[0] = cd_previous.end - 1

    if cd_next is None and codon_pos[2] > cd_current.end:
        raise IndexError(f'Codon at site {record.CHROM}:{record.POS} '
                         f'ends after current CDS while no following CDS was given.')

    # use initial positions from subsequent coding sequence if current codon
    # end before end position of current coding sequence
    if codon_pos[2] == cd_current.end + 1:
        codon_pos[2] = cd_next.start
    elif codon_pos[1] == cd_current.end + 1:
        codon_pos[1] = cd_next.start
        codon_pos[2] = cd_next.start + 1

    # seq uses 0-based positions
    codon = ''.join(seq[pos - 1] for pos in codon_pos).upper()

    return codon, codon_pos, codon_start, pos_codon, pos_rel


# parse reverse coding sequence
def parse_codon_reverse(record, seq, cd_current, cd_previous, cd_next):
    # position relative to end of coding sequence
    pos_rel = (cd_current.end - cd_current.phase) - record.POS

    # position relative to codon end
    pos_codon = pos_rel % 3

    # inclusive codon start, 1-based
    codon_start = record.POS + pos_codon

    # the codon positions
    codon_pos = [codon_start, codon_start - 1, codon_start - 2]

    if cd_previous is None and codon_pos[2] < cd_current.start:
        raise IndexError(f'Codon at site {record.CHROM}:{record.POS} '
                         f'starts before current CDS while no previous CDS was given.')

    # we assume here that cd_previous and cd_next have the same orientation
    # use final positions from previous coding sequence if current codon
    # ends before start position of current coding sequence
    if codon_pos[1] == cd_current.start:
        codon_pos[2] = cd_previous.end
    elif codon_pos[0] == cd_current.start:
        codon_pos[1] = cd_previous.end
        codon_pos[2] = cd_previous.end - 1

    if cd_next is None and codon_pos[0] > cd_current.end:
        raise IndexError(f'Codon at site {record.CHROM}:{record.POS} '
                         f'ends after current CDS while no following CDS was given.')

    # use initial positions from subsequent coding sequence if current codon
    # starts before end position of current coding sequence
    if codon_pos[0] == cd_current.end + 1:
        codon_pos[0] = cd_next.start
    elif codon_pos[1] == cd_current.end + 1:
        codon_pos[0] = cd_next.start + 1
        codon_pos[1] = cd_next.start

    # we use 0-based positions here
    codon = ''.join(seq[pos - 1] for pos in codon_pos)

    # take complement and convert to uppercase ('n' might be lowercase)
    codon = str(Seq(codon).complement()).upper()

    return codon, codon_pos, codon_start, pos_codon, pos_rel


# parse the default and alternative codons from the VEP annotation if present
def parse_codons_vep(record):
    codon_pairs = []
    if 'CSQ' in record.INFO:

        # iterate over fields
        # we might have several codon pairs when SNPs are adjacent
        # as there are several possibles of combining them in this case
        for field in record.INFO['CSQ']:
            if match := re.search("(.{3})/(.{3})", field):
                codon_pairs.append([m.upper() for m in [match[1], match[2]]])

    return codon_pairs


# obtain alternative allele from given record
def get_alt_allele(record, strand):
    alt = None
    for allele in record.ALT:

        # assume there is at most one alternative allele
        if str(allele) in bases:

            if strand == '-':
                alt = Seq(str(allele)).complement().__str__()
            else:
                alt = str(allele)

    return alt


# translate codon into amino acid
def get_degeneracy(codon, pos):
    amino_acid = codon_table[codon]
    codon = list(codon)
    alt = np.array([])

    for b in bases[bases != codon[pos]]:
        codon[pos] = b
        alt = np.append(alt, codon_table[''.join(codon)])

    return unique_to_degeneracy[sum(amino_acid == alt)]


# create codon degeneracy table
def get_degeneracy_table():
    codon_degeneracy = {}
    for a in bases:
        for b in bases:
            for c in bases:
                codon = ''.join([a, b, c])
                codon_degeneracy[codon] = ''.join([str(get_degeneracy(codon, pos)) for pos in range(0, 3)])

    return codon_degeneracy


# mutate the given codon to 'alt' at position 'pos'
def mutate(codon, alt, pos):
    return codon[0:pos] + alt + codon[pos + 1:]


# whether the given codons code for the same amino acid
def is_synonymous(codon1, codon2):
    # handle case where there are stop codons
    if codon1 in stop_codons or codon2 in stop_codons:
        return codon1 in stop_codons and codon2 in stop_codons

    return codon_table[codon1] == codon_table[codon2]


# apply the given callback on every coding site
# callback signature: (record, cd, codon, codon_pos, codon_start, pos_codon, contig)
# n_errors_tol: maximum number of assertion errors until raising an exception
# we pass a lot of codon information here
# the returned record is added to vcf_writer
def map_sites(vcf_reader, vcf_writer, gff_file, ref_reader, callback, n_errors_tol=3, info_fields={}):
    # the erroneous records
    records_err = []

    cds = load_cds(gff_file)
    genes = load_genes(gff_file)

    try:
        record = init_record(next(vcf_reader), info_fields)
        contig, seq = get_contig(record, ref_reader)
        gene, cd = None, None

        while True:

            # fetch gene if not up to date
            if gene is None or record.CHROM != gene.seqname or not (gene.start <= record.POS <= gene.end):

                try:
                    gene = get_gene(record, genes)
                except LookupError as err:
                    record = write_and_init(record, vcf_writer, vcf_reader, info_fields)
                    # print(f"No gene found, skipping record {record.CHROM}:{record.POS}")
                    continue

                print(f'Found gene {gene.attribute}')

            # fetch coding sequence if not up to date
            if cd is None or cd.seqname != record.CHROM or not (cd.start <= record.POS <= cd.end):

                try:
                    cd_previous, cd, cd_next = get_cd_context(record, cds)
                except LookupError as err:
                    record = write_and_init(record, vcf_writer, vcf_reader, info_fields)
                    # print(f"No coding sequence found, skipping record {record.CHROM}:{record.POS}")
                    continue

                print(f'Found coding sequence: {cd.seqname}:{cd.start}-{cd.end}, '
                      f'reminder: {(cd.end - cd.start + 1) % 3}, '
                      f'phase: {cd.phase}, orientation: {cd.strand}, '
                      f'current position: {record.CHROM}:{record.POS}')

            # fetch contig if not up to date
            if contig.id != record.CHROM:
                contig, seq = get_contig(record, ref_reader)

                print(f'Fetching contig: {contig.id}')

            # annotate if record is in coding sequence
            if cd.start <= record.POS <= cd.end:

                try:
                    # parse forward strand
                    if cd.strand == '+':
                        codon, codon_pos, codon_start, pos_codon, pos_rel = \
                            parse_codon_forward(record, seq, cd, cd_previous, cd_next)
                    else:
                        # parse reverse strand
                        codon, codon_pos, codon_start, pos_codon, pos_rel = \
                            parse_codon_reverse(record, seq, cd, cd_previous, cd_next)

                # skip site on IndexError
                except IndexError as err:
                    print(err)
                    record = write_and_init(record, vcf_writer, vcf_reader, info_fields)
                    continue

                # make sure the reference allele matches
                # with the position on the reference genome
                if contig[record.POS - 1] != record.REF[0]:
                    records_err.append(record)

                if len(records_err) > n_errors_tol:
                    raise AssertionError(f"Reference allele of site {record.CHROM}:{record.POS} does "
                                         f"not match with position in reference genome")

                # call callback function
                record = callback(record, cd, codon, codon_pos, codon_start, pos_codon, contig)

            record = write_and_init(record, vcf_writer, vcf_reader, info_fields)

    except StopIteration:
        vcf_writer.close()
        print(f"In total, {len(records_err)} reference allele assertion "
              f"errors occurred at {str([str(r) for r in records_err])}.")
