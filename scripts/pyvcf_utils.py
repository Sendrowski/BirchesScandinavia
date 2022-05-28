import re

import numpy as np
import vcf


# apply callback to sites and save result to file_out
# the record is omitted if 'callback' returns None
def map_sites(file_in, file_out, callback):
    reader = vcf.Reader(filename=file_in)
    writer = vcf.Writer(open(file_out, 'w'), reader)

    i = 0
    for record in reader:
        record = callback(record, i)

        if record is not None:
            writer.write_record(record)

        i += 1

    writer.close()


# count number of sites where the given callback returned true
def count_sites(file_in, callback):
    reader = vcf.Reader(filename=file_in)

    i = 0
    for j, record in enumerate(reader):
        count_site = callback(record, i)

        if count_site is True:
            i += 1

        if j % 1000 == 0:
            print(f"{j} sites processed")

    return i


# obtain a subsample from the given set of haplotypes
def subsample(haplotypes, size):
    # return an empty array when there are no haplotypes
    # this can happen for an uncalled outgroup site
    if len(haplotypes) == 0:
        return []

    replace = False

    # sample with replacement if size is greater than the number of haplotypes
    # this shouldn't happen very often
    if size > len(haplotypes):
        replace = True
        print(f'Only {len(haplotypes)} haplotypes found: subsampling {size} samples with replacement.')

    return np.array(haplotypes)[np.random.choice(len(haplotypes), size=size, replace=replace)]


# returns dict of counts indexed by A, C, G and T
def count(haplotypes):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for key, value in zip(*np.unique(haplotypes, return_counts=True)):
        counts[key] = value

    return counts


# returns 'A,C,G,T'
def base_dict_to_string(d):
    return ','.join(map(str, [d['A'], d['C'], d['G'], d['T']]))


# get a list of all called haplotypes
def haplotypes(calls):
    haplotypes = []

    for call in calls:
        if call.gt_bases:
            bases = re.split("/|\|", call.gt_bases)

            for base in bases:
                if base in ['A', 'C', 'G', 'T']:
                    haplotypes.append(base)

    return haplotypes


# restrict the calls to the sample set
def restrict(calls, names, exclude=False):
    if exclude:
        return list(filter((lambda call: call.sample not in names), calls))

    return list(filter((lambda call: call.sample in names), calls))

# get the alternative allele assuming that the record has at most one
def get_alt_allele(record):
    for allele in record.alleles:
        if allele and str(allele) != record.REF:
            return str(allele)
