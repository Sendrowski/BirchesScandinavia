import os

import numpy as np
import pandas as pd

import utils

try:
    targets_file = snakemake.input.targets
    gff = snakemake.input.gff
    n_files = snakemake.params.n_files
    out = snakemake.output
except NameError:
    # testing
    targets_file = "resources/reference/targeted_positions.txt"
    gff = "resources/reference/genome.corrected.gff.gz"
    n_files = 1
    out = ["scratch/intervals1.bed"]

pd.options.mode.chained_assignment = None

# the list of targeted genes whose positions are relative to the gene
targets = pd.read_csv(targets_file, sep='\t')

# load the genes from the annotation file
genes = utils.load_gff(gff)
genes = genes[genes.feature == 'gene']

# only keep rows where TanjaKeep is true
targets = targets.loc[targets.TanjaKeep]
targets = targets.reset_index(drop=True)

# initialize array
rows = []

# loop through genes and fetch their positions relative to the contigs
for i, target in targets.iterrows():
    gene = genes[genes.attribute.str.contains('ID=' + target.CHROMOSOME + ';')].iloc[0]

    contig = gene.seqname

    # We use 0-based coordinates for start and 1-based for end positions
    # as is done BED files. The GFF file uses 1-based coordinates so we need
    # to offset start by 1
    start = (gene.start - 1) + target.START

    # The stop position in the BED file is non-inclusive so we
    # just add the length to the start value
    stop = start + target.LENGTH

    rows.append([contig, start, stop])

    if i % 100 == 0: print(f"Processed genes: {i}/{len(targets)}", flush=True)

# randomize order and split arrays into chunks
# we do this to obtain chunks of similar sequence length
np.random.seed(0)
np.random.shuffle(rows)

# partition set of intervals
rows = np.array_split(rows, n_files)

# write contents to files
for i, chunk in enumerate(rows):
    # sort chunk in ascending contig order
    chunk = sorted(chunk, key=lambda x: (len(x[0]), x[0]))

    with open(out[i], 'w') as f:
        f.write(os.linesep.join(['\t'.join(row) for row in chunk]) + os.linesep)
