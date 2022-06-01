"""
Trim reads using TRIMMOMATIC.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

from snakemake.shell import shell

input = ' '.join(snakemake.input)
output = ' '.join(snakemake.output)
adapter = "$CONDA_PREFIX/share/trimmomatic-0.36-5/adapters/TruSeq3-PE-2.fa"
threads = snakemake.threads

shell(f"trimmomatic PE -threads {threads} {input} {output} \
    ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36")
