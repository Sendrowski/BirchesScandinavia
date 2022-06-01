"""
Create pairwise lists of linkage disequilibrium using PLINK.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2022-05-31"

import os.path

from snakemake.shell import shell

try:
    bed = snakemake.input.bed
    window_size = snakemake.params.get('window_size', 50)
    step_size = snakemake.params.get('step_size', 10)
    R2 = snakemake.params.get('R2', 0.1)
    out = snakemake.output[0]
except NameError:
    # testing
    input = 'output/default/snps/pendula/biallelic/snps.bed'
    window_size = 50
    step_size = 10
    R2 = 0.1
    out = 'scratch/'

output_prefix = os.path.dirname(out)
input_file_prefix = os.path.basename(bed).replace('.bed', '')

shell(f"""
    cd {output_prefix}
    plink --bfile {input_file_prefix} --indep-pairwise {window_size} {step_size} {R2} --allow-extra-chr
""")
