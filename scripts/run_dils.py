import os

from snakemake.shell import shell

wdir = os.getcwd()
config = "resources/dils/sim/config.yaml"
binpath = "resources/dils/bin"
snakefile = "resources/dils/bin/Snakefile_2pop"

input = snakemake.input[0]
sample_set = snakemake.params.sample_set
sample_class = snakemake.params.sample_class
outdir = snakemake.output[0]

# os.environ['R_LIBS_USER'] = os.environ['CONDA_PREFIX'] + '/bin'

# There were difficulties installing the required R packages into
# the corresponding Conda environment. That's why the module
# system is used here.
shell(f"""
    module load bioinfo-tools snakemake/6.9.1
    module load R/3.6.0
    module load R_packages/3.6.0
    snakemake --snakefile {snakefile} --configfile '{wdir}/{config}' --config \
    infile='{wdir}/{input}' config_yaml='{wdir}/{config}' binpath='{wdir}/{binpath}' \
    timeStamp={outdir} --profile slurm/dils
""")

# move Rplots.pdf which ended up in the root directory
# shell(f"mv Rplots.pdf {outdir}")
