try:
    data_files = snakemake.input.data
    out_data = snakemake.output.data
    out_seed = snakemake.output.seed
except NameError:
    # testing
    data_files = [f"output/default/est-sfs/data/{n}.txt" for n in range(1, 10)]
    out_data = "scratch/data.combined.txt"
    out_seed = "scratch/seed.txt"

# create seed file
open(out_seed, 'w').write('0')

# combine contents of all data files
with open(out_data, 'w') as out:
    for file in data_files:
        out.writelines(open(file, 'r').readlines())
