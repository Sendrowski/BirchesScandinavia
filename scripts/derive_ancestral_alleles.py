from snakemake.shell import shell

data = snakemake.input.data
seed = snakemake.input.seed
config = snakemake.input.config
bin = snakemake.input.bin
out_sfs = snakemake.output.sfs
out_probs = snakemake.output.probs
tmp_dir = snakemake.resources.tmpdir

# execute command
shell(f"{bin} {config} {data} {seed} {out_sfs} {out_probs}")
