import dadi

try:
    vcf = snakemake.input.vcf
    pops = snakemake.input.pops
    out = snakemake.output[0]
except NameError:
    # testing
    vcf = "output/default/snps/pendula_pubescens/biallelic/snps.vcf.gz"
    pops = "output/default/sample_sets/subpopulations/pendula_pubescens/2_pops.txt"
    out = "scratch/pendula_pubescens.2_pops.dadi"

# project down to a smaller sample size
n_proj = 10

# create data dict from VCF file
dd = dadi.Misc.make_data_dict_vcf(vcf, pops, subsample={'pop0': n_proj, 'pop1': n_proj})

# convert to dadi SNP format as described at
# https://dadi.readthedocs.io/en/latest/user-guide/importing-data/
sep = '\t'
with open(out, 'w') as f:
    # write header
    f.write(sep.join(['birch', 'elder', 'Allele1', 'pop0', 'pop1', 'Allele2', 'pop0', 'pop1', 'contig', 'positions']) + '\n')

    # write sites
    for name, site in dd.items():
        f.write(sep.join([site['context'],
                          site['outgroup_context'],
                          site['segregating'][0],
                          str(site['calls']['pop0'][0]),
                          str(site['calls']['pop1'][0]),
                          site['segregating'][1],
                          str(site['calls']['pop0'][1]),
                          str(site['calls']['pop1'][1]),
                          *name.split('_')]) + '\n')
