
if (exists('snakemake')) {
    input = snakemake@input[[1]]
    n = snakemake@params[['n']]
    prefix = snakemake@params[['prefix']]
    postprocessing_source = snakemake@config$polydfe_postprocessing_source
} else {
    input = "output/default/polydfe/pendula.1_pops.C.full_anc.out"
    n = 10
    prefix = "scratch/polydfe"
    postprocessing_source = "resources/polydfe/postprocessing/script.R"
}

source(postprocessing_source)

bootstrapData(input, outputfile = prefix, rep = n)
