###
# Create init files for bootstraps using polyDFE.
###

if (exists('snakemake')) {
    input = snakemake@input[[1]]
    prefix = snakemake@params[['prefix']]
    type = snakemake@params[['type']]
    postprocessing_source = snakemake@config$polydfe_postprocessing_source
    
    source(postprocessing_source)
} else {
    input = "output/default/polydfe/pendula/out/C.deleterious.out"
    prefix = "scratch/init_C.deleterious"
    type = "deleterious"
    postprocessing_source = "https://github.com/paula-tataru/polyDFE/raw/master/postprocessing.R"
    
    library(devtools)
    source_url(postprocessing_source)
}

# For some reason we need to pass the parameters to be fixed explicitly.
# We only support model C here
fix = c("eps_cont")

if (type == 'deleterious_anc' | type == 'deleterious') {
    fix = c(fix, c("p_b", "S_b"))
}

if (type == 'deleterious' | type == 'full') {
    fix = c(fix, "eps_an")
}

createInitLines(input, prefix, fix = fix)
