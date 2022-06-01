###
# Parse polyDFE output.
###

if (exists('snakemake')) {
    input = snakemake@input[[1]]
    output = snakemake@output[[1]]
    postprocessing_source = snakemake@config$polydfe_postprocessing_source
} else {
    # testing
    input = "output/default/polydfe/pendula/output/C.full_anc.out"
    output = "scratch/dfe_parsed_output.txt"
    postprocessing_source = "resources/polydfe/postprocessing/script.R"
}


source(postprocessing_source)

intervals = c(-100, -10, -1, 0, 1)

est = c(sapply(input, parseOutput))
dfe = t(sapply(est, getDiscretizedDFE, intervals))

# write intervals and discretized DFE values
write(paste(intervals, collapse = " "), file = output, append = TRUE)
write(paste(dfe, collapse = " "), file = output, append = TRUE)

# write gradient
write(est[[1]]$criteria, file = output, append = TRUE)