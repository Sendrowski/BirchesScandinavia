
if (exists('snakemake')) {
    full_anc = snakemake@input[['full_anc']]
    full = snakemake@input[['full']]
    deleterious_anc = snakemake@input[['deleterious_anc']]
    deleterious = snakemake@input[['deleterious']]
    output = snakemake@output[[1]]
    postprocessing_source = snakemake@config$polydfe_postprocessing_source
} else {
    # testing
    full_anc = "output/default/polydfe/pendula/output/C.full_anc.out"
    full = "output/default/polydfe/pendula/output/C.full.out"
    deleterious_anc = "output/default/polydfe/pendula/output/C.deleterious_anc.out"
    deleterious = "output/default/polydfe/pendula/output/C.deleterious.out"
    output = "scratch/model_comparison.txt"
    postprocessing_source = "resources/polydfe/postprocessing/script.R"
}

source(postprocessing_source)

tt = rbind(compareModels(full_anc, full)$LRT,
           compareModels(deleterious_anc, deleterious)$LRT,
           compareModels(full, deleterious)$LRT,
           compareModels(full_anc, deleterious_anc)$LRT)

tt = data.frame(type1=c('full_anc', 'deleterious_anc', 'full', 'full_anc'),tt)
tt = data.frame(type2=c('full', 'deleterious', 'deleterious', 'deleterious_anc'),tt)

write.table(tt, output)