#!/usr/bin/env Rscript
# patch_gene_names.R
# Patches existing mashr RDS files to add external_gene_name from var table.
# Priority: human ortholog (one2one) > external_gene_name > ensembl_id

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

# load ortholog map
orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[orth$Human.homology.type == 'ortholog_one2one' & orth$Human.gene.name != '',
            c('Gene.stable.ID','Human.gene.name')]
orth = orth[!duplicated(orth$Gene.stable.ID),]
rownames(orth) = orth$Gene.stable.ID

for (ct in cell.type.levels) {
    rds.file  = file.path(opt$checkpoints, paste0(ct, '_dglm_mashr_results.rds'))
    name.file = file.path(opt$checkpoints, paste0(ct, '_gene_names.csv'))

    if (!file.exists(rds.file)) {
        message('Skipping ', ct, ' — no RDS file')
        next
    }
    if (!file.exists(name.file)) {
        message('Skipping ', ct, ' — no gene names CSV yet')
        next
    }

    message('Patching ', ct)
    obj       = readRDS(rds.file)
    gene.ids  = dimnames(obj$dglm_results)[[1]]
    gene.names = read.csv(name.file, row.names=1)

    # build symbol vector: human ortholog > external_gene_name > ensembl_id
    symbols = sapply(gene.ids, function(g) {
        if (g %in% rownames(orth) && orth[g,'Human.gene.name'] != '') {
            orth[g,'Human.gene.name']
        } else if (g %in% rownames(gene.names) &&
                   gene.names[g,'external_gene_name'] != g &&
                   gene.names[g,'external_gene_name'] != '') {
            as.character(gene.names[g,'external_gene_name'])
        } else {
            g
        }
    })
    names(symbols) = gene.ids

    obj$human_symbols = symbols
    saveRDS(obj, file=rds.file)
    message('  patched: ', sum(symbols != gene.ids), ' genes renamed')
}

message('done.')
