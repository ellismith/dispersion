#!/usr/bin/env Rscript

# Number of cores (0 = auto-detect, capped at 16)
n.cores = 16

# Predictor of interest
predictor       = 'age'
predictor.label = 'Age'
predictor.units = 'Years'

# Sex variable
sex.variable = 'sex'

# Age variable
age.variable = 'age'

# Covariates in mean model — following Chiou et al. approach
# age + sex + sequencing batch + library size proxies
model.covariates = c('age', 'sex', 'mean_n_umi', 'n_cells')

# Region levels
region.levels = c('ACC', 'CN', 'dlPFC', 'EC', 'HIP', 'IPP', 'lCb', 'M1', 'MB', 'mdTN', 'NAc')

# Cell types
cell.type.levels = c(
    'astrocytes',
    'basket_cells',
    'cerebellar_neurons',
    'ependymal_cells',
    'GABAergic_neurons',
    'glutamatergic_neurons',
    'medium_spiny_neurons',
    'microglia',
    'midbrain_neurons',
    'opc',
    'oligodendrocytes',
    'vascular_cells'
)

# h5ad map
h5ad.map = list(
    astrocytes           = 'Res1_astrocytes_update.h5ad',
    basket_cells         = 'Res1_basket-cells_update.h5ad',
    cerebellar_neurons   = 'Res1_cerebellar-neurons_subset.h5ad',
    ependymal_cells      = 'Res1_ependymal-cells_new.h5ad',
    GABAergic_neurons    = 'Res1_GABAergic-neurons_subset.h5ad',
    glutamatergic_neurons= 'Res1_glutamatergic-neurons_update.h5ad',
    medium_spiny_neurons = 'Res1_medium-spiny-neurons_subset.h5ad',
    microglia            = 'Res1_microglia_new.h5ad',
    midbrain_neurons     = 'Res1_midbrain-neurons_update.h5ad',
    opc                  = 'Res1_opc-olig_subset.h5ad',
    oligodendrocytes     = 'Res1_opc-olig_subset.h5ad',
    vascular_cells       = 'Res1_vascular-cells_subset.h5ad'
)

# h5ad base directory
h5ad.dir = '/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update'

# OPC louvain clusters
opc.louvain.clusters = c('12', '13')

# Minimum age filter
min.age = 1.0

# Minimum animals expressing gene
min.animals = 10

# False sign rate cutoff (mashr)
fsr.cutoff = 0.2

# Strong-subset q-value cutoff (mashr) -- used ONLY to select which genes are
# "strong" enough to drive data-driven covariance estimation (cov_pca/cov_ed).
# Matches Chiou et al.'s dglm_mashr.R, which hardcodes q<0.05 here (a
# different, stricter number than their final LFSR significance cutoff above).
strong.subset.qval.cutoff = 0.05

# Fraction of regions considered shared
fraction.shared.cutoff = 0.75

# Fraction of regions considered unique
fraction.unique.cutoff = 0.45

# Random seed
seed = 42

# Region colors
region.colors = c(
    ACC  = '#1b9e77',
    CN   = '#d95f02',
    dlPFC= '#7570b3',
    EC   = '#e7298a',
    HIP  = '#66a61e',
    IPP  = '#e6ab02',
    lCb  = '#a6761d',
    M1   = '#666666',
    MB   = '#ed1c24',
    mdTN = '#00aeef',
    NAc  = '#86328c'
)

# Auto-detect cores, cap at n.cores
n.cores = min(n.cores, parallel::detectCores(logical=FALSE))
message('Using ', n.cores, ' cores')

# Cell type acronym labels (for plotting)
cell.type.labels = c(
    astrocytes            = 'AST',
    basket_cells          = 'BC',
    cerebellar_neurons    = 'CER',
    ependymal_cells       = 'EPEN',
    GABAergic_neurons     = 'INH',
    glutamatergic_neurons = 'EXC',
    medium_spiny_neurons  = 'MSN',
    microglia             = 'MGL',
    midbrain_neurons      = 'MBN',
    opc                   = 'OPC',
    oligodendrocytes      = 'OLIG',
    vascular_cells        = 'VASC'
)
