"""
plot_style.py
Shared style constants for summary heatmaps across all analyses.
"""

# Canonical internal names (underscores) — used by DGLM and gene_variance
CELL_TYPES = [
    'astrocytes', 'basket_cells', 'cerebellar_neurons', 'ependymal_cells',
    'GABAergic_neurons', 'glutamatergic_neurons', 'medium_spiny_neurons',
    'microglia', 'midbrain_neurons', 'opc', 'oligodendrocytes', 'vascular_cells',
]

CT_LABELS = {
    'astrocytes':             'AST',
    'basket_cells':           'BC',
    'cerebellar_neurons':     'CER',
    'ependymal_cells':        'EPEN',
    'GABAergic_neurons':      'INH',
    'glutamatergic_neurons':  'EXC',
    'medium_spiny_neurons':   'MSN',
    'microglia':              'MGL',
    'midbrain_neurons':       'MBN',
    'opc':                    'OPC',
    'oligodendrocytes':       'OLIG',
    'vascular_cells':         'VASC',
}

# Centroid summary CSVs use hyphens; opc is 'opc-olig'
CENTROID_CT_ORDER = [
    'astrocytes', 'basket-cells', 'cerebellar-neurons', 'ependymal-cells',
    'GABAergic-neurons', 'glutamatergic-neurons', 'medium-spiny-neurons',
    'microglia', 'midbrain-neurons', 'opc-olig', 'oligodendrocytes', 'vascular-cells',
]

CENTROID_CT_LABELS = {
    'astrocytes':            'AST',
    'basket-cells':          'BC',
    'cerebellar-neurons':    'CER',
    'ependymal-cells':       'EPEN',
    'GABAergic-neurons':     'INH',
    'glutamatergic-neurons': 'EXC',
    'medium-spiny-neurons':  'MSN',
    'microglia':             'MGL',
    'midbrain-neurons':      'MBN',
    'opc-olig':              'OPC',
    'oligodendrocytes':      'OLIG',
    'vascular-cells':        'VASC',
}

# Diversity filenames use full underscore names; opc is 'oligodendrocyte_precursor_cells'
DIVERSITY_CT_ORDER = [
    'astrocytes', 'basket_cells', 'cerebellar_neurons', 'ependymal_cells',
    'GABAergic_neurons', 'glutamatergic_neurons', 'medium_spiny_neurons',
    'microglia', 'midbrain_neurons', 'oligodendrocyte_precursor_cells',
    'oligodendrocytes', 'vascular_cells',
]

DIVERSITY_CT_LABELS = {
    'astrocytes':                      'AST',
    'basket_cells':                    'BC',
    'cerebellar_neurons':              'CER',
    'ependymal_cells':                 'EPEN',
    'GABAergic_neurons':               'INH',
    'glutamatergic_neurons':           'EXC',
    'medium_spiny_neurons':            'MSN',
    'microglia':                       'MGL',
    'midbrain_neurons':                'MBN',
    'oligodendrocyte_precursor_cells': 'OPC',
    'oligodendrocytes':                'OLIG',
    'vascular_cells':                  'VASC',
}

REGIONS = ['ACC', 'CN', 'dlPFC', 'EC', 'HIP', 'IPP', 'lCb', 'M1', 'MB', 'mdTN', 'NAc']

TICK_FS  = 13
LABEL_FS = 15
TITLE_FS = 14
CBAR_FS  = 13

# PC-residual analysis (factor_analysis/) uses hyphens + splits OPC/OLIG
PC_RESIDUAL_CT_ORDER = ["AST","BC","CER","EPEN","INH","EXC","MSN","MGL","MBN","OLIG","OPC","VASC"]
PC_RESIDUAL_CT_LABELS = {
    "astrocytes":"AST","basket-cells":"BC","cerebellar-neurons":"CER",
    "ependymal-cells":"EPEN","GABAergic-neurons":"INH","glutamatergic-neurons":"EXC",
    "medium-spiny-neurons":"MSN","microglia":"MGL","midbrain-neurons":"MBN",
    "oligodendrocytes":"OLIG","OPCs":"OPC","vascular-cells":"VASC",
}
