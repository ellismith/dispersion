import pandas as pd
import os

prefix_map = {
    'GABA': 'GABAergic_neurons', 'Glutamatergic': 'glutamatergic_neurons',
    'Cerebellar': 'cerebellar_neurons', 'MSN': 'medium_spiny_neurons',
    'Midbrain': 'midbrain_neurons', 'Basket': 'basket_cells',
    'Microglia': 'microglia', 'Astrocyte': 'astrocytes',
    'OPC': 'opc', 'Oligo': 'oligodendrocytes',
    'Vascular': 'vascular_cells', 'Ependymal': 'ependymal_cells'
}
disp_to_prefix = {v: k for k, v in prefix_map.items()}
region_map = {'ACC':'acc','CN':'cn','dlPFC':'dlpfc','EC':'ec','HIP':'hip',
              'IPP':'ipp','lCb':'lcb','M1':'m1','MB':'mb','mdTN':'mdtn','NAc':'nac'}

CCC_BASE = '/scratch/easmit31/cell_cell/results/within_region_analysis_corrected/hypergeometric_all_regions/subtype_pair_summary'
ccc = pd.concat([
    pd.read_csv(f'{CCC_BASE}/strengthening_subtype_pairs_sig_enriched.csv'),
    pd.read_csv(f'{CCC_BASE}/weakening_subtype_pairs_sig_enriched.csv')
], ignore_index=True)

def parse_sub(s):
    parts = s.rsplit('_', 1)
    return parts[0], int(parts[1])

ccc[['sender_prefix','sender_num']] = pd.DataFrame(ccc['sender'].map(parse_sub).tolist(), index=ccc.index)
ccc[['receiver_prefix','receiver_num']] = pd.DataFrame(ccc['receiver'].map(parse_sub).tolist(), index=ccc.index)
ccc['sender_broad'] = ccc['sender_prefix'].map(prefix_map)
ccc['receiver_broad'] = ccc['receiver_prefix'].map(prefix_map)

lcb = {'opc','microglia','astrocytes','vascular_cells','oligodendrocytes'}
cell_types = ['GABAergic_neurons','glutamatergic_neurons','cerebellar_neurons',
              'medium_spiny_neurons','midbrain_neurons','basket_cells',
              'microglia','astrocytes','opc','oligodendrocytes','vascular_cells','ependymal_cells']

DGLM_BASE = '/scratch/easmit31/dispersion/dglm'
disp_rows = []
for ct in cell_types:
    base = f"disp_age__FINAL_autosomeX_300c_{'noLcb_' if ct in lcb else ''}QCfiltered"
    f = f"{DGLM_BASE}/{base}/{ct}/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv"
    if not os.path.exists(f):
        print(f"SKIP: {ct}"); continue
    df = pd.read_csv(f, sep='\t')
    sig = df[df['mash_lfsr'] < 0.05].copy()
    for (subcl, reg), grp in sig.groupby(['cell_type','region']):
        lnum = int(subcl.rsplit('_',1)[1])
        n_inc = grp[grp['mash_beta']>0]['ensembl_id'].nunique()
        n_dec = grp[grp['mash_beta']<=0]['ensembl_id'].nunique()
        disp_rows.append({
            'cell_type_broad': ct, 'subcluster': subcl, 'louvain_num': lnum,
            'region': reg, 'region_ccc': region_map.get(reg, reg.lower()),
            'n_sig_inc': n_inc, 'n_sig_dec': n_dec,
            'net_direction': 'increasing' if n_inc > n_dec else 'decreasing'
        })

disp = pd.DataFrame(disp_rows)
print(f"Dispersion: {len(disp)} subtype x region combos")

results = []
for _, row in disp.iterrows():
    s = ccc[(ccc['sender_broad']==row['cell_type_broad']) & (ccc['sender_num']==row['louvain_num']) & (ccc['region']==row['region_ccc'])]
    r = ccc[(ccc['receiver_broad']==row['cell_type_broad']) & (ccc['receiver_num']==row['louvain_num']) & (ccc['region']==row['region_ccc'])]
    results.append({**row.to_dict(),
        'n_str_sender': len(s[s['direction']=='strengthening']),
        'n_wk_sender': len(s[s['direction']=='weakening']),
        'n_str_receiver': len(r[r['direction']=='strengthening']),
        'n_wk_receiver': len(r[r['direction']=='weakening'])})

out = pd.DataFrame(results)
out['total_ccc'] = out['n_str_sender']+out['n_wk_sender']+out['n_str_receiver']+out['n_wk_receiver']
out['in_ccc'] = out['total_ccc'] > 0
out.to_csv(f'{DGLM_BASE}/dispersion_ccc_overlap.csv', index=False)
print(f"Saved. Combos with CCC overlap: {out['in_ccc'].sum()} of {len(out)}")
print("\nTop 20 by total CCC interactions:")
print(out[out['in_ccc']].sort_values('total_ccc',ascending=False)[
    ['cell_type_broad','subcluster','region','net_direction',
     'n_sig_inc','n_sig_dec','n_str_sender','n_wk_sender','n_str_receiver','n_wk_receiver']
].head(20).to_string(index=False))
