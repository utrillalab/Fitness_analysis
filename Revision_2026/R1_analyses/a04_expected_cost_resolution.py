"""R1-M3: quantitative justification for comparing fitness categories with the ME-model.

Question: could a single-gene RB-TnSeq assay detect the proteome cost of one protein at all?
Linear cost of unneeded protein (growth laws, Scott et al. 2010; Dekel & Alon 2005; Bienick et al. 2014):
    dmu/mu = phi / phi_max        (phi = mass fraction of the protein, phi_max = fraction that would stop growth)
Bienick et al. 2014: unused protein at 15% of mass halves growth  ->  phi_max ~ 0.30 (linear extrapolation).
Scott et al. 2010 extrapolate to ~0.48; both are used as bounds.
Pooled fitness is a log2 ratio after N generations:  f ~ N * dmu/mu  (N = 4-8 population doublings,
Price et al. 2018 Methods: "typically 4-8 population doublings").
So the smallest detectable gain (f = +0.5) needs  phi >= 0.5 * phi_max / N.
"""
import os
import numpy as np, pandas as pd
from common import *

df = long_table()
flux = load_flux()
tot = df.groupby('condition').mass_fg.transform('sum')
df['phi'] = df.mass_fg / tot

scen = {'optimistic (phi_max 0.30, 8 gen)': (0.30, 8), 'central (phi_max 0.30, 6 gen)': (0.30, 6),
        'conservative (phi_max 0.48, 4 gen)': (0.48, 4)}
rows = []
for name, (pmax, N) in scen.items():
    phi_min = 0.5 * pmax / N
    df['f_pred'] = N * df.phi / pmax
    for c, d in df.groupby('condition'):
        me = d[d.category == 'mean-effect']
        rows.append({'scenario': name, 'phi_needed_for_f=0.5_%': 100 * phi_min, 'condition': c,
                     'proteins_above_threshold': int((d.phi >= phi_min).sum()),
                     'of_which_mean_effect': int((me.phi >= phi_min).sum()),
                     'median_phi_mean_effect_%': 100 * me.phi.median(),
                     'median_predicted_f_mean_effect': me.f_pred.median(),
                     '%_mean_effect_with_predicted_f>=0.1': 100 * (me.f_pred >= 0.1).mean()})
res = pd.DataFrame(rows)
res.round(5).to_csv(f'{OUT}/a04_single_gene_detectability.csv', index=False)
print(res.groupby('scenario')[['phi_needed_for_f=0.5_%', 'proteins_above_threshold', 'of_which_mean_effect',
                               'median_phi_mean_effect_%', 'median_predicted_f_mean_effect',
                               '%_mean_effect_with_predicted_f>=0.1']].agg(['min', 'max']).round(4).T.to_string())

# aggregate: what the ME-model says could be spared, and the growth gain the linear law would predict
agg = []
for c in ME_CONDS:
    d = df[df.condition == c].assign(flux=lambda x: x.bnum.map(flux[c]))
    me = d[d.category == 'mean-effect']
    spare = me.loc[(me.flux == 0) | me.flux.isna(), 'phi'].sum()
    zero = me.loc[me.flux == 0, 'phi'].sum()
    n_zero = int((me.flux == 0).sum())
    agg.append({'condition': c, 'mean_effect_zero_flux_%': 100 * zero, 'n_genes_zero_flux': n_zero,
                'mean_per_gene_%': 100 * zero / max(n_zero, 1),
                'mean_effect_zero_flux_or_absent_%': 100 * spare,
                'predicted_gain_zero_flux_%_(phi_max 0.30)': 100 * zero / 0.30,
                'predicted_gain_zero_flux_%_(phi_max 0.48)': 100 * zero / 0.48})
agg = pd.DataFrame(agg).set_index('condition')
agg.round(3).to_csv(f'{OUT}/a04_aggregate_unused_proteome.csv')
print('\nAggregate unused (zero-flux) mean-effect proteome and the growth gain the linear cost law predicts:')
print(agg.round(2).to_string())
print('\nmedians:', agg.median().round(2).to_dict())

# ---- exact generations per experiment (Price et al. 2018 expsUsed; see a05_t_statistic.py) ------
gen_file = f'{OUT}/a05_consistency_and_generations.csv'
if os.path.exists(gen_file):
    gen = pd.read_csv(gen_file).groupby('condition').generations.mean().dropna()
    ex_rows = []
    for c, N in gen.items():
        d = df[df.condition == c]; me = d[d.category == 'mean-effect']
        for pmax in (0.30, 0.48):
            phi_min = 0.5 * pmax / N
            ex_rows.append({'condition': c, 'generations': round(N, 2), 'phi_max': pmax,
                            'phi_needed_for_f=0.5_%': 100 * phi_min, 'proteins_above': int((d.phi >= phi_min).sum()),
                            'median_predicted_f_mean_effect': (N * me.phi / pmax).median()})
    ex_rows = pd.DataFrame(ex_rows)
    ex_rows.round(5).to_csv(f'{OUT}/a04_detectability_exact_generations.csv', index=False)
    print('\nWith exact generations per condition:')
    print(ex_rows.groupby('phi_max')[['generations', 'phi_needed_for_f=0.5_%', 'proteins_above',
                                       'median_predicted_f_mean_effect']].agg(['min', 'max']).round(4).T.to_string())
