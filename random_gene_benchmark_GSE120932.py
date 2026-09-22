import pandas as pd
import numpy as np
from scipy.stats import hypergeom

TOP_TABLE = 'GSE120932.top.table.tsv'
N_RANDOM = 100000
SEED = 20260910

b56 = {'PPP2R5A','PPP2R5B','PPP2R5C','PPP2R5D','PPP2R5E'}

df = pd.read_csv(TOP_TABLE, sep='\t')
df = df.dropna(subset=['Gene.symbol']).copy()

# Preserve the probe-level B56 input used in the manuscript.
b56_probe = df[df['Gene.symbol'].isin(b56)].copy()
b56_cols = [c for c in ['ID','Gene.symbol','logFC','adj.P.Val','P.Value','B'] if c in b56_probe.columns]
b56_probe[b56_cols].sort_values('adj.P.Val').to_csv('GSE120932_B56_probe_results.csv', index=False)

rows = []
for gene, x in df.groupby('Gene.symbol'):
    rows.append({
        'gene': gene,
        'n_probes': len(x),
        'qualifying_downregulated_gene': bool((x['logFC'] < 0).all() and (x['adj.P.Val'] < 0.05).all())
    })
genes = pd.DataFrame(rows)

one = genes[genes.n_probes == 1]
three = genes[genes.n_probes == 3]

M1, K1 = len(one), int(one.qualifying_downregulated_gene.sum())
M3, K3 = len(three), int(three.qualifying_downregulated_gene.sum())

rng = np.random.default_rng(SEED)
one_q = one.qualifying_downregulated_gene.to_numpy(dtype=int)
three_q = three.qualifying_downregulated_gene.to_numpy(dtype=int)
counts = np.empty(N_RANDOM, dtype=int)
for i in range(N_RANDOM):
    counts[i] = rng.choice(one_q, 3, replace=False).sum() + rng.choice(three_q, 2, replace=False).sum()

empirical_p = ((counts >= 3).sum() + 1) / (N_RANDOM + 1)

exact_p = 0.0
for k1 in range(4):
    for k3 in range(3):
        if k1 + k3 >= 3:
            exact_p += hypergeom.pmf(k1, M1, K1, 3) * hypergeom.pmf(k3, M3, K3, 2)

# Save compact distribution and summary so the distributed outputs are recreated directly.
vals, freqs = np.unique(counts, return_counts=True)
dist = pd.DataFrame({
    'qualifying_gene_count': vals.astype(int),
    'frequency': freqs.astype(int),
    'proportion': freqs / N_RANDOM
})
dist.to_csv('GSE120932_random_gene_benchmark_distribution.csv', index=False)

b56_gene = genes[genes.gene.isin(b56)]
b56_qualifying = int(b56_gene.qualifying_downregulated_gene.sum())
summary = pd.DataFrame({
    'metric': [
        'B56 qualifying genes','B56 panel size','matched random draws','empirical P value',
        'exact matched probability','one probe eligible genes','one probe qualifying genes',
        'three probe eligible genes','three probe qualifying genes'
    ],
    'value': [b56_qualifying, len(b56), N_RANDOM, empirical_p, exact_p, M1, K1, M3, K3]
})
summary.to_csv('GSE120932_random_gene_benchmark_summary.csv', index=False)

print(f'B56 qualifying genes = {b56_qualifying}/{len(b56)}')
print(f'Empirical P = {empirical_p:.6f}')
print(f'Exact matched probability = {exact_p:.6f}')
