"""Exploratory GSE120932 B56 benchmark matched for probes, expression and variance.

Input: public GSE120932 series matrix and the study's complete GEO2R top table.
This script does not refit limma; use run_limma_GSE120932_scripted.R for that.
"""

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


PARENTAL = ("GSM3421746", "GSM3421747", "GSM3421748")
RESISTANT = ("GSM3421752", "GSM3421753", "GSM3421754")
B56 = ("PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E")
DEFAULT_POOL_SIZE = 200


def load_matrix(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as source:
        metadata = {}
        for line in source:
            if line.startswith("!Sample_title\t") or line.startswith("!Sample_geo_accession\t"):
                metadata[line.split("\t", 1)[0]] = next(csv.reader([line], delimiter="\t"))[1:]
            if line.startswith("!series_matrix_table_begin"):
                break
        matrix = pd.read_csv(source, sep="\t", low_memory=False)
    matrix = matrix[matrix.ID_REF.astype(str).str.startswith("ILMN_")].copy()
    matrix = matrix.set_index("ID_REF")
    accessions = metadata["!Sample_geo_accession"]
    titles = metadata["!Sample_title"]
    title_map = dict(zip(accessions, titles))
    expected = list(PARENTAL + RESISTANT)
    assert all(sample in matrix.columns for sample in expected)
    assert all("K562, replicate" in title_map[sample] for sample in PARENTAL)
    assert all("K562-IR w/o imatinib" in title_map[sample] for sample in RESISTANT)
    assert matrix.index.is_unique
    return matrix[expected].apply(pd.to_numeric), title_map


def gene_summary(table, matrix):
    required = {"ID", "Gene.symbol", "logFC", "adj.P.Val"}
    assert required.issubset(table.columns), required - set(table.columns)
    assert table.ID.is_unique and matrix.index.is_unique
    selected = matrix.loc[table.ID, list(PARENTAL + RESISTANT)]
    reconstructed = selected[list(RESISTANT)].mean(axis=1).to_numpy() - selected[list(PARENTAL)].mean(axis=1).to_numpy()
    error = np.abs(reconstructed - table.logFC.to_numpy())
    assert np.max(error) < 1e-5, f"GEO2R contrast does not match selected arrays: {max(error)}"

    annotated = table.copy()
    annotated["Gene.symbol"] = annotated["Gene.symbol"].fillna("").astype(str).str.strip()
    annotated = annotated[annotated["Gene.symbol"] != ""].copy()
    genes = []
    for gene, probes in annotated.groupby("Gene.symbol", sort=True):
        sample_values = matrix.loc[probes.ID, list(PARENTAL + RESISTANT)].mean(axis=0).to_numpy(dtype=float)
        variance = float(np.var(sample_values, ddof=1))
        # Pooled within-group variance excludes the difference between the
        # parental and resistant means from the matching covariate.
        residual_sum_squares = (((sample_values[:3] - sample_values[:3].mean()) ** 2).sum()
                                + ((sample_values[3:] - sample_values[3:].mean()) ** 2).sum())
        within_group_variance = float(residual_sum_squares / 4)
        genes.append({
            "gene": gene,
            "n_probes": len(probes),
            "baseline_log2_expression": float(sample_values[:3].mean()),
            "log2_variance_all_six": float(np.log2(max(variance, 1e-12))),
            "log2_variance_within_group": float(np.log2(max(within_group_variance, 1e-12))),
            "qualifies": bool(((probes.logFC < 0) & (probes["adj.P.Val"] < 0.05)).all()),
        })
    return pd.DataFrame(genes).set_index("gene"), float(error.max())


def build_pools(genes, pool_size, variance_feature="log2_variance_all_six"):
    rows = []
    pools = {}
    for b56 in B56:
        target = genes.loc[b56]
        peers = genes[(genes.n_probes == target.n_probes) & (~genes.index.isin(B56))].copy()
        assert len(peers) >= pool_size, (b56, len(peers))
        features = ["baseline_log2_expression", variance_feature]
        std = peers[features].std(ddof=1)
        assert (std > 0).all()
        distance = np.sqrt((((peers[features] - target[features].astype(float)) / std) ** 2).sum(axis=1))
        peers["distance"] = distance
        peers = peers.sort_values(["distance"], kind="stable").head(pool_size)
        assert (peers.n_probes == target.n_probes).all()
        pools[b56] = peers
        for gene, row in peers.iterrows():
            rows.append({"target": b56, "candidate": gene, **row.to_dict()})
    return pools, pd.DataFrame(rows)


def simulate(pools, draws, seed):
    rng = np.random.default_rng(seed)
    names = [pools[target].index.to_numpy() for target in B56]
    successes = [pools[target].qualifies.to_numpy(dtype=int) for target in B56]
    pool_size = len(pools[B56[0]])
    sample_indices = rng.integers(0, pool_size, size=(draws, len(B56)))
    sampled_names = np.column_stack([names[j][sample_indices[:, j]] for j in range(len(B56))])
    # A panel must contain five distinct genes even if gene-specific pools overlap.
    for j in range(1, len(B56)):
        collisions = np.any(sampled_names[:, j, None] == sampled_names[:, :j], axis=1)
        while np.any(collisions):
            sample_indices[collisions, j] = rng.integers(0, pool_size, size=collisions.sum())
            sampled_names[collisions, j] = names[j][sample_indices[collisions, j]]
            collisions = np.any(sampled_names[:, j, None] == sampled_names[:, :j], axis=1)
    hit_counts = np.zeros(draws, dtype=int)
    for j in range(len(B56)):
        hit_counts += successes[j][sample_indices[:, j]]
    return hit_counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series_matrix", type=Path)
    parser.add_argument("geo2r_top_table", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--draws", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260922)
    parser.add_argument("--pool-size", type=int, default=DEFAULT_POOL_SIZE)
    parser.add_argument("--variance-mode", choices=("all-six", "within-group"), default="all-six",
                        help="Use all six arrays or pooled residual variance within the two groups")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix, titles = load_matrix(args.series_matrix)
    table = pd.read_csv(args.geo2r_top_table, sep="\t")
    genes, max_fc_difference = gene_summary(table, matrix)
    variance_feature = ("log2_variance_all_six" if args.variance_mode == "all-six"
                        else "log2_variance_within_group")
    pools, pool_table = build_pools(genes, args.pool_size, variance_feature)
    observed = int(genes.loc[list(B56), "qualifies"].sum())
    assert observed == 3, observed
    counts = simulate(pools, args.draws, args.seed)
    pvalue = (int((counts >= observed).sum()) + 1) / (args.draws + 1)

    prefix = ("GSE120932_expression_variance_matched" if args.variance_mode == "all-six"
              else "GSE120932_within_group_variance_matched")
    pool_table.to_csv(args.output_dir / f"{prefix}_pools.csv", index=False)
    pd.DataFrame({"qualifying_gene_count": np.arange(6),
                  "frequency": np.bincount(counts, minlength=6)}).to_csv(
        args.output_dir / f"{prefix}_distribution.csv", index=False)
    summary = pd.DataFrame([
        {"target": target, "observed_qualifies": bool(genes.loc[target, "qualifies"]),
         "target_probe_count": int(genes.loc[target, "n_probes"]),
         "candidate_pool_size": len(pools[target]),
         "pool_qualifying_rate": float(pools[target].qualifies.mean()),
         "baseline_log2_expression": float(genes.loc[target, "baseline_log2_expression"]),
         variance_feature: float(genes.loc[target, variance_feature])} for target in B56
    ])
    summary.to_csv(args.output_dir / f"{prefix}_targets.csv", index=False)
    with open(args.output_dir / f"{prefix}_summary.txt", "w") as out:
        out.write(f"selected_parental={','.join(PARENTAL)}\n")
        out.write(f"selected_resistant_drug_withdrawn={','.join(RESISTANT)}\n")
        out.write(f"matrix_probes={len(matrix)}\n")
        out.write(f"max_absolute_difference_from_GEO2R_logFC={max_fc_difference:.12g}\n")
        out.write("gene_expression=mean_log2_signal_over_mapped_probes_per_array\n")
        out.write(f"variance_mode={args.variance_mode}\n")
        out.write("variance=sample_variance_across_six_arrays_or_pooled_within_group_residual_variance\n")
        out.write("distance=Euclidean_over_baseline_mean_and_log2_variance_scaled_by_stratum_SD\n")
        out.write(f"pool_size_per_B56_gene={args.pool_size}\n")
        out.write("exclude_the_five_B56_genes_from_each_pool=true\n")
        out.write("draw_five_distinct_genes=true\n")
        out.write(f"seed={args.seed}\n")
        out.write(f"random_sets={args.draws}\n")
        out.write(f"observed_qualifying_B56_genes={observed}\n")
        out.write(f"random_sets_with_at_least_three={(counts >= observed).sum()}\n")
        out.write(f"empirical_P_one_count_correction={pvalue:.9f}\n")
    print(f"GSE120932: {len(matrix)} probes, GEO2R logFC max difference {max_fc_difference:.2g}")
    print(f"Expression/variance matched empirical P={pvalue:.6f} ({args.draws} panels; seed {args.seed})")
    print(summary[["target", "target_probe_count", "pool_qualifying_rate"]].to_string(index=False))


if __name__ == "__main__":
    main()
