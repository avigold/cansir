"""
Build ML-ready feature matrix from FinnGen + Open Targets data.

Labels come from FinnGen (real cohort data):
- Positive: genes with significant protective variants (negative beta)
  across multiple cancer types
- Negative: genes with significant risk variants (positive beta)
  across multiple cancer types
- Unlabeled: genes with pathway/function data but no strong GWAS signal

Features come from Open Targets (pathway, function, association data).
"""

from pathlib import Path

import numpy as np
import pandas as pd


def load_finngen_protective(path: str = "data/processed/finngen_protective_genes.csv") -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{path} not found. Run finngen.py first.")
    return pd.read_csv(p)


def load_finngen_risk(path: str = "data/processed/finngen_risk_genes.csv") -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{path} not found. Run finngen.py first.")
    return pd.read_csv(p)


def load_open_targets_features(path: str = "data/raw/open_targets_features.csv") -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        print(f"Warning: {path} not found. Run open_targets.py first.")
        return pd.DataFrame()
    return pd.read_csv(p)


def build_feature_matrix(
    protective_df: pd.DataFrame,
    risk_df: pd.DataFrame,
    ot_df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Build per-gene feature matrix with labels from FinnGen and features
    from both FinnGen summary stats and Open Targets.

    FinnGen provides the labels (based on real genotype-phenotype data):
    - Protective genes: significant negative beta across cancer types
    - Risk genes: significant positive beta across cancer types

    FinnGen also provides features:
    - n_protective/risk_variants, n_cancer_types, mean_beta, etc.

    Open Targets provides biological features:
    - Pathway memberships, function annotations, disease associations
    """
    # Prepare FinnGen protective gene features
    prot = protective_df.rename(columns={
        "n_protective_variants": "fg_n_protective_variants",
        "n_cancer_types": "fg_n_cancer_types_protective",
        "mean_beta": "fg_mean_beta_protective",
        "min_beta": "fg_min_beta",
        "best_pval": "fg_best_pval_protective",
        "mean_af_diff": "fg_mean_af_diff_protective",
        "mean_af_cases": "fg_mean_af_cases",
        "mean_af_controls": "fg_mean_af_controls",
        "protection_score": "fg_protection_score",
    })[["gene", "fg_n_protective_variants", "fg_n_cancer_types_protective",
        "fg_mean_beta_protective", "fg_min_beta", "fg_best_pval_protective",
        "fg_mean_af_diff_protective", "fg_mean_af_cases", "fg_mean_af_controls",
        "fg_protection_score"]]

    # Prepare FinnGen risk gene features
    risk_cols = {
        "n_risk_variants": "fg_n_risk_variants",
        "n_cancer_types": "fg_n_cancer_types_risk",
        "mean_beta": "fg_mean_beta_risk",
        "max_beta": "fg_max_beta",
        "best_pval": "fg_best_pval_risk",
        "mean_af_diff": "fg_mean_af_diff_risk",
        "risk_score": "fg_risk_score",
    }
    available_risk_cols = {k: v for k, v in risk_cols.items() if k in risk_df.columns}
    rsk = risk_df.rename(columns=available_risk_cols)
    rsk_keep = ["gene"] + list(available_risk_cols.values())
    rsk = rsk[[c for c in rsk_keep if c in rsk.columns]]

    # Merge protective + risk
    merged = pd.merge(prot, rsk, on="gene", how="outer")

    # Merge with Open Targets features
    if not ot_df.empty:
        merged = pd.merge(merged, ot_df, on="gene", how="outer")

    # --- Create labels from FinnGen data ---
    # This is the key improvement: labels are based on real cohort data
    labels = pd.Series(-1, index=merged.index, name="label")

    # Positive (protective): significant protective variants across 2+ cancer types
    # OR very strong protective effect in any single cancer type
    if "fg_n_cancer_types_protective" in merged.columns:
        is_protective = (
            (merged["fg_n_cancer_types_protective"] >= 2)
            & (merged["fg_mean_beta_protective"] < -0.05)
        ) | (
            merged["fg_min_beta"].notna()
            & (merged["fg_min_beta"] < -0.3)
        )
        labels[is_protective] = 1

    # Negative (risk): significant risk variants across 2+ cancer types
    if "fg_n_cancer_types_risk" in merged.columns:
        is_risk = (
            (merged["fg_n_cancer_types_risk"] >= 2)
            & (merged.get("fg_mean_beta_risk", pd.Series(0.0, index=merged.index)) > 0.05)
        ) | (
            merged.get("fg_max_beta", pd.Series(0.0, index=merged.index)) > 0.3
        )
        labels[is_risk] = 0

    # Genes that appear in both lists: label by dominant signal
    both = (labels == 1) & merged["fg_n_cancer_types_risk"].notna()
    for idx in merged[both].index:
        prot_score = merged.loc[idx].get("fg_protection_score", 0) or 0
        risk_score = merged.loc[idx].get("fg_risk_score", 0) or 0
        if risk_score > prot_score:
            labels[idx] = 0

    merged["label"] = labels

    # --- Define feature columns ---
    feature_cols = [
        # FinnGen features (from real cohort data)
        "fg_n_protective_variants",
        "fg_n_cancer_types_protective",
        "fg_mean_beta_protective",
        "fg_min_beta",
        "fg_best_pval_protective",
        "fg_mean_af_diff_protective",
        "fg_mean_af_cases",
        "fg_mean_af_controls",
        "fg_protection_score",
        "fg_n_risk_variants",
        "fg_n_cancer_types_risk",
        "fg_mean_beta_risk",
        "fg_max_beta",
        "fg_best_pval_risk",
        "fg_mean_af_diff_risk",
        "fg_risk_score",
        # Open Targets features (biological context)
        "n_cancer_associations",
        "max_cancer_assoc_score",
        "mean_cancer_assoc_score",
        "n_pathways",
        "in_dna_repair",
        "in_apoptosis",
        "in_cell_cycle",
        "in_immune",
        "in_signal_transduction",
        "func_mentions_repair",
        "func_mentions_tumor_suppressor",
        "func_mentions_apoptosis",
        "func_mentions_immune",
    ]

    available_features = [c for c in feature_cols if c in merged.columns]

    # Fill NaN appropriately
    binary_cols = [c for c in available_features if c.startswith(("in_", "func_"))]
    continuous_cols = [c for c in available_features if c not in binary_cols]

    for col in binary_cols:
        merged[col] = merged[col].fillna(0).astype(int)

    for col in continuous_cols:
        merged[col] = merged[col].fillna(0)

    return merged, available_features


def save_feature_matrix(
    df: pd.DataFrame,
    features: list[str],
    output_dir: str = "data/processed",
):
    """Save the feature matrix for model training."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path / "feature_matrix.csv", index=False)

    labels = df["label"]
    n_pos = (labels == 1).sum()
    n_neg = (labels == 0).sum()
    n_unlabeled = (labels == -1).sum()

    print(f"\nFeature matrix: {len(df)} genes, {len(features)} features")
    print(f"Labels from FinnGen cohort data:")
    print(f"  {n_pos} protective (negative beta, multi-cancer)")
    print(f"  {n_neg} risk (positive beta, multi-cancer)")
    print(f"  {n_unlabeled} unlabeled (to be predicted)")


if __name__ == "__main__":
    print("=== Building Feature Matrix (FinnGen + Open Targets) ===\n")

    protective_df = load_finngen_protective()
    risk_df = load_finngen_risk()
    ot_df = load_open_targets_features()

    merged, features = build_feature_matrix(protective_df, risk_df, ot_df)
    save_feature_matrix(merged, features)
