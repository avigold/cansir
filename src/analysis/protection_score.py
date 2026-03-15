"""
Polygenic Protection Score (PPS) construction and analysis.

The PPS captures genetic predisposition to cancer PROTECTION, as distinct
from simply inverting a standard Polygenic Risk Score (PRS).

Key hypothesis: PPS and inverted PRS are NOT equivalent — they capture
different genetic architectures (detoxification/adhesion/apoptosis vs
the absence of transcription/replication/signaling risk variants).

Using summary statistics only (no individual-level data), we can:
1. Construct PPS weights from protective variant effects
2. Compare the variant composition of PPS vs inverted PRS
3. Show they tag different pathways
4. Estimate their correlation (expected to be moderate, not perfect)
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def build_pps_weights(
    protective_variants: pd.DataFrame,
    pval_threshold: float = 5e-6,
) -> pd.DataFrame:
    """
    Build PPS weights from protective variants.

    Each variant gets a weight = abs(beta) * -sign(beta)
    (positive weight = protective contribution).
    """
    df = protective_variants.copy()

    # Filter
    if "pval" in df.columns:
        df = df[df["pval"] <= pval_threshold]
    if "effect" in df.columns:
        df = df[df["effect"] == "protective"]

    if df.empty:
        return pd.DataFrame()

    weights = pd.DataFrame({
        "rsid": df["rsids"].values if "rsids" in df.columns else range(len(df)),
        "chrom": df["#chrom"].values if "#chrom" in df.columns else "",
        "pos": df["pos"].values if "pos" in df.columns else 0,
        "gene": df["nearest_genes"].values if "nearest_genes" in df.columns else "",
        "beta": df["beta"].values,
        "pval": df["pval"].values if "pval" in df.columns else np.nan,
        "af": df["af_alt"].values if "af_alt" in df.columns else np.nan,
        "pps_weight": np.abs(df["beta"].values),  # All protective, so weight = |beta|
    })

    return weights.sort_values("pps_weight", ascending=False)


def build_prs_weights(
    risk_variants: pd.DataFrame,
    pval_threshold: float = 5e-6,
) -> pd.DataFrame:
    """Build standard PRS weights from risk variants."""
    df = risk_variants.copy()

    if "pval" in df.columns:
        df = df[df["pval"] <= pval_threshold]
    if "effect" in df.columns:
        df = df[df["effect"] == "risk"]

    if df.empty:
        return pd.DataFrame()

    weights = pd.DataFrame({
        "rsid": df["rsids"].values if "rsids" in df.columns else range(len(df)),
        "chrom": df["#chrom"].values if "#chrom" in df.columns else "",
        "pos": df["pos"].values if "pos" in df.columns else 0,
        "gene": df["nearest_genes"].values if "nearest_genes" in df.columns else "",
        "beta": df["beta"].values,
        "pval": df["pval"].values if "pval" in df.columns else np.nan,
        "af": df["af_alt"].values if "af_alt" in df.columns else np.nan,
        "prs_weight": df["beta"].values,  # Positive beta = risk contribution
    })

    return weights.sort_values("prs_weight", ascending=False)


def compare_pps_vs_inverted_prs(
    pps_weights: pd.DataFrame,
    prs_weights: pd.DataFrame,
    output_dir: str = "data/processed",
) -> dict:
    """
    Compare PPS variant composition against inverted PRS.

    If PPS == inverted PRS, they should share the same variants
    (just with flipped signs). If PPS captures distinct biology,
    the variant sets should be largely different.
    """
    pps_variants = set(pps_weights["rsid"].values)
    prs_variants = set(prs_weights["rsid"].values)

    shared = pps_variants & prs_variants
    pps_only = pps_variants - prs_variants
    prs_only = prs_variants - pps_variants

    # Gene-level comparison
    pps_genes = set(pps_weights["gene"].dropna().unique())
    prs_genes = set(prs_weights["gene"].dropna().unique())
    shared_genes = pps_genes & prs_genes
    pps_only_genes = pps_genes - prs_genes
    prs_only_genes = prs_genes - pps_genes

    gene_jaccard = (
        len(shared_genes) / len(pps_genes | prs_genes)
        if (pps_genes | prs_genes) else 0
    )

    results = {
        "n_pps_variants": len(pps_variants),
        "n_prs_variants": len(prs_variants),
        "n_shared_variants": len(shared),
        "n_pps_only_variants": len(pps_only),
        "n_prs_only_variants": len(prs_only),
        "variant_jaccard": (
            len(shared) / len(pps_variants | prs_variants)
            if (pps_variants | prs_variants) else 0
        ),
        "n_pps_genes": len(pps_genes),
        "n_prs_genes": len(prs_genes),
        "n_shared_genes": len(shared_genes),
        "n_pps_only_genes": len(pps_only_genes),
        "n_prs_only_genes": len(prs_only_genes),
        "gene_jaccard": gene_jaccard,
    }

    print("\n=== PPS vs Inverted PRS Comparison ===")
    print(f"\nVariant level:")
    print(f"  PPS variants:     {results['n_pps_variants']}")
    print(f"  PRS variants:     {results['n_prs_variants']}")
    print(f"  Shared:           {results['n_shared_variants']}")
    print(f"  Variant Jaccard:  {results['variant_jaccard']:.3f}")

    print(f"\nGene level:")
    print(f"  PPS genes:        {results['n_pps_genes']}")
    print(f"  PRS genes:        {results['n_prs_genes']}")
    print(f"  Shared genes:     {results['n_shared_genes']}")
    print(f"  Gene Jaccard:     {gene_jaccard:.3f}")

    if gene_jaccard < 0.3:
        results["interpretation"] = (
            "PPS and PRS tag LARGELY DIFFERENT gene sets. "
            "Protection is not the inverse of risk."
        )
    elif gene_jaccard > 0.7:
        results["interpretation"] = (
            "PPS and PRS share most genes. "
            "Protection may be the inverse of risk."
        )
    else:
        results["interpretation"] = (
            "PPS and PRS show PARTIAL overlap. "
            "Protection has both shared and unique genetic components."
        )

    print(f"\n  → {results['interpretation']}")

    # Save
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pps_vs_prs_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    return results


def pathway_composition_comparison(
    pps_weights: pd.DataFrame,
    prs_weights: pd.DataFrame,
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    output_dir: str = "data/processed",
) -> pd.DataFrame:
    """
    Compare which pathways are represented in PPS vs PRS.

    Uses the enrichment results to annotate whether PPS-unique genes
    fall in protective-unique pathways.
    """
    pps_genes = set(pps_weights["gene"].dropna().unique())
    prs_genes = set(prs_weights["gene"].dropna().unique())

    pps_only_genes = pps_genes - prs_genes
    prs_only_genes = prs_genes - pps_genes

    # Get pathway annotations
    prot_terms = set(prot_enrichment["native"]) if not prot_enrichment.empty else set()
    risk_terms = set(risk_enrichment["native"]) if not risk_enrichment.empty else set()

    prot_only_terms = prot_terms - risk_terms
    risk_only_terms = risk_terms - prot_terms

    summary = {
        "pps_only_genes": len(pps_only_genes),
        "prs_only_genes": len(prs_only_genes),
        "prot_only_pathways": len(prot_only_terms),
        "risk_only_pathways": len(risk_only_terms),
    }

    print(f"\n=== PPS/PRS Pathway Composition ===")
    print(f"  PPS-unique genes: {len(pps_only_genes)}")
    print(f"  PRS-unique genes: {len(prs_only_genes)}")
    print(f"  Protective-unique pathways: {len(prot_only_terms)}")
    print(f"  Risk-unique pathways:       {len(risk_only_terms)}")

    return pd.DataFrame([summary])


def estimate_pps_prs_correlation(
    all_variants: pd.DataFrame,
    output_dir: str = "data/processed",
) -> dict:
    """
    Estimate the expected correlation between PPS and inverted PRS
    using variant-level effect sizes.

    If PPS == -PRS, correlation should be ~1.0.
    If they're independent, correlation should be ~0.
    Moderate correlation means partially overlapping but distinct.
    """
    if all_variants.empty:
        return {"correlation": None, "note": "No data"}

    # For each variant, compute its contribution to PPS and to -PRS
    df = all_variants.copy()

    # PPS contribution: |beta| for protective variants, 0 for risk
    df["pps_contrib"] = np.where(df["beta"] < 0, np.abs(df["beta"]), 0)

    # Inverted PRS contribution: |beta| for risk variants (inverted = protective),
    # 0 for protective
    df["inv_prs_contrib"] = np.where(df["beta"] > 0, np.abs(df["beta"]), 0)

    # At the variant level, these are by definition non-overlapping
    # (a variant is either protective or risk, not both)
    # So the question is: at the GENE level, do the same genes contribute to both?

    if "nearest_genes" not in df.columns:
        return {"correlation": None, "note": "No gene column"}

    gene_level = df.groupby("nearest_genes").agg(
        pps_total=("pps_contrib", "sum"),
        inv_prs_total=("inv_prs_contrib", "sum"),
    ).reset_index()

    # Only genes with nonzero contribution to at least one score
    gene_level = gene_level[
        (gene_level["pps_total"] > 0) | (gene_level["inv_prs_total"] > 0)
    ]

    if len(gene_level) < 10:
        return {"correlation": None, "note": "Too few genes"}

    rho, p_val = stats.spearmanr(gene_level["pps_total"], gene_level["inv_prs_total"])

    results = {
        "gene_level_spearman_rho": rho,
        "p_value": p_val,
        "n_genes": len(gene_level),
        "n_pps_only_genes": (
            (gene_level["pps_total"] > 0) & (gene_level["inv_prs_total"] == 0)
        ).sum(),
        "n_inv_prs_only_genes": (
            (gene_level["pps_total"] == 0) & (gene_level["inv_prs_total"] > 0)
        ).sum(),
        "n_both": (
            (gene_level["pps_total"] > 0) & (gene_level["inv_prs_total"] > 0)
        ).sum(),
    }

    print(f"\n=== PPS vs Inverted PRS: Gene-Level Correlation ===")
    print(f"  Spearman rho: {rho:.3f} (p={p_val:.2e})")
    print(f"  Genes in PPS only:      {results['n_pps_only_genes']}")
    print(f"  Genes in inv-PRS only:  {results['n_inv_prs_only_genes']}")
    print(f"  Genes in both:          {results['n_both']}")

    if abs(rho) < 0.3:
        print(f"  → PPS and inverted PRS are largely INDEPENDENT")
        results["interpretation"] = "Independent — protection ≠ inverse of risk"
    elif abs(rho) > 0.7:
        print(f"  → PPS and inverted PRS are strongly CORRELATED")
        results["interpretation"] = "Correlated — protection ≈ inverse of risk"
    else:
        print(f"  → PPS and inverted PRS are PARTIALLY correlated")
        results["interpretation"] = "Partial overlap — distinct but related"

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pps_prs_correlation.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    return results


if __name__ == "__main__":
    print("=== Polygenic Protection Score Analysis ===\n")

    # Load FinnGen data
    variants = pd.read_csv("data/raw/finngen_significant_variants.csv",
                           dtype={"#chrom": str})

    protective = variants[variants["effect"] == "protective"]
    risk = variants[variants["effect"] == "risk"]

    # Build weights
    pps = build_pps_weights(protective)
    prs = build_prs_weights(risk)

    print(f"PPS: {len(pps)} variants")
    print(f"PRS: {len(prs)} variants")

    # Compare
    compare_pps_vs_inverted_prs(pps, prs)

    # Correlation estimate
    estimate_pps_prs_correlation(variants)
