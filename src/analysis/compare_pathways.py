"""
Statistical comparison of protective vs risk pathway enrichment.

Tests the central hypothesis: the pathway architecture of cancer
protection is distinct from cancer risk.

Methods:
- Jaccard similarity of enriched term sets
- Fisher's exact test on 2x2 contingency table
- Spearman rank correlation of enrichment p-values
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def compute_pathway_overlap(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
) -> dict:
    """
    Compute overlap statistics between protective and risk enriched term sets.

    Returns dict with Jaccard, Fisher's exact results, term breakdowns.
    """
    prot_terms = set(prot_enrichment["native"].values) if not prot_enrichment.empty else set()
    risk_terms = set(risk_enrichment["native"].values) if not risk_enrichment.empty else set()

    shared = prot_terms & risk_terms
    prot_only = prot_terms - risk_terms
    risk_only = risk_terms - prot_terms
    union = prot_terms | risk_terms

    jaccard = len(shared) / len(union) if union else 0.0

    return {
        "n_protective_terms": len(prot_terms),
        "n_risk_terms": len(risk_terms),
        "n_shared": len(shared),
        "n_protective_only": len(prot_only),
        "n_risk_only": len(risk_only),
        "jaccard_similarity": jaccard,
        "shared_terms": shared,
        "protective_only_terms": prot_only,
        "risk_only_terms": risk_only,
    }


def fishers_exact_test(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    all_tested_terms: set[str] | None = None,
) -> dict:
    """
    Fisher's exact test on pathway membership.

    2x2 table:
                        Enriched in risk    Not enriched in risk
    Enriched in prot        a                   b
    Not enriched in prot    c                   d

    H0: pathway membership is independent of protective/risk status.
    A significant p-value means the two gene sets enrich DIFFERENT pathways.
    """
    prot_terms = set(prot_enrichment["native"].values) if not prot_enrichment.empty else set()
    risk_terms = set(risk_enrichment["native"].values) if not risk_enrichment.empty else set()

    if all_tested_terms is None:
        all_tested_terms = prot_terms | risk_terms

    a = len(prot_terms & risk_terms)                         # both
    b = len(prot_terms - risk_terms)                         # prot only
    c = len(risk_terms - prot_terms)                         # risk only
    d = len(all_tested_terms - prot_terms - risk_terms)      # neither

    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")

    return {
        "contingency_table": table,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "interpretation": (
            "Protective and risk gene sets enrich DISTINCT pathways"
            if p_value < 0.05
            else "No significant difference in pathway enrichment"
        ),
    }


def rank_correlation(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
) -> dict:
    """
    Spearman rank correlation of enrichment p-values for shared terms.

    If protection is simply the inverse of risk, we'd expect a strong
    positive correlation (same pathways, similar significance).
    Low or negative correlation supports distinctness.
    """
    if prot_enrichment.empty or risk_enrichment.empty:
        return {"spearman_rho": None, "p_value": None, "n_shared": 0}

    # Merge on term ID
    merged = pd.merge(
        prot_enrichment[["native", "p_value"]],
        risk_enrichment[["native", "p_value"]],
        on="native",
        suffixes=("_prot", "_risk"),
    )

    if len(merged) < 3:
        return {
            "spearman_rho": None,
            "p_value": None,
            "n_shared": len(merged),
            "note": "Too few shared terms for correlation",
        }

    # Use -log10(p) for ranking
    merged["mlog_prot"] = -np.log10(merged["p_value_prot"].clip(lower=1e-300))
    merged["mlog_risk"] = -np.log10(merged["p_value_risk"].clip(lower=1e-300))

    rho, p_val = stats.spearmanr(merged["mlog_prot"], merged["mlog_risk"])

    return {
        "spearman_rho": rho,
        "p_value": p_val,
        "n_shared": len(merged),
        "interpretation": (
            "Strong positive correlation — pathways are similar"
            if rho > 0.5
            else "Weak/no correlation — pathways are distinct"
            if rho < 0.3
            else "Moderate correlation"
        ),
    }


def per_source_breakdown(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
) -> pd.DataFrame:
    """
    Break down overlap statistics by annotation source (GO:BP, Reactome, etc).
    """
    rows = []

    sources = set()
    if not prot_enrichment.empty:
        sources.update(prot_enrichment["source"].unique())
    if not risk_enrichment.empty:
        sources.update(risk_enrichment["source"].unique())

    for source in sorted(sources):
        prot_src = prot_enrichment[prot_enrichment["source"] == source] if not prot_enrichment.empty else pd.DataFrame()
        risk_src = risk_enrichment[risk_enrichment["source"] == source] if not risk_enrichment.empty else pd.DataFrame()

        overlap = compute_pathway_overlap(prot_src, risk_src)
        rows.append({
            "source": source,
            "n_protective": overlap["n_protective_terms"],
            "n_risk": overlap["n_risk_terms"],
            "n_shared": overlap["n_shared"],
            "n_protective_only": overlap["n_protective_only"],
            "n_risk_only": overlap["n_risk_only"],
            "jaccard": overlap["jaccard_similarity"],
        })

    return pd.DataFrame(rows)


def full_comparison(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    output_dir: str = "data/processed",
) -> dict:
    """
    Run all comparison analyses and save results.
    """
    print("\n=== Pathway Comparison: Protective vs Risk ===")

    # 1. Overlap
    overlap = compute_pathway_overlap(prot_enrichment, risk_enrichment)
    print(f"\nOverlap:")
    print(f"  Protective-only pathways: {overlap['n_protective_only']}")
    print(f"  Risk-only pathways:       {overlap['n_risk_only']}")
    print(f"  Shared pathways:          {overlap['n_shared']}")
    print(f"  Jaccard similarity:       {overlap['jaccard_similarity']:.3f}")

    # 2. Fisher's exact test
    fisher = fishers_exact_test(prot_enrichment, risk_enrichment)
    print(f"\nFisher's exact test:")
    print(f"  Odds ratio: {fisher['odds_ratio']:.3f}")
    print(f"  P-value:    {fisher['p_value']:.2e}")
    print(f"  → {fisher['interpretation']}")

    # 3. Rank correlation
    rank_corr = rank_correlation(prot_enrichment, risk_enrichment)
    print(f"\nSpearman rank correlation (shared terms):")
    if rank_corr["spearman_rho"] is not None:
        print(f"  rho:     {rank_corr['spearman_rho']:.3f}")
        print(f"  P-value: {rank_corr['p_value']:.2e}")
        print(f"  → {rank_corr['interpretation']}")
    else:
        print(f"  {rank_corr.get('note', 'Insufficient data')}")

    # 4. Per-source breakdown
    breakdown = per_source_breakdown(prot_enrichment, risk_enrichment)
    print(f"\nPer-source breakdown:")
    print(breakdown.to_string(index=False))

    # Save
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    breakdown.to_csv(out_path / "pathway_comparison_by_source.csv", index=False)

    # Save detailed term lists
    _save_term_details(prot_enrichment, risk_enrichment, overlap, out_path)

    return {
        "overlap": overlap,
        "fisher": fisher,
        "rank_correlation": rank_corr,
        "breakdown": breakdown,
    }


def _save_term_details(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    overlap: dict,
    out_path: Path,
):
    """Save lists of protective-only, risk-only, and shared pathways."""
    # Protective-only
    if not prot_enrichment.empty:
        prot_only_df = prot_enrichment[
            prot_enrichment["native"].isin(overlap["protective_only_terms"])
        ].sort_values("p_value")
        prot_only_df.to_csv(out_path / "pathways_protective_only.csv", index=False)

    # Risk-only
    if not risk_enrichment.empty:
        risk_only_df = risk_enrichment[
            risk_enrichment["native"].isin(overlap["risk_only_terms"])
        ].sort_values("p_value")
        risk_only_df.to_csv(out_path / "pathways_risk_only.csv", index=False)

    # Shared — merge to show p-values from both sides
    if not prot_enrichment.empty and not risk_enrichment.empty:
        shared_prot = prot_enrichment[
            prot_enrichment["native"].isin(overlap["shared_terms"])
        ][["native", "name", "source", "p_value", "intersection_size"]].rename(
            columns={"p_value": "p_value_protective", "intersection_size": "genes_protective"}
        )
        shared_risk = risk_enrichment[
            risk_enrichment["native"].isin(overlap["shared_terms"])
        ][["native", "p_value", "intersection_size"]].rename(
            columns={"p_value": "p_value_risk", "intersection_size": "genes_risk"}
        )
        shared = pd.merge(shared_prot, shared_risk, on="native")
        shared.to_csv(out_path / "pathways_shared.csv", index=False)
