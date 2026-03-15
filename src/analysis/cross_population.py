"""
Cross-population validation of the protective vs risk pathway distinction.

Tests whether the pathway architecture difference found in FinnGen
(European/Finnish) replicates in an independent East Asian population (BBJ).

This is critical for the paper: if the finding is population-specific,
it may reflect population structure rather than biology. If it replicates,
the finding is likely universal.
"""

from pathlib import Path

import pandas as pd

from src.analysis.compare_pathways import (
    compute_pathway_overlap,
    fishers_exact_test,
    per_source_breakdown,
    rank_correlation,
)
from src.analysis.enrichment import run_enrichment_both


def run_cross_population_enrichment(
    bbj_protective_genes: list[str],
    bbj_risk_genes: list[str],
    background: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run pathway enrichment on BBJ gene sets."""
    print(f"\n  BBJ enrichment: {len(bbj_protective_genes)} protective, "
          f"{len(bbj_risk_genes)} risk genes")

    prot_enrichment, risk_enrichment = run_enrichment_both(
        bbj_protective_genes,
        bbj_risk_genes,
        background=background,
    )

    return prot_enrichment, risk_enrichment


def compare_enrichment_patterns(
    finngen_prot_enrichment: pd.DataFrame,
    finngen_risk_enrichment: pd.DataFrame,
    bbj_prot_enrichment: pd.DataFrame,
    bbj_risk_enrichment: pd.DataFrame,
    output_dir: str = "data/processed",
) -> dict:
    """
    Compare the pathway enrichment patterns between FinnGen and BBJ.

    Key questions:
    1. Does BBJ also show distinct protective vs risk pathways?
    2. Are the protective pathways similar between FinnGen and BBJ?
    3. Are the risk pathways similar between FinnGen and BBJ?
    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("  CROSS-POPULATION PATHWAY COMPARISON")
    print("=" * 60)

    results = {}

    # 1. Does BBJ show the same protective/risk distinction?
    print("\n--- BBJ internal: protective vs risk ---")
    bbj_overlap = compute_pathway_overlap(bbj_prot_enrichment, bbj_risk_enrichment)
    bbj_fisher = fishers_exact_test(bbj_prot_enrichment, bbj_risk_enrichment)

    results["bbj_jaccard"] = bbj_overlap["jaccard_similarity"]
    results["bbj_fisher_p"] = bbj_fisher["p_value"]

    print(f"  BBJ Jaccard (prot vs risk): {bbj_overlap['jaccard_similarity']:.3f}")
    print(f"  BBJ Fisher's p-value:       {bbj_fisher['p_value']:.2e}")
    print(f"  → {bbj_fisher['interpretation']}")

    # 2. Are protective pathways consistent across populations?
    print("\n--- Cross-population: FinnGen protective vs BBJ protective ---")
    prot_cross = compute_pathway_overlap(finngen_prot_enrichment, bbj_prot_enrichment)
    prot_rank = rank_correlation(finngen_prot_enrichment, bbj_prot_enrichment)

    results["protective_cross_jaccard"] = prot_cross["jaccard_similarity"]
    results["protective_cross_spearman"] = prot_rank.get("spearman_rho")

    print(f"  Jaccard (FinnGen prot ∩ BBJ prot): {prot_cross['jaccard_similarity']:.3f}")
    if prot_rank.get("spearman_rho") is not None:
        print(f"  Spearman rho: {prot_rank['spearman_rho']:.3f} "
              f"(p={prot_rank['p_value']:.2e})")

    # 3. Are risk pathways consistent across populations?
    print("\n--- Cross-population: FinnGen risk vs BBJ risk ---")
    risk_cross = compute_pathway_overlap(finngen_risk_enrichment, bbj_risk_enrichment)
    risk_rank = rank_correlation(finngen_risk_enrichment, bbj_risk_enrichment)

    results["risk_cross_jaccard"] = risk_cross["jaccard_similarity"]
    results["risk_cross_spearman"] = risk_rank.get("spearman_rho")

    print(f"  Jaccard (FinnGen risk ∩ BBJ risk): {risk_cross['jaccard_similarity']:.3f}")
    if risk_rank.get("spearman_rho") is not None:
        print(f"  Spearman rho: {risk_rank['spearman_rho']:.3f} "
              f"(p={risk_rank['p_value']:.2e})")

    # 4. Verify: are cross-set comparisons LESS similar?
    print("\n--- Control: FinnGen protective vs BBJ risk (should be LOW) ---")
    cross_set = compute_pathway_overlap(finngen_prot_enrichment, bbj_risk_enrichment)
    results["control_cross_jaccard"] = cross_set["jaccard_similarity"]
    print(f"  Jaccard (FinnGen prot ∩ BBJ risk): {cross_set['jaccard_similarity']:.3f}")

    # Summary
    print("\n" + "=" * 60)
    print("  REPLICATION SUMMARY")
    print("=" * 60)

    replicates = (
        results.get("bbj_jaccard", 1.0) < 0.5
        and results.get("bbj_fisher_p", 1.0) < 0.05
    )
    if replicates:
        print("  The protective/risk pathway distinction REPLICATES in BBJ")
        print("  → Finding is likely universal, not population-specific")
    else:
        print("  Replication is partial or unclear")
        print("  → Further investigation needed")

    same_direction = (
        results.get("protective_cross_jaccard", 0) > results.get("control_cross_jaccard", 1)
        and results.get("risk_cross_jaccard", 0) > results.get("control_cross_jaccard", 1)
    )
    if same_direction:
        print("  Same-type pathways are more similar across populations")
        print("  than cross-type pathways → consistent biology")

    # Save
    import json
    with open(out_path / "cross_population_results.json", "w") as f:
        json.dump({k: v for k, v in results.items()
                   if not isinstance(v, set)}, f, indent=2, default=str)

    return results
