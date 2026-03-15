#!/usr/bin/env python3
"""
CANSIR: Cancer-Protective Gene Discovery Pipeline

Tests the hypothesis: "Genetic protection against cancer is biologically
distinct from genetic risk — not simply its inverse."

Data source: FinnGen R12 (500,348 individuals, linked health registry)
Method: Pathway enrichment comparison of protective vs risk gene sets

Pipeline:
1. Download FinnGen GWAS summary stats for 19 cancer endpoints
2. LD-clump and separate protective (beta<0) from risk (beta>0) variants
3. Map to gene sets with proper background
4. Run pathway enrichment (GO, Reactome, KEGG) on each set via g:Profiler
5. Statistically compare the two pathway profiles
6. Generate publication figures
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def step1_finngen_data():
    """Download FinnGen summary stats, clump, and extract gene sets."""
    print("\n" + "=" * 60)
    print("STEP 1: FinnGen R12 data")
    print("  Downloading 19 cancer endpoint GWAS summary statistics")
    print("  Collecting background gene set for enrichment analysis")
    print("=" * 60)

    from src.data.finngen import (
        clump_variants,
        expand_gene_names,
        extract_protective_genes,
        extract_risk_genes,
        fetch_all_cancer_sumstats,
    )

    # Download and collect background genes
    df, background_genes = fetch_all_cancer_sumstats(collect_background=True)
    if df.empty:
        print("ERROR: No FinnGen data retrieved.")
        sys.exit(1)

    # LD clump to avoid counting correlated variants multiple times
    df = clump_variants(df)

    # Expand multi-gene entries
    df = expand_gene_names(df)

    # Extract gene sets
    protective = extract_protective_genes(df)
    protective.to_csv("data/processed/finngen_protective_genes.csv", index=False)

    risk = extract_risk_genes(df)
    risk.to_csv("data/processed/finngen_risk_genes.csv", index=False)

    print(f"\n  Protective genes: {len(protective)}")
    print(f"  Risk genes:       {len(risk)}")
    print(f"  Background genes: {len(background_genes)}")

    return protective, risk, background_genes


def step2_enrichment(protective_genes, risk_genes, background_genes):
    """Run pathway enrichment on both gene sets."""
    print("\n" + "=" * 60)
    print("STEP 2: Pathway enrichment analysis")
    print("  Sources: GO (BP, MF, CC), Reactome, KEGG, WikiPathways")
    print("  Method: g:Profiler with g:SCS multiple testing correction")
    print("  Background: all genes in FinnGen summary statistics")
    print("=" * 60)

    from src.analysis.enrichment import run_enrichment_both

    prot_genes = protective_genes["gene"].tolist()
    risk_genes_list = risk_genes["gene"].tolist()
    bg = list(background_genes) if background_genes else None

    prot_enrichment, risk_enrichment = run_enrichment_both(
        prot_genes, risk_genes_list, background=bg
    )

    # Save enrichment results
    Path("data/processed").mkdir(parents=True, exist_ok=True)
    prot_enrichment.to_csv("data/processed/enrichment_protective.csv", index=False)
    risk_enrichment.to_csv("data/processed/enrichment_risk.csv", index=False)

    return prot_enrichment, risk_enrichment


def step3_comparison(prot_enrichment, risk_enrichment):
    """Statistical comparison of pathway profiles."""
    print("\n" + "=" * 60)
    print("STEP 3: Statistical comparison")
    print("  Testing: are protective and risk pathways distinct?")
    print("  Methods: Jaccard similarity, Fisher's exact, Spearman rank")
    print("=" * 60)

    from src.analysis.compare_pathways import full_comparison

    results = full_comparison(prot_enrichment, risk_enrichment)

    # Save summary
    summary = {
        "jaccard_similarity": results["overlap"]["jaccard_similarity"],
        "n_protective_only": results["overlap"]["n_protective_only"],
        "n_risk_only": results["overlap"]["n_risk_only"],
        "n_shared": results["overlap"]["n_shared"],
        "fisher_odds_ratio": results["fisher"]["odds_ratio"],
        "fisher_p_value": results["fisher"]["p_value"],
        "fisher_interpretation": results["fisher"]["interpretation"],
        "spearman_rho": results["rank_correlation"]["spearman_rho"],
        "spearman_p_value": results["rank_correlation"]["p_value"],
    }
    with open("data/processed/comparison_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    return results


def step4_figures(prot_enrichment, risk_enrichment, comparison_results):
    """Generate publication-quality figures."""
    print("\n" + "=" * 60)
    print("STEP 4: Generating figures")
    print("=" * 60)

    from src.analysis.visualization import generate_all_figures

    generate_all_figures(
        prot_enrichment,
        risk_enrichment,
        comparison_results["overlap"],
        comparison_results["breakdown"],
    )


def step5_sensitivity(background_genes):
    """Run enrichment at multiple significance thresholds."""
    print("\n" + "=" * 60)
    print("STEP 5: Sensitivity analysis")
    print("  Testing robustness across p-value thresholds")
    print("=" * 60)

    import pandas as pd

    from src.analysis.compare_pathways import compute_pathway_overlap
    from src.analysis.enrichment import run_enrichment_both

    # Load the raw significant variants (already downloaded)
    raw_path = Path("data/raw/finngen_significant_variants.csv")
    if not raw_path.exists():
        print("  Skipping: raw variants file not found")
        return

    df = pd.read_csv(raw_path, dtype={"#chrom": str})

    bg = list(background_genes) if background_genes else None

    thresholds = {"5e-8": 5e-8, "5e-7": 5e-7, "5e-6": 5e-6}
    sensitivity_results = []

    for label, threshold in thresholds.items():
        filtered = df[df["pval"] <= threshold]
        prot = filtered[filtered["effect"] == "protective"]["nearest_genes"].dropna().unique().tolist()
        risk = filtered[filtered["effect"] == "risk"]["nearest_genes"].dropna().unique().tolist()

        if len(prot) < 5 or len(risk) < 5:
            print(f"  {label}: too few genes (prot={len(prot)}, risk={len(risk)})")
            continue

        print(f"\n  Threshold {label}: {len(prot)} protective, {len(risk)} risk genes")
        prot_e, risk_e = run_enrichment_both(prot, risk, background=bg)
        overlap = compute_pathway_overlap(prot_e, risk_e)

        sensitivity_results.append({
            "threshold": label,
            "n_protective_genes": len(prot),
            "n_risk_genes": len(risk),
            "n_protective_pathways": overlap["n_protective_terms"],
            "n_risk_pathways": overlap["n_risk_terms"],
            "n_shared": overlap["n_shared"],
            "jaccard": overlap["jaccard_similarity"],
        })

    if sensitivity_results:
        sens_df = pd.DataFrame(sensitivity_results)
        sens_df.to_csv("data/processed/sensitivity_analysis.csv", index=False)
        print(f"\n  Sensitivity results:")
        print(sens_df.to_string(index=False))


def step6_bbj_validation(finngen_prot_enrichment, finngen_risk_enrichment):
    """Cross-population validation using BioBank Japan."""
    print("\n" + "=" * 60)
    print("STEP 6: Cross-population validation (BioBank Japan)")
    print("  Testing if pathway distinction replicates in East Asians")
    print("=" * 60)

    from src.data.bbj import (
        extract_bbj_gene_sets,
        fetch_bbj_cancer_sumstats,
    )

    bbj_df = fetch_bbj_cancer_sumstats()

    if bbj_df.empty:
        print("\n  BBJ data unavailable (NBDC server may be down).")
        print("  Skipping cross-population validation.")
        print("  Re-run later when https://humandbs.dbcls.jp is accessible.")
        return

    bbj_prot_genes, bbj_risk_genes = extract_bbj_gene_sets(bbj_df)

    if not bbj_prot_genes or not bbj_risk_genes:
        print("  Insufficient BBJ gene sets for enrichment analysis.")
        return

    from src.analysis.cross_population import (
        compare_enrichment_patterns,
        run_cross_population_enrichment,
    )

    bbj_prot_enrich, bbj_risk_enrich = run_cross_population_enrichment(
        bbj_prot_genes, bbj_risk_genes
    )

    compare_enrichment_patterns(
        finngen_prot_enrichment, finngen_risk_enrichment,
        bbj_prot_enrich, bbj_risk_enrich,
    )


def main():
    print("=" * 60)
    print("  CANSIR: Cancer-Protective Gene Discovery")
    print()
    print("  Hypothesis: Genetic protection against cancer is")
    print("  biologically distinct from genetic risk.")
    print()
    print("  Data: FinnGen R12 (N=500,348)")
    print("  Method: Pathway enrichment comparison")
    print("=" * 60)

    # Core analysis
    protective, risk, background = step1_finngen_data()
    prot_enrichment, risk_enrichment = step2_enrichment(
        protective, risk, background
    )
    comparison = step3_comparison(prot_enrichment, risk_enrichment)
    step4_figures(prot_enrichment, risk_enrichment, comparison)

    # Robustness check
    step5_sensitivity(background)

    # Cross-population validation
    step6_bbj_validation(prot_enrichment, risk_enrichment)

    # Final summary
    print("\n" + "=" * 60)
    print("  ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\n  Key result:")

    jaccard = comparison["overlap"]["jaccard_similarity"]
    fisher_p = comparison["fisher"]["p_value"]
    n_prot_only = comparison["overlap"]["n_protective_only"]
    n_risk_only = comparison["overlap"]["n_risk_only"]
    n_shared = comparison["overlap"]["n_shared"]

    print(f"    Jaccard similarity: {jaccard:.3f}")
    print(f"    Fisher's exact p-value: {fisher_p:.2e}")
    print(f"    Protective-only pathways: {n_prot_only}")
    print(f"    Risk-only pathways: {n_risk_only}")
    print(f"    Shared pathways: {n_shared}")

    if fisher_p < 0.05 and jaccard < 0.5:
        print(f"\n  → SUPPORTS HYPOTHESIS: Protective and risk gene sets")
        print(f"    enrich distinct biological pathways (p={fisher_p:.2e}).")
    elif jaccard > 0.7:
        print(f"\n  → REJECTS HYPOTHESIS: Protective and risk pathways")
        print(f"    substantially overlap (Jaccard={jaccard:.2f}).")
    else:
        print(f"\n  → PARTIAL SUPPORT: Significant distinction with moderate overlap.")

    print(f"\n  Output files:")
    for f in sorted(Path("data/processed").glob("fig*")):
        print(f"    {f}")
    print(f"    data/processed/comparison_summary.json")
    print(f"    data/processed/enrichment_protective.csv")
    print(f"    data/processed/enrichment_risk.csv")
    print("=" * 60)


if __name__ == "__main__":
    main()
