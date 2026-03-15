"""
Pathway enrichment analysis using g:Profiler.

Runs Gene Ontology (BP, MF, CC), Reactome, KEGG, and WikiPathways
enrichment on protective and risk gene sets separately.

The background gene set is all genes appearing in FinnGen summary
statistics (before p-value filtering), ensuring unbiased enrichment.
"""

from gprofiler import GProfiler

import pandas as pd


SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"]

gp = GProfiler(return_dataframe=True)


def run_enrichment(
    gene_list: list[str],
    background: list[str] | None = None,
    sources: list[str] = SOURCES,
    significance_threshold: float = 0.05,
    organism: str = "hsapiens",
) -> pd.DataFrame:
    """
    Run pathway/GO enrichment for a gene list via g:Profiler.

    Args:
        gene_list: Gene symbols to test for enrichment.
        background: Background gene set (all testable genes).
            If None, g:Profiler uses all annotated genes.
        sources: Annotation sources to query.
        significance_threshold: Adjusted p-value cutoff.
        organism: Organism code.

    Returns:
        DataFrame of significantly enriched terms.
    """
    kwargs = {
        "organism": organism,
        "query": gene_list,
        "sources": sources,
        "significance_threshold_method": "g_SCS",  # g:Profiler's correction
        "user_threshold": significance_threshold,
        "no_evidences": False,
    }
    if background:
        kwargs["background"] = background

    result = gp.profile(**kwargs)

    if result.empty:
        return result

    # Clean up columns
    keep_cols = [
        "source", "native", "name", "p_value", "significant",
        "description", "term_size", "query_size", "intersection_size",
        "precision", "recall", "effective_domain_size", "intersections",
    ]
    available = [c for c in keep_cols if c in result.columns]
    result = result[available].copy()

    # Add gene ratio
    if "intersection_size" in result.columns and "query_size" in result.columns:
        result["gene_ratio"] = result["intersection_size"] / result["query_size"]

    return result.sort_values("p_value")


def run_enrichment_both(
    protective_genes: list[str],
    risk_genes: list[str],
    background: list[str] | None = None,
    sources: list[str] = SOURCES,
    significance_threshold: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run enrichment for both protective and risk gene sets.

    Returns:
        (protective_enrichment, risk_enrichment) DataFrames.
    """
    print(f"  Running enrichment for {len(protective_genes)} protective genes...")
    prot_enrichment = run_enrichment(
        protective_genes, background, sources, significance_threshold
    )
    print(f"    {len(prot_enrichment)} enriched terms")

    print(f"  Running enrichment for {len(risk_genes)} risk genes...")
    risk_enrichment = run_enrichment(
        risk_genes, background, sources, significance_threshold
    )
    print(f"    {len(risk_enrichment)} enriched terms")

    return prot_enrichment, risk_enrichment


def sensitivity_analysis(
    protective_genes_by_threshold: dict[str, list[str]],
    risk_genes_by_threshold: dict[str, list[str]],
    background: list[str] | None = None,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """
    Run enrichment at multiple p-value thresholds to test robustness.

    Args:
        protective_genes_by_threshold: {threshold_label: gene_list}
        risk_genes_by_threshold: {threshold_label: gene_list}

    Returns:
        {threshold_label: (prot_enrichment, risk_enrichment)}
    """
    results = {}
    for label in protective_genes_by_threshold:
        print(f"\n  Threshold: {label}")
        prot_e, risk_e = run_enrichment_both(
            protective_genes_by_threshold[label],
            risk_genes_by_threshold[label],
            background,
        )
        results[label] = (prot_e, risk_e)
    return results
