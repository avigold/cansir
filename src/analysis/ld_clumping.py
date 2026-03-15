"""
LD-aware clumping of GWAS variants.

Ensures that one physical locus contributes one gene to the enrichment
analysis, preventing inflated gene counts from correlated variants.

Uses distance-based clumping (no reference panel needed) as a practical
approximation sufficient for pathway enrichment analysis.
"""

import pandas as pd


def clump_by_distance(
    df: pd.DataFrame,
    window_kb: int = 500,
    pval_col: str = "pval",
    chrom_col: str = "#chrom",
    pos_col: str = "pos",
) -> pd.DataFrame:
    """
    Distance-based LD clumping: within each window, keep the variant
    with the smallest p-value.

    Args:
        df: DataFrame with at minimum chrom, pos, and pval columns.
        window_kb: Clumping window in kilobases.
        pval_col: Column name for p-values.
        chrom_col: Column name for chromosome.
        pos_col: Column name for position.

    Returns:
        DataFrame with one representative variant per locus.
    """
    if df.empty:
        return df

    window_bp = window_kb * 1000
    clumped_rows = []

    for chrom, chrom_df in df.groupby(chrom_col):
        sorted_df = chrom_df.sort_values(pval_col).copy()
        used = set()

        for idx, row in sorted_df.iterrows():
            if idx in used:
                continue

            # This variant is the lead SNP for its locus
            clumped_rows.append(row)

            # Mark all variants within the window as used
            pos = row[pos_col]
            nearby = sorted_df[
                (sorted_df[pos_col] >= pos - window_bp)
                & (sorted_df[pos_col] <= pos + window_bp)
            ].index
            used.update(nearby)

    return pd.DataFrame(clumped_rows).reset_index(drop=True)
