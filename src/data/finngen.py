"""
Fetch and process FinnGen R12 GWAS summary statistics for cancer endpoints.

FinnGen provides summary statistics from ~500k Finnish individuals with linked
health registry data. Each row is a variant with:
- beta: effect size (negative = protective against cancer)
- af_alt_cases vs af_alt_controls: allele frequency difference
- pval: statistical significance

This gives us REAL effect sizes from a REAL cohort — not just curated
literature results like the GWAS Catalog.
"""

import gzip
import io
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from tqdm import tqdm

BASE_URL = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"

# Major cancer endpoints with sufficient case counts (>1000)
# Using EXALLC variants (controls exclude other cancer patients) for cleaner signal
CANCER_ENDPOINTS = {
    "breast": "C3_BREAST_EXALLC",
    "prostate": "C3_PROSTATE_EXALLC",
    "colorectal": "C3_COLORECTAL_EXALLC",
    "lung": "C3_BRONCHUS_LUNG_EXALLC",
    "melanoma": "C3_MELANOMA_EXALLC",
    "bladder": "C3_BLADDER_EXALLC",
    "kidney": "C3_KIDNEY_NOTRENALPELVIS_EXALLC",
    "pancreas": "C3_PANCREAS_EXALLC",
    "stomach": "C3_STOMACH_EXALLC",
    "thyroid": "C3_THYROID_GLAND_EXALLC",
    "brain": "C3_BRAIN_EXALLC",
    "corpus_uteri": "C3_CORPUS_UTERI_EXALLC",
    "ovary": "C3_OVARY_EXALLC",
    "head_neck": "C3_HEAD_AND_NECK_EXALLC",
    "liver": "C3_HEPATOCELLU_CARC_EXALLC",
    "any_cancer": "C3_CANCER_EXALLC",
    "skin": "C3_SKIN_EXALLC",
    "colon": "C3_COLON_EXALLC",
    "non_small_cell_lung": "C3_LUNG_NONSMALL_EXALLC",
}

# Significance thresholds
PVAL_THRESHOLD = 5e-6  # suggestive significance (genome-wide is 5e-8)
BETA_PROTECTIVE_THRESHOLD = 0  # negative beta = protective


def download_and_filter(
    endpoint_id: str,
    pval_threshold: float = PVAL_THRESHOLD,
    collect_background_genes: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, set[str]]:
    """
    Download a FinnGen summary stats file and filter for significant variants.

    Args:
        endpoint_id: FinnGen endpoint ID.
        pval_threshold: P-value threshold for significance.
        collect_background_genes: If True, also return the set of ALL gene
            symbols in the file (before filtering). Needed for proper
            enrichment analysis background.

    Returns:
        DataFrame of significant variants. If collect_background_genes=True,
        returns (DataFrame, set_of_all_genes).
    """
    url = f"{BASE_URL}/finngen_R12_{endpoint_id}.gz"
    cache_dir = Path("data/raw/finngen_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached_file = cache_dir / f"finngen_R12_{endpoint_id}.gz"

    # Download to cache (skip if already cached)
    if cached_file.exists():
        print(f"    Using cached file: {cached_file}")
    else:
        result = subprocess.run(
            ["curl", "-s", "-o", str(cached_file), url],
            capture_output=True,
            timeout=300,
        )
        if result.returncode != 0:
            print(f"  Error downloading {endpoint_id}: {result.stderr.decode()}")
            cached_file.unlink(missing_ok=True)
            if collect_background_genes:
                return pd.DataFrame(), set()
            return pd.DataFrame()

    # Read the cached file
    try:
        df = pd.read_csv(
            cached_file,
            sep="\t",
            compression="gzip",
            dtype={
                "#chrom": str,
                "pos": int,
                "ref": str,
                "alt": str,
                "rsids": str,
                "nearest_genes": str,
                "pval": float,
                "mlogp": float,
                "beta": float,
                "sebeta": float,
                "af_alt": float,
                "af_alt_cases": float,
                "af_alt_controls": float,
            },
        )
    except Exception as e:
        print(f"  Error reading {endpoint_id}: {e}")
        if collect_background_genes:
            return pd.DataFrame(), set()
        return pd.DataFrame()

    # Extract background genes before filtering
    background_genes = set()
    if collect_background_genes:
        all_genes = df["nearest_genes"].dropna().str.split(",")
        for gene_list in all_genes:
            background_genes.update(g.strip() for g in gene_list if g.strip())
        print(f"    Background genes from this endpoint: {len(background_genes)}")

    # Filter for significance
    significant = df[df["pval"] <= pval_threshold].copy()

    if collect_background_genes:
        return significant, background_genes
    return significant


def classify_variant_effect(row: pd.Series) -> str:
    """Classify a variant as protective, risk, or ambiguous."""
    beta = row["beta"]
    af_cases = row.get("af_alt_cases", None)
    af_controls = row.get("af_alt_controls", None)

    # Primary signal: beta direction
    if beta < 0:
        return "protective"
    elif beta > 0:
        return "risk"
    return "ambiguous"


def fetch_all_cancer_sumstats(
    pval_threshold: float = PVAL_THRESHOLD,
    output_dir: str = "data/raw",
    collect_background: bool = True,
) -> tuple[pd.DataFrame, set[str]]:
    """
    Download and filter FinnGen summary stats for all cancer endpoints.

    Args:
        pval_threshold: P-value cutoff.
        output_dir: Directory for raw output.
        collect_background: If True, collect all gene symbols across all
            endpoints (before filtering) for enrichment background.

    Returns:
        (combined_significant_variants, background_gene_set)
    """
    all_results = []
    all_background_genes = set()

    for cancer_name, endpoint_id in tqdm(CANCER_ENDPOINTS.items(), desc="FinnGen endpoints"):
        print(f"\n  Downloading {cancer_name} ({endpoint_id})...")

        if collect_background:
            df, bg_genes = download_and_filter(
                endpoint_id, pval_threshold, collect_background_genes=True
            )
            all_background_genes.update(bg_genes)
        else:
            df = download_and_filter(endpoint_id, pval_threshold)

        if df.empty:
            print(f"  No significant variants for {cancer_name}")
            continue

        df["cancer_type"] = cancer_name
        df["endpoint_id"] = endpoint_id
        df["effect"] = df.apply(classify_variant_effect, axis=1)

        n_protective = (df["effect"] == "protective").sum()
        n_risk = (df["effect"] == "risk").sum()
        print(f"  {len(df)} significant variants: "
              f"{n_protective} protective, {n_risk} risk")

        all_results.append(df)

    if not all_results:
        print("No significant variants found across any cancer type!")
        return pd.DataFrame(), all_background_genes

    combined = pd.concat(all_results, ignore_index=True)

    # Save raw results
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path / "finngen_significant_variants.csv", index=False)

    if collect_background:
        # Save background gene list
        with open(out_path / "background_genes.txt", "w") as f:
            f.write("\n".join(sorted(all_background_genes)))
        print(f"\nBackground gene set: {len(all_background_genes)} genes")

    print(f"Total: {len(combined)} significant variants across "
          f"{combined['cancer_type'].nunique()} cancer types")

    return combined, all_background_genes


def extract_protective_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate protective variants by gene.

    For each gene, compute:
    - How many cancer types it's protective against
    - Strength of protective effect (mean beta)
    - Allele frequency difference between cases and controls
    - Consistency of the protective signal
    """
    protective = df[df["effect"] == "protective"].copy()

    if protective.empty:
        return pd.DataFrame()

    # Compute case-control frequency difference
    protective["af_diff"] = (
        protective["af_alt_controls"] - protective["af_alt_cases"]
    )

    # Gene-level aggregation
    gene_stats = (
        protective.groupby("nearest_genes")
        .agg(
            n_protective_variants=("rsids", "nunique"),
            n_cancer_types=("cancer_type", "nunique"),
            cancer_types=("cancer_type", lambda x: ",".join(sorted(set(x)))),
            mean_beta=("beta", "mean"),
            min_beta=("beta", "min"),
            best_pval=("pval", "min"),
            mean_af_diff=("af_diff", "mean"),
            mean_af_cases=("af_alt_cases", "mean"),
            mean_af_controls=("af_alt_controls", "mean"),
            variants=("rsids", lambda x: ",".join(sorted(set(str(v) for v in x)))),
        )
        .reset_index()
        .rename(columns={"nearest_genes": "gene"})
    )

    # Composite protection score:
    # - Stronger negative beta = more protective
    # - More cancer types = more broadly protective
    # - Larger AF difference = more enriched in cancer-free population
    gene_stats["protection_score"] = (
        gene_stats["n_cancer_types"]
        * abs(gene_stats["mean_beta"])
        * (1 + gene_stats["mean_af_diff"].clip(lower=0) * 100)
    )

    gene_stats = gene_stats.sort_values("protection_score", ascending=False)

    return gene_stats


def extract_risk_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Same aggregation but for risk variants — these serve as negative
    examples for the ML model.
    """
    risk = df[df["effect"] == "risk"].copy()
    if risk.empty:
        return pd.DataFrame()

    risk["af_diff"] = risk["af_alt_cases"] - risk["af_alt_controls"]

    gene_stats = (
        risk.groupby("nearest_genes")
        .agg(
            n_risk_variants=("rsids", "nunique"),
            n_cancer_types=("cancer_type", "nunique"),
            cancer_types=("cancer_type", lambda x: ",".join(sorted(set(x)))),
            mean_beta=("beta", "mean"),
            max_beta=("beta", "max"),
            best_pval=("pval", "min"),
            mean_af_diff=("af_diff", "mean"),
        )
        .reset_index()
        .rename(columns={"nearest_genes": "gene"})
    )

    gene_stats["risk_score"] = (
        gene_stats["n_cancer_types"]
        * gene_stats["mean_beta"]
        * (1 + gene_stats["mean_af_diff"].clip(lower=0) * 100)
    )

    return gene_stats.sort_values("risk_score", ascending=False)


def clump_variants(df: pd.DataFrame, window_kb: int = 500) -> pd.DataFrame:
    """Apply LD clumping to significant variants before gene extraction."""
    from src.analysis.ld_clumping import clump_by_distance

    clumped_parts = []
    for cancer_type, group in df.groupby("cancer_type"):
        clumped = clump_by_distance(group, window_kb=window_kb)
        clumped_parts.append(clumped)

    result = pd.concat(clumped_parts, ignore_index=True)
    print(f"  LD clumping: {len(df)} → {len(result)} variants "
          f"(window={window_kb}kb)")
    return result


def expand_gene_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expand rows where nearest_genes contains multiple comma-separated genes.
    Each gene gets its own row with the same variant data.
    """
    if "nearest_genes" not in df.columns:
        return df

    rows = []
    for _, row in df.iterrows():
        genes = str(row["nearest_genes"]).split(",")
        for gene in genes:
            gene = gene.strip()
            if gene and gene != "nan":
                new_row = row.copy()
                new_row["nearest_genes"] = gene
                rows.append(new_row)

    return pd.DataFrame(rows).reset_index(drop=True)


if __name__ == "__main__":
    print("=" * 60)
    print("  FinnGen R12 Cancer-Protective Variant Discovery")
    print(f"  Significance threshold: p < {PVAL_THRESHOLD}")
    print("=" * 60)

    # Step 1: Download and filter
    df, background = fetch_all_cancer_sumstats()
    if df.empty:
        exit(1)

    # Step 2: LD clump
    df = clump_variants(df)

    # Step 3: Expand multi-gene entries
    df = expand_gene_names(df)

    # Step 4: Extract protective genes
    protective_genes = extract_protective_genes(df)
    protective_genes.to_csv("data/processed/finngen_protective_genes.csv", index=False)

    print(f"\nTop 30 protective genes (across multiple cancer types):")
    display_cols = ["gene", "n_cancer_types", "mean_beta", "best_pval",
                    "protection_score", "cancer_types"]
    print(protective_genes[display_cols].head(30).to_string(index=False))

    # Step 5: Extract risk genes
    risk_genes = extract_risk_genes(df)
    risk_genes.to_csv("data/processed/finngen_risk_genes.csv", index=False)

    print(f"\nTop 10 risk genes:")
    print(risk_genes.head(10).to_string(index=False))
