"""
Two-sample Mendelian Randomization for cancer protection.

Uses protective variants as genetic instruments to test causal claims:
"Does genetically-determined higher activity of pathway X CAUSE lower
cancer incidence?"

This strengthens the enrichment finding from correlation to causation.

Methods implemented:
- Inverse-Variance Weighted (IVW) — primary estimate
- MR-Egger — tests for directional pleiotropy
- Weighted Median — robust to up to 50% invalid instruments

Exposure: gene expression (eQTL data from GTEx or FinnGen itself)
Outcome: cancer incidence (FinnGen summary stats)
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def ivw_estimate(
    beta_exposure: np.ndarray,
    se_exposure: np.ndarray,
    beta_outcome: np.ndarray,
    se_outcome: np.ndarray,
) -> dict:
    """
    Inverse-Variance Weighted MR estimate.

    The Wald ratio for each instrument is beta_outcome / beta_exposure.
    IVW combines these using inverse-variance weighting.
    """
    # Wald ratios
    ratio = beta_outcome / beta_exposure
    ratio_se = se_outcome / np.abs(beta_exposure)

    # IVW: weighted average of ratios
    weights = 1.0 / (ratio_se ** 2)
    ivw_beta = np.sum(weights * ratio) / np.sum(weights)
    ivw_se = np.sqrt(1.0 / np.sum(weights))

    # Test statistic
    z = ivw_beta / ivw_se
    p_value = 2 * stats.norm.sf(np.abs(z))

    return {
        "method": "IVW",
        "beta": ivw_beta,
        "se": ivw_se,
        "z": z,
        "p_value": p_value,
        "n_instruments": len(beta_exposure),
        "causal_direction": "protective" if ivw_beta < 0 else "risk",
    }


def mr_egger(
    beta_exposure: np.ndarray,
    se_exposure: np.ndarray,
    beta_outcome: np.ndarray,
    se_outcome: np.ndarray,
) -> dict:
    """
    MR-Egger regression.

    Unlike IVW, includes an intercept term. A non-zero intercept indicates
    directional pleiotropy (instruments affect outcome through paths
    other than the exposure).
    """
    if len(beta_exposure) < 3:
        return {
            "method": "MR-Egger",
            "beta": None,
            "intercept": None,
            "p_value": None,
            "note": "Need >= 3 instruments for MR-Egger",
        }

    # Orientation: ensure majority of exposure effects are positive
    sign = np.sign(beta_exposure)
    beta_exp_oriented = beta_exposure * sign
    beta_out_oriented = beta_outcome * sign

    # Weighted regression: beta_outcome ~ intercept + slope * beta_exposure
    # Weights = 1 / se_outcome^2
    weights = 1.0 / (se_outcome ** 2)
    W = np.diag(weights)

    X = np.column_stack([np.ones(len(beta_exp_oriented)), beta_exp_oriented])
    y = beta_out_oriented

    # Weighted least squares
    XtWX = X.T @ W @ X
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return {"method": "MR-Egger", "beta": None, "note": "Singular matrix"}

    coeffs = XtWX_inv @ X.T @ W @ y

    # Residuals for SE estimation
    residuals = y - X @ coeffs
    sigma2 = (residuals.T @ W @ residuals) / (len(y) - 2)
    cov = sigma2 * XtWX_inv
    se_coeffs = np.sqrt(np.diag(cov))

    intercept = coeffs[0]
    slope = coeffs[1]
    se_intercept = se_coeffs[0]
    se_slope = se_coeffs[1]

    # P-values
    z_intercept = intercept / se_intercept if se_intercept > 0 else 0
    z_slope = slope / se_slope if se_slope > 0 else 0
    p_intercept = 2 * stats.norm.sf(np.abs(z_intercept))
    p_slope = 2 * stats.norm.sf(np.abs(z_slope))

    return {
        "method": "MR-Egger",
        "beta": slope,
        "se": se_slope,
        "p_value": p_slope,
        "intercept": intercept,
        "intercept_se": se_intercept,
        "intercept_p": p_intercept,
        "pleiotropy_detected": p_intercept < 0.05,
        "n_instruments": len(beta_exposure),
    }


def weighted_median(
    beta_exposure: np.ndarray,
    se_exposure: np.ndarray,
    beta_outcome: np.ndarray,
    se_outcome: np.ndarray,
    n_bootstrap: int = 1000,
) -> dict:
    """
    Weighted median MR estimate.

    Consistent even if up to 50% of instruments are invalid.
    Uses bootstrapping for SE estimation.
    """
    ratio = beta_outcome / beta_exposure
    ratio_se = se_outcome / np.abs(beta_exposure)
    weights = 1.0 / (ratio_se ** 2)
    weights = weights / np.sum(weights)

    # Weighted median
    order = np.argsort(ratio)
    sorted_ratio = ratio[order]
    sorted_weights = weights[order]
    cumsum = np.cumsum(sorted_weights)
    median_idx = np.searchsorted(cumsum, 0.5)
    wm_beta = sorted_ratio[min(median_idx, len(sorted_ratio) - 1)]

    # Bootstrap SE
    rng = np.random.default_rng(42)
    boot_estimates = []
    for _ in range(n_bootstrap):
        boot_ratio = ratio + rng.normal(0, ratio_se)
        order_b = np.argsort(boot_ratio)
        cumsum_b = np.cumsum(weights[order_b])
        idx_b = np.searchsorted(cumsum_b, 0.5)
        boot_estimates.append(boot_ratio[order_b[min(idx_b, len(boot_ratio) - 1)]])

    wm_se = np.std(boot_estimates)
    z = wm_beta / wm_se if wm_se > 0 else 0
    p_value = 2 * stats.norm.sf(np.abs(z))

    return {
        "method": "Weighted Median",
        "beta": wm_beta,
        "se": wm_se,
        "z": z,
        "p_value": p_value,
        "n_instruments": len(beta_exposure),
    }


def run_mr_for_pathway(
    instruments: pd.DataFrame,
    exposure_col: str = "beta_exposure",
    se_exposure_col: str = "se_exposure",
    outcome_col: str = "beta_outcome",
    se_outcome_col: str = "se_outcome",
) -> dict:
    """
    Run all three MR methods for a set of instruments.

    Args:
        instruments: DataFrame with columns for exposure and outcome
            beta/SE values. Each row is one genetic instrument.

    Returns:
        Dict with results from all three methods.
    """
    beta_exp = instruments[exposure_col].values
    se_exp = instruments[se_exposure_col].values
    beta_out = instruments[outcome_col].values
    se_out = instruments[se_outcome_col].values

    # Remove invalid instruments
    valid = (
        np.isfinite(beta_exp) & np.isfinite(se_exp)
        & np.isfinite(beta_out) & np.isfinite(se_out)
        & (se_exp > 0) & (se_out > 0)
    )
    beta_exp = beta_exp[valid]
    se_exp = se_exp[valid]
    beta_out = beta_out[valid]
    se_out = se_out[valid]

    if len(beta_exp) < 2:
        return {"error": f"Only {len(beta_exp)} valid instruments (need >= 2)"}

    results = {
        "ivw": ivw_estimate(beta_exp, se_exp, beta_out, se_out),
        "egger": mr_egger(beta_exp, se_exp, beta_out, se_out),
        "weighted_median": weighted_median(beta_exp, se_exp, beta_out, se_out),
    }

    return results


def prepare_instruments_from_finngen(
    protective_variants: pd.DataFrame,
    cancer_type: str = "any_cancer",
    min_pval: float = 5e-6,
) -> pd.DataFrame:
    """
    Prepare MR instruments from FinnGen data.

    For a simple MR where exposure = carrying the protective allele
    and outcome = cancer risk, the instrument betas are straightforward:
    - beta_exposure: allele frequency effect (proxy for "exposure" to the allele)
    - beta_outcome: effect on cancer (the FinnGen beta)

    For a more rigorous approach, we'd use eQTL data as the exposure.
    This function prepares the simpler version first.
    """
    if cancer_type != "all":
        df = protective_variants[
            protective_variants["cancer_type"] == cancer_type
        ].copy()
    else:
        df = protective_variants.copy()

    if df.empty:
        return pd.DataFrame()

    # For this simplified MR:
    # Exposure = carrying the alt allele (beta_exposure = 1 for all, se arbitrary)
    # Outcome = cancer risk (beta from FinnGen)
    # This tests: "Does the allele causally reduce cancer risk?"
    # A more rigorous version would use eQTL betas as exposure.
    instruments = pd.DataFrame({
        "rsid": df["rsids"].values,
        "gene": df["nearest_genes"].values if "nearest_genes" in df.columns else "",
        "beta_exposure": np.ones(len(df)),  # Simplified: unit exposure
        "se_exposure": np.full(len(df), 0.1),  # Placeholder
        "beta_outcome": df["beta"].values,
        "se_outcome": df["sebeta"].values,
        "pval": df["pval"].values,
        "af": df["af_alt"].values if "af_alt" in df.columns else np.nan,
    })

    # Filter by significance
    instruments = instruments[instruments["pval"] <= min_pval]

    return instruments


def run_pathway_level_mr(
    variants_df: pd.DataFrame,
    gene_pathway_map: dict[str, list[str]],
    output_dir: str = "data/processed",
) -> pd.DataFrame:
    """
    Run MR grouped by pathway.

    For each protective pathway, use the variants in that pathway's genes
    as instruments and test whether they causally reduce cancer risk.

    Args:
        variants_df: Significant protective variants with beta, sebeta, etc.
        gene_pathway_map: {pathway_name: [gene1, gene2, ...]}

    Returns:
        DataFrame of MR results per pathway.
    """
    results = []

    for pathway, genes in gene_pathway_map.items():
        gene_set = set(genes)
        pathway_variants = variants_df[
            variants_df["nearest_genes"].isin(gene_set)
        ] if "nearest_genes" in variants_df.columns else pd.DataFrame()

        if len(pathway_variants) < 3:
            continue

        instruments = pd.DataFrame({
            "beta_exposure": np.ones(len(pathway_variants)),
            "se_exposure": np.full(len(pathway_variants), 0.1),
            "beta_outcome": pathway_variants["beta"].values,
            "se_outcome": pathway_variants["sebeta"].values,
        })

        mr_results = run_mr_for_pathway(instruments)

        if "error" not in mr_results:
            ivw = mr_results["ivw"]
            egger = mr_results["egger"]

            results.append({
                "pathway": pathway,
                "n_instruments": ivw["n_instruments"],
                "ivw_beta": ivw["beta"],
                "ivw_p": ivw["p_value"],
                "egger_beta": egger.get("beta"),
                "egger_p": egger.get("p_value"),
                "egger_intercept_p": egger.get("intercept_p"),
                "pleiotropy": egger.get("pleiotropy_detected", False),
            })

    results_df = pd.DataFrame(results)

    if not results_df.empty:
        results_df = results_df.sort_values("ivw_p")
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(out / "mr_pathway_results.csv", index=False)

    return results_df


if __name__ == "__main__":
    print("=== Mendelian Randomization: Protective Pathways ===")

    # Load FinnGen protective variants
    variants = pd.read_csv("data/raw/finngen_significant_variants.csv",
                           dtype={"#chrom": str})
    protective = variants[variants["effect"] == "protective"]

    instruments = prepare_instruments_from_finngen(protective, cancer_type="any_cancer")
    print(f"Instruments: {len(instruments)}")

    if len(instruments) >= 3:
        results = run_mr_for_pathway(instruments)
        for method, res in results.items():
            print(f"\n{method}:")
            for k, v in res.items():
                print(f"  {k}: {v}")
