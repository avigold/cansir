"""
Fetch cancer-related associations from the GWAS Catalog REST API.

Identifies variants with protective effects (OR < 1) across cancer phenotypes,
mapped to their associated genes.
"""

import json
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"

# EFO terms for major cancer types
CANCER_EFOS = {
    "cancer": "EFO_0000311",
    "breast_cancer": "EFO_0000305",
    "lung_cancer": "EFO_0001071",
    "colorectal_cancer": "EFO_0005842",
    "prostate_cancer": "EFO_0001663",
    "melanoma": "EFO_0000756",
    "pancreatic_cancer": "EFO_0002618",
    "ovarian_cancer": "EFO_0001075",
    "liver_cancer": "EFO_0002890",
    "kidney_cancer": "EFO_0003086",
    "bladder_cancer": "EFO_0003032",
    "stomach_cancer": "EFO_0000503",
    "leukemia": "EFO_0000574",
    "lymphoma": "EFO_0000574",
    "thyroid_cancer": "EFO_0002892",
    "brain_cancer": "EFO_0000519",
}


def fetch_associations_for_trait(efo_id: str, page_size: int = 500) -> list[dict]:
    """Fetch all associations for a given EFO trait ID."""
    associations = []
    page = 0

    while True:
        url = (
            f"{BASE_URL}/efoTraits/{efo_id}/associations"
            f"?page={page}&size={page_size}"
        )
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  Error fetching page {page} for {efo_id}: {e}")
            break

        data = resp.json()
        embedded = data.get("_embedded", {})
        assocs = embedded.get("associations", [])

        if not assocs:
            break

        associations.extend(assocs)

        # Check if there are more pages
        links = data.get("_links", {})
        if "next" not in links:
            break

        page += 1
        time.sleep(0.2)  # Be polite to the API

    return associations


def parse_association(assoc: dict) -> list[dict]:
    """Extract relevant fields from a GWAS Catalog association record."""
    rows = []

    risk_allele_freq = assoc.get("riskFrequency")
    p_value = assoc.get("pvalue")
    p_value_mlog = assoc.get("pvalueMantissa")
    p_value_exponent = assoc.get("pvalueExponent")

    # Extract OR or beta
    or_value = assoc.get("orPerCopyNum")
    beta = assoc.get("betaNum")
    beta_direction = assoc.get("betaDirection")
    ci = assoc.get("range")

    # Extract SNPs and genes from loci
    loci = assoc.get("loci", [])
    for locus in loci:
        for risk_allele in locus.get("strongestRiskAlleles", []):
            snp_name = risk_allele.get("riskAlleleName", "")
            # Format is typically "rs12345-A"
            rsid = snp_name.split("-")[0] if "-" in snp_name else snp_name

        for gene_entry in locus.get("authorReportedGenes", []):
            gene_name = gene_entry.get("geneName", "")

            rows.append({
                "rsid": rsid,
                "gene": gene_name,
                "or_value": or_value,
                "beta": beta,
                "beta_direction": beta_direction,
                "p_value": p_value,
                "p_value_mlog": p_value_mlog,
                "p_value_exponent": p_value_exponent,
                "risk_allele_freq": risk_allele_freq,
                "ci": ci,
            })

    return rows


def fetch_all_cancer_associations(output_dir: str = "data/raw") -> pd.DataFrame:
    """Fetch associations for all cancer types and combine into a DataFrame."""
    all_rows = []

    for cancer_name, efo_id in tqdm(CANCER_EFOS.items(), desc="Cancer types"):
        print(f"\nFetching {cancer_name} ({efo_id})...")
        associations = fetch_associations_for_trait(efo_id)
        print(f"  Found {len(associations)} associations")

        for assoc in associations:
            rows = parse_association(assoc)
            for row in rows:
                row["cancer_type"] = cancer_name
            all_rows.extend(rows)

        time.sleep(1)  # Pause between cancer types

    df = pd.DataFrame(all_rows)

    if df.empty:
        print("No associations found!")
        return df

    # Save raw data
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path / "gwas_cancer_associations.csv", index=False)
    print(f"\nSaved {len(df)} total association records")

    return df


def filter_protective_variants(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for variants with protective effects.

    Protective = OR < 1 (the 'risk' allele actually reduces risk)
    or negative beta with 'decrease' direction.
    """
    protective = df[
        (df["or_value"].notna() & (df["or_value"] < 1.0))
        | (
            df["beta"].notna()
            & df["beta_direction"].notna()
            & (df["beta_direction"].str.lower() == "decrease")
        )
    ].copy()

    # Add a protection_score: how far below 1.0 the OR is
    protective["protection_score"] = protective["or_value"].apply(
        lambda x: 1.0 - x if pd.notna(x) else None
    )

    return protective.sort_values("protection_score", ascending=False)


def summarize_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate protective variants by gene, counting evidence across cancer types."""
    if df.empty:
        return df

    gene_summary = (
        df.groupby("gene")
        .agg(
            n_variants=("rsid", "nunique"),
            n_cancer_types=("cancer_type", "nunique"),
            cancer_types=("cancer_type", lambda x: ", ".join(sorted(set(x)))),
            min_or=("or_value", "min"),
            mean_or=("or_value", "mean"),
            best_pvalue=("p_value", "min"),
            variants=("rsid", lambda x: ", ".join(sorted(set(x)))),
        )
        .reset_index()
    )

    # Score: genes protective across MORE cancer types rank higher
    gene_summary["multi_cancer_score"] = (
        gene_summary["n_cancer_types"] * (1.0 - gene_summary["mean_or"].fillna(0.5))
    )

    return gene_summary.sort_values("multi_cancer_score", ascending=False)


if __name__ == "__main__":
    print("=== GWAS Catalog Cancer-Protective Variant Finder ===\n")

    # Step 1: Fetch all cancer associations
    df = fetch_all_cancer_associations()

    if df.empty:
        print("No data retrieved. Check network connection.")
        exit(1)

    print(f"\nTotal associations: {len(df)}")
    print(f"Unique variants: {df['rsid'].nunique()}")
    print(f"Unique genes: {df['gene'].nunique()}")

    # Step 2: Filter for protective variants
    protective = filter_protective_variants(df)
    print(f"\nProtective variants (OR < 1): {len(protective)}")

    # Save protective variants
    protective.to_csv("data/processed/protective_variants.csv", index=False)

    # Step 3: Summarize by gene
    gene_summary = summarize_by_gene(protective)
    gene_summary.to_csv("data/processed/protective_genes_summary.csv", index=False)

    print(f"\nTop 20 genes with protective variants across multiple cancer types:")
    print(gene_summary.head(20).to_string(index=False))
