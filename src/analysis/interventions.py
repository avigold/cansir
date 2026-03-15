"""
Map cancer-protective pathways to actionable interventions.

Uses DGIdb to find drugs/compounds that interact with protective genes,
then categorizes by pathway and interaction type (activator, inhibitor, etc.).

The goal: identify supplements, drugs, or lifestyle interventions that
could enhance the cancer-protective biological pathways we identified.
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests


DGIDB_GRAPHQL = "https://dgidb.org/api/graphql"

# Batch size for DGIdb queries (too many genes at once can timeout)
BATCH_SIZE = 25


def query_dgidb(gene_names: list[str]) -> list[dict]:
    """Query DGIdb for drug-gene interactions."""
    all_interactions = []

    for i in range(0, len(gene_names), BATCH_SIZE):
        batch = gene_names[i:i + BATCH_SIZE]
        gene_list = '["' + '", "'.join(batch) + '"]'

        query = f"""
        {{
          interactions(geneNames: {gene_list}) {{
            nodes {{
              gene {{ name longName }}
              drug {{ name approved conceptId }}
              interactionScore
              interactionTypes {{ type directionality }}
              publications {{ pmid }}
            }}
          }}
        }}
        """

        try:
            resp = requests.post(
                DGIDB_GRAPHQL,
                json={"query": query},
                timeout=30,
            )
            data = resp.json()
            nodes = data.get("data", {}).get("interactions", {}).get("nodes", [])
            all_interactions.extend(nodes)
        except Exception as e:
            print(f"  DGIdb error for batch {i}: {e}")

        time.sleep(0.3)

    return all_interactions


def parse_interactions(interactions: list[dict]) -> pd.DataFrame:
    """Parse DGIdb interactions into a clean DataFrame."""
    rows = []
    for ix in interactions:
        gene = ix.get("gene", {}).get("name", "")
        gene_long = ix.get("gene", {}).get("longName", "")
        drug = ix.get("drug", {}).get("name", "")
        approved = ix.get("drug", {}).get("approved", False)
        score = ix.get("interactionScore")

        types = ix.get("interactionTypes", [])
        interaction_type = ", ".join(t.get("type", "") or "" for t in types) if types else ""
        directionality = ", ".join(t.get("directionality", "") or "" for t in types) if types else ""

        pmids = [p.get("pmid", "") for p in ix.get("publications", [])]

        rows.append({
            "gene": gene,
            "gene_full_name": gene_long,
            "drug": drug,
            "approved": approved,
            "interaction_type": interaction_type,
            "directionality": directionality,
            "score": score,
            "n_publications": len(pmids),
            "pmids": ";".join(str(p) for p in pmids[:5]),
        })

    return pd.DataFrame(rows)


def categorize_interventions(
    interactions_df: pd.DataFrame,
    pathway_genes: dict[str, list[str]],
) -> dict[str, pd.DataFrame]:
    """
    Group drug interactions by protective pathway.

    Args:
        interactions_df: All drug-gene interactions.
        pathway_genes: {pathway_name: [gene1, gene2, ...]}

    Returns:
        {pathway_name: DataFrame of relevant drug interactions}
    """
    pathway_drugs = {}

    for pathway, genes in pathway_genes.items():
        gene_set = set(genes)
        pathway_ix = interactions_df[interactions_df["gene"].isin(gene_set)]

        if not pathway_ix.empty:
            pathway_drugs[pathway] = pathway_ix.sort_values(
                ["approved", "n_publications"], ascending=[False, False]
            )

    return pathway_drugs


def identify_activators(interactions_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for compounds that ACTIVATE or INDUCE protective genes.

    These are the most actionable — they enhance protective pathways.
    """
    activating_terms = ["activator", "agonist", "inducer", "stimulator",
                        "positive modulator", "potentiator"]
    activating_dirs = ["ACTIVATING"]

    mask = (
        interactions_df["interaction_type"].str.lower().apply(
            lambda x: any(t in str(x) for t in activating_terms)
        )
        | interactions_df["directionality"].isin(activating_dirs)
    )

    return interactions_df[mask].copy()


def identify_natural_compounds(interactions_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for natural compounds, supplements, and nutraceuticals.
    """
    natural_terms = [
        "rutin", "quercetin", "curcumin", "resveratrol", "sulforaphane",
        "indole", "chrysin", "kaempferol", "genistein", "epigallocatechin",
        "lycopene", "lutein", "selenium", "zinc", "vitamin", "folate",
        "folic acid", "retinol", "tocopherol", "ascorbic", "melatonin",
        "berberine", "apigenin", "luteolin", "naringenin", "hesperidin",
        "catechin", "ellagic acid", "phytosterol", "beta-carotene",
        "diindolylmethane", "calcium d-glucarate",
    ]

    mask = interactions_df["drug"].str.lower().apply(
        lambda x: any(t in str(x).lower() for t in natural_terms)
    )

    # Also include approved supplements
    supplement_mask = interactions_df["drug"].str.lower().apply(
        lambda x: any(t in str(x).lower() for t in [
            "vitamin", "mineral", "omega", "probiotic", "fiber",
        ])
    )

    return interactions_df[mask | supplement_mask].copy()


def run_intervention_analysis(
    protective_genes: list[str],
    pathway_gene_map: dict[str, list[str]] | None = None,
    output_dir: str = "data/processed",
) -> dict:
    """
    Full intervention analysis pipeline.

    1. Query DGIdb for all protective gene drug interactions
    2. Categorize by pathway
    3. Identify activators (most actionable)
    4. Identify natural compounds / supplements
    5. Rank by evidence strength
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  Intervention Analysis: Protective Pathway Targets")
    print("=" * 60)

    # Step 1: Query DGIdb
    print(f"\n  Querying DGIdb for {len(protective_genes)} protective genes...")
    raw_interactions = query_dgidb(protective_genes)
    print(f"  Found {len(raw_interactions)} drug-gene interactions")

    if not raw_interactions:
        print("  No interactions found.")
        return {}

    # Step 2: Parse
    all_ix = parse_interactions(raw_interactions)
    all_ix.to_csv(out / "all_drug_gene_interactions.csv", index=False)

    n_genes = all_ix["gene"].nunique()
    n_drugs = all_ix["drug"].nunique()
    n_approved = all_ix[all_ix["approved"]]["drug"].nunique()
    print(f"  {n_genes} genes with interactions")
    print(f"  {n_drugs} unique compounds ({n_approved} approved)")

    # Step 3: Identify activators
    activators = identify_activators(all_ix)
    activators.to_csv(out / "protective_activators.csv", index=False)
    print(f"\n  Activators of protective genes: {len(activators)}")
    if not activators.empty:
        print("  Top activators:")
        top_act = activators.drop_duplicates("drug").head(15)
        for _, row in top_act.iterrows():
            app = " [APPROVED]" if row["approved"] else ""
            print(f"    {row['drug']:35s} → {row['gene']:10s} "
                  f"({row['interaction_type']}){app}")

    # Step 4: Natural compounds
    natural = identify_natural_compounds(all_ix)
    natural.to_csv(out / "natural_compounds.csv", index=False)
    print(f"\n  Natural compounds / supplements: {len(natural)}")
    if not natural.empty:
        print("  Supplement interactions with protective genes:")
        for _, row in natural.drop_duplicates("drug").iterrows():
            app = " [APPROVED]" if row["approved"] else ""
            print(f"    {row['drug']:35s} → {row['gene']:10s} "
                  f"({row['interaction_type']}){app}")

    # Step 5: Pathway-level summary
    if pathway_gene_map:
        print(f"\n  === Pathway-Level Drug Mapping ===")
        pathway_drugs = categorize_interventions(all_ix, pathway_gene_map)
        for pathway, drugs in pathway_drugs.items():
            n = drugs["drug"].nunique()
            n_app = drugs[drugs["approved"]]["drug"].nunique()
            print(f"\n  {pathway}: {n} compounds ({n_app} approved)")
            top = drugs.drop_duplicates("drug").head(5)
            for _, row in top.iterrows():
                app = " [APPROVED]" if row["approved"] else ""
                print(f"    {row['drug']:30s} → {row['gene']}{app}")

    # Step 6: Summary report
    summary = {
        "total_interactions": len(all_ix),
        "genes_with_interactions": n_genes,
        "unique_compounds": n_drugs,
        "approved_compounds": n_approved,
        "activators": len(activators),
        "natural_compounds": len(natural),
    }
    with open(out / "intervention_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    return summary


if __name__ == "__main__":
    # Load protective genes
    genes_df = pd.read_csv("data/processed/finngen_protective_genes.csv")
    protective_genes = genes_df["gene"].tolist()

    run_intervention_analysis(protective_genes)
