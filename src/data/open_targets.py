"""
Fetch gene-level evidence from Open Targets Platform via GraphQL API.

Collects disease associations, pathway memberships, and tractability data
for genes of interest.
"""

import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Key cancer-related pathways and biological processes
PATHWAY_QUERIES = {
    "DNA_repair": "DNA repair",
    "apoptosis": "apoptosis",
    "tumor_suppressor": "tumor suppressor",
    "immune_surveillance": "immune surveillance",
    "cell_cycle_control": "cell cycle checkpoint",
    "telomere_maintenance": "telomere maintenance",
    "autophagy": "autophagy",
}


def query_disease_associations(gene_id: str) -> dict:
    """Get disease association scores for a gene from Open Targets."""
    query = """
    query GeneAssociations($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        approvedName
        associatedDiseases(page: {size: 100, index: 0}) {
          rows {
            disease {
              id
              name
            }
            score
            datatypeScores {
              componentId: id
              score
            }
          }
        }
        pathways {
          pathway
          pathwayId
          topLevelTerm
        }
        functionDescriptions
      }
    }
    """
    try:
        resp = requests.post(
            OT_GRAPHQL,
            json={"query": query, "variables": {"ensemblId": gene_id}},
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json().get("data", {}).get("target", {})
    except requests.RequestException as e:
        print(f"  Error querying {gene_id}: {e}")
        return {}


def search_gene_by_symbol(symbol: str) -> str | None:
    """Look up Ensembl ID for a gene symbol via Open Targets search."""
    query = """
    query GeneSearch($queryString: String!) {
      search(queryString: $queryString, entityNames: ["target"], page: {size: 1, index: 0}) {
        hits {
          id
          name
          entity
        }
      }
    }
    """
    try:
        resp = requests.post(
            OT_GRAPHQL,
            json={"query": query, "variables": {"queryString": symbol}},
            timeout=30,
        )
        resp.raise_for_status()
        hits = resp.json().get("data", {}).get("search", {}).get("hits", [])
        if hits:
            return hits[0]["id"]
    except requests.RequestException:
        pass
    return None


def get_gene_features(gene_symbol: str) -> dict:
    """Build a feature dict for a gene from Open Targets data."""
    ensembl_id = search_gene_by_symbol(gene_symbol)
    if not ensembl_id:
        return {"gene": gene_symbol, "ensembl_id": None}

    data = query_disease_associations(ensembl_id)
    if not data:
        return {"gene": gene_symbol, "ensembl_id": ensembl_id}

    features = {
        "gene": gene_symbol,
        "ensembl_id": ensembl_id,
        "approved_name": data.get("approvedName", ""),
    }

    # Extract cancer-related association scores
    assoc_rows = (
        data.get("associatedDiseases", {}).get("rows", [])
    )
    cancer_scores = []
    for row in assoc_rows:
        disease_name = row.get("disease", {}).get("name", "").lower()
        if any(
            term in disease_name
            for term in ["cancer", "carcinoma", "tumor", "tumour", "neoplasm",
                         "leukemia", "lymphoma", "melanoma", "sarcoma"]
        ):
            cancer_scores.append(row.get("score", 0))

    features["n_cancer_associations"] = len(cancer_scores)
    features["max_cancer_assoc_score"] = max(cancer_scores) if cancer_scores else 0
    features["mean_cancer_assoc_score"] = (
        sum(cancer_scores) / len(cancer_scores) if cancer_scores else 0
    )

    # Extract pathway information
    pathways = data.get("pathways", [])
    pathway_names = [p.get("pathway", "").lower() for p in pathways]
    top_level_terms = [p.get("topLevelTerm", "").lower() for p in pathways]
    all_pathway_text = " ".join(pathway_names + top_level_terms)

    features["n_pathways"] = len(pathways)
    features["in_dna_repair"] = int(
        any("repair" in p for p in pathway_names)
    )
    features["in_apoptosis"] = int(
        any("apoptosis" in p or "death" in p for p in pathway_names)
    )
    features["in_cell_cycle"] = int(
        any("cell cycle" in p or "checkpoint" in p for p in pathway_names)
    )
    features["in_immune"] = int(
        any("immune" in p or "interferon" in p for p in pathway_names)
    )
    features["in_signal_transduction"] = int(
        "signal transduction" in all_pathway_text
    )

    # Function descriptions
    func_desc = " ".join(data.get("functionDescriptions", [])).lower()
    features["func_mentions_repair"] = int("repair" in func_desc)
    features["func_mentions_tumor_suppressor"] = int(
        "tumor suppressor" in func_desc or "tumour suppressor" in func_desc
    )
    features["func_mentions_apoptosis"] = int("apoptosis" in func_desc)
    features["func_mentions_immune"] = int("immune" in func_desc)

    return features


def fetch_features_for_genes(
    gene_symbols: list[str], output_dir: str = "data/raw"
) -> pd.DataFrame:
    """Fetch Open Targets features for a list of gene symbols."""
    all_features = []

    for symbol in tqdm(gene_symbols, desc="Open Targets lookup"):
        features = get_gene_features(symbol)
        all_features.append(features)
        time.sleep(0.3)  # Rate limiting

    df = pd.DataFrame(all_features)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path / "open_targets_features.csv", index=False)

    return df


if __name__ == "__main__":
    # Example: fetch features for some well-known cancer-related genes
    test_genes = ["TP53", "BRCA1", "BRCA2", "APC", "MLH1", "ATM", "CHEK2", "PTEN"]
    print("Testing Open Targets data collection with known cancer genes...\n")
    df = fetch_features_for_genes(test_genes)
    print(df.to_string(index=False))
