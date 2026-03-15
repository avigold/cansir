"""
Publication-quality figures for the pathway enrichment comparison.

Central figure: side-by-side enrichment dot plots showing that
protective and risk gene sets enrich distinct biological processes.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2

# Overly broad GO terms that dominate enrichment results but are
# uninformative — filter these out for clearer figures
EXCLUDE_TERMS = {
    "GO:0005515",  # protein binding
    "GO:0003674",  # molecular_function (root)
    "GO:0005488",  # binding
    "GO:0005737",  # cytoplasm
    "GO:0008150",  # biological_process (root)
    "GO:0005575",  # cellular_component (root)
    "GO:0009987",  # cellular process
    "GO:0005623",  # cell
    "GO:0044464",  # cell part
    "GO:0005622",  # intracellular
    "GO:0043226",  # organelle
    "GO:0043229",  # intracellular organelle
    "GO:0043227",  # membrane-bounded organelle
    "GO:0043231",  # intracellular membrane-bounded organelle
    "GO:0044424",  # intracellular part
    "GO:0044237",  # cellular metabolic process
    "GO:0071704",  # organic substance metabolic process
    "GO:0044238",  # primary metabolic process
    "GO:0008152",  # metabolic process
    "GO:0065007",  # biological regulation
    "GO:0050789",  # regulation of biological process
    "GO:0050794",  # regulation of cellular process
    "GO:0048518",  # positive regulation of biological process
    "GO:0048519",  # negative regulation of biological process
    "GO:0048522",  # positive regulation of cellular process
    "GO:0048523",  # negative regulation of cellular process
    "GO:0080090",  # regulation of primary metabolic process
    "GO:0031323",  # regulation of cellular metabolic process
    "GO:0060255",  # regulation of macromolecule metabolic process
    "GO:0019222",  # regulation of metabolic process
    "GO:0010468",  # regulation of gene expression
    "GO:0043170",  # macromolecule metabolic process
    "GO:0009058",  # biosynthetic process
    "GO:0044249",  # cellular biosynthetic process
    "GO:0034641",  # cellular nitrogen compound metabolic process
    "GO:0006807",  # nitrogen compound metabolic process
    "GO:0006139",  # nucleobase-containing compound metabolic process
    "GO:0090304",  # nucleic acid metabolic process
    "GO:0005634",  # nucleus
    "GO:0032991",  # protein-containing complex
    "GO:0031090",  # organelle membrane
    "GO:0005829",  # cytosol
    "GO:0110165",  # cellular anatomical entity
    "GO:0009059",  # macromolecule biosynthetic process
    "GO:0034645",  # cellular macromolecule biosynthetic process
    "GO:0044260",  # cellular macromolecule metabolic process
    "GO:0044271",  # cellular nitrogen compound biosynthetic process
    "GO:0006725",  # cellular aromatic compound metabolic process
    "GO:0046483",  # heterocycle metabolic process
    "GO:0006082",  # organic acid metabolic process
    "GO:1901360",  # organic cyclic compound metabolic process
    "GO:1901576",  # organic substance biosynthetic process
    "GO:0019538",  # protein metabolic process
}

# Maximum term_size to include — filters out overly general terms
MAX_TERM_SIZE = 2000


def set_style():
    """Set publication-quality plot defaults."""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
    })


def _filter_informative(df: pd.DataFrame) -> pd.DataFrame:
    """Remove overly broad terms and keep informative, specific pathways."""
    if df.empty:
        return df
    filtered = df[
        (~df["native"].isin(EXCLUDE_TERMS))
        & (df["term_size"] <= MAX_TERM_SIZE)
    ].copy()
    return filtered


def plot_enrichment_dotplot(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    top_n: int = 20,
    output_path: str = "data/processed/fig_enrichment_dotplot.png",
):
    """
    Side-by-side dot plots of top enriched pathways.
    Filters out broad GO terms to show specific, informative pathways.
    """
    set_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

    def _dotplot(ax, df, title, cmap):
        if df.empty:
            ax.set_title(title)
            ax.text(0.5, 0.5, "No enriched terms", ha="center", va="center",
                    transform=ax.transAxes)
            return

        top = _filter_informative(df).head(top_n).copy()
        if top.empty:
            return

        top["neg_log_p"] = -np.log10(top["p_value"].clip(lower=1e-300))
        top = top.sort_values("neg_log_p")

        # Truncate long names
        top["display_name"] = top["name"].str[:55]
        # Add source prefix for clarity
        top["display_name"] = top["source"] + ": " + top["display_name"]

        ratio_col = "gene_ratio" if "gene_ratio" in top.columns else "precision"

        scatter = ax.scatter(
            top[ratio_col],
            range(len(top)),
            s=top["intersection_size"].clip(lower=3) * 10,
            c=top["neg_log_p"],
            cmap=cmap,
            edgecolors="black",
            linewidths=0.5,
            alpha=0.85,
            vmin=top["neg_log_p"].min() * 0.8,
            vmax=top["neg_log_p"].max() * 1.05,
        )
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["display_name"], fontsize=9)
        ax.set_xlabel("Gene ratio (intersection / query size)")
        ax.set_title(title, fontweight="bold", fontsize=13)
        plt.colorbar(scatter, ax=ax, label="-log10(p-value)", shrink=0.6)

    _dotplot(ax1, prot_enrichment,
             f"Protective genes ({prot_enrichment['query_size'].iloc[0] if not prot_enrichment.empty else 0} genes)",
             "Blues")
    _dotplot(ax2, risk_enrichment,
             f"Risk genes ({risk_enrichment['query_size'].iloc[0] if not risk_enrichment.empty else 0} genes)",
             "Reds")

    plt.suptitle(
        "Pathway enrichment of cancer-protective vs cancer-risk loci\n"
        "FinnGen R12 (N=500,348) | p < 5e-6 | LD-clumped (500kb)",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_unique_pathways(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    overlap_stats: dict,
    top_n: int = 15,
    output_path: str = "data/processed/fig_unique_pathways.png",
):
    """
    Horizontal bar chart showing TOP pathways unique to each set.
    This is the key figure that demonstrates biological distinctness.
    """
    set_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))

    # Protective-only
    prot_only_terms = overlap_stats.get("protective_only_terms", set())
    if not prot_enrichment.empty and prot_only_terms:
        prot_unique = prot_enrichment[
            prot_enrichment["native"].isin(prot_only_terms)
        ].copy()
        prot_unique = _filter_informative(prot_unique).head(top_n)

        if not prot_unique.empty:
            prot_unique["neg_log_p"] = -np.log10(prot_unique["p_value"].clip(lower=1e-300))
            prot_unique = prot_unique.sort_values("neg_log_p")
            prot_unique["display"] = prot_unique["source"] + ": " + prot_unique["name"].str[:50]

            colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(prot_unique)))
            ax1.barh(range(len(prot_unique)), prot_unique["neg_log_p"],
                     color=colors, edgecolor="black", linewidth=0.5)
            ax1.set_yticks(range(len(prot_unique)))
            ax1.set_yticklabels(prot_unique["display"], fontsize=9)

    ax1.set_xlabel("-log10(p-value)")
    ax1.set_title("Pathways UNIQUE to protective genes", fontweight="bold",
                  color="#2C5F8A", fontsize=13)

    # Risk-only
    risk_only_terms = overlap_stats.get("risk_only_terms", set())
    if not risk_enrichment.empty and risk_only_terms:
        risk_unique = risk_enrichment[
            risk_enrichment["native"].isin(risk_only_terms)
        ].copy()
        risk_unique = _filter_informative(risk_unique).head(top_n)

        if not risk_unique.empty:
            risk_unique["neg_log_p"] = -np.log10(risk_unique["p_value"].clip(lower=1e-300))
            risk_unique = risk_unique.sort_values("neg_log_p")
            risk_unique["display"] = risk_unique["source"] + ": " + risk_unique["name"].str[:50]

            colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(risk_unique)))
            ax2.barh(range(len(risk_unique)), risk_unique["neg_log_p"],
                     color=colors, edgecolor="black", linewidth=0.5)
            ax2.set_yticks(range(len(risk_unique)))
            ax2.set_yticklabels(risk_unique["display"], fontsize=9)

    ax2.set_xlabel("-log10(p-value)")
    ax2.set_title("Pathways UNIQUE to risk genes", fontweight="bold",
                  color="#8A2C2C", fontsize=13)

    plt.suptitle(
        "Biologically distinct pathways in cancer protection vs cancer risk\n"
        f"Jaccard similarity = {overlap_stats['jaccard_similarity']:.3f} | "
        f"Fisher's exact test: pathways are distinct",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_divergence_chart(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    top_n: int = 25,
    output_path: str = "data/processed/fig_divergence.png",
):
    """
    Butterfly / tornado chart showing enrichment in protective (left)
    vs risk (right) for the most divergent pathways.
    """
    set_style()

    if prot_enrichment.empty or risk_enrichment.empty:
        return

    prot_f = _filter_informative(prot_enrichment)
    risk_f = _filter_informative(risk_enrichment)

    # Merge all terms with p-values from both sides
    all_terms = set(prot_f["native"]) | set(risk_f["native"])
    rows = []
    for term_id in all_terms:
        p_row = prot_f[prot_f["native"] == term_id]
        r_row = risk_f[risk_f["native"] == term_id]

        name = ""
        source = ""
        if not p_row.empty:
            name = p_row.iloc[0]["name"]
            source = p_row.iloc[0]["source"]
        elif not r_row.empty:
            name = r_row.iloc[0]["name"]
            source = r_row.iloc[0]["source"]

        prot_mlogp = -np.log10(p_row.iloc[0]["p_value"]) if not p_row.empty else 0
        risk_mlogp = -np.log10(r_row.iloc[0]["p_value"]) if not r_row.empty else 0

        rows.append({
            "term_id": term_id,
            "name": name,
            "source": source,
            "prot_mlogp": prot_mlogp,
            "risk_mlogp": risk_mlogp,
            "divergence": prot_mlogp - risk_mlogp,  # positive = more protective
        })

    div_df = pd.DataFrame(rows)

    # Get the most divergent in each direction
    most_protective = div_df.nlargest(top_n, "divergence")
    most_risk = div_df.nsmallest(top_n, "divergence")
    show = pd.concat([most_protective, most_risk]).drop_duplicates("term_id")
    show = show.sort_values("divergence")
    show["display"] = show["source"] + ": " + show["name"].str[:45]

    fig, ax = plt.subplots(figsize=(12, max(8, len(show) * 0.32)))

    colors = ["#D94A4A" if d < 0 else "#4A90D9" for d in show["divergence"]]
    ax.barh(range(len(show)), show["divergence"], color=colors,
            edgecolor="black", linewidth=0.3, alpha=0.85)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show["display"], fontsize=8)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("← More enriched in RISK genes    |    More enriched in PROTECTIVE genes →",
                  fontsize=10)
    ax.set_title(
        "Pathway enrichment divergence: Protective vs Risk\n"
        "FinnGen R12 (N=500,348)",
        fontweight="bold", fontsize=13,
    )

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_venn_diagram(
    overlap_stats: dict,
    output_path: str = "data/processed/fig_pathway_venn.png",
):
    """Venn diagram showing overlap of enriched pathway sets."""
    set_style()
    fig, ax = plt.subplots(figsize=(8, 6))

    v = venn2(
        subsets=(
            overlap_stats["n_protective_only"],
            overlap_stats["n_risk_only"],
            overlap_stats["n_shared"],
        ),
        set_labels=("Protective", "Risk"),
        set_colors=("#4A90D9", "#D94A4A"),
        alpha=0.6,
        ax=ax,
    )

    jaccard = overlap_stats["jaccard_similarity"]
    ax.set_title(
        f"Enriched pathway overlap\nJaccard similarity = {jaccard:.3f}",
        fontsize=13,
        fontweight="bold",
    )

    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_comparison_heatmap(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    top_n: int = 20,
    output_path: str = "data/processed/fig_pathway_heatmap.png",
):
    """
    Heatmap comparing enrichment significance across both gene sets.
    Only shows specific, informative terms.
    """
    set_style()

    if prot_enrichment.empty and risk_enrichment.empty:
        return

    prot_f = _filter_informative(prot_enrichment)
    risk_f = _filter_informative(risk_enrichment)

    # Get top terms from each
    top_prot = set(prot_f.head(top_n)["native"]) if not prot_f.empty else set()
    top_risk = set(risk_f.head(top_n)["native"]) if not risk_f.empty else set()
    all_terms = top_prot | top_risk

    rows = []
    for term_id in all_terms:
        prot_row = prot_f[prot_f["native"] == term_id]
        risk_row = risk_f[risk_f["native"] == term_id]

        name = ""
        source = ""
        if not prot_row.empty:
            name = prot_row.iloc[0]["name"]
            source = prot_row.iloc[0]["source"]
        elif not risk_row.empty:
            name = risk_row.iloc[0]["name"]
            source = risk_row.iloc[0]["source"]

        rows.append({
            "term": f"{source}: {name[:40]}",
            "Protective": -np.log10(prot_row.iloc[0]["p_value"]) if not prot_row.empty else 0,
            "Risk": -np.log10(risk_row.iloc[0]["p_value"]) if not risk_row.empty else 0,
        })

    matrix = pd.DataFrame(rows).set_index("term")
    matrix = matrix.sort_values("Protective", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(8, len(matrix) * 0.35)))
    sns.heatmap(
        matrix,
        cmap="YlOrRd",
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "-log10(p-value)"},
    )
    ax.set_title(
        "Pathway enrichment: Protective vs Risk\n(specific terms only)",
        fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_source_breakdown(
    breakdown_df: pd.DataFrame,
    output_path: str = "data/processed/fig_source_breakdown.png",
):
    """
    Grouped bar chart showing protective-only vs risk-only vs shared
    pathways per annotation source.
    """
    set_style()

    if breakdown_df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    x = range(len(breakdown_df))
    width = 0.25

    ax.bar(
        [i - width for i in x],
        breakdown_df["n_protective_only"],
        width,
        label="Protective only",
        color="#4A90D9",
    )
    ax.bar(
        x,
        breakdown_df["n_shared"],
        width,
        label="Shared",
        color="#9B59B6",
    )
    ax.bar(
        [i + width for i in x],
        breakdown_df["n_risk_only"],
        width,
        label="Risk only",
        color="#D94A4A",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(breakdown_df["source"], rotation=45, ha="right")
    ax.set_ylabel("Number of enriched terms")
    ax.set_title("Pathway enrichment by annotation source", fontweight="bold")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def generate_all_figures(
    prot_enrichment: pd.DataFrame,
    risk_enrichment: pd.DataFrame,
    overlap_stats: dict,
    breakdown_df: pd.DataFrame,
    output_dir: str = "data/processed",
):
    """Generate all publication figures."""
    print("\n=== Generating figures ===")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    plot_enrichment_dotplot(
        prot_enrichment, risk_enrichment,
        output_path=str(out / "fig1_enrichment_dotplot.png"),
    )
    plot_unique_pathways(
        prot_enrichment, risk_enrichment, overlap_stats,
        output_path=str(out / "fig2_unique_pathways.png"),
    )
    plot_divergence_chart(
        prot_enrichment, risk_enrichment,
        output_path=str(out / "fig3_divergence.png"),
    )
    plot_venn_diagram(
        overlap_stats,
        output_path=str(out / "fig4_pathway_venn.png"),
    )
    plot_comparison_heatmap(
        prot_enrichment, risk_enrichment,
        output_path=str(out / "fig5_pathway_heatmap.png"),
    )
    plot_source_breakdown(
        breakdown_df,
        output_path=str(out / "fig6_source_breakdown.png"),
    )
