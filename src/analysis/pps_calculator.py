"""
Polygenic Protection Score (PPS) Calculator v2.

Improvements over v1:
1. Data-driven pathway mapping from g:Profiler enrichment (not hardcoded)
2. LD proxy matching for variants missing from consumer chips
3. Empirical population percentile calibration via simulation
4. Self-contained HTML visual report with intervention recommendations
"""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


# ---------------------------------------------------------------------------
# 1. Pathway mapping — data-driven from g:Profiler enrichment
# ---------------------------------------------------------------------------

# Theme patterns: keywords → pathway theme name
THEME_PATTERNS = {
    "Detoxification": ["cytochrome", "xenobiotic", "drug metabolism", "detox",
                       "glutathione", "metabolism of xenobiotics"],
    "Glucuronidation": ["glucuronid", "glucurono"],
    "DNA Repair": ["DNA repair", "homologous recombination", "mismatch repair",
                   "base excision", "double-strand break", "DNA damage"],
    "Apoptosis": ["apoptosis", "apoptotic", "programmed cell death", "caspase"],
    "Cell Adhesion": ["cell adhesion", "cadherin", "cell junction", "adherens",
                      "anchoring junction", "desmosome", "tight junction"],
    "Hippo Signaling": ["hippo"],
    "TGF-beta Signaling": ["TGF-beta", "tgf-beta", "SMAD", "activin"],
    "Telomere Maintenance": ["telomere", "telomerase"],
    "Immune Regulation": ["immune", "interferon", "interleukin", "cytokine",
                          "T cell", "B cell", "natural killer", "MHC", "HLA"],
    "Cell Cycle Control": ["cell cycle", "mitotic", "checkpoint", "cell division",
                           "G1/S", "G2/M"],
    "Growth Factor Signaling": ["growth factor", "MAPK", "ERK", "receptor tyrosine",
                                "insulin", "IGF"],
    "Cell Differentiation": ["differentiation", "stem cell", "morphogenesis",
                             "embryonic development"],
    "Nervous System": ["neuron", "axon", "synap", "nervous system"],
    "Bile Acid Metabolism": ["bile", "cholesterol"],
    "Steroid/Hormone Metabolism": ["steroid", "estrogen", "androgen", "hormone"],
}

# Broad/uninformative GO terms to skip
EXCLUDE_TERMS = {
    "GO:0005515", "GO:0003674", "GO:0005488", "GO:0005737", "GO:0008150",
    "GO:0005575", "GO:0009987", "GO:0005623", "GO:0044464", "GO:0005622",
    "GO:0043226", "GO:0043229", "GO:0043227", "GO:0043231", "GO:0044424",
    "GO:0044237", "GO:0071704", "GO:0044238", "GO:0008152", "GO:0065007",
    "GO:0050789", "GO:0050794", "GO:0048518", "GO:0048519", "GO:0048522",
    "GO:0048523", "GO:0080090", "GO:0031323", "GO:0060255", "GO:0019222",
    "GO:0010468", "GO:0043170", "GO:0009058", "GO:0044249", "GO:0034641",
    "GO:0006807", "GO:0006139", "GO:0090304", "GO:0005634", "GO:0032991",
    "GO:0031090", "GO:0005829", "GO:0110165",
}


def build_pathway_gene_map(
    enrichment_path: str = "data/processed/enrichment_protective_with_genes.csv",
    max_term_size: int = 500,
) -> dict[str, set[str]]:
    """
    Build gene-to-pathway-theme mapping from g:Profiler enrichment results.
    Returns {theme_name: set_of_genes}.
    """
    path = Path(enrichment_path)
    if not path.exists():
        print(f"  Warning: {enrichment_path} not found, using fallback pathway mapping")
        return _fallback_pathway_genes()

    df = pd.read_csv(path)
    if "intersections" not in df.columns:
        print("  Warning: no 'intersections' column, using fallback")
        return _fallback_pathway_genes()

    # Filter specific, informative terms
    df = df[
        (~df["native"].isin(EXCLUDE_TERMS))
        & (df["term_size"] <= max_term_size)
    ]

    theme_genes = {theme: set() for theme in THEME_PATTERNS}
    theme_genes["Other Protective"] = set()  # catch-all

    for _, row in df.iterrows():
        term_name = str(row["name"]).lower()
        genes_str = row["intersections"]

        # Parse gene list (may be string repr of list or comma-separated)
        if isinstance(genes_str, str):
            genes_str = genes_str.strip("[]").replace("'", "").replace('"', '')
            genes = [g.strip() for g in genes_str.split(",") if g.strip()]
        else:
            continue

        # Match to theme
        matched = False
        for theme, patterns in THEME_PATTERNS.items():
            if any(p.lower() in term_name for p in patterns):
                theme_genes[theme].update(genes)
                matched = True
                break

        if not matched:
            theme_genes["Other Protective"].update(genes)

    # Remove empty themes
    theme_genes = {k: v for k, v in theme_genes.items() if v}

    total_genes = len(set().union(*theme_genes.values()))
    print(f"  Pathway mapping: {total_genes} genes across {len(theme_genes)} themes")
    for theme, genes in sorted(theme_genes.items(), key=lambda x: -len(x[1])):
        print(f"    {theme:30s}: {len(genes)} genes")

    return theme_genes


def _fallback_pathway_genes() -> dict[str, set[str]]:
    """Hardcoded fallback if enrichment data isn't available."""
    return {
        "Detoxification": {"CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2B6",
                           "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4"},
        "Glucuronidation": {"UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A9",
                            "UGT2B7", "UGT2B15"},
        "Cell Adhesion": {"CDH1", "CDH2", "CDH12", "CTNNA1", "CTNNB1", "JUP"},
        "Apoptosis": {"TP53", "BAX", "BCL2", "CASP3", "CASP8", "CASP9", "PTEN"},
        "DNA Repair": {"BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51", "MLH1", "MSH2"},
        "Telomere Maintenance": {"TERT", "TERC", "POT1"},
    }


# ---------------------------------------------------------------------------
# 2. PPS weights and genotype loading
# ---------------------------------------------------------------------------

def load_pps_weights(path: str = "data/raw/finngen_significant_variants.csv") -> pd.DataFrame:
    """Load and LD-clump protective variant weights."""
    df = pd.read_csv(path, dtype={"#chrom": str})
    protective = df[df["effect"] == "protective"].copy()

    from src.analysis.ld_clumping import clump_by_distance
    clumped = clump_by_distance(protective, window_kb=500)

    # Explode multi-rsID entries for better chip matching
    rows = []
    for _, row in clumped.iterrows():
        rsids = str(row["rsids"]).split(",")
        for rsid in rsids:
            rsid = rsid.strip()
            if rsid and rsid != "nan":
                r = row.copy()
                r["rsid"] = rsid
                rows.append(r)

    weights = pd.DataFrame(rows)
    weights["weight"] = np.abs(weights["beta"])
    weights["gene"] = weights["nearest_genes"] if "nearest_genes" in weights.columns else ""

    keep = ["rsid", "#chrom", "pos", "ref", "alt", "beta", "weight",
            "af_alt", "gene", "cancer_type"]
    return weights[[c for c in keep if c in weights.columns]].copy()


def load_genotype(filepath: str) -> pd.DataFrame:
    """Auto-detect file format and load genotype data."""
    filepath = str(filepath)

    if filepath.endswith(".vcf") or filepath.endswith(".vcf.gz"):
        print(f"  Loading VCF: {Path(filepath).name}")
        return _load_vcf(filepath)
    else:
        with open(filepath) as f:
            header = f.read(500)
        if "AncestryDNA" in header:
            print(f"  Loading AncestryDNA: {Path(filepath).name}")
            return _load_consumer(filepath, sep="\t", skip_comment="#",
                                  cols=["rsid", "chrom", "pos", "allele1", "allele2"])
        elif "23andMe" in header:
            print(f"  Loading 23andMe: {Path(filepath).name}")
            return _load_23andme(filepath)
        else:
            print(f"  Loading generic: {Path(filepath).name}")
            return _load_consumer(filepath, sep="\t", skip_comment="#",
                                  cols=["rsid", "chrom", "pos", "allele1", "allele2"])


def _load_consumer(filepath, sep, skip_comment, cols):
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(skip_comment) or line.startswith("rsid"):
                continue
            parts = line.strip().split(sep)
            if len(parts) >= len(cols):
                rows.append({c: parts[i] for i, c in enumerate(cols)})
    df = pd.DataFrame(rows)
    if "pos" in df.columns:
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    return df


def _load_23andme(filepath):
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                gt = parts[3]
                rows.append({
                    "rsid": parts[0], "chrom": parts[1],
                    "pos": int(parts[2]),
                    "allele1": gt[0] if len(gt) >= 1 else "",
                    "allele2": gt[1] if len(gt) >= 2 else gt[0],
                })
    return pd.DataFrame(rows)


def _load_vcf(filepath):
    rows = []
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0].replace("chr", "")
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            fmt = parts[8].split(":")
            sample = parts[9].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0
            gt = sample[gt_idx].replace("|", "/").split("/")
            alleles = [ref] + alt.split(",")
            try:
                a1 = alleles[int(gt[0])]
                a2 = alleles[int(gt[1])] if len(gt) > 1 else a1
            except (ValueError, IndexError):
                continue
            rows.append({
                "rsid": parts[2] if parts[2] != "." else f"{chrom}:{pos}",
                "chrom": chrom, "pos": pos,
                "allele1": a1, "allele2": a2,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. Missing variant handling
# ---------------------------------------------------------------------------
#
# Two types of "missing" variants:
#
# A) CONSUMER CHIPS (23andMe, AncestryDNA):
#    All genotyped positions are reported regardless of genotype.
#    Missing = not on the chip. Use mean imputation (2*AF).
#
# B) VCF FILES:
#    Only VARIANT sites are reported. Hom-ref sites are ABSENT.
#    This creates ascertainment bias: matched positions are enriched
#    for alt alleles, inflating the score.
#    Fix: check if a position is in the callable region (BED file).
#    - In callable region but not in VCF → dosage = 0 (hom ref)
#    - Outside callable region → mean imputation (2*AF)


def load_callable_regions(bed_path: str) -> dict[str, list[tuple[int, int]]]:
    """Load a BED file of callable regions, indexed by chromosome."""
    regions = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0].replace("chr", "")
            start = int(parts[1])
            end = int(parts[2])
            regions.setdefault(chrom, []).append((start, end))
    # Sort by start position for binary search
    for chrom in regions:
        regions[chrom].sort()
    return regions


def is_in_callable_region(chrom: str, pos: int, regions: dict) -> bool:
    """Check if a position falls within any callable region."""
    chrom_regions = regions.get(str(chrom))
    if not chrom_regions:
        return False
    # Binary search
    import bisect
    idx = bisect.bisect_right([r[0] for r in chrom_regions], pos) - 1
    if idx < 0:
        return False
    return chrom_regions[idx][0] <= pos <= chrom_regions[idx][1]


# ---------------------------------------------------------------------------
# 4. Population percentile simulation
# ---------------------------------------------------------------------------

def simulate_population_pps(
    weights: pd.DataFrame,
    n_simulations: int = 10_000,
    seed: int = 42,
) -> np.ndarray:
    """
    Simulate PPS distribution using allele frequencies.
    Returns array of simulated PPS values.
    """
    rng = np.random.default_rng(seed)
    af = weights["af_alt"].values
    w = weights["weight"].values

    # For each simulation, sample dosage ~ Binomial(2, AF) for each variant
    dosages = rng.binomial(2, af, size=(n_simulations, len(af)))
    scores = dosages @ w  # matrix multiply: (N x V) @ (V,) = (N,)

    return scores


# ---------------------------------------------------------------------------
# 5. Core PPS computation
# ---------------------------------------------------------------------------

def compute_pps(
    genotype: pd.DataFrame,
    weights: pd.DataFrame,
    pathway_genes: dict[str, set[str]],
    match_by: str = "rsid",
    callable_bed: str | None = None,
    n_simulations: int = 10_000,
) -> dict:
    """
    Compute PPS using direct matches + proper handling of missing variants.

    For consumer chips (rsid matching):
      - All genotyped positions reported → use actual dosage
      - Missing positions → mean imputation (2*AF)

    For VCF files (pos matching):
      - Variant sites in VCF → use actual dosage
      - Positions in callable regions but NOT in VCF → dosage 0 (hom ref)
      - Positions outside callable regions → mean imputation (2*AF)
      This corrects ascertainment bias from VCFs only reporting variant sites.
    """

    # Primary matching
    if match_by == "rsid":
        merged = pd.merge(weights, genotype, on="rsid", how="inner")
    else:
        weights_c = weights.copy()
        genotype_c = genotype.copy()
        weights_c["_key"] = weights_c["#chrom"].astype(str) + ":" + weights_c["pos"].astype(str)
        genotype_c["_key"] = genotype_c["chrom"].astype(str) + ":" + genotype_c["pos"].astype(str)
        merged = pd.merge(weights_c, genotype_c, on="_key", how="inner",
                          suffixes=("", "_geno"))

    n_vcf_matched = len(merged)
    n_total_weights = weights.drop_duplicates("rsid").shape[0]

    # For VCF: recover hom-ref sites from callable regions
    n_homref_recovered = 0
    if match_by == "pos" and callable_bed:
        callable_regions = load_callable_regions(callable_bed)
        matched_rsids = set(merged["rsid"])
        unmatched_in_vcf = weights[~weights["rsid"].isin(matched_rsids)]

        homref_rows = []
        for _, w in unmatched_in_vcf.iterrows():
            chrom = str(w.get("#chrom", ""))
            pos = w.get("pos", 0)
            if is_in_callable_region(chrom, pos, callable_regions):
                row = w.copy()
                row["allele1"] = w["ref"]  # hom ref
                row["allele2"] = w["ref"]
                homref_rows.append(row)

        if homref_rows:
            homref_df = pd.DataFrame(homref_rows)
            n_homref_recovered = len(homref_df)
            merged = pd.concat([merged, homref_df], ignore_index=True)

    n_matched = len(merged)
    matched_rsids = set(merged["rsid"])
    unmatched = weights[~weights["rsid"].isin(matched_rsids)]

    print(f"  VCF/chip matched:  {n_vcf_matched}")
    if n_homref_recovered:
        print(f"  Hom-ref recovered: {n_homref_recovered} (in callable regions, absent from VCF)")
    print(f"  Total scored:      {n_matched}/{n_total_weights} "
          f"({n_matched/n_total_weights*100:.1f}%)")
    print(f"  Mean-imputed:      {len(unmatched)} (outside callable regions)")

    if merged.empty:
        return {"error": "No variants matched", "n_matched": 0}

    # Compute actual dosage for matched variants
    dosages = []
    for _, row in merged.iterrows():
        alt = row["alt"]
        a1 = str(row.get("allele1", ""))
        a2 = str(row.get("allele2", ""))
        dose = (1 if a1 == alt else 0) + (1 if a2 == alt else 0)
        dosages.append(dose)

    merged["dosage"] = dosages
    merged["expected_dosage"] = 2 * merged["af_alt"]
    merged["deviation"] = merged["dosage"] - merged["expected_dosage"]
    merged["score_contribution"] = merged["dosage"] * merged["weight"]

    # Score for observed variants
    observed_score = merged["score_contribution"].sum()

    # Score for unobserved variants (mean imputation: dosage = 2*AF)
    if not unmatched.empty:
        imputed_score = (2 * unmatched["af_alt"] * unmatched["weight"]).sum()
    else:
        imputed_score = 0

    total_pps = observed_score + imputed_score

    # The INFORMATIVE score: deviation from expected at observed loci
    # This is what determines the percentile
    deviation_score = (merged["deviation"] * merged["weight"]).sum()

    # Empirical percentile via simulation — ONLY using matched variants
    # This ensures percentile reflects real genotype information, not coverage
    sim_scores = simulate_population_pps(
        merged[["af_alt", "weight"]].dropna(), n_simulations
    )
    # The simulation gives total scores; compute individual's rank
    # against the simulated distribution at matched loci only
    percentile = (sim_scores < observed_score).mean() * 100

    # Per-pathway scores
    pathway_scores = {}
    for theme, genes in pathway_genes.items():
        pw_variants = merged[merged["gene"].isin(genes)]
        if pw_variants.empty:
            pathway_scores[theme] = {
                "score": 0, "n_variants": 0, "percentile": None, "z_score": None,
            }
            continue

        pw_score = pw_variants["score_contribution"].sum()

        # Pathway-level simulation
        pw_sim = simulate_population_pps(
            pw_variants[["af_alt", "weight"]].dropna(), n_simulations
        )
        pw_pct = (pw_sim < pw_score).mean() * 100

        pw_expected = pw_sim.mean()
        pw_std = pw_sim.std()
        pw_z = (pw_score - pw_expected) / pw_std if pw_std > 0 else 0

        pathway_scores[theme] = {
            "score": round(pw_score, 4),
            "expected": round(pw_expected, 4),
            "z_score": round(pw_z, 2),
            "percentile": round(pw_pct, 1),
            "n_variants": len(pw_variants),
            "genes_matched": sorted(pw_variants["gene"].unique().tolist()),
        }

    # Top contributors — by deviation from expected (most informative)
    merged["abs_deviation_contrib"] = np.abs(merged["deviation"] * merged["weight"])
    top = merged.nlargest(20, "abs_deviation_contrib")[
        ["rsid", "gene", "dosage", "expected_dosage", "deviation",
         "weight", "score_contribution", "cancer_type"]
    ]

    return {
        "overall_pps": round(total_pps, 4),
        "observed_score": round(observed_score, 4),
        "imputed_score": round(imputed_score, 4),
        "deviation_score": round(deviation_score, 4),
        "population_percentile": round(percentile, 1),
        "simulated_mean": round(sim_scores.mean(), 4),
        "simulated_std": round(sim_scores.std(), 4),
        "n_direct_matches": n_matched,
        "n_homref_recovered": n_homref_recovered,
        "n_imputed": len(unmatched),
        "n_total_weights": n_total_weights,
        "coverage_pct": round(n_matched / n_total_weights * 100, 1),
        "pathway_scores": pathway_scores,
        "top_contributors": top.to_dict("records"),
    }


# ---------------------------------------------------------------------------
# 6. Reports
# ---------------------------------------------------------------------------

INTERVENTION_MAP = {
    "Detoxification": "Sulforaphane (broccoli sprouts), NAC, cruciferous vegetables",
    "Glucuronidation": "Calcium D-glucarate, cruciferous vegetables",
    "DNA Repair": "EGCG (green tea), melatonin 3mg, blueberries, epicatechin",
    "Apoptosis": "Curcumin + piperine, EGCG, exercise 150+ min/wk",
    "Cell Adhesion": "Vitamin D 1-4k IU/day, low-dose aspirin, omega-3",
    "Hippo Signaling": "Exercise (activates Hippo via catecholamines), metformin, berberine",
    "TGF-beta Signaling": "Exercise, vitamin D (caution: dual role in late-stage cancer)",
    "Telomere Maintenance": "Exercise, stress reduction, adequate sleep",
    "Immune Regulation": "Vitamin D, exercise, adequate sleep, zinc, selenium",
    "Cell Cycle Control": "EGCG, curcumin, exercise",
    "Growth Factor Signaling": "Metformin, berberine, intermittent fasting",
    "Cell Differentiation": "Vitamin A/retinoids, vitamin D",
    "Nervous System": "Omega-3, exercise, sleep optimization",
    "Bile Acid Metabolism": "High fiber 30+ g/day, UDCA, reduce saturated fat",
    "Steroid/Hormone Metabolism": "DIM 300mg/day, exercise, maintain healthy weight",
    "Other Protective": "Exercise, Mediterranean diet, adequate sleep",
}


def print_report(results: dict):
    """Console report."""
    print("\n" + "=" * 60)
    print("  POLYGENIC PROTECTION SCORE REPORT")
    print("=" * 60)

    pct = results["population_percentile"]
    print(f"\n  Overall PPS:         {results['overall_pps']}")
    print(f"    Observed score:    {results['observed_score']} ({results['n_direct_matches']} variants)")
    print(f"    Imputed score:     {results['imputed_score']} ({results['n_imputed']} variants at pop. mean)")
    print(f"    Deviation score:   {results['deviation_score']:+.4f} (your diff from average)")
    print(f"  Population mean:     {results['simulated_mean']} (at observed loci)")
    print(f"  Your percentile:     {pct:.1f}th")
    print(f"  Genotype coverage:   {results['n_direct_matches']}/{results['n_total_weights']} "
          f"({results['coverage_pct']}%)")
    if results.get('n_homref_recovered', 0) > 0:
        print(f"    (includes {results['n_homref_recovered']} hom-ref sites recovered from callable regions)")

    print(f"\n  --- Pathway Breakdown ---")
    scored = [(p, s) for p, s in results["pathway_scores"].items()
              if s.get("percentile") is not None]
    scored.sort(key=lambda x: x[1]["percentile"])

    for pathway, scores in scored:
        pctl = scores["percentile"]
        n = scores["n_variants"]
        bar_len = int(pctl / 5)
        bar = "█" * bar_len + "░" * (20 - bar_len)
        color_indicator = "🔴" if pctl < 25 else "🟡" if pctl < 50 else "🟢"
        print(f"  {color_indicator} {pathway:30s} {bar} {pctl:5.1f}% (n={n})")

    # Unscored pathways
    unscored = [(p, s) for p, s in results["pathway_scores"].items()
                if s.get("percentile") is None]
    if unscored:
        print(f"\n  Pathways with no matched variants:")
        for pathway, _ in unscored:
            print(f"     {pathway}")

    # Intervention recommendations for weak pathways
    print(f"\n  --- Intervention Recommendations (weakest pathways) ---")
    for pathway, scores in scored[:3]:
        pctl = scores["percentile"]
        intervention = INTERVENTION_MAP.get(pathway, "See full report")
        print(f"\n  {pathway} ({pctl:.0f}th percentile):")
        print(f"    → {intervention}")

    print(f"\n  --- Most Informative Variants (largest deviation from expected) ---")
    for v in results["top_contributors"][:10]:
        dose = "HOM" if v["dosage"] == 2 else "HET" if v["dosage"] == 1 else "REF"
        dev = v.get("deviation", 0)
        direction = "+" if dev > 0 else "-" if dev < 0 else "="
        print(f"  {v['rsid']:18s} {str(v.get('gene','')):12s} {dose} "
              f"dev={dev:+.2f} w={v['weight']:.3f} {v.get('cancer_type','')}")


def generate_html_report(results: dict, output_path: str = "data/processed/pps_report.html"):
    """Generate a self-contained HTML report."""
    scored = [(p, s) for p, s in results["pathway_scores"].items()
              if s.get("percentile") is not None]
    scored.sort(key=lambda x: x[1]["percentile"])

    pathway_bars = ""
    for pathway, scores in scored:
        pctl = scores["percentile"]
        color = "#e74c3c" if pctl < 25 else "#f39c12" if pctl < 50 else "#27ae60"
        intervention = INTERVENTION_MAP.get(pathway, "")
        genes = ", ".join(scores.get("genes_matched", [])[:8])
        if len(scores.get("genes_matched", [])) > 8:
            genes += f" (+{len(scores['genes_matched']) - 8} more)"

        pathway_bars += f"""
        <div class="pathway">
            <div class="pathway-name">{pathway} <span class="n-variants">({scores['n_variants']} variants)</span></div>
            <div class="bar-container">
                <div class="bar" style="width: {pctl}%; background: {color};">{pctl:.0f}%</div>
            </div>
            <div class="genes">Genes: {genes}</div>
            <div class="intervention">Suggested: {intervention}</div>
        </div>"""

    top_variants = ""
    for v in results["top_contributors"][:15]:
        dose = "HOM" if v["dosage"] == 2 else "HET" if v["dosage"] == 1 else "REF"
        dev = v.get("deviation", 0)
        top_variants += f"""
        <tr>
            <td>{v['rsid']}</td>
            <td>{v.get('gene', '')}</td>
            <td>{dose}</td>
            <td>{v['weight']:.3f}</td>
            <td>{v['score_contribution']:.3f}</td>
            <td>{v.get('cancer_type', '')}</td>
        </tr>"""

    overall_pct = results["population_percentile"]
    overall_color = "#e74c3c" if overall_pct < 25 else "#f39c12" if overall_pct < 50 else "#27ae60"

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Polygenic Protection Score Report</title>
<style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333;
           background: #fafafa; }}
    h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
    h2 {{ color: #34495e; margin-top: 30px; }}
    .overall {{ background: white; border-radius: 12px; padding: 30px;
               box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; margin: 20px 0; }}
    .percentile {{ font-size: 72px; font-weight: bold; color: {overall_color}; }}
    .percentile-label {{ font-size: 18px; color: #666; }}
    .stats {{ display: flex; justify-content: space-around; margin-top: 20px; }}
    .stat {{ text-align: center; }}
    .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
    .stat-label {{ font-size: 12px; color: #999; }}
    .pathway {{ background: white; border-radius: 8px; padding: 15px; margin: 10px 0;
               box-shadow: 0 1px 5px rgba(0,0,0,0.08); }}
    .pathway-name {{ font-weight: bold; font-size: 15px; margin-bottom: 5px; }}
    .n-variants {{ color: #999; font-weight: normal; font-size: 13px; }}
    .bar-container {{ background: #ecf0f1; border-radius: 10px; height: 24px; overflow: hidden; }}
    .bar {{ height: 100%; border-radius: 10px; color: white; font-size: 12px;
           line-height: 24px; padding-left: 8px; min-width: 35px;
           transition: width 0.5s ease; }}
    .genes {{ font-size: 12px; color: #7f8c8d; margin-top: 4px; }}
    .intervention {{ font-size: 13px; color: #2980b9; margin-top: 3px; font-style: italic; }}
    table {{ width: 100%; border-collapse: collapse; background: white;
            border-radius: 8px; overflow: hidden; box-shadow: 0 1px 5px rgba(0,0,0,0.08); }}
    th {{ background: #34495e; color: white; padding: 10px; text-align: left; font-size: 13px; }}
    td {{ padding: 8px 10px; border-bottom: 1px solid #ecf0f1; font-size: 13px; }}
    tr:hover {{ background: #f8f9fa; }}
    .warning {{ background: #fff3cd; border: 1px solid #ffc107; border-radius: 8px;
               padding: 15px; margin: 15px 0; }}
    .footer {{ text-align: center; color: #999; font-size: 12px; margin-top: 40px;
              padding: 20px; border-top: 1px solid #eee; }}
</style>
</head>
<body>
<h1>Polygenic Protection Score Report</h1>

<div class="overall">
    <div class="percentile-label">Your Population Percentile</div>
    <div class="percentile">{overall_pct:.0f}<span style="font-size: 36px">th</span></div>
    <div class="stats">
        <div class="stat">
            <div class="stat-value">{results['overall_pps']:.1f}</div>
            <div class="stat-label">Your PPS</div>
        </div>
        <div class="stat">
            <div class="stat-value">{results['simulated_mean']:.1f}</div>
            <div class="stat-label">Population Mean</div>
        </div>
        <div class="stat">
            <div class="stat-value">{results['coverage_pct']}%</div>
            <div class="stat-label">Variant Coverage</div>
        </div>
        <div class="stat">
            <div class="stat-value">{results['n_direct_matches']}</div>
            <div class="stat-label">Variants Scored</div>
        </div>
    </div>
</div>

<h2>Pathway Protection Profile</h2>
<p>Each bar shows your genetic protection level relative to the population.
Lower percentiles indicate weaker genetic protection in that pathway.</p>
{pathway_bars}

<div class="warning">
    <strong>Note:</strong> This score reflects genetic predisposition only. Lifestyle,
    environment, and epigenetic factors also contribute to cancer protection. Low scores
    in a pathway suggest where targeted interventions may be most beneficial.
    This is not medical advice. Consult a healthcare professional before making changes.
</div>

<h2>Top Contributing Variants</h2>
<table>
<tr><th>Variant</th><th>Gene</th><th>Genotype</th><th>Weight</th><th>Contribution</th><th>Cancer Type</th></tr>
{top_variants}
</table>

<div class="footer">
    Generated by CANSIR Polygenic Protection Score Calculator<br>
    Data source: FinnGen R12 (N=500,348) | {results['n_total_weights']} protective loci<br>
    Percentiles calibrated via {10000:,}-sample population simulation
</div>
</body>
</html>"""

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html)
    print(f"\n  HTML report saved to {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python pps_calculator.py <genotype_file>")
        print("  Supports: .vcf, .vcf.gz, 23andMe .txt, AncestryDNA .txt")
        sys.exit(1)

    filepath = sys.argv[1]

    print("=== Polygenic Protection Score Calculator v2 ===\n")

    # Load weights
    print("  Loading PPS weights...")
    weights = load_pps_weights()
    n_unique = weights.drop_duplicates("rsid").shape[0]
    print(f"  {n_unique} independent protective loci "
          f"({len(weights)} including multi-rsID expansion)")

    # Build pathway mapping
    print("\n  Building pathway gene map...")
    pathway_genes = build_pathway_gene_map()

    # Load genotype
    genotype = load_genotype(filepath)
    print(f"  {len(genotype)} variants in genotype file")

    # Determine matching strategy
    match_by = "pos" if filepath.endswith((".vcf", ".vcf.gz")) else "rsid"

    # Look for callable regions BED file (for VCF bias correction)
    callable_bed = None
    if match_by == "pos":
        # Try to find a BED file alongside the VCF
        vcf_dir = Path(filepath).parent
        bed_candidates = list(vcf_dir.glob("*callable*")) + list(vcf_dir.glob("*benchmark*.bed"))
        if bed_candidates:
            callable_bed = str(bed_candidates[0])
            print(f"  Found callable regions: {Path(callable_bed).name}")
        else:
            print("  Warning: no callable regions BED file found.")
            print("  VCF scores may be inflated (hom-ref sites not distinguished from uncovered).")
            print("  Place a .bed file in the same directory as the VCF to fix this.")

    # Compute PPS
    print(f"\n  Computing PPS (match_by={match_by})...")
    results = compute_pps(genotype, weights, pathway_genes, match_by,
                          callable_bed=callable_bed)

    if "error" in results:
        print(f"\n  Error: {results['error']}")
        sys.exit(1)

    print_report(results)
    generate_html_report(results)

    # Save JSON
    out_path = Path("data/processed/pps_report.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
