"""
Cross-population validation using BioBank Japan (BBJ) summary statistics.

BBJ provides GWAS summary stats from ~200k East Asian individuals.
Primary download source: JENGER (jenger.riken.jp) — RIKEN's own server.
Fallback: NBDC (humandbs.dbcls.jp)

Genome build: GRCh37 (hg19)
"""

import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# JENGER download URLs (RIKEN's server — reliable)
# Each URL downloads a tar.gz or .gz file directly
JENGER_CANCER_ENDPOINTS = {
    "breast": {"url": "http://jenger.riken.jp/102/", "cases": 5552, "ext": "tar.gz"},
    "prostate": {"url": "http://jenger.riken.jp/158/", "cases": 8645, "ext": "gz"},
    "colorectal": {"url": "http://jenger.riken.jp/49/", "cases": 7062, "ext": "tar.gz"},
    "lung": {"url": "http://jenger.riken.jp/75/", "cases": 4050, "ext": "tar.gz"},
    "gastric": {"url": "http://jenger.riken.jp/61/", "cases": 6563, "ext": "tar.gz"},
    "liver": {"url": "http://jenger.riken.jp/100/", "cases": 1866, "ext": "tar.gz"},
    "renal": {"url": "http://jenger.riken.jp/156/", "cases": 523, "ext": "gz"},
}

PVAL_THRESHOLD = 5e-6


def download_from_jenger(
    cancer_name: str,
    endpoint: dict,
    cache_dir: str = "data/raw/bbj_cache",
) -> Path | None:
    """Download BBJ summary stats from JENGER."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    # Check for cached extracted file
    cached_tsvs = list(cache.glob(f"bbj_{cancer_name}*.tsv")) + \
                  list(cache.glob(f"bbj_{cancer_name}*.txt"))
    if cached_tsvs:
        print(f"    Using cached: {cached_tsvs[0].name}")
        return cached_tsvs[0]

    url = endpoint["url"]
    ext = endpoint["ext"]
    download_path = cache / f"bbj_{cancer_name}.{ext}"

    print(f"    Downloading from JENGER: {url}")
    result = subprocess.run(
        ["curl", "-sL", "--max-time", "600", "-o", str(download_path), url],
        capture_output=True,
        timeout=660,
    )
    if result.returncode != 0:
        print(f"    Download failed: {result.stderr.decode()[:200]}")
        download_path.unlink(missing_ok=True)
        return None

    # Check file size
    if download_path.stat().st_size < 10000:
        print(f"    File too small — likely an error page")
        download_path.unlink(missing_ok=True)
        return None

    # Extract
    output_path = cache / f"bbj_{cancer_name}.tsv"
    try:
        if ext == "tar.gz":
            with tarfile.open(download_path, "r:gz") as tar:
                members = tar.getmembers()
                # Find the largest data file (skip READMEs etc)
                data_members = sorted(
                    [m for m in members if m.size > 1000],
                    key=lambda m: m.size,
                    reverse=True,
                )
                if not data_members:
                    data_members = members

                target = data_members[0]
                f = tar.extractfile(target)
                if f:
                    if target.name.endswith(".gz"):
                        # Write inner .gz to disk, then decompress
                        # (in-memory decompress fails CRC on large files)
                        import gzip
                        import shutil
                        inner_gz = cache / f"bbj_{cancer_name}_inner.gz"
                        inner_gz.write_bytes(f.read())
                        with gzip.open(inner_gz, "rb") as gz_in:
                            with open(output_path, "wb") as out:
                                shutil.copyfileobj(gz_in, out)
                        inner_gz.unlink()
                    else:
                        output_path.write_bytes(f.read())
                    print(f"    Extracted: {target.name} → {output_path.name}")

            download_path.unlink()
        elif ext == "gz":
            import gzip
            import shutil
            with gzip.open(download_path, "rb") as f_in:
                with open(output_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            download_path.unlink()
            print(f"    Extracted: {output_path.name}")
        return output_path

    except Exception as e:
        print(f"    Extraction error: {e}")
        download_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)
        return None


def read_bbj_sumstats(file_path: Path) -> pd.DataFrame:
    """
    Read a BBJ summary statistics file, handling multiple format variants.
    """
    # Detect separator and read
    with open(file_path) as f:
        header = f.readline()

    sep = "\t" if "\t" in header else " "
    df = pd.read_csv(file_path, sep=sep, low_memory=False)

    # Harmonize column names across BBJ format variants
    col_map = {
        # Format 1: newer standardized (prostate v2)
        "chromosome": "#chrom", "base_pair_location": "pos",
        "effect_allele": "alt", "other_allele": "ref",
        "beta": "beta", "standard_error": "sebeta",
        "effect_allele_frequency": "af_alt", "p_value": "pval",
        "variant_id": "rsids",
        # Format 2: older BBJ (renal, etc.)
        "Chromosome": "#chrom", "Position": "pos",
        "Rare allele": "alt", "Common allele": "ref",
        "Beta": "beta", "Standard error": "sebeta",
        "GWAS association P-value": "pval",
        "Variant": "rsids",
        "MAF, cases": "af_alt_cases", "MAF, non-Cases": "af_alt_controls",
        # Format 3: SAIGE output (BBJ large-scale GWAS)
        "p.value": "pval", "p.value.NA": "pval_na",
        "SNPID": "rsids",
        "Allele1": "ref", "Allele2": "alt",
        "AF_Allele2": "af_alt",
        "AF.Cases": "af_alt_cases", "AF.Controls": "af_alt_controls",
        "Tstat": "tstat",
        # Format 4: generic
        "CHR": "#chrom", "POS": "pos", "SNP": "rsids",
        "A1": "alt", "A2": "ref",
        "BETA": "beta", "SE": "sebeta", "P": "pval",
        "FRQ": "af_alt",
    }

    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Ensure required columns exist
    required = ["beta", "pval"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"    Warning: missing columns {missing}")
        print(f"    Available: {list(df.columns)}")
        return pd.DataFrame()

    # Convert types
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    if "sebeta" in df.columns:
        df["sebeta"] = pd.to_numeric(df["sebeta"], errors="coerce")
    if "#chrom" in df.columns:
        df["#chrom"] = df["#chrom"].astype(str)

    # Try to get gene names — BBJ files don't have nearest_genes
    # We'll map by rsid later if needed

    return df


def classify_effect(row):
    if row["beta"] < 0:
        return "protective"
    elif row["beta"] > 0:
        return "risk"
    return "ambiguous"


def fetch_bbj_cancer_sumstats(
    pval_threshold: float = PVAL_THRESHOLD,
    output_dir: str = "data/raw",
) -> pd.DataFrame:
    """Download and filter BBJ summary stats from JENGER."""
    all_results = []

    for cancer, endpoint in tqdm(JENGER_CANCER_ENDPOINTS.items(), desc="BBJ/JENGER"):
        print(f"\n  {cancer} ({endpoint['cases']} cases)...")
        file_path = download_from_jenger(cancer, endpoint)

        if file_path is None:
            continue

        df = read_bbj_sumstats(file_path)
        if df.empty:
            continue

        # Filter significant
        significant = df[df["pval"] <= pval_threshold].copy()
        if significant.empty:
            print(f"  No significant variants for {cancer}")
            continue

        significant["cancer_type"] = cancer
        significant["effect"] = significant.apply(classify_effect, axis=1)

        n_prot = (significant["effect"] == "protective").sum()
        n_risk = (significant["effect"] == "risk").sum()
        print(f"  {len(significant)} significant: {n_prot} protective, {n_risk} risk")

        all_results.append(significant)

    if not all_results:
        print("No BBJ data retrieved.")
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out / "bbj_significant_variants.csv", index=False)

    print(f"\nBBJ total: {len(combined)} significant variants")
    return combined


def extract_bbj_gene_sets(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """
    Extract gene lists from BBJ data.

    BBJ files typically don't have nearest_genes. We use rsid-based
    gene mapping via the GWAS Catalog or Ensembl, or fall back to
    matching rsids with FinnGen data (which does have gene annotations).
    """
    # BBJ uses GRCh37, so map variants to genes using GRCh37 gene annotations
    gene_annot_path = Path("data/raw/grch37_genes.csv")
    if gene_annot_path.exists() and "#chrom" in df.columns and "pos" in df.columns:
        genes = pd.read_csv(gene_annot_path, dtype={"chrom": str})

        # Build a lookup: for each chrom, find the nearest gene to each position
        print("  Mapping BBJ variants to nearest genes (GRCh37)...")
        gene_map = {}
        for chrom, chrom_genes in genes.groupby("chrom"):
            starts = chrom_genes["start"].values
            ends = chrom_genes["end"].values
            names = chrom_genes["gene"].values
            # Midpoints for distance calculation
            mids = (starts + ends) / 2
            gene_map[str(chrom)] = (starts, ends, mids, names)

        mapped_genes = []
        for _, row in df.iterrows():
            chrom = str(row["#chrom"])
            pos = row["pos"]
            if chrom in gene_map:
                starts, ends, mids, names = gene_map[chrom]
                # Check if inside a gene first
                inside = (pos >= starts) & (pos <= ends)
                if inside.any():
                    mapped_genes.append(names[inside][0])
                else:
                    # Nearest gene
                    dists = np.abs(mids - pos)
                    nearest_idx = np.argmin(dists)
                    if dists[nearest_idx] < 500_000:  # within 500kb
                        mapped_genes.append(names[nearest_idx])
                    else:
                        mapped_genes.append(None)
            else:
                mapped_genes.append(None)

        df["gene"] = mapped_genes
        mapped = df["gene"].notna().sum()
        print(f"  Mapped {mapped}/{len(df)} BBJ variants to genes via GRCh37 annotations")
    else:
        if "gene" not in df.columns:
            print("  Warning: no gene mapping available. Download GRCh37 annotations first.")
            return [], []

    protective = df[(df["effect"] == "protective") & df["gene"].notna()]
    risk = df[(df["effect"] == "risk") & df["gene"].notna()]

    prot_genes = protective["gene"].unique().tolist()
    risk_genes = risk["gene"].unique().tolist()

    print(f"  BBJ protective genes: {len(prot_genes)}")
    print(f"  BBJ risk genes:       {len(risk_genes)}")

    return prot_genes, risk_genes


if __name__ == "__main__":
    print("=" * 60)
    print("  BioBank Japan Cross-Population Validation")
    print("  Source: JENGER (jenger.riken.jp)")
    print("=" * 60)

    df = fetch_bbj_cancer_sumstats()
    if not df.empty:
        prot_genes, risk_genes = extract_bbj_gene_sets(df)
