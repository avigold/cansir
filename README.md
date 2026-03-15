# CANSIR: Cancer-Protective Gene Discovery

Analysis pipeline for the paper: **"The genetic architecture of cancer protection is biologically distinct from cancer risk: a pan-cancer pathway enrichment analysis"**

## Key finding

Cancer-protective genetic variants (GWAS β < 0) enrich for biologically distinct pathways compared to cancer-risk variants (β > 0). This holds across 19 cancer types in FinnGen (N=500,348) and replicates in BioBank Japan (N~200,000).

- **Fisher's exact P = 8.57×10⁻²⁴** (pathway profiles are distinct)
- **Jaccard similarity = 0.320** (68% of pathways are unique to one set)
- **PPS ≠ inverted PRS** (gene Jaccard = 0.217, Spearman ρ = −0.306)

## Pipeline

```bash
# Install dependencies
pip install -r requirements.txt

# Run the full analysis
python run_pipeline.py
```

The pipeline:
1. Downloads FinnGen R12 GWAS summary statistics (19 cancer endpoints)
2. LD-clumps and separates protective from risk variants
3. Runs pathway enrichment via g:Profiler (GO, Reactome, KEGG)
4. Statistically compares pathway profiles (Fisher's exact, Jaccard, Spearman)
5. Runs sensitivity analysis across P-value thresholds
6. Cross-population validation using BioBank Japan (via JENGER)

## Polygenic Protection Score (PPS) Calculator

```bash
# Score an individual's genotype
PYTHONPATH=. python src/analysis/pps_calculator.py <genotype_file>
```

Accepts 23andMe, AncestryDNA, or VCF files. Outputs an HTML report with per-pathway protection percentiles and intervention recommendations.

## Project structure

```
src/
├── data/
│   ├── finngen.py              # FinnGen R12 data download and processing
│   ├── bbj.py                  # BioBank Japan cross-population validation
│   ├── open_targets.py         # Gene pathway/function features
│   └── gwas_catalog.py         # GWAS Catalog associations
├── analysis/
│   ├── enrichment.py           # g:Profiler pathway enrichment
│   ├── compare_pathways.py     # Jaccard, Fisher's exact, Spearman comparison
│   ├── visualization.py        # Publication figures
│   ├── ld_clumping.py          # Distance-based LD clumping
│   ├── pps_calculator.py       # Polygenic Protection Score calculator
│   ├── protection_score.py     # PPS vs PRS comparison
│   ├── mendelian_randomization.py  # MR analysis
│   ├── cross_population.py     # FinnGen-BBJ comparison
│   ├── interventions.py        # Drug-gene interaction mapping (DGIdb)
│   └── imputation.py           # Genotype imputation integration
├── features/
│   └── build_features.py       # ML feature matrix
└── models/
    └── classifier.py           # Semi-supervised gene classifier
paper/
└── preprint_draft.md           # Manuscript
```

## Data sources

- **FinnGen R12**: https://r12.finngen.fi/ (public, no application required)
- **BioBank Japan**: http://jenger.riken.jp/ (public, no application required)
- **Drug-gene interactions**: DGIdb GraphQL API (https://dgidb.org/)

## Citation

Preprint forthcoming on bioRxiv.

## Licence

MIT
