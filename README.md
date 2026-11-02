# CRISPR-targeted macrophage bulk RNA-seq analysis

This repository-ready script reproduces the transcriptomic analyses reported in
the manuscript: count import and metadata construction, donor-paired DESeq2
differential expression, global PCA, donor-adjusted pairwise PCA, Hallmark GSEA,
and PROGENy pathway inference followed by donor-paired limma testing.

Exploratory analyses that were not used in the manuscript are intentionally not
included (quanTIseq, edgeR differential-expression cross-checks, custom TAM
modules, Reactome GSEA, contamination checks, CAMERA, ROAST, and FRY).

## Repository layout

```text
.
├── analysis.R
├── README.md
├── data/
│   └── gene_count.xls
├── results/          # created by the script
└── figures/          # created by the script
```

`gene_count.xls` is the original tab-delimited count table despite its `.xls`
extension. It must contain `gene_id`, `gene_name`, and the 28 raw-count sample
columns described below. The script never modifies this input.

## Samples and naming

Four donors are expected: `HD164`, `HD180`, `HD183`, and `HD181`. Each donor
has seven samples in this order:

```text
NOEP, NTC, CLEC4E, CLEC6A, CD209, SELL, SIG9
```

The raw matrix column suffix is `SIG9`, for example `HD164_SIG9`. The biological
condition is consistently reported as `SIGLEC9` in metadata, contrasts, tables,
and figures. The script performs this mapping explicitly.

## Package requirements

Use a current R release and install the following packages.

Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "limma", "edgeR", "decoupleR", "clusterProfiler"
))
```

CRAN packages:

```r
install.packages(c(
  "msigdbr", "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel",
  "stringr", "readr"
))
```

`edgeR` is used only to filter/TMM-normalize counts before `limma::voom` for the
PROGENy input matrix. It is not used as a second differential-expression test.
The first PROGENy run may require internet access to obtain the human PROGENy
network through `decoupleR`; the retrieved network is exported for provenance.

## Running the analysis

From the repository root:

```sh
Rscript analysis.R
```

The script stops early if required columns, packages, integer counts, or sample
order are missing or inconsistent. Randomized steps use seed `123`. Package and
R versions are recorded in `results/sessionInfo.txt`.

## Analysis details

- DESeq2 uses the paired design `~ donor + condition`. Genes are retained when
  they have at least 10 counts in at least four samples. The five targeted
  conditions are contrasted with NTC; NTC is also contrasted with NOEP.
- The global PCA uses the 500 most variable VST-transformed genes.
- Each targeted-condition-versus-NTC PCA uses its own 500 most variable genes.
  Donor-associated variation is removed with `limma::removeBatchEffect` while
  preserving the condition term. Plot polygons are 95% normal ellipses.
- Hallmark GSEA ranks genes by the signed DESeq2 Wald statistic. The combined
  manuscript plot retains pathways significant in at least two targeted
  conditions or FDR-significant pathways in the top decile of absolute NES.
- PROGENy uses the top 500 human model genes per pathway and `run_mlm`. Inferred
  activities are tested with the paired model `~ donor + condition` in limma.
  An asterisk in the heatmap denotes nominal `p < 0.05`; full nominal and
  multiple-testing-adjusted values are retained in the result table.

## Main outputs

| Manuscript analysis | Output |
|---|---|
| Paired differential expression | `results/deseq2/` |
| Global PCA (Supplementary Fig. S5A) | `figures/pca/global_PCA.*` |
| CLEC4E pairwise PCA (Fig. 5A) | `figures/pca/CLEC4E_vs_NTC_donor_adjusted_PCA.*` |
| Other pairwise PCA panels | `figures/pca/*_vs_NTC_donor_adjusted_PCA.*` |
| Hallmark GSEA tables | `results/hallmark_gsea/` |
| Hallmark manuscript dot plot | `figures/hallmark_gsea/Hallmark_GSEA_combined_dotplot.*` |
| PROGENy/limma results | `results/progeny/limma_results.csv` |
| PROGENy heatmap (Supplementary Fig. S5G) | `figures/progeny/PROGENy_pathway_activity.*` |

Both PDF and 600-dpi PNG versions are written for manuscript figures. PCA
coordinates and figure-level plotting data are exported as CSV files so that
the displayed results can be checked without rerunning the full analysis.

## Reproducibility and reuse

The script uses project-relative paths only and contains no private filesystem
locations or interactive calls. `results/analysis_objects.rds` stores the fitted
DESeq2 object, transformed expression data, and final Hallmark and PROGENy result
tables. For a permanent public release, archive the repository version together
with its `sessionInfo.txt` and cite the resulting release DOI.
