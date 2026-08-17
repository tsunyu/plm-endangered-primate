# Setup

## 1. Analysis root

Scripts write under an **analysis root** that contains at least:

```text
data/     # genotypes, reference, phenotypes
output/   # phase products (can be populated by re-running pipelines)
works/    # publication products: figures/, tables/, results/ (created by figure scripts)
```

Set paths once:

```bash
bash configure_base_dir.sh /path/to/analysis_root
```

This writes `analysis/base_dir.env` (gitignored). Scripts read `PLM_BASE_DIR` from that file or from the environment.

Alternatively, export before each session:

```bash
export PLM_BASE_DIR=/path/to/analysis_root
```

## 2. Interpreters

| Software | Version | Notes |
| --- | --- | --- |
| Python | 3.8+ | `python3` on `PATH`. Use a venv (below). |
| R | 4.x | Needed only for phase 3a plotting (`Rscript`). |
| Java | 11+ | Needed only if SnpEff is invoked as `java -jar …`. |

Add the analysis utilities directory to `PYTHONPATH` when calling phase scripts that `import utils`:

```bash
export PYTHONPATH="$(pwd)/analysis:${PYTHONPATH}"
```

## 3. What you need by task

| Task | Python | R | External tools |
| --- | --- | --- | --- |
| `bash pipeline.sh --list` | — | — | — |
| Redraw figures / tables from precomputed `output/` | core packages + `openpyxl` (Excel) | — | analysis root with `output/` (S16 also needs `works/enrichment/`) |
| Phase 2 annotation | core + `pysam` | — | minimap2, Java/SnpEff, bcftools, bgzip, tabix |
| Phase 3a population genomics | core | ggplot2 stack | PLINK 1.9, bcftools, vcftools, KING, bgzip, tabix |
| Phase 3b demography | core | — | fastsimcoal2 (`fsc28` / `fsc27` / `fsc26`); GNU parallel recommended |
| Phase 4 ESM-2 scoring | core + `torch` + `transformers` | — | GPU recommended; optional `ESM2_LOCAL_PATH` |
| Phase 5 genetic load | core | — | — |
| Phase 8 phenotype / GWAS | core | — | PLINK 1.9, **plink2**, GCTA, GEMMA |
| Phase 9 validation / DFE sims | core + `scikit-learn` | — | (optional) SLiM |

## 4. Python packages

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r analysis/requirements.txt
```

Pinned minima live in [`analysis/requirements.txt`](analysis/requirements.txt). Summary:

| Package | Used for |
| --- | --- |
| `numpy`, `pandas`, `scipy` | Numeric work in phases and figure redraws |
| `matplotlib`, `seaborn` | Phase plots and publication redraws |
| `pyyaml` | `analysis/utils.py` / `config.yaml` |
| `statsmodels` | Phenotype models (phase 8) |
| `scikit-learn` | ClinVar / fitness validation (phase 9) |
| `openpyxl` | `figures_tables/export_supplementary_tables_xlsx.py` |
| `pysam` | Ancestral FASTA access (phase 2) |
| `torch`, `transformers` | ESM-2 (phase 4 only) |

ESM-2 scoring needs a GPU-capable PyTorch build; install the wheel for your CUDA/CPU stack **before** `transformers`. Optional: `export ESM2_LOCAL_PATH=/path/to/esm2_t33_650M_UR50D` to skip Hugging Face download (`facebook/esm2_t33_650M_UR50D`).

## 5. R packages (phase 3a plotting)

```bash
Rscript analysis/install_r_packages.R
```

List: [`analysis/r_requirements.txt`](analysis/r_requirements.txt). These are the packages actually imported:

| Package | Used in |
| --- | --- |
| `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork` | `phase3a_generate_plots.R` and inline plots in phase 3a shell scripts |
| `yaml`, `data.table` | `analysis/utils.R` |

Gene-set enrichment for Supplementary Fig. S16 uses **precomputed g:Profiler tables** on the analysis root (`$PLM_BASE_DIR/works/enrichment/`); no Bioconductor enrichment packages are required. Those tables are not stored in this repository.

## 6. External bioinformatics tools

Install binaries on `PATH` unless an override is listed. Commands in parentheses are what the scripts call.

| Tool | Command | Used in | Override |
| --- | --- | --- | --- |
| PLINK 1.9 | `plink` | Phase 3a (QC helpers, ROH, PCA); phase 8 | — |
| PLINK 2 | `plink2` | Phase 8 GRM / genotype extraction | — |
| bcftools | `bcftools` | VCF handling; `bcftools roh` in phase 3a | — |
| htslib | `bgzip`, `tabix` | VCF compression and indexing | — |
| VCFtools | `vcftools` | Nucleotide diversity / Tajima’s *D* (phase 3a) | — |
| KING | `king` | Kinship (phase 3a structure) | — |
| minimap2 | `minimap2` | Outgroup alignments for ancestral alleles (phase 2) | — |
| SnpEff | `snpEff` or `java -jar snpEff.jar` | Variant annotation (phase 2) | `SNPEFF_BIN`, `SNPEFF_JAR` |
| fastsimcoal2 | `fsc28` (fallback `fsc27`, `fsc26`) | Demographic inference (phase 3b) | — |
| GCTA | `gcta64` or `gcta` | Mixed-model GWAS (phase 8) | `GCTA_BIN` |
| GEMMA | `gemma` | Mixed-model GWAS (phase 8) | `GEMMA_BIN` |
| GNU parallel | `parallel` | Recommended for fastsimcoal2 parallel/bootstrap | — |
| SLiM | `slim` | Optional DFE cross-check (phase 9) | — |

SnpEff expects a custom database named `Rhinopithecus_roxellana_ASM756505v1` (see `analysis/phase2_annotation/README.md`).

## 7. Environment variables

| Variable | Purpose |
| --- | --- |
| `PLM_BASE_DIR` | Analysis root (`data/`, `output/`, publication-products `works/`) |
| `PYTHONPATH` | Include `analysis/` so `import utils` works |
| `ESM2_LOCAL_PATH` | Local ESM-2 checkpoint (phase 4) |
| `GCTA_BIN` | GCTA executable if not on `PATH` |
| `GEMMA_BIN` | GEMMA executable if not on `PATH` |
| `SNPEFF_JAR` | Path to `snpEff.jar` |
| `SNPEFF_BIN` | Full SnpEff command (e.g. `java -jar /path/to/snpEff.jar`) |

## 8. Memory

Large VCF scans and ESM-2 inference are memory-intensive. See `analysis/memory_config.sh` for suggested limits on shared servers.

## 9. Figure outputs

Publication redraws write to the analysis-root publication-products directory:

```text
$PLM_BASE_DIR/works/figures/
$PLM_BASE_DIR/works/figures/supplementary/
$PLM_BASE_DIR/works/tables/
$PLM_BASE_DIR/works/results/          # e.g. polarized DAF cache
```
