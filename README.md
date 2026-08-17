# Protein language modelling reveals a latent drift load underlying health risk in an endangered primate

Analysis and figure code accompanying the study of genetic load, inbreeding architecture, and morbidity in the Shennongjia (SNJ) population of the golden snub-nosed monkey (*Rhinopithecus roxellana*).

**Paper:** *Protein language modelling reveals a latent drift load underlying health risk in an endangered primate*

**Repository:** [https://github.com/tsunyu/plm-endangered-primate](https://github.com/tsunyu/plm-endangered-primate)

This repository contains **code only** (analysis and figure scripts). It does not ship genotypes, phenotypes, gene lists, or other analysis tables. Full re-runs and figure redraws read local data under `$PLM_BASE_DIR` (see [DATA.md](DATA.md) and [SETUP.md](SETUP.md)).

---

## Highlights

- Zero-shot **ESM-2** scoring of 9,725 missense variants on the monkey’s own protein sequences
- ClinVar-calibrated pathogenicity weights and dominance-aware **genetic load**
- Polarized derived-allele frequency spectra for LoF and top-decile missense variants
- **fastsimcoal2** demography and forward DFE simulations under the inferred *N*<sub>e</sub> history
- Morbidity associations (CHS) beyond *F*<sub>ROH</sub>, plus exploratory GWAS and burden-gene enrichment

## Repository layout

```text
.
├── README.md                 # this file
├── SETUP.md                  # environment, BASE_DIR, dependencies
├── DATA.md                   # inputs / outputs expected on disk
├── LICENSE
├── MANIFEST.csv              # step ↔ display-item map
├── pipeline.sh               # list / optionally run paper-order steps
├── configure_base_dir.sh     # writes analysis/base_dir.env (gitignored)
├── analysis/                 # phase pipelines (annotation → validation)
│   ├── config.yaml           # paths use /path/to/analysis_root placeholders
│   ├── config.yaml.example
│   ├── base_dir.env.example
│   ├── requirements.txt
│   ├── utils.py / utils.R
│   └── phase*/               # analysis modules (directory names follow the paper pipeline)
└── figures_tables/           # publication redraws (Figs 2–6, Supp. S1–S16)
    ├── run_all_figures.py
    ├── plot_figure0*.py
    ├── plot_supplementary.py
    ├── regenerate_tables.py
    ├── export_supplementary_tables_xlsx.py
    ├── enrichment/README.md  # where S16 tables live on the analysis root
    └── lib/                  # shared figure utilities (paths, DAF, style)
```

Figure 1 (study schematic) is **manual artwork** and has no generator script here.

## Quick start

```bash
# 1. Point every script at your analysis root (directory with data/ and output/)
bash configure_base_dir.sh /path/to/your/analysis_root

# 2. Install Python deps (R + command-line tools: SETUP.md)
pip install -r analysis/requirements.txt

# 3. List paper-order steps (does not execute by default)
bash pipeline.sh --list

# 4. Redraw publication figures/tables only (needs precomputed output/ + works/)
bash pipeline.sh --run-figures
```

Expected layout of the analysis root after configuration:

```text
analysis_root/
├── data/          # VCF / PLINK / reference (not in this repo)
├── output/        # precomputed phase products (not in this repo)
└── works/         # publication products: figures/, tables/, results/, enrichment/
```

## Dependencies

Full inventory (Python, R, and command-line tools, mapped to each phase): **[SETUP.md](SETUP.md)**.

- Python: `pip install -r analysis/requirements.txt`
- R (phase 3a plots only): `Rscript analysis/install_r_packages.R`
- External binaries (PLINK, bcftools, SnpEff, fastsimcoal2, GCTA, GEMMA, …): [SETUP.md §6](SETUP.md)

## Analysis order

| Step | Module | Main display items |
| --- | --- | --- |
| 1 | Annotation (ancestral states, SnpEff) | Supp. Fig. S2; Table S12 |
| 2 | Population genomics (diversity, ROH, structure) | Figs 2b–d; Supp. Figs S1, S3; Tables S1–S2 |
| 3 | Demography (fastsimcoal2) | Fig. 6a–c; Supp. Fig. S4; Tables S8–S9 |
| 4 | ESM-2 missense scoring | Fig. 3a; Supp. Fig. S5; Table S3 |
| 5 | Genetic load | Supp. Figs S6–S7; Table S4 |
| 6 | Phenotype / GWAS | Fig. 4a,c; Supp. Figs S8, S15; Tables S6, S11 |
| 7 | Method validation / DFE sims | Figs 3b–c, 4b,d, 6d; Supp. Figs S9–S14; Tables S5, S10 |
| 8–10 | Publication figures / tables / Excel | Figs 2–6; Supp. Figs S1–S16; Tables S1–S12 |

Directory names under `analysis/` keep the original module numbers (phase2–5, 8–9). Unused modules are not shipped. Run order is the table above / `pipeline.sh --list`, not the folder numbers.

Full mapping: [`MANIFEST.csv`](MANIFEST.csv).

## Citation

If you use this code, please cite the paper:

> Protein language modelling reveals a latent drift load underlying health risk in an endangered primate

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). Underlying genomic and phenotype data remain subject to the data-use terms of the study and collaborating institutions (see [DATA.md](DATA.md)).

## Contact

Questions about the code package: open an issue on GitHub or contact the corresponding authors listed in the manuscript.
