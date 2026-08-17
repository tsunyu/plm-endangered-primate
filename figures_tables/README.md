# Publication figure and table scripts

Redraws Figures 2–6 and Supplementary Figures S1–S16 from precomputed `output/` products.

```bash
# from the package root
export PLM_BASE_DIR=/path/to/analysis_root   # directory with data/ and output/
cd figures_tables
python3 run_all_figures.py          # Figs 2–6 + Supp. S1–S16
python3 plot_supplementary.py       # Supp. Figs S1–S16 only
python3 regenerate_tables.py
python3 export_supplementary_tables_xlsx.py
```

Figure 1 is not generated here. Polarized DAF logic lives in `lib/daf.py`.

## Outputs

Written under `$PLM_BASE_DIR/works/` (publication-products directory on the analysis root):

- `figures/` — main Figures 2–6
- `figures/supplementary/` — Supp. Figs S1–S16
- `tables/` — supplementary tables / Excel workbook
- `results/` — caches (e.g. polarized DAF)

## Enrichment (Supp. Fig. S16)

Inputs are **not** in this repository. The redraw reads
`$PLM_BASE_DIR/works/enrichment/` (see `enrichment/README.md`).
