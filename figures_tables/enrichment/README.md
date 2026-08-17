# Enrichment tables for Supplementary Fig. S16

This code repository does **not** ship gene lists or g:Profiler result tables.

Supp. Fig. S16 is redrawn by `plot_supplementary.py` from precomputed files on the analysis root:

```text
$PLM_BASE_DIR/works/enrichment/
```

Place at least:

| File | Used for |
| --- | --- |
| `gprofiler_top10pct_missense_plot_terms.csv` | Panel a (required) |
| `gprofiler_lof_plot_terms.csv` | Panel b (preferred) |
| `gprofiler_lof_informative_terms.csv` | Panel b fallback if the plot-terms file is empty |

Those CSVs are g:Profiler FDR tables for top-decile missense (*P* ≥ 0.612) and LoF gene sets. They belong on the analysis root, not in this GitHub checkout.
