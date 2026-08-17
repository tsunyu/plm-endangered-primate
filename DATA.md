# Data layout

This repository ships **code only**. It does not contain genotypes, phenotypes, gene lists, enrichment tables, figures, or other analysis products.

Place restricted inputs and precomputed products on a local **analysis root** (`$PLM_BASE_DIR`), not in this GitHub checkout.

## Restricted inputs (not in GitHub)

Place these under `$PLM_BASE_DIR/data/` (names match the analysis scripts):

| Path (typical) | Description |
| --- | --- |
| `data/monkey_snp_sex_qc.vcf(.gz)` | Cohort VCF after sex/QC filters (phases 2, 3a, 8) |
| `data/monkey_snp_sex_qc.{bed,bim,fam}` | PLINK files for the same set |
| `data/hardfilted.snp.pass.autosomes.vcf.gz` | Autosome-only VCF used by phase 3b SFS |
| `data/reference/` | *R. roxellana* genome, GFF/GTF, proteins |
| `data/disease_record_noredundancy.csv` | Field morbidity records used to build CHS |

`filter_autosomes_only.sh` can produce `data/monkey_snp_autosomes_only.vcf.gz` from the sex/QC VCF; phase 3b step scripts currently read `hardfilted.snp.pass.autosomes.vcf.gz`.

Access to individual-level genomic and health data is governed by the study’s ethics approvals and data-use agreements. Contact the corresponding authors for collaboration or controlled access.

## Precomputed analysis products

Re-drawing manuscript figures expects products under `$PLM_BASE_DIR/output/`, including (non-exhaustive):

```text
output/phase2_annotation/
output/phase3a_population_genomics/
output/phase3b_fastsimcoal2/
output/phase4_plm_predictions/
output/phase5_genetic_load/
output/phenotype_genotype_analysis/
output/method_validation/
output/allele_frequency_spectrum/
```

These directories are produced by the phase pipelines in `analysis/` and are the precomputed analysis products used to redraw the paper figures.

Supplementary Fig. S16 additionally reads g:Profiler tables from `$PLM_BASE_DIR/works/enrichment/` (the publication-products directory on the analysis root; not from this repository). See [figures_tables/enrichment/README.md](figures_tables/enrichment/README.md).

## What this repository does contain

| Path | Contents |
| --- | --- |
| `analysis/` | Phase pipelines and helpers |
| `analysis/phase3b_fastsimcoal2/models/` | Demographic model templates (`.tpl` / `.est`) and fitted `.par` files from this study’s SFS |
| `analysis/phase3b_fastsimcoal2/seed.txt` | RNG seed used for a fastsimcoal2 run |
| `analysis/config.yaml.example` | Path template |
| `figures_tables/` | Publication redraw scripts |

The `.par` files are inferred demographic parameters (not genotypes). They can be used to reproduce forward simulations under the fitted history (Fig. 6d).

## Key methodological constants (manuscript)

| Parameter | Value |
| --- | --- |
| Missense pathogenicity sigmoid | *P* = 1 / [1 + exp{0.5287(LLR + 6.8920)}] |
| Top-decile deleterious missense | *P* ≥ 0.612 |
| LoF weight | *P* = 0.95 |
| Dominance coefficient *h* | 0.25 |
| Generation time (demography) | 10 years |
| Mutation rate (fastsimcoal2) | 1.36 × 10⁻⁸ (as in Methods) |
