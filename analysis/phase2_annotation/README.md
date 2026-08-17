# Phase 2: Ancestral State and Variant Annotation

## Overview

Phase 2 performs comprehensive variant annotation by inferring ancestral alleles and adding functional annotations to genetic variants. This phase integrates multiple annotation sources to provide a complete picture of variant effects on genes and proteins.

## Pipeline Structure

The annotation pipeline consists of the following steps, executed in sequence:

```
phase2_step0_annotation_pipeline.sh (Main Pipeline)
├── Step 1: Ancestral State Inference
│   └── phase2_step1_ancestral_states.sh
├── Step 2: SnpEff Functional Annotation
│   └── phase2_step2_snpeff_annotation.sh
├── Step 3: Functional Annotation Extraction
│   └── phase2_step3_functional_annotation.sh
└── Step 4: Visualization
    ├── phase2_step4_create_figures.sh
    ├── phase2_step4.1_visualize_snpeff.py
    └── phase2_step4.2_visualize_functional.py
```

## Scripts

### Main Pipeline

#### `phase2_step0_annotation_pipeline.sh`
**Purpose**: Orchestrates the complete annotation pipeline with checkpoint management.

**Features**:
- Checkpoint-based execution (resume from interruptions)
- Memory monitoring
- Automatic step skipping for completed analyses
- Progress logging

**Usage**:
```bash
bash phase2_step0_annotation_pipeline.sh [OPTIONS]
```

**Options**:
- `--force`: Clear all checkpoints and force rerun all steps
- `--clear-checkpoints`: Clear all checkpoints before starting
- `--resume-from STEP`: Resume from a specific step

**Example**:
```bash
# Normal run with checkpointing
bash phase2_step0_annotation_pipeline.sh

# Force rerun all steps
bash phase2_step0_annotation_pipeline.sh --force

# Resume from SnpEff annotation
bash phase2_step0_annotation_pipeline.sh --resume-from "SnpEff Annotation"
```

### Step 1: Ancestral State Inference

#### `phase2_step1_ancestral_states.sh`
**Purpose**: Infers ancestral alleles using phylogenetic outgroups.

**Method**:
- Uses minimap2 for fast genome-to-genome alignment
- Parsimony-based inference with multiple outgroup species
- Maps reference positions to outgroup coordinates using PAF CIGAR strings

**Input**:
- Reference genome: `data/reference/Rhinopithecus_roxellana.fna`
- Outgroup genomes:
  - `data/reference/Pygathrix_nemaeus.fna` (Primary: Douc langur)
  - `data/reference/Rhinopithecus_bieti.fna` (Secondary: Black snub-nosed monkey)
  - `data/reference/Macaca_mulatta.fna` (Tertiary: Rhesus macaque)
- VCF file: `data/monkey_snp_sex_qc.vcf`

**Output**:
- `output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz`: VCF with ancestral allele annotations (AA, AA_CONF, AA_METHOD)
- `output/phase2_annotation/ancestral_states/ancestral_summary.txt`: Summary statistics
- Alignment files: `alignment_outgroup{1,2,3}.paf`

**Usage**:
```bash
bash phase2_step1_ancestral_states.sh [--clear-checkpoints]
```

**Key Parameters**:
- `MIN_OUTGROUPS_AGREE=2`: Minimum outgroups with same allele for confidence
- `MIN_ALIGNMENT_QUALITY=20`: Minimum alignment quality threshold
- `THREADS=32`: Number of parallel threads

### Step 2: SnpEff Functional Annotation

#### `phase2_step2_snpeff_annotation.sh`
**Purpose**: Annotates variants with functional consequences using SnpEff.

**Features**:
- Gene region identification
- Coding effect prediction (missense, nonsense, frameshift, etc.)
- Loss-of-function (LOF) variant identification
- Impact categorization (HIGH, MODERATE, LOW, MODIFIER)

**Input**:
- VCF with ancestral states: `output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz`
- Falls back to: `data/monkey_snp_sex_qc.vcf` if ancestral VCF not available

**Output**:
- `output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz`: SnpEff annotated VCF
- `output/phase2_annotation/snpeff_annotation/snpeff_stats.html`: SnpEff statistics report
- `output/phase2_annotation/snpeff_annotation/annotation_summary.txt`: Summary of annotations
- `output/phase2_annotation/snpeff_annotation/high_impact_variants.csv`: HIGH impact variants
- `output/phase2_annotation/snpeff_annotation/lof_variants.csv`: Loss-of-function variants

**Usage**:
```bash
bash phase2_step2_snpeff_annotation.sh
```

**Prerequisites**:
- SnpEff database: `Rhinopithecus_roxellana_ASM756505v1` (custom built)
- SnpEff on PATH, or `SNPEFF_JAR=/path/to/snpEff.jar` / `SNPEFF_BIN="java -jar /path/to/snpEff.jar"`

**Key Parameters**:
- `SNPEFF_DB="Rhinopithecus_roxellana_ASM756505v1"`: Database name
- `THREADS=8`: Number of threads
- `JAVA_OPTS="-Xmx80g -Xms4g"`: Java memory settings

### Step 3: Functional Annotation Extraction

#### `phase2_step3_functional_annotation.sh`
**Purpose**: Extracts and classifies functional annotations from SnpEff results.

**Outputs**:
- Variant-level summaries:
  - `v2_variant_counts_by_impact.tsv`: Counts by impact category
  - `v2_variant_counts_by_effect.tsv`: Counts by effect type
  - `v2_variant_counts_by_category.tsv`: Counts by functional category
  - `v2_high_impact_variants.csv`: HIGH impact variants
  - `v2_lof_variants.csv`: Loss-of-function variants

- Gene-level summaries:
  - `v2_gene_summary.tsv`: Comprehensive gene-level statistics
  - `v2_genes_{high,moderate,low,modifier}_{all,known,predicted}.txt`: Gene lists by impact
  - `v2_biotype_counts.tsv`: Gene counts by biotype
  - `v2_gene_category_counts.tsv`: Gene counts by functional category

- Optional visualizations:
  - `v2_impact_distribution.png`
  - `v2_top_effects.png`
  - `v2_category_distribution.png`
  - `v2_gene_counts_by_impact.png`
  - `v2_biotype_vs_impact_genes.png`
  - `v2_top_genes_by_variant_count.png`
  - `v2_top_genes_by_high_impact.png`

**Usage**:
```bash
bash phase2_step3_functional_annotation.sh
```

**Input**:
- SnpEff annotated VCF: `output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz`

**Output Directory**:
- `output/phase2_annotation/functional_annotation/`

### Step 4: Visualization

#### `phase2_step4_create_figures.sh`
**Purpose**: Unified script to generate all annotation visualization figures.

**Usage**:
```bash
bash phase2_step4_create_figures.sh [OPTIONS]
```

**Options**:
- `--snpeff-only`: Generate only SnpEff figures
- `--functional-only`: Generate only functional annotation figures
- `--dpi DPI`: Figure DPI (default: 300)
- `--format FORMAT`: Figure format: png, pdf, svg (default: png)

**Example**:
```bash
# Generate all figures
bash phase2_step4_create_figures.sh

# Generate only SnpEff figures with high DPI
bash phase2_step4_create_figures.sh --snpeff-only --dpi 600

# Generate PDF format
bash phase2_step4_create_figures.sh --format pdf
```

#### `phase2_step4.1_visualize_snpeff.py`
**Purpose**: Creates publication-quality figures from SnpEff annotation results.

**Features**:
- Impact distribution plots (bar and pie charts)
- Top effect types visualization
- High-impact genes analysis
- Summary statistics panels

**Usage**:
```bash
python3 phase2_step4.1_visualize_snpeff.py [OPTIONS]
```

**Options**:
- `--input-vcf PATH`: Path to annotated VCF file
- `--output-dir PATH`: Output directory for figures
- `--dpi INT`: Figure DPI (default: 300)
- `--format FORMAT`: Figure format: png, pdf, svg (default: png)
- `--top-effects INT`: Number of top effect types to show (default: 20)
- `--top-genes INT`: Number of top genes to show (default: 30)

#### `phase2_step4.2_visualize_functional.py`
**Purpose**: Creates publication-quality figures from functional annotation gene lists.

**Features**:
- Gene counts by impact category
- Comparison across gene classification levels (all, genes, known)
- Gene category distribution
- Biotype distribution (top N biotypes)

**Usage**:
```bash
python3 phase2_step4.2_visualize_functional.py [OPTIONS]
```

**Options**:
- `--input-dir PATH`: Directory with gene list files
- `--output-dir PATH`: Output directory for figures (default: same as input)
- `--gene-type TYPE`: Gene type to visualize: all, genes, known (default: genes)
- `--dpi INT`: Figure DPI (default: 300)
- `--format FORMAT`: Figure format: png, pdf, svg (default: png)
- `--all-types`: Generate plots for all gene types
- `--comparison`: Generate comparison plot across gene types
- `--biotype-top N`: Show top N biotypes (default: 15)

### Preparation Script

#### `phase2_prep_extract_gene_protein_mapping.py`
**Purpose**: Extracts gene and protein ID/name mappings from GFF annotation file.

**Outputs**:
- `gene_info.csv/tsv`: Gene ID, Gene Name, Description
- `transcript_protein_mapping.csv/tsv`: Transcript ID → Protein ID
- `gene_transcript_protein.csv/tsv`: Complete Gene → Transcript → Protein mapping
- `protein_info.csv/tsv`: Protein ID, Protein Name, Description

**Usage**:
```bash
python3 phase2_prep_extract_gene_protein_mapping.py [OPTIONS]
```

**Options**:
- `--gff FILE`: GFF file path (default: `data/reference/Rhinopithecus_roxellana.gff`)
- `--output DIR`: Output directory (default: `output/gene_protein_mapping`)
- `--format FORMAT`: Output format: csv, tsv, or both (default: both)

**Example**:
```bash
# Extract with defaults
python3 phase2_prep_extract_gene_protein_mapping.py

# Specify custom GFF file
python3 phase2_prep_extract_gene_protein_mapping.py --gff /path/to/annotation.gff

# Output only CSV format
python3 phase2_prep_extract_gene_protein_mapping.py --format csv
```

## Input Requirements

### Required Files
1. **VCF file**: `data/monkey_snp_sex_qc.vcf` (or gzipped)
2. **Reference genome**: `data/reference/Rhinopithecus_roxellana.fna`
3. **GFF annotation**: `data/reference/Rhinopithecus_roxellana.gff`
4. **Outgroup genomes** (for ancestral state inference):
   - `data/reference/Pygathrix_nemaeus.fna`
   - `data/reference/Rhinopithecus_bieti.fna`
   - `data/reference/Macaca_mulatta.fna`

### Software Dependencies
- **minimap2**: For genome alignment (Step 1)
- **SnpEff**: For functional annotation (Step 2)
- **Python 3.8+**: With packages: pandas, numpy, matplotlib, seaborn, pysam
- **bcftools/tabix**: For VCF indexing and manipulation

## Output Structure

```
output/phase2_annotation/
├── .checkpoints/                    # Checkpoint files for resuming
│   ├── ancestral_state_inference.done
│   ├── snpeff_annotation.done
│   └── additional_functional_annotation.done
├── ancestral_states/
│   ├── variants_with_ancestral.vcf.gz
│   ├── ancestral_summary.txt
│   └── alignment_outgroup{1,2,3}.paf
├── snpeff_annotation/
│   ├── annotated_variants.vcf.gz
│   ├── snpeff_stats.html
│   ├── annotation_summary.txt
│   ├── high_impact_variants.csv
│   └── lof_variants.csv
├── functional_annotation/
│   ├── v2_variant_counts_by_impact.tsv
│   ├── v2_variant_counts_by_effect.tsv
│   ├── v2_gene_summary.tsv
│   ├── v2_genes_*.txt
│   └── *.png (optional visualizations)
└── phase2_pipeline.log
```

## Workflow

### Complete Pipeline Run
```bash
cd analysis/phase2_annotation
bash phase2_step0_annotation_pipeline.sh
```

### Individual Step Execution

#### Step 1: Ancestral States
```bash
bash phase2_step1_ancestral_states.sh
```

#### Step 2: SnpEff Annotation
```bash
bash phase2_step2_snpeff_annotation.sh
```

#### Step 3: Functional Annotation
```bash
bash phase2_step3_functional_annotation.sh
```

#### Step 4: Visualization
```bash
bash phase2_step4_create_figures.sh
```

## Checkpoint System

The pipeline uses a checkpoint system to allow resuming interrupted runs:

- **Checkpoint location**: `output/phase2_annotation/.checkpoints/`
- **Checkpoint files**: Each step creates a `.done` file with completion timestamp
- **Automatic resumption**: Re-running the pipeline automatically skips completed steps
- **Force rerun**: Use `--force` flag to clear checkpoints and rerun all steps

## Memory Management

- **Memory monitoring**: Pipeline checks memory usage every 60 seconds
- **Memory limit**: 120GB maximum (reserves 8GB for system)
- **SnpEff memory**: 80GB Java heap space
- **Warnings**: Pipeline will warn if memory usage exceeds thresholds

## Troubleshooting

### SnpEff Database Not Found
**Error**: `Database Rhinopithecus_roxellana_ASM756505v1 not found`

**Solution**:
1. Ensure SnpEff database is built and available
2. Check database location: `snpeff config` shows data directory
3. Verify database name matches `SNPEFF_DB` in script

### Outgroup Genomes Missing
**Error**: `Outgroup genome not found`

**Solution**:
1. Download required outgroup genomes
2. Update paths in `phase2_step1_ancestral_states.sh`:
   - `OUTGROUP1`, `OUTGROUP2`, `OUTGROUP3` variables
3. Ensure genomes are in FASTA format

### Memory Issues
**Error**: `High memory usage detected`

**Solution**:
1. Reduce parallel threads (`THREADS` parameter)
2. Increase system swap space
3. Process chromosomes/scaffolds in batches
4. Close other memory-intensive applications

### Checkpoint Issues
**Problem**: Pipeline not resuming from checkpoint

**Solution**:
1. Check checkpoint files exist: `ls output/phase2_annotation/.checkpoints/`
2. Clear checkpoints: `bash phase2_step0_annotation_pipeline.sh --clear-checkpoints`
3. Force rerun: `bash phase2_step0_annotation_pipeline.sh --force`

## Integration with Other Phases

### Input
- Quality-controlled VCF file (`data/monkey_snp_sex_qc.vcf`)

### Output for subsequent steps
- **ESM-2 scoring (phase 4)**: Uses annotated variants for protein sequence preparation
- **Genetic load (phase 5)**: Uses functional annotations to identify deleterious variants

## Key Output Files

### For Downstream Analysis
1. **Annotated VCF**: `snpeff_annotation/annotated_variants.vcf.gz`
   - Contains: SnpEff annotations, ancestral alleles, functional impacts

2. **High Impact Variants**: `snpeff_annotation/high_impact_variants.csv`
   - List of variants with HIGH functional impact

3. **LOF Variants**: `snpeff_annotation/lof_variants.csv`
   - Loss-of-function variants for genetic load calculation

4. **Gene Summary**: `functional_annotation/v2_gene_summary.tsv`
   - Comprehensive gene-level annotation statistics

## Performance Notes

- **Step 1 (Ancestral States)**: ~2-4 hours (depends on genome size and alignment)
- **Step 2 (SnpEff)**: ~1-2 hours (depends on variant count)
- **Step 3 (Functional Extraction)**: ~30 minutes
- **Step 4 (Visualization)**: ~10 minutes

**Total runtime**: ~4-7 hours for complete pipeline

## Citation

If using this pipeline, please cite:
- SnpEff: Cingolani et al. (2012) Fly 6:80-92
- minimap2: Li (2018) Bioinformatics 34:3094-3100

## Support

For issues or questions:
1. Check log files in `output/phase2_annotation/`
2. Review checkpoint status
3. Verify input file paths and formats
4. Check software dependencies and versions
