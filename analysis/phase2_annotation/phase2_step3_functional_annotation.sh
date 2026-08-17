#!/bin/bash
#
# Phase 3: Functional Annotation (v2) — refined extraction, classification, statistics
#
# Inputs:
#   - ${BASE_DIR}/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz
#
# Outputs (into functional_annotation directory, with v2_ prefix):
#   - v2_variant_counts_by_impact.tsv
#   - v2_variant_counts_by_effect.tsv
#   - v2_variant_counts_by_category.tsv
#   - v2_high_impact_variants.csv
#   - v2_lof_variants.csv
#   - v2_gene_summary.tsv
#   - v2_genes_<IMPACT>_{all,known,predicted}.txt
#   - v2_biotype_counts.tsv
#   - v2_gene_category_counts.tsv
#   - Optional figures: v2_impact_distribution.png, v2_top_effects.png, v2_category_distribution.png, v2_gene_counts_by_impact.png
#
# Usage:
#   bash analysis/phase2_annotation/phase2_step3_functional_annotation.sh
#

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
SNPEFF_VCF="${BASE_DIR}/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
OUT_DIR="${BASE_DIR}/output/phase2_annotation/functional_annotation"

mkdir -p "${OUT_DIR}"

LOGFILE="${OUT_DIR}/functional_annotation_v2.log"
exec > >(tee -a "$LOGFILE") 2>&1

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_message "====================================================================="
log_message "FUNCTIONAL ANNOTATION V2 — EXTRACTION / CLASSIFICATION / STATISTICS"
log_message "====================================================================="
log_message "Input VCF: ${SNPEFF_VCF}"
log_message "Output Dir: ${OUT_DIR}"

if [ ! -f "$SNPEFF_VCF" ]; then
    log_message "ERROR: SnpEff annotated VCF not found: $SNPEFF_VCF"
    exit 1
fi

python3 - << 'PYCODE'
import gzip
import sys
from collections import Counter, defaultdict
import csv
import os

vcf_file = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
out_dir = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/functional_annotation"

# SnpEff ANN indices (reference format):
#  0 Allele
#  1 Annotation (effect)
#  2 Annotation_Impact
#  3 Gene_Name
#  4 Gene_ID
#  5 Feature_Type
#  6 Feature_ID (e.g., transcript id)
#  7 Transcript_BioType
#  8 Rank/Total
#  9 HGVS.c
# 10 HGVS.p
# 11 cDNA.pos/cDNA.length
# 12 CDS.pos/CDS.length
# 13 AA.pos/AA.length
# 14 Distance
# 15 ERRORS/WARNINGS/INFO

IDX_EFFECT = 1
IDX_IMPACT = 2
IDX_GENE_NAME = 3
IDX_GENE_ID = 4
IDX_FEATURE_ID = 6
IDX_BIOTYPE = 7

impacts_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

# Mapping of SnpEff effects to broader categories
CODING_EFFECTS = {
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "start_gained",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "protein_altering_variant",
}
SPLICE_EFFECTS = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
}
UTR5_EFFECTS = {"5_prime_UTR_variant"}
UTR3_EFFECTS = {"3_prime_UTR_variant"}
INTRONIC_EFFECTS = {"intron_variant"}
UPSTREAM_EFFECTS = {"upstream_gene_variant"}
DOWNSTREAM_EFFECTS = {"downstream_gene_variant"}
INTERGENIC_EFFECTS = {"intergenic_region"}
REGULATORY_EFFECTS = {
    "regulatory_region_variant",
    "TF_binding_site_variant",
    "regulatory_region_amplification",
    "regulatory_region_ablation",
}

def effect_category(effect: str) -> str:
    if effect in CODING_EFFECTS:
        return "coding"
    if effect in SPLICE_EFFECTS:
        return "splice"
    if effect in UTR5_EFFECTS:
        return "utr5"
    if effect in UTR3_EFFECTS:
        return "utr3"
    if effect in INTRONIC_EFFECTS:
        return "intronic"
    if effect in UPSTREAM_EFFECTS:
        return "upstream"
    if effect in DOWNSTREAM_EFFECTS:
        return "downstream"
    if effect in INTERGENIC_EFFECTS:
        return "intergenic"
    if effect in REGULATORY_EFFECTS:
        return "regulatory"
    return "other"

def clean_gene_name(name: str) -> str:
    if not name or name == ".":
        return ""
    s = name.strip()
    for pref in ("GENE_", "gene-", "rna-", "CHR_START-", "CHR_END-"):
        if s.startswith(pref):
            s = s[len(pref):]
    for tok in ("-exon-", "-intron-"):
        if tok in s:
            s = s.split(tok, 1)[0]
    return s.strip()

def is_predicted_gene(name: str) -> bool:
    if not name:
        return False
    prefixes = ("LOC", "XLOC_", "MSTRG.", "TCONS_", "novel", "novelGene")
    return name.startswith(prefixes)

variant_impact_counts = Counter()
variant_effect_counts = Counter()
variant_category_counts = Counter()

lof_variants = []
high_impact_variants = []

# Gene-level aggregation
gene_impacts = defaultdict(Counter)  # gene -> impact counts
gene_effects = defaultdict(Counter)  # gene -> effect counts
gene_categories = defaultdict(Counter)  # gene -> category counts
gene_biotypes = defaultdict(set)
gene_variant_counts = Counter()  # number of variants that hit the gene
gene_info = {}  # gene -> dict(meta)

# Biotype-level aggregation
biotype_genes = defaultdict(set)
biotype_impact_genes = defaultdict(lambda: defaultdict(set))

def update_gene(gene_key: str, impact: str, effect: str, category: str, biotype: str):
    if not gene_key:
        return
    gene_impacts[gene_key][impact] += 1
    gene_effects[gene_key][effect] += 1
    gene_categories[gene_key][category] += 1
    if biotype:
        gene_biotypes[gene_key].add(biotype)
        biotype_genes[biotype].add(gene_key)
        biotype_impact_genes[biotype][impact].add(gene_key)

# Track per variant the set of genes affected (to count gene_variant_counts once per variant per gene)

opener = gzip.open if vcf_file.endswith('.gz') else open
with opener(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 8:
            continue
        chrom, pos, _vid, ref, alt, _qual, _filt, info = fields[:8]

        info_map = {}
        for item in info.split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                info_map[k] = v

        # Parse ANN
        ann_items = info_map.get('ANN')
        per_variant_genes = set()
        if ann_items:
            for ann in ann_items.split(','):
                parts = ann.split('|')
                if len(parts) <= IDX_IMPACT:
                    continue
                effect = parts[IDX_EFFECT] if len(parts) > IDX_EFFECT else ''
                impact = parts[IDX_IMPACT] if len(parts) > IDX_IMPACT else ''
                gene_raw = parts[IDX_GENE_NAME] if len(parts) > IDX_GENE_NAME else ''
                gene_id = parts[IDX_GENE_ID] if len(parts) > IDX_GENE_ID else ''
                biotype = parts[IDX_BIOTYPE] if len(parts) > IDX_BIOTYPE else ''
                tx_id = parts[IDX_FEATURE_ID] if len(parts) > IDX_FEATURE_ID else ''

                effect = effect or ''
                impact = impact or 'MODIFIER'

                category = effect_category(effect)
                variant_impact_counts[impact] += 1
                variant_effect_counts[effect] += 1
                variant_category_counts[category] += 1

                gene_clean = clean_gene_name(gene_raw)
                gene_key = gene_clean or gene_raw or ''
                # Update gene-level structures
                update_gene(gene_key, impact, effect, category, biotype)
                if gene_key:
                    per_variant_genes.add(gene_key)
                    if gene_key not in gene_info:
                        gene_info[gene_key] = {
                            'gene_id': gene_id,
                            'predicted': is_predicted_gene(gene_clean),
                        }

                if impact == 'HIGH':
                    high_impact_variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gene': gene_key or '.',
                        'effect': effect,
                        'biotype': biotype,
                        'transcript': tx_id,
                    })

        # Increment gene_variant_counts once per gene per variant
        for g in per_variant_genes:
            gene_variant_counts[g] += 1

        # LOF parsing from INFO: LOF=<entries>
        # SnpEff LOF format example: LOF=(GENE|GENEID|TR|EXON#|PERC|...)
        lof_raw = info_map.get('LOF')
        if lof_raw:
            # Normalize: remove parentheses if any
            lof_items = []
            cur = ''
            paren = 0
            for ch in lof_raw:
                if ch == '(':
                    paren += 1
                elif ch == ')':
                    paren -= 1
                if ch == ',' and paren == 0:
                    lof_items.append(cur)
                    cur = ''
                else:
                    cur += ch
            if cur:
                lof_items.append(cur)

            for entry in lof_items:
                e = entry.strip()
                if e.startswith('(') and e.endswith(')'):
                    e = e[1:-1]
                parts = e.split('|')
                gene = parts[0] if len(parts) > 0 else ''
                gene_clean = clean_gene_name(gene)
                lof_variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gene': gene_clean or gene or '.',
                    'lof': e,
                })

# Write variant-level summaries
with open(os.path.join(out_dir, 'v2_variant_counts_by_impact.tsv'), 'w') as f:
    f.write("Impact\tCount\n")
    total = sum(variant_impact_counts.values()) or 1
    for imp in impacts_order:
        f.write(f"{imp}\t{variant_impact_counts.get(imp, 0)}\n")

with open(os.path.join(out_dir, 'v2_variant_counts_by_effect.tsv'), 'w') as f:
    f.write("Effect\tCount\n")
    for effect, count in variant_effect_counts.most_common():
        f.write(f"{effect}\t{count}\n")

with open(os.path.join(out_dir, 'v2_variant_counts_by_category.tsv'), 'w') as f:
    f.write("Category\tCount\n")
    for cat, count in variant_category_counts.most_common():
        f.write(f"{cat}\t{count}\n")

# Write variants of interest
if high_impact_variants:
    with open(os.path.join(out_dir, 'v2_high_impact_variants.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['chrom','pos','ref','alt','gene','effect','biotype','transcript'])
        w.writeheader()
        w.writerows(high_impact_variants)

if lof_variants:
    with open(os.path.join(out_dir, 'v2_lof_variants.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['chrom','pos','ref','alt','gene','lof'])
        w.writeheader()
        w.writerows(lof_variants)

# Gene lists by impact and category (known/predicted/all)
def highest_impact(imp_counter: Counter) -> str:
    for imp in impacts_order:
        if imp_counter.get(imp, 0) > 0:
            return imp
    return 'MODIFIER'

genes_by_impact_all = defaultdict(set)
genes_by_impact_known = defaultdict(set)
genes_by_impact_pred = defaultdict(set)

for gene, imp_counter in gene_impacts.items():
    hi = highest_impact(imp_counter)
    genes_by_impact_all[hi].add(gene)
    pred = gene_info.get(gene, {}).get('predicted', False)
    if pred:
        genes_by_impact_pred[hi].add(gene)
    else:
        genes_by_impact_known[hi].add(gene)

for imp in impacts_order:
    if genes_by_impact_all[imp]:
        with open(os.path.join(out_dir, f'v2_genes_{imp.lower()}_all.txt'), 'w') as f:
            for g in sorted(genes_by_impact_all[imp]):
                f.write(g + '\n')
    if genes_by_impact_known[imp]:
        with open(os.path.join(out_dir, f'v2_genes_{imp.lower()}_known.txt'), 'w') as f:
            for g in sorted(genes_by_impact_known[imp]):
                f.write(g + '\n')
    if genes_by_impact_pred[imp]:
        with open(os.path.join(out_dir, f'v2_genes_{imp.lower()}_predicted.txt'), 'w') as f:
            for g in sorted(genes_by_impact_pred[imp]):
                f.write(g + '\n')

# Biotype counts summary
with open(os.path.join(out_dir, 'v2_biotype_counts.tsv'), 'w') as f:
    f.write('Biotype\tTotal\tHIGH\tMODERATE\tLOW\tMODIFIER\n')
    for bt in sorted(biotype_genes.keys()):
        total = len(biotype_genes[bt])
        row = [
            bt,
            str(total),
            str(len(biotype_impact_genes[bt].get('HIGH', set()))),
            str(len(biotype_impact_genes[bt].get('MODERATE', set()))),
            str(len(biotype_impact_genes[bt].get('LOW', set()))),
            str(len(biotype_impact_genes[bt].get('MODIFIER', set()))),
        ]
        f.write('\t'.join(row) + '\n')

# Gene category counts summary (per gene, aggregate categories across its effects)
category_gene_sets = defaultdict(set)
for gene, cat_counter in gene_categories.items():
    for cat in cat_counter.keys():
        category_gene_sets[cat].add(gene)

with open(os.path.join(out_dir, 'v2_gene_category_counts.tsv'), 'w') as f:
    f.write('Category\tGenes\n')
    for cat, gset in sorted(category_gene_sets.items(), key=lambda x: (-len(x[1]), x[0])):
        f.write(f"{cat}\t{len(gset)}\n")

# Comprehensive gene summary table
with open(os.path.join(out_dir, 'v2_gene_summary.tsv'), 'w') as f:
    f.write('Gene\tGeneID\tPredicted\tVariants\tHighestImpact\tImpacts\tTopEffects\tCategories\tBiotypes\n')
    for gene in sorted(gene_impacts.keys()):
        gi = gene_info.get(gene, {})
        imp_counts = gene_impacts[gene]
        eff_counts = gene_effects[gene]
        cat_counts = gene_categories[gene]
        hi = highest_impact(imp_counts)
        effects_sorted = ','.join([f"{e}:{c}" for e, c in eff_counts.most_common(10)])
        cats_sorted = ','.join([f"{c}:{n}" for c, n in sorted(cat_counts.items(), key=lambda x: (-x[1], x[0]))])
        imps_sorted = ','.join([f"{imp}:{imp_counts.get(imp,0)}" for imp in impacts_order])
        biotypes_sorted = ','.join(sorted(gene_biotypes.get(gene, [])))
        f.write('\t'.join([
            gene,
            gi.get('gene_id', '') or '.',
            'Yes' if gi.get('predicted', False) else 'No',
            str(gene_variant_counts.get(gene, 0)),
            hi,
            imps_sorted,
            effects_sorted,
            cats_sorted,
            biotypes_sorted or '.',
        ]) + '\n')

# Optional: basic figures if matplotlib is available
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Match phase2_step4.2_visualize_functional.py style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)

    # Impact distribution
    counts = [variant_impact_counts.get(i, 0) for i in impacts_order]
    plt.figure(figsize=(6,4))
    bars = plt.bar(impacts_order, counts,
                   color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"],
                   alpha=0.85, edgecolor='black', linewidth=1.2)
    plt.title('Variant Impact Distribution', fontsize=13, fontweight='bold')
    plt.ylabel('Count', fontsize=12, fontweight='bold')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'v2_impact_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Top effects
    top_effects = variant_effect_counts.most_common(20)
    if top_effects:
        labels = [e for e,_ in top_effects]
        values = [c for _,c in top_effects]
        plt.figure(figsize=(8,6))
        plt.barh(labels[::-1], values[::-1], color="#4e79a7", alpha=0.9, edgecolor='black', linewidth=1.0)
        plt.xlabel('Count', fontsize=12, fontweight='bold')
        plt.title('Top Consequences', fontsize=13, fontweight='bold')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'v2_top_effects.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Category distribution
    cats, cat_vals = zip(*sorted(variant_category_counts.items(), key=lambda x: -x[1])) if variant_category_counts else ([], [])
    if cats:
        plt.figure(figsize=(8,5))
        plt.barh(list(cats)[::-1], list(cat_vals)[::-1], color="#59a14f", alpha=0.9, edgecolor='black', linewidth=1.0)
        plt.xlabel('Count', fontsize=12, fontweight='bold')
        plt.title('Variant Category Distribution', fontsize=13, fontweight='bold')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'v2_category_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Gene counts by highest impact
    g_counts = [len(genes_by_impact_all[i]) for i in impacts_order]
    plt.figure(figsize=(6,4))
    bars = plt.bar(impacts_order, g_counts,
                   color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"],
                   alpha=0.85, edgecolor='black', linewidth=1.2)
    plt.title('Genes by Highest Impact', fontsize=13, fontweight='bold')
    plt.ylabel('Unique genes', fontsize=12, fontweight='bold')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'v2_gene_counts_by_impact.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Biotype vs impact (top 15 biotypes by total genes) — stacked horizontal bars
    if biotype_genes:
        top_biotypes = sorted(biotype_genes.keys(), key=lambda bt: len(biotype_genes[bt]), reverse=True)[:15]
        data = {imp: [len(biotype_impact_genes[bt].get(imp, set())) for bt in top_biotypes] for imp in impacts_order}
        y = list(range(len(top_biotypes)))
        left = [0]*len(top_biotypes)
        colors = ["#d62728", "#ff7f0e", "#1f77b4", "#2ca02c"]
        plt.figure(figsize=(9, max(4, 0.4*len(top_biotypes)+2)))
        for idx, imp in enumerate(impacts_order):
            vals = data[imp]
            plt.barh(y, vals, left=left, color=colors[idx], label=imp, alpha=0.9, edgecolor='black', linewidth=0.8)
            left = [l+v for l, v in zip(left, vals)]
        plt.yticks(y, top_biotypes)
        plt.xlabel('Gene count', fontsize=12, fontweight='bold')
        plt.title('Biotype vs Impact (top biotypes)', fontsize=13, fontweight='bold')
        plt.legend(loc='lower right', frameon=True, fancybox=True, shadow=True)
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'v2_biotype_vs_impact_genes.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Top genes by number of variants
    if gene_variant_counts:
        top_genes = gene_variant_counts.most_common(25)
        labels = [g for g, _ in top_genes]
        values = [c for _, c in top_genes]
        plt.figure(figsize=(8, max(4, 0.35*len(labels)+2)))
        plt.barh(labels[::-1], values[::-1], color="#6b6ecf", alpha=0.9, edgecolor='black', linewidth=1.0)
        plt.xlabel('Variants affecting gene', fontsize=12, fontweight='bold')
        plt.title('Top genes by variant count', fontsize=13, fontweight='bold')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'v2_top_genes_by_variant_count.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # Top genes by HIGH-impact annotation counts
    high_counts = [(g, gene_impacts[g].get('HIGH', 0)) for g in gene_impacts]
    high_counts = [(g, c) for g, c in high_counts if c > 0]
    if high_counts:
        high_counts.sort(key=lambda x: x[1], reverse=True)
        top_high = high_counts[:25]
        labels = [g for g, _ in top_high]
        values = [c for _, c in top_high]
        plt.figure(figsize=(8, max(4, 0.35*len(labels)+2)))
        plt.barh(labels[::-1], values[::-1], color="#d62728", alpha=0.9, edgecolor='black', linewidth=1.0)
        plt.xlabel('High-impact annotations', fontsize=12, fontweight='bold')
        plt.title('Top genes by High impact', fontsize=13, fontweight='bold')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'v2_top_genes_by_high_impact.png'), dpi=300, bbox_inches='tight')
        plt.close()
except Exception as e:
    print(f"WARNING: Plotting skipped ({e})", file=sys.stderr)

print("\n" + "="*70)
print("FUNCTIONAL ANNOTATION V2 SUMMARY")
print("="*70)
total_var = sum(variant_impact_counts.values())
print(f"Total annotated records: {total_var}")
print("Impact distribution:")
for imp in impacts_order:
    c = variant_impact_counts.get(imp, 0)
    pct = 100.0 * c / total_var if total_var else 0.0
    print(f"  {imp}: {c} ({pct:.1f}%)")
print("Top 10 effects:")
for eff, cnt in variant_effect_counts.most_common(10):
    print(f"  {eff}: {cnt}")
print(f"Genes summarized: {len(gene_impacts)}")
print("Outputs written with prefix 'v2_' in:")
print(out_dir)
print("="*70)
PYCODE

log_message "Functional annotation v2 complete"

exit 0


