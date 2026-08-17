#!/usr/bin/env python3
"""
Update Parameter Bounds for fastsimcoal2 Models
================================================

This script creates updated .est files with more appropriate parameter bounds
based on biological knowledge and expected ranges.

Usage:
    python3 update_parameter_bounds.py [--conservative|--expanded|--custom]

Options:
    --conservative  : Keep current bounds (default)
    --expanded      : Use wider bounds (recommended)
    --custom        : Interactive mode to set custom bounds

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import shutil

# Configuration
BASE_DIR = get_base_dir()
MODEL_DIR = Path(__file__).resolve().parent / "models"


def create_expanded_bounds():
    """Return expanded parameter bounds."""
    return {
        'NCUR': {'min': 50, 'max': 200000, 'note': 'Current Ne - very wide range'},
        'NBOT': {'min': 20, 'max': 100000, 'note': 'Bottleneck Ne - allow severe bottlenecks'},
        'NANC': {'min': 500, 'max': 500000, 'note': 'Ancestral Ne - allow large ancient populations'},
        'NINTER': {'min': 100, 'max': 200000, 'note': 'Intermediate Ne between two bottlenecks'},
        'NRECOVER': {'min': 100, 'max': 200000, 'note': 'Short recovery Ne before recent contraction'},
        'N1': {'min': 100, 'max': 200000, 'note': 'Complex model intermediate Ne 1'},
        'N2': {'min': 100, 'max': 200000, 'note': 'Complex model intermediate Ne 2'},
        'N3': {'min': 100, 'max': 200000, 'note': 'Complex model intermediate Ne 3'},
        'RMID': {'min': 1.01, 'max': 10.0, 'note': 'NMID/NCUR ratio for decline model'},
        'RBOT': {'min': 1.01, 'max': 10.0, 'note': 'NBOT/NMID ratio for decline model'},
        'RANC': {'min': 1.01, 'max': 10.0, 'note': 'NANC/NBOT ratio for decline model'},
        'DTBOT': {'min': 1, 'max': 20000, 'note': 'Positive time gap to enforce TBOT > recent time'},
        'DTANC': {'min': 1, 'max': 30000, 'note': 'Positive time gap to enforce TANC > TBOT'},
        'DTRECENT_BOT': {'min': 1, 'max': 10000, 'note': 'Positive gap: TRECENT_BOT > TRECENT_RECOVERY'},
        'DTOLD_RECOVERY': {'min': 1, 'max': 20000, 'note': 'Positive gap: TOLD_RECOVERY > TRECENT_BOT'},
        'DTOLD_BOT': {'min': 1, 'max': 40000, 'note': 'Positive gap: TOLD_BOT > TOLD_RECOVERY'},
        'DTRECOVERY_OLD': {'min': 1, 'max': 10000, 'note': 'Positive gap: TRECOVERY_OLD > TRECENT'},
        'DTBOT_OLD': {'min': 1, 'max': 30000, 'note': 'Positive gap: TBOT_OLD > TRECOVERY_OLD'},
        'DT2': {'min': 1, 'max': 5000, 'note': 'Positive gap: T2 > T1'},
        'DT3': {'min': 1, 'max': 15000, 'note': 'Positive gap: T3 > T2'},
        'DT4': {'min': 1, 'max': 50000, 'note': 'Positive gap: T4 > T3'},
        'TBOT': {'min': 5, 'max': 10000, 'note': 'Bottleneck time - wider temporal range'},
        'TRECOVER': {'min': 50, 'max': 50000, 'note': 'Recovery time - wider range'},
        'TRECOVERY': {'min': 5, 'max': 10000, 'note': 'Single-bottleneck recovery time'},
        'TRECENT': {'min': 20, 'max': 2000, 'note': 'Recent contraction/decline time'},
        'TRECOVERY_OLD': {'min': 100, 'max': 10000, 'note': 'Ancient recovery after bottleneck'},
        'TBOT_OLD': {'min': 200, 'max': 30000, 'note': 'Ancient bottleneck onset'},
        'TRECENT_RECOVERY': {'min': 5, 'max': 1500, 'note': 'Recent bottleneck recovery time'},
        'TRECENT_BOT': {'min': 20, 'max': 6000, 'note': 'Recent bottleneck onset'},
        'TOLD_RECOVERY': {'min': 50, 'max': 15000, 'note': 'Old bottleneck recovery time'},
        'TOLD_BOT': {'min': 100, 'max': 50000, 'note': 'Old bottleneck onset'},
        'T1': {'min': 5, 'max': 1000, 'note': 'Complex model time 1'},
        'T2': {'min': 20, 'max': 5000, 'note': 'Complex model time 2'},
        'T3': {'min': 100, 'max': 15000, 'note': 'Complex model time 3'},
        'T4': {'min': 300, 'max': 50000, 'note': 'Complex model time 4'},
        'TCHANGE': {'min': 5, 'max': 20000, 'note': 'Size change time'},
        'TGROWTH': {'min': 5, 'max': 10000, 'note': 'Growth start time'},
        'TANC': {'min': 500, 'max': 100000, 'note': 'Ancient time'},
        'GROWTH': {'min': -0.02, 'max': 0.005, 'note': 'Growth rate (negative = expansion)'}
    }


def create_conservative_bounds():
    """Return conservative (current) parameter bounds."""
    return {
        'NCUR': {'min': 500, 'max': 100000, 'note': 'Current Ne'},
        'NBOT': {'min': 100, 'max': 50000, 'note': 'Bottleneck Ne'},
        'NANC': {'min': 500, 'max': 100000, 'note': 'Ancestral Ne'},
        'NINTER': {'min': 100, 'max': 100000, 'note': 'Intermediate Ne between bottlenecks'},
        'NRECOVER': {'min': 300, 'max': 100000, 'note': 'Recovery Ne before recent contraction'},
        'N1': {'min': 100, 'max': 100000, 'note': 'Complex model Ne 1'},
        'N2': {'min': 100, 'max': 100000, 'note': 'Complex model Ne 2'},
        'N3': {'min': 100, 'max': 100000, 'note': 'Complex model Ne 3'},
        'RMID': {'min': 1.01, 'max': 5.0, 'note': 'NMID/NCUR ratio for decline model'},
        'RBOT': {'min': 1.01, 'max': 5.0, 'note': 'NBOT/NMID ratio for decline model'},
        'RANC': {'min': 1.01, 'max': 5.0, 'note': 'NANC/NBOT ratio for decline model'},
        'DTBOT': {'min': 1, 'max': 8000, 'note': 'Positive time gap to enforce TBOT > recent time'},
        'DTANC': {'min': 1, 'max': 20000, 'note': 'Positive time gap to enforce TANC > TBOT'},
        'DTRECENT_BOT': {'min': 1, 'max': 3000, 'note': 'Positive gap: TRECENT_BOT > TRECENT_RECOVERY'},
        'DTOLD_RECOVERY': {'min': 1, 'max': 8000, 'note': 'Positive gap: TOLD_RECOVERY > TRECENT_BOT'},
        'DTOLD_BOT': {'min': 1, 'max': 20000, 'note': 'Positive gap: TOLD_BOT > TOLD_RECOVERY'},
        'DTRECOVERY_OLD': {'min': 1, 'max': 3800, 'note': 'Positive gap: TRECOVERY_OLD > TRECENT'},
        'DTBOT_OLD': {'min': 1, 'max': 16000, 'note': 'Positive gap: TBOT_OLD > TRECOVERY_OLD'},
        'DT2': {'min': 1, 'max': 2500, 'note': 'Positive gap: T2 > T1'},
        'DT3': {'min': 1, 'max': 7000, 'note': 'Positive gap: T3 > T2'},
        'DT4': {'min': 1, 'max': 20000, 'note': 'Positive gap: T4 > T3'},
        'TBOT': {'min': 10, 'max': 5000, 'note': 'Bottleneck time'},
        'TRECOVER': {'min': 100, 'max': 20000, 'note': 'Recovery time'},
        'TRECOVERY': {'min': 10, 'max': 3000, 'note': 'Single-bottleneck recovery time'},
        'TRECENT': {'min': 20, 'max': 1000, 'note': 'Recent contraction/decline time'},
        'TRECOVERY_OLD': {'min': 200, 'max': 4000, 'note': 'Ancient recovery time'},
        'TBOT_OLD': {'min': 400, 'max': 20000, 'note': 'Ancient bottleneck onset time'},
        'TRECENT_RECOVERY': {'min': 10, 'max': 1000, 'note': 'Recent bottleneck recovery time'},
        'TRECENT_BOT': {'min': 30, 'max': 4000, 'note': 'Recent bottleneck onset time'},
        'TOLD_RECOVERY': {'min': 100, 'max': 10000, 'note': 'Old bottleneck recovery time'},
        'TOLD_BOT': {'min': 300, 'max': 30000, 'note': 'Old bottleneck onset time'},
        'T1': {'min': 10, 'max': 500, 'note': 'Complex model time 1'},
        'T2': {'min': 50, 'max': 3000, 'note': 'Complex model time 2'},
        'T3': {'min': 200, 'max': 10000, 'note': 'Complex model time 3'},
        'T4': {'min': 500, 'max': 30000, 'note': 'Complex model time 4'},
        'TCHANGE': {'min': 10, 'max': 10000, 'note': 'Size change time'},
        'TGROWTH': {'min': 10, 'max': 5000, 'note': 'Growth start time'},
        'TANC': {'min': 500, 'max': 10000, 'note': 'Ancient time'},
        'GROWTH': {'min': -0.01, 'max': 0.001, 'note': 'Growth rate'}
    }


def update_est_file(model_name, bounds_dict):
    """Update .est file with new bounds."""
    est_file = MODEL_DIR / f"{model_name}.est"
    
    if not est_file.exists():
        print(f"  ⚠️  File not found: {est_file}")
        return False
    
    # Backup original
    backup_file = est_file.with_suffix('.est.original')
    if not backup_file.exists():
        shutil.copy(est_file, backup_file)
        print(f"  📄 Backup created: {backup_file.name}")
    
    # Read file
    with open(est_file, 'r') as f:
        lines = f.readlines()
    
    # Update lines
    updated_lines = []
    params_updated = []
    
    for line in lines:
        if line.strip() and not line.startswith('//') and not line.startswith('['):
            # Check if this is a parameter line
            parts = line.split()
            if len(parts) >= 5:
                param_name = parts[1]
                if param_name in bounds_dict:
                    # Update bounds
                    parts[3] = str(bounds_dict[param_name]['min'])
                    parts[4] = str(bounds_dict[param_name]['max'])
                    line = '  '.join(parts) + '\n'
                    params_updated.append(param_name)
        
        updated_lines.append(line)
    
    # Write updated file
    with open(est_file, 'w') as f:
        f.writelines(updated_lines)
    
    if params_updated:
        print(f"  ✅ Updated {len(params_updated)} parameters: {', '.join(params_updated)}")
    else:
        print(f"  ℹ️  No parameters updated")
    
    return True


def update_all_models(bounds_type='expanded'):
    """Update all model .est files."""
    print("=" * 80)
    print("UPDATING PARAMETER BOUNDS FOR FASTSIMCOAL2 MODELS")
    print("=" * 80)
    print()
    
    if bounds_type == 'expanded':
        bounds = create_expanded_bounds()
        print("Using EXPANDED bounds (recommended)")
    else:
        bounds = create_conservative_bounds()
        print("Using CONSERVATIVE bounds (original)")
    
    print()
    print("Bounds to apply:")
    print("-" * 80)
    for param, info in bounds.items():
        print(f"  {param:<12}: {info['min']:>8} - {info['max']:<10}  ({info['note']})")
    print()
    
    models = [
        'constant_ne',
        'single_bottleneck',
        'two_consecutive_bottlenecks',
        'bottleneck_continuous_decline',
        'bottleneck_recent_contraction',
        'complex_multi_event'
    ]
    
    print("Updating models:")
    print("-" * 80)
    
    for model in models:
        print(f"\n{model}:")
        update_est_file(model, bounds)
    
    print()
    print("=" * 80)
    print("UPDATE COMPLETE")
    print("=" * 80)
    print()
    print("Original files backed up with .est.original extension")
    print()
    print("To revert:")
    print("  for f in models/*.est.original; do")
    print("    mv \"$f\" \"${f%.original}\"")
    print("  done")
    print()


def show_comparison():
    """Show comparison of conservative vs expanded bounds."""
    print("=" * 80)
    print("PARAMETER BOUNDS COMPARISON")
    print("=" * 80)
    print()
    
    conservative = create_conservative_bounds()
    expanded = create_expanded_bounds()
    
    print(f"{'Parameter':<12} {'Conservative':<25} {'Expanded':<25} {'Change':<20}")
    print("-" * 80)
    
    for param in sorted(conservative.keys()):
        cons_range = f"{conservative[param]['min']} - {conservative[param]['max']}"
        exp_range = f"{expanded[param]['min']} - {expanded[param]['max']}"
        
        # Calculate change
        cons_width = conservative[param]['max'] - conservative[param]['min']
        exp_width = expanded[param]['max'] - expanded[param]['min']
        change = f"{exp_width/cons_width:.1f}× wider" if cons_width > 0 else "N/A"
        
        print(f"{param:<12} {cons_range:<25} {exp_range:<25} {change:<20}")
    
    print()


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Update parameter bounds for fastsimcoal2 models')
    parser.add_argument('--type', choices=['conservative', 'expanded', 'compare'],
                       default='compare', help='Bounds type')
    parser.add_argument('--apply', action='store_true', help='Apply the changes')
    
    args = parser.parse_args()
    
    if args.type == 'compare':
        show_comparison()
        print()
        print("To apply expanded bounds:")
        print("  python3 update_parameter_bounds.py --type expanded --apply")
        print()
        
    else:
        if args.apply:
            update_all_models(args.type)
        else:
            print("Dry run mode. Add --apply to make changes.")
            print()
            if args.type == 'expanded':
                bounds = create_expanded_bounds()
            else:
                bounds = create_conservative_bounds()
            
            print("Bounds that would be applied:")
            for param, info in bounds.items():
                print(f"  {param}: {info['min']} - {info['max']}")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        # Default: show comparison
        show_comparison()
        print()
        print("To apply expanded bounds (RECOMMENDED):")
        print("  python3 update_parameter_bounds.py --type expanded --apply")
        print()
        print("To keep conservative bounds:")
        print("  (No action needed - current bounds are conservative)")
        print()
    else:
        main()
