#!/usr/bin/env python3
"""Plot validation results for visual verification."""

import argparse
import csv
from pathlib import Path

# Use non-interactive backend for headless/WSL environments
import matplotlib
matplotlib.use('Agg')  # Must be before importing pyplot
import matplotlib.pyplot as plt
import numpy as np

# Plot style (with fallback for compatibility)
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    try:
        plt.style.use('seaborn-darkgrid')
    except:
        pass

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

def read_csv(filepath):
    """Read CSV file and return data as dict of arrays."""
    data = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, val in row.items():
                if key not in data:
                    data[key] = []
                try:
                    data[key].append(float(val))
                except ValueError:
                    data[key].append(val)
    return {k: np.array(v) for k, v in data.items()}

def plot_ortho(test_root, output_dir, compare=False):
    """Plot all 6 orthotropic tests on ONE panel."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    tests = [
        ('x_tension', 'X-Tension', 'strain_xx', 'stress_xx', COLORS[0]),
        ('y_tension', 'Y-Tension', 'strain_yy', 'stress_yy', COLORS[1]),
        ('z_tension', 'Z-Tension', 'strain_zz', 'stress_zz', COLORS[2]),
        ('xy_shear', 'XY-Shear', 'strain_xy', 'stress_xy', COLORS[3]),
        ('xz_shear', 'XZ-Shear', 'strain_xz', 'stress_xz', COLORS[4]),
        ('yz_shear', 'YZ-Shear', 'strain_yz', 'stress_yz', COLORS[5]),
    ]
    
    for test_name, label, strain_col, stress_col, color in tests:
        csv_file = test_root / 'ortho' / test_name / f'{test_name}_out.csv'
        if csv_file.exists():
            data = read_csv(csv_file)
            if strain_col in data and stress_col in data:
                ax.plot(data[strain_col] * 100, data[stress_col] / 1000,
                       color=color, linewidth=2, label=label)
    
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.set_ylabel('Stress (GPa)', fontsize=12)
    ax.set_title('Orthotropic Validation (All Directions)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'ortho.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  → ortho.png")

def plot_asymmetry(test_root, output_dir, compare=False):
    """Plot asymmetry tests."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # X-compression
    ax = axes[0]
    csv_file = test_root / 'asymmetry' / 'x_compression' / 'x_compression_out.csv'
    if csv_file.exists():
        data = read_csv(csv_file)
        if 'strain_xx' in data and 'stress_xx' in data:
            ax.plot(data['strain_xx'] * 100, data['stress_xx'] / 1000,
                   'b-', linewidth=2, label='Compression')
            
            if compare:
                gold_file = test_root / 'asymmetry' / 'x_compression' / 'gold' / 'x_compression_out.csv'
                if gold_file.exists():
                    gold = read_csv(gold_file)
                    if 'strain_xx' in gold and 'stress_xx' in gold:
                        ax.plot(gold['strain_xx'] * 100, gold['stress_xx'] / 1000,
                               'r--', linewidth=1.5, alpha=0.7, label='Gold')
    
    ax.axhline(y=-3.9, color='gray', linestyle='--', alpha=0.5, label='σ_c = -3900 MPa')
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.set_ylabel('Stress (GPa)', fontsize=12)
    ax.set_title('X-Direction Compression', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Yield surface shift
    ax = axes[1]
    csv_file = test_root / 'asymmetry' / 'yield_surface_shift' / 'yield_surface_shift_out.csv'
    if csv_file.exists():
        data = read_csv(csv_file)
        if 'strain_xx' in data and 'stress_xx' in data and 'stress_yy' in data:
            ax.plot(data['strain_xx'] * 100, data['stress_xx'] / 1000,
                   'b-', linewidth=2, label='σ_xx')
            ax.plot(data['strain_xx'] * 100, data['stress_yy'] / 1000,
                   'g-', linewidth=2, label='σ_yy')
    
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.set_ylabel('Stress (GPa)', fontsize=12)
    ax.set_title('Biaxial Loading', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'asymmetry.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  → asymmetry.png")

def plot_postyield(test_root, output_dir, compare=False):
    """Plot all 5 post-yield modes on ONE panel."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    tests = [
        ('perfect_plasticity', 'Perfect Plasticity', COLORS[0]),
        ('exp_hardening', 'Exponential Hardening', COLORS[1]),
        ('linear_hardening', 'Linear Hardening', COLORS[2]),
        ('exp_softening', 'Exponential Softening', COLORS[3]),
        ('piecewise_softening', 'Piecewise Softening', COLORS[4]),
    ]
    
    for test_name, label, color in tests:
        csv_file = test_root / 'postyield' / test_name / f'{test_name}_out.csv'
        if csv_file.exists():
            data = read_csv(csv_file)
            if 'strain_xx' in data and 'stress_xx' in data:
                ax.plot(data['strain_xx'] * 100, data['stress_xx'] / 1000,
                       color=color, linewidth=2, label=label)
    
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.set_ylabel('Stress (GPa)', fontsize=12)
    ax.set_title('Post-Yield Evolution (All Modes)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'postyield.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  → postyield.png")

def extract_angle_from_name(test_name):
    """Extract misorientation angle from test directory name."""
    import re
    
    # Check for specific patterns
    if 'twin' in test_name.lower():
        return 'Twin', None
    
    # Try to find angle in name (e.g., "bicrystal_45deg", "45_deg", "bicrystal_90")
    match = re.search(r'(\d+\.?\d*)\s*deg', test_name, re.IGNORECASE)
    if match:
        angle = float(match.group(1))
        return f'{angle:.1f}°' if angle != int(angle) else f'{int(angle)}°', angle
    
    # Check for number followed by underscore or end (e.g., "bicrystal_90")
    match = re.search(r'_(\d+)(?:_|$)', test_name)
    if match:
        angle = int(match.group(1))
        return f'{angle}°', angle
    
    # Fallback
    return test_name.replace('_', ' ').title(), None

def plot_grains(test_root, output_dir, compare=False):
    """Plot all grain tests - one color per angle, solid/dashed for grain 1/2."""
    grains_dir = test_root / 'grains'
    if not grains_dir.exists():
        return
    
    test_dirs = sorted([d for d in grains_dir.iterdir() if d.is_dir()])
    if len(test_dirs) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(11, 8))
    
    # Assign colors based on unique angles
    angle_color_map = {}
    color_idx = 0
    
    for test_dir in test_dirs:
        test_name = test_dir.name
        label_base, angle = extract_angle_from_name(test_name)
        
        # Assign color
        if label_base not in angle_color_map:
            angle_color_map[label_base] = COLORS[color_idx % len(COLORS)]
            color_idx += 1
        
        color = angle_color_map[label_base]
        
        csv_file = test_dir / f'{test_name}_out.csv'
        if not csv_file.exists():
            continue
        
        data = read_csv(csv_file)
        
        # Plot grain-specific curves
        if 'strain_xx_grain1' in data and 'stress_xx_grain1' in data:
            ax.plot(data['strain_xx_grain1'] * 100, data['stress_xx_grain1'] / 1000,
                   color=color, linewidth=2.5, linestyle='-', 
                   label=f'{label_base} Grain 1')
        
        if 'strain_xx_grain2' in data and 'stress_xx_grain2' in data:
            ax.plot(data['strain_xx_grain2'] * 100, data['stress_xx_grain2'] / 1000,
                   color=color, linewidth=2.5, linestyle='--',
                   label=f'{label_base} Grain 2')
    
    # Add annotation about strain redistribution
    ax.text(0.98, 0.02, 
            'Note: Strain decrease in one grain indicates unloading due to load redistribution\n'
            '(fixed grain (Grain 2) unloads as pulled grain (Grain 2) continue to soften)',
            transform=ax.transAxes, fontsize=9, 
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.set_ylabel('Stress (GPa)', fontsize=12)
    ax.set_title('Multi-Grain Validation (Orientation Effects)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=9, ncol=2)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'grains.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  → grains.png")

def plot_czm(test_root, output_dir, compare=False):
    """Plot all 4 CZM tests on ONE panel with correct column names."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Correct column names based on actual CSV headers
    tests = [
        ('mode_i', 'Mode I (Opening)', 'avg_delta_eff', 'avg_traction', COLORS[0]),
        ('mode_ii', 'Mode II (Shear)', 'avg_delta_eff', 'avg_tangent_traction', COLORS[1]),  # FIXED!
        ('mixed_mode', 'Mixed Mode', 'avg_delta_eff', 'T_mag', COLORS[2]),
        ('complete_separation', 'Complete Separation', 'avg_delta_eff', 'avg_traction', COLORS[3]),
    ]
    
    for test_name, label, sep_col, trac_col, color in tests:
        csv_file = test_root / 'czm' / test_name / f'{test_name}_out.csv'
        if csv_file.exists():
            data = read_csv(csv_file)
            if sep_col in data and trac_col in data:
                ax.plot(data[sep_col], data[trac_col] / 1000,
                       color=color, linewidth=2, label=label)
            else:
                print(f"  Warning: {test_name} missing columns: {sep_col} or {trac_col}")
    
    ax.set_xlabel('Effective Separation (nm)', fontsize=12)
    ax.set_ylabel('Traction (GPa)', fontsize=12)
    ax.set_title('Cohesive Zone Model Validation', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'czm.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  → czm.png")

def main():
    parser = argparse.ArgumentParser(description='Plot validation results')
    parser.add_argument('--category', required=True)
    parser.add_argument('--test-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--compare', action='store_true')
    args = parser.parse_args()
    
    # Call appropriate plotting function
    if args.category == 'ortho':
        plot_ortho(args.test_root, args.output_dir, args.compare)
    elif args.category == 'asymmetry':
        plot_asymmetry(args.test_root, args.output_dir, args.compare)
    elif args.category == 'postyield':
        plot_postyield(args.test_root, args.output_dir, args.compare)
    elif args.category == 'grains':
        plot_grains(args.test_root, args.output_dir, args.compare)
    elif args.category == 'czm':
        plot_czm(args.test_root, args.output_dir, args.compare)
    else:
        print(f"Unknown category: {args.category}")

if __name__ == '__main__':
    main()