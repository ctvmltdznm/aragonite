#!/usr/bin/env python3
"""
Extract per-grain stress and strain statistics from MOOSE exodus output
Reads needle_id element variable and computes averages per needle
"""

import netCDF4
import numpy as np
import pandas as pd
import sys

def extract_per_needle_stats(exodus_file, output_csv='per_needle_stats.csv'):
    """
    Extract per-needle statistics from exodus file
    
    Parameters:
    -----------
    exodus_file : str
        Path to exodus output file (e.g., 'tiny_aragonite_out.e')
    output_csv : str
        Path to output CSV file
    """
    
    print(f"Reading {exodus_file}...")
    
    # Open exodus file
    nc = netCDF4.Dataset(exodus_file, 'r')
    
    # Get time values
    times = nc.variables['time_whole'][:]
    n_times = len(times)
    print(f"Found {n_times} time steps")
    
    # Get number of elements
    n_elem = nc.dimensions['num_elem'].size
    print(f"Found {n_elem} elements")
    
    # Find variable indices (handle MaskedConstant from netCDF4)
    var_names = nc.variables['name_elem_var'][:]
    var_names_str = []
    for name in var_names:
        # Convert to string, handling masked values and bytes
        chars = []
        for c in name:
            # Skip masked constants (padding)
            if hasattr(c, 'mask') or str(type(c)) == "<class 'numpy.ma.core.MaskedConstant'>":
                continue
            # Decode bytes to string
            if isinstance(c, bytes):
                chars.append(c.decode('utf-8'))
            else:
                chars.append(str(c))
        var_str = ''.join(chars).strip()
        var_names_str.append(var_str)
    
    print(f"Available element variables: {var_names_str}")
    
    # Find indices of variables we need
    try:
        needle_id_idx = var_names_str.index('needle_id')
        stress_xx_idx = var_names_str.index('stress_xx')
        strain_xx_idx = var_names_str.index('strain_xx')
        plastic_strain_idx = var_names_str.index('plastic_strain_equiv')
        vonmises_idx = var_names_str.index('vonmises_stress')
    except ValueError as e:
        print(f"Error: Required variable not found: {e}")
        print(f"Available variables: {var_names_str}")
        nc.close()
        return
    
    # Read needle_id (constant in time)
    # Variable name format: vals_elem_varXeb1 where X is index+1
    needle_id_var = f'vals_elem_var{needle_id_idx+1}eb1'
    needle_ids = nc.variables[needle_id_var][0, :]  # First timestep
    
    unique_needles = np.unique(needle_ids)
    n_needles = len(unique_needles)
    print(f"Found {n_needles} unique needles: {unique_needles}")
    
    # Prepare output data
    results = []
    
    # Loop over time steps
    for t_idx, time in enumerate(times):
        print(f"Processing time step {t_idx+1}/{n_times} (t={time:.2f})...", end='\r')
        
        # Read variables at this time step
        stress_xx_var = f'vals_elem_var{stress_xx_idx+1}eb1'
        strain_xx_var = f'vals_elem_var{strain_xx_idx+1}eb1'
        plastic_var = f'vals_elem_var{plastic_strain_idx+1}eb1'
        vonmises_var = f'vals_elem_var{vonmises_idx+1}eb1'
        
        stress_xx = nc.variables[stress_xx_var][t_idx, :]
        strain_xx = nc.variables[strain_xx_var][t_idx, :]
        plastic_strain = nc.variables[plastic_var][t_idx, :]
        vonmises = nc.variables[vonmises_var][t_idx, :]
        
        # Compute statistics per needle
        for needle in unique_needles:
            mask = needle_ids == needle
            
            row = {
                'time': time,
                'needle_id': int(needle),
                'n_elements': np.sum(mask),
                'stress_xx_avg': np.mean(stress_xx[mask]),
                'stress_xx_max': np.max(stress_xx[mask]),
                'stress_xx_min': np.min(stress_xx[mask]),
                'stress_xx_std': np.std(stress_xx[mask]),
                'strain_xx_avg': np.mean(strain_xx[mask]),
                'strain_xx_max': np.max(strain_xx[mask]),
                'strain_xx_min': np.min(strain_xx[mask]),
                'plastic_strain_avg': np.mean(plastic_strain[mask]),
                'plastic_strain_max': np.max(plastic_strain[mask]),
                'vonmises_avg': np.mean(vonmises[mask]),
                'vonmises_max': np.max(vonmises[mask]),
            }
            
            results.append(row)
    
    print()  # New line after progress
    
    # Close exodus file
    nc.close()
    
    # Create DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    
    print(f"\n✓ Saved per-needle statistics to {output_csv}")
    print(f"  Total rows: {len(df)}")
    print(f"  Time steps: {n_times}")
    print(f"  needles: {n_needles}")
    
    # Print summary
    print("\n=== Summary Statistics ===")
    for needle in unique_needles:
        needle_data = df[df['needle_id'] == needle]
        max_stress = needle_data['stress_xx_avg'].max()
        max_plastic = needle_data['plastic_strain_max'].max()
        print(f"needle {int(needle)}: max_stress={max_stress:.1f} MPa, max_plastic={max_plastic:.4f}")
    
    return df

def plot_per_needle(csv_file='per_needle_stats.csv'):
    """
    Create plots of per-needle stress-strain curves
    """
    import matplotlib.pyplot as plt
    
    df = pd.read_csv(csv_file)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Stress vs Time per needle
    ax = axes[0, 0]
    for needle in df['needle_id'].unique():
        needle_data = df[df['needle_id'] == needle]
        ax.plot(needle_data['time'], needle_data['stress_xx_avg'], 
                label=f'needle {int(needle)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Stress XX (MPa)')
    ax.set_title('Stress vs Time per needle')
#    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Stress-Strain per needle
    ax = axes[0, 1]
    for needle in df['needle_id'].unique():
        needle_data = df[df['needle_id'] == needle]
        ax.plot(needle_data['strain_xx_avg'], needle_data['stress_xx_avg'],
                label=f'needle {int(needle)}', linewidth=2)
    ax.set_xlabel('Average Strain XX')
    ax.set_ylabel('Average Stress XX (MPa)')
    ax.set_title('Stress-Strain per needle')
#    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Plastic strain vs Time
    ax = axes[1, 0]
    for needle in df['needle_id'].unique():
        needle_data = df[df['needle_id'] == needle]
        ax.plot(needle_data['time'], needle_data['plastic_strain_avg'],
                label=f'needle {int(needle)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Plastic Strain')
    ax.set_title('Plastic Strain vs Time per needle')
#    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Von Mises vs Time
    ax = axes[1, 1]
    for needle in df['needle_id'].unique():
        needle_data = df[df['needle_id'] == needle]
        ax.plot(needle_data['time'], needle_data['vonmises_avg'],
                label=f'needle {int(needle)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Von Mises Stress (MPa)')
    ax.set_title('Von Mises Stress vs Time per needle')
#    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('per_needle_plots.png', dpi=150)
    print(f"\n✓ Saved plots to per_needle_plots.png")
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python extract_per_needle_stats.py <exodus_file> [output_csv]")
        print("Example: python extract_per_needle_stats.py tiny_aragonite_out.e")
        sys.exit(1)
    
    exodus_file = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else 'per_needle_stats.csv'
    
    # Extract statistics
    df = extract_per_needle_stats(exodus_file, output_csv)
    
    # Create plots if matplotlib available
    try:
        plot_per_needle(output_csv)
    except ImportError:
        print("\nNote: Install matplotlib to generate plots:")
        print("  pip install matplotlib")