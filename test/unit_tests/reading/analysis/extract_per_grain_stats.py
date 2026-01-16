#!/usr/bin/env python3
"""
Extract per-grain stress and strain statistics from MOOSE exodus output
Reads grain_id element variable and computes averages per grain
"""

import netCDF4
import numpy as np
import pandas as pd
import sys

def extract_per_grain_stats(exodus_file, output_csv='per_grain_stats.csv'):
    """
    Extract per-grain statistics from exodus file
    
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
        grain_id_idx = var_names_str.index('grain_id')
        stress_xx_idx = var_names_str.index('stress_xx')
        strain_xx_idx = var_names_str.index('strain_xx')
        plastic_strain_idx = var_names_str.index('plastic_strain_equiv')
        vonmises_idx = var_names_str.index('vonmises_stress')
    except ValueError as e:
        print(f"Error: Required variable not found: {e}")
        print(f"Available variables: {var_names_str}")
        nc.close()
        return
    
    # Read grain_id (constant in time)
    # Variable name format: vals_elem_varXeb1 where X is index+1
    grain_id_var = f'vals_elem_var{grain_id_idx+1}eb1'
    grain_ids = nc.variables[grain_id_var][0, :]  # First timestep
    
    unique_grains = np.unique(grain_ids)
    n_grains = len(unique_grains)
    print(f"Found {n_grains} unique grains: {unique_grains}")
    
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
        
        # Compute statistics per grain
        for grain in unique_grains:
            mask = grain_ids == grain
            
            row = {
                'time': time,
                'grain_id': int(grain),
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
    
    print(f"\n✓ Saved per-grain statistics to {output_csv}")
    print(f"  Total rows: {len(df)}")
    print(f"  Time steps: {n_times}")
    print(f"  Grains: {n_grains}")
    
    # Print summary
    print("\n=== Summary Statistics ===")
    for grain in unique_grains:
        grain_data = df[df['grain_id'] == grain]
        max_stress = grain_data['stress_xx_avg'].max()
        max_plastic = grain_data['plastic_strain_max'].max()
        print(f"Grain {int(grain)}: max_stress={max_stress:.1f} MPa, max_plastic={max_plastic:.4f}")
    
    return df

def plot_per_grain(csv_file='per_grain_stats.csv'):
    """
    Create plots of per-grain stress-strain curves
    """
    import matplotlib.pyplot as plt
    
    df = pd.read_csv(csv_file)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Stress vs Time per grain
    ax = axes[0, 0]
    for grain in df['grain_id'].unique():
        grain_data = df[df['grain_id'] == grain]
        ax.plot(grain_data['time'], grain_data['stress_xx_avg'], 
                label=f'Grain {int(grain)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Stress XX (MPa)')
    ax.set_title('Stress vs Time per Grain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Stress-Strain per grain
    ax = axes[0, 1]
    for grain in df['grain_id'].unique():
        grain_data = df[df['grain_id'] == grain]
        ax.plot(grain_data['strain_xx_avg'], grain_data['stress_xx_avg'],
                label=f'Grain {int(grain)}', linewidth=2)
    ax.set_xlabel('Average Strain XX')
    ax.set_ylabel('Average Stress XX (MPa)')
    ax.set_title('Stress-Strain per Grain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Plastic strain vs Time
    ax = axes[1, 0]
    for grain in df['grain_id'].unique():
        grain_data = df[df['grain_id'] == grain]
        ax.plot(grain_data['time'], grain_data['plastic_strain_avg'],
                label=f'Grain {int(grain)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Plastic Strain')
    ax.set_title('Plastic Strain vs Time per Grain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Von Mises vs Time
    ax = axes[1, 1]
    for grain in df['grain_id'].unique():
        grain_data = df[df['grain_id'] == grain]
        ax.plot(grain_data['time'], grain_data['vonmises_avg'],
                label=f'Grain {int(grain)}', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average Von Mises Stress (MPa)')
    ax.set_title('Von Mises Stress vs Time per Grain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('per_grain_plots.png', dpi=150)
    print(f"\n✓ Saved plots to per_grain_plots.png")
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python extract_per_grain_stats.py <exodus_file> [output_csv]")
        print("Example: python extract_per_grain_stats.py tiny_aragonite_out.e")
        sys.exit(1)
    
    exodus_file = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else 'per_grain_stats.csv'
    
    # Extract statistics
    df = extract_per_grain_stats(exodus_file, output_csv)
    
    # Create plots if matplotlib available
    try:
        plot_per_grain(output_csv)
    except ImportError:
        print("\nNote: Install matplotlib to generate plots:")
        print("  pip install matplotlib")