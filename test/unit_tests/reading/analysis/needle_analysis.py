#!/usr/bin/env python3
"""
Simplified per-needle analysis - minimal dependencies
Only requires: netCDF4 (no numpy, no pandas needed)
"""

import netCDF4
import sys

def simple_needle_stats(exodus_file):
    """Extract per-needle stats with minimal dependencies"""
    
    print(f"Reading {exodus_file}...")
    
    # Open exodus file
    nc = netCDF4.Dataset(exodus_file, 'r')
    
    # Get variable names (handle MaskedConstant from netCDF4)
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
    
    print(f"Available variables: {var_names_str}")
    
    # Find indices
    try:
        needle_id_idx = var_names_str.index('needle_id')
        stress_idx = var_names_str.index('stress_xx')
        strain_idx = var_names_str.index('strain_xx')
        plastic_idx = var_names_str.index('plastic_strain_equiv')
    except ValueError as e:
        print(f"Error: Required variable not found: {e}")
        nc.close()
        return
    
    # Read needle_id (time=0)
    needle_var = f'vals_elem_var{needle_id_idx+1}eb1'
    needle_ids = list(nc.variables[needle_var][0, :])
    
    # Find unique needles (manual, no numpy)
    unique_needles = sorted(list(set(needle_ids)))
    n_needles = len(unique_needles)
    
    print(f"\nFound {n_needles} unique needles: {unique_needles}")
    
    # Get times
    times = list(nc.variables['time_whole'][:])
    
    # Analyze a few time steps
    print("\n" + "="*70)
    print("PER-needle ANALYSIS")
    print("="*70)
    
    # Time indices to check
    check_times = [0, len(times)//4, len(times)//2, 3*len(times)//4, -1]
    
    for t_idx in check_times:
        if t_idx >= len(times):
            continue
            
        time = times[t_idx]
        
        print(f"\nTime = {time:.1f}")
        print("-" * 70)
        
        # Read variables
        stress_var = f'vals_elem_var{stress_idx+1}eb1'
        strain_var = f'vals_elem_var{strain_idx+1}eb1'
        plastic_var = f'vals_elem_var{plastic_idx+1}eb1'
        
        stress = list(nc.variables[stress_var][t_idx, :])
        strain = list(nc.variables[strain_var][t_idx, :])
        plastic = list(nc.variables[plastic_var][t_idx, :])
        
        # Compute stats per needle (manual loops, no numpy)
        for needle in unique_needles:
            # Find elements in this needle
            stress_needle = []
            strain_needle = []
            plastic_needle = []
            
            for i, gid in enumerate(needle_ids):
                if gid == needle:
                    stress_needle.append(stress[i])
                    strain_needle.append(strain[i])
                    plastic_needle.append(plastic[i])
            
            n_elem = len(stress_needle)
            
            # Compute averages
            stress_avg = sum(stress_needle) / n_elem if n_elem > 0 else 0
            strain_avg = sum(strain_needle) / n_elem if n_elem > 0 else 0
            plastic_max = max(plastic_needle) if plastic_needle else 0
            plastic_avg = sum(plastic_needle) / n_elem if n_elem > 0 else 0
            
            # Compute E_effective
            E_eff = stress_avg / strain_avg if strain_avg > 1e-10 else 0
            
            print(f"  needle {int(needle)}: "
                  f"σ={stress_avg:7.1f} MPa, "
                  f"ε={strain_avg:.6f}, "
                  f"E={E_eff:6.0f} MPa, "
                  f"εₚ_max={plastic_max:.6f}")
    
    # Find first yield per needle
    print("\n" + "="*70)
    print("YIELDING SEQUENCE")
    print("="*70)
    
    yield_times = {}
    
    for needle in unique_needles:
        for t_idx, time in enumerate(times):
            plastic_var = f'vals_elem_var{plastic_idx+1}eb1'
            plastic = list(nc.variables[plastic_var][t_idx, :])
            
            # Get plastic strain for this needle
            plastic_needle = [plastic[i] for i, gid in enumerate(needle_ids) if gid == needle]
            plastic_max = max(plastic_needle) if plastic_needle else 0
            
            if plastic_max > 1e-6:
                # First yield detected
                if needle not in yield_times:
                    stress_var = f'vals_elem_var{stress_idx+1}eb1'
                    stress = list(nc.variables[stress_var][t_idx, :])
                    stress_needle = [stress[i] for i, gid in enumerate(needle_ids) if gid == needle]
                    stress_avg = sum(stress_needle) / len(stress_needle)
                    
                    yield_times[needle] = (time, stress_avg)
                    print(f"needle {int(needle)} yields at time={time:6.1f}, σ={stress_avg:7.1f} MPa")
                break
    
    # Report needles that haven't yielded
    for needle in unique_needles:
        if needle not in yield_times:
            print(f"needle {int(needle)} has not yielded yet")
    
    nc.close()
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python simple_needle_analysis.py <exodus_file>")
        print("Example: python simple_needle_analysis.py tiny_aragonite_out.e")
        sys.exit(1)
    
    exodus_file = sys.argv[1]
    simple_needle_stats(exodus_file)