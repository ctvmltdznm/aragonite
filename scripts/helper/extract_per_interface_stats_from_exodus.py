#!/usr/bin/env python3
"""
extract_per_interface_stats.py  (PATCHED - side_set version)
=============================================================

CHANGE from original: reads interfaces from side_sets instead of
an interface_id element variable.

Your exodus side set names:
    block0_block1, block0_block2, block0_block3,
    block1_block2, block1_block3, block2_block3

Usage:
    python extract_per_interface_stats.py faster_load.e
    → per_interface_stats.csv   (feeds directly to plot_polycrystal_czm_rve.py)
"""

import netCDF4
import numpy as np
import pandas as pd
import argparse
import sys


def get_variable_names(nc, var_type='elem'):
    """Extract variable names from exodus, handling masked arrays."""
    key_map = {'elem': 'name_elem_var', 'nodal': 'name_nod_var'}
    var_key = key_map.get(var_type, 'name_elem_var')
    if var_key not in nc.variables:
        return []
    var_names_str = []
    for name in nc.variables[var_key][:]:
        chars = []
        for c in name:
            if hasattr(c, 'mask') or 'MaskedConstant' in str(type(c)):
                continue
            chars.append(c.decode('utf-8') if isinstance(c, bytes) else str(c))
        var_str = ''.join(chars).strip()
        if var_str:
            var_names_str.append(var_str)
    return var_names_str


def get_sideset_names(nc):
    """
    Return dict {sideset_name: 0-based_index}.
    PATCHED: this replaces the interface_id lookup from the original script.
    """
    if 'ss_names' not in nc.variables:
        print("  Warning: 'ss_names' not found. Trying numeric IDs only.")
        return {}
    names = {}
    for i, name_arr in enumerate(nc.variables['ss_names'][:]):
        chars = []
        for c in name_arr:
            if hasattr(c, 'mask') or 'MaskedConstant' in str(type(c)):
                continue
            chars.append(c.decode('utf-8') if isinstance(c, bytes) else str(c))
        name = ''.join(chars).strip()
        if name:
            names[name] = i
    return names


def get_sideset_elements(nc, ss_index):
    """
    Return 0-based element indices for side set at ss_index.
    PATCHED: replaces per-element interface_id masking from original script.
    """
    var = f'elem_ss{ss_index + 1}'
    if var not in nc.variables:
        print(f"  Warning: '{var}' not in file. Available ss vars: "
              f"{[k for k in nc.variables if k.startswith('elem_ss')]}")
        return np.array([], dtype=int)
    return nc.variables[var][:].astype(int) - 1   # 1-based → 0-based


def extract_per_interface_stats(exodus_file, output_csv='per_interface_stats.csv'):
    """
    Extract per-interface CZM statistics using side_set names.
    """
    print(f"Reading {exodus_file}...")
    nc = netCDF4.Dataset(exodus_file, 'r')

    times = nc.variables['time_whole'][:]
    n_times = len(times)
    print(f"  {n_times} time steps")

    # ── Element variables ──────────────────────────────────────────────────
    elem_vars = get_variable_names(nc, 'elem')
    print(f"\nElement variables ({len(elem_vars)}):")
    for i, v in enumerate(elem_vars):
        print(f"  {i+1:2d}. {v}")

    # ── Side sets (PATCHED key section) ───────────────────────────────────
    all_sidesets = get_sideset_names(nc)
    print(f"\nAll side sets: {list(all_sidesets.keys())}")

    # Keep only block*_block* names (your interface convention)
    interface_sidesets = {
        name: idx for name, idx in all_sidesets.items()
        if name.startswith('block') and '_block' in name[len('block0'):]
    }
    if not interface_sidesets:
        # Fallback: any sideset not in ['left','right','top','bottom','back','front']
        skip = {'left', 'right', 'top', 'bottom', 'back', 'front'}
        interface_sidesets = {
            name: idx for name, idx in all_sidesets.items()
            if name not in skip
        }

    if not interface_sidesets:
        print("❌ No interface side sets identified. Exiting.")
        nc.close()
        return None

    print(f"\nInterface side sets selected:")
    sorted_names = sorted(interface_sidesets.keys())
    interface_id_map = {name: i+1 for i, name in enumerate(sorted_names)}
    interface_elements = {}
    for name in sorted_names:
        elems = get_sideset_elements(nc, interface_sidesets[name])
        interface_elements[name] = elems
        print(f"  ID {interface_id_map[name]}  '{name}'  → {len(elems)} elements")

    # ── Build block offset map for global→local element index ─────────────
    # Exodus stores element variables per-block (eb1, eb2, ...).
    # Side set elem IDs are GLOBAL (1-based across all blocks).
    # We must concatenate all blocks to build a global array.
    n_elem_blk = nc.dimensions['num_el_blk'].size if 'num_el_blk' in nc.dimensions else 1
    block_sizes = []
    for b in range(n_elem_blk):
        dim = f'num_el_in_blk{b+1}'
        block_sizes.append(nc.dimensions[dim].size if dim in nc.dimensions else 0)
    total_elems = sum(block_sizes)
    print(f"\n  Element blocks: {n_elem_blk}, sizes: {block_sizes}, total: {total_elems}")

    def read_all_blocks(nc, var_idx, t_idx):
        """Read element variable across ALL blocks, return global concatenated array."""
        arrays = []
        for b, bsize in enumerate(block_sizes):
            vname = f'vals_elem_var{var_idx+1}eb{b+1}'
            if vname in nc.variables and bsize > 0:
                arrays.append(np.array(nc.variables[vname][t_idx, :], dtype=float))
            elif bsize > 0:
                arrays.append(np.zeros(bsize))
        return np.concatenate(arrays) if arrays else np.array([])

    # ── Variable index lookup ──────────────────────────────────────────────
    czm_var_candidates = {
        'damage':              ['damage_variable', 'damage', 'interface_damage'],
        'normal_jump':         ['interface_normal_jump', 'normal_jump'],
        'tangent_jump':     ['interface_tangent_jump', 'tangent_jump'],
        'normal_traction':     ['interface_normal_traction', 'normal_traction'],
        'tangent_traction': ['interface_tangent_traction', 'tangent_traction'],
        'mode_ratio':          ['mode_ratio', 'phi'],
        'n_contacts':          ['n_contacts', 'num_contacts'],
    }
    found_vars = {}
    for key, candidates in czm_var_candidates.items():
        for c in candidates:
            if c in elem_vars:
                found_vars[key] = elem_vars.index(c)
                print(f"  ✓ {key} → '{c}'")
                break

    if not found_vars:
        print("❌ No CZM element variables found. Available:", elem_vars)
        nc.close()
        return None

    # Warn about missing critical variables
    print("\n  Missing CZM variables:")
    for key in ['tangent_jump', 'tangent_traction']:
        if key not in found_vars:
            print(f"    ⚠️  {key} — mode ratio calculation will fail; panel (d) will be empty")
            print(f"        Add to MOOSE: [AuxVariables] [{key}] + appropriate AuxKernel")

    # ── Main extraction loop ───────────────────────────────────────────────
    results = []
    print(f"\nExtracting data...")

    for t_idx, time in enumerate(times):
        print(f"  timestep {t_idx+1}/{n_times}  t={time:.4f}", end='\r')

        # Load all relevant element variables for this timestep — ALL BLOCKS
        timestep_data = {}
        for key, var_idx in found_vars.items():
            timestep_data[key] = read_all_blocks(nc, var_idx, t_idx)

        for name in sorted_names:
            elem_ids = interface_elements[name]
            if len(elem_ids) == 0:
                continue

            row = {
                'time':           float(time),
                'interface_id':   interface_id_map[name],
                'interface_name': name,
                'n_elements':     len(elem_ids),
            }

            # Per-variable stats
            for key, data in timestep_data.items():
                vals = data[elem_ids]
                row[f'{key}_avg'] = float(np.mean(vals))
                row[f'{key}_max'] = float(np.max(vals))
                row[f'{key}_min'] = float(np.min(vals))
                row[f'{key}_std'] = float(np.std(vals))

            # Derived: δ_eff
            if 'normal_jump' in timestep_data and 'tangent_jump' in timestep_data:
                dn = timestep_data['normal_jump'][elem_ids]
                dt = timestep_data['tangent_jump'][elem_ids]
                d_eff = np.sqrt(dn**2 + dt**2)
                row['delta_eff_avg'] = float(np.mean(d_eff))
                row['delta_eff_max'] = float(np.max(d_eff))

            # Derived: T_eff
            if 'normal_traction' in timestep_data and 'tangent_traction' in timestep_data:
                Tn = timestep_data['normal_traction'][elem_ids]
                Tt = timestep_data['tangent_traction'][elem_ids]
                T_eff = np.sqrt(Tn**2 + Tt**2)
                row['traction_eff_avg'] = float(np.mean(T_eff))
                row['traction_eff_max'] = float(np.max(T_eff))

            # Derived: mode ratio φ
            if 'normal_jump' in timestep_data and 'tangent_jump' in timestep_data:
                dn_avg = row['normal_jump_avg']
                dt_avg = row['tangent_jump_avg']
                row['phi_avg'] = abs(dt_avg) / (abs(dn_avg) + 1e-30)

            results.append(row)

    print()
    nc.close()

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)

    print(f"\n✓ Saved: {output_csv}")
    print(f"  {len(df)} rows  |  {len(sorted_names)} interfaces  |  {n_times} timesteps")
    print(f"  Columns: {list(df.columns)}")

    # Final state summary
    print("\n=== Final State ===")
    t_last = df['time'].max()
    for row in df[df['time'] == t_last].itertuples():
        line = f"  {row.interface_name:20s}"
        if hasattr(row, 'damage_avg'):     line += f"  D={row.damage_avg:.3f}"
        if hasattr(row, 'traction_eff_avg'): line += f"  T={row.traction_eff_avg:.1f} MPa"
        if hasattr(row, 'phi_avg'):        line += f"  φ={row.phi_avg:.2f}"
        print(line)

    return df


def plot_per_interface(csv_file='per_interface_stats.csv'):
    """Quick 4-panel preview."""
    import matplotlib.pyplot as plt

    df = pd.read_csv(csv_file)
    interfaces = sorted(df['interface_id'].unique())
    cmap = plt.cm.tab10
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    def label(df, iface):
        """Convert interface name to readable format"""
        row = df[df['interface_id'] == iface]
        if 'interface_name' in df.columns:
            name = row['interface_name'].iloc[0]
            # Convert "Block0_Block1" → "Interface 0-1"
            if name.startswith('Block'):
                parts = name.replace('Block', '').split('_')
                #return f'Interface {parts[0]}'
                if len(parts) == 2:
                    return f'Interface {parts[0]}-{parts[1]}'
            return name
        return f'Interface {iface}'

    # (a) Damage
    ax = axes[0, 0]
    if 'damage_avg' in df.columns:
        for i, iface in enumerate(interfaces):
            d = df[df['interface_id'] == iface]
            ax.plot(d['time'], d['damage_avg'], lw=2, label=label(df, iface), color=cmap(i%10))
            #ax.set_yscale('log')

    ax.set(xlabel='Time [s]', ylabel='Damage D', title='(a) Interface Damage Evolution',
           ylim=[0, 0.52])
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # (b) Traction
    ax = axes[0, 1]
    tcol = 'traction_eff_avg' if 'traction_eff_avg' in df.columns else 'normal_traction_avg'
    if tcol in df.columns:
        for i, iface in enumerate(interfaces):
            d = df[df['interface_id'] == iface]
            ax.plot(d['time'], d[tcol], lw=2, label=label(df, iface), color=cmap(i%10))
    ax.set(xlabel='Time [s]', ylabel='Traction [MPa]', title='(b) Traction Evolution')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # (c) Traction-separation
    ax = axes[1, 0]
    jcol = 'delta_eff_avg' if 'delta_eff_avg' in df.columns else 'normal_jump_avg'
    if jcol in df.columns and tcol in df.columns:
        for i, iface in enumerate(interfaces):
            d = df[df['interface_id'] == iface]
            # Data is in mm, convert to nm (×1e6)
            ax.plot(d[jcol]*1e3, d[tcol], lw=2, label=label(df, iface), color=cmap(i%10))
    ax.set(xlabel='Separation [nm]', ylabel='Traction [MPa]',
           title='(c) Traction–Separation')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # (d) Mode ratio - SKIP t=0 initialization
    ax = axes[1, 1]
    if 'phi_avg' in df.columns:
        for i, iface in enumerate(interfaces):
            d = df[df['interface_id'] == iface]
            # Filter out t=0
            d_filtered = d[d['time'] > 0]
            if len(d_filtered) > 0:
                ax.plot(d_filtered['time'], d_filtered['phi_avg'], lw=2, 
                       label=label(df, iface), color=cmap(i%10))
        ax.set(xlabel='Time [s]', ylabel='Mode Ratio φ', title='(d) Mode Mixity')
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'Mode ratio not available\n(add tangent_jump & tangent_traction to exodus)', 
                ha='center', va='center', transform=ax.transAxes, fontsize=10)

    plt.tight_layout()
    out = csv_file.replace('.csv', '_plots.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"✓ Saved plots: {out}")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract per-interface CZM stats (side_set version)')
    parser.add_argument('exodus_file')
    parser.add_argument('--output', '-o', default='per_interface_stats.csv')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()

    df = extract_per_interface_stats(args.exodus_file, args.output)
    if df is not None and not args.no_plot:
        try:
            plot_per_interface(args.output)
        except Exception as e:
            print(f"Plot failed: {e}")

