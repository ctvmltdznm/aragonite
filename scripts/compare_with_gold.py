#!/usr/bin/env python3
"""Compare MOOSE CSV output with gold file."""

import sys
import csv
import argparse
from pathlib import Path

def compare_csv(test_file, gold_file, verbose=False, max_values_only=True):
    """Compare two CSV files with variable-specific tolerances.
    
    Args:
        test_file: Path to test output CSV
        gold_file: Path to gold reference CSV
        verbose: Print detailed comparison info
        max_values_only: Compare only maximum values (default: True)
                        If False, compares every timestep
    """
    
    # Tolerances
    STRESS_REL_TOL = 1e-3  # 0.1%
    STRESS_ABS_TOL = 1.0   # 1 MPa
    STRAIN_REL_TOL = 1e-3  # 0.1%
    STRAIN_ABS_TOL = 1e-8
    DEFAULT_REL_TOL = 1e-3
    DEFAULT_ABS_TOL = 1e-6
    
    # Read files
    with open(test_file) as f:
        test_data = list(csv.DictReader(f))
    with open(gold_file) as f:
        gold_data = list(csv.DictReader(f))
    
    # Check columns match
    test_cols = set(test_data[0].keys())
    gold_cols = set(gold_data[0].keys())
    if test_cols != gold_cols:
        if verbose:
            print(f"Column mismatch:")
            print(f"  Test: {sorted(test_cols)}")
            print(f"  Gold: {sorted(gold_cols)}")
        return False
    
    # Compare based on mode
    if max_values_only:
        return compare_max_values(test_data, gold_data, test_cols, verbose)
    else:
        return compare_all_timesteps(test_data, gold_data, test_cols, verbose)

def get_tolerance(col):
    """Get tolerance based on variable name."""
    STRESS_REL_TOL = 1e-3
    STRESS_ABS_TOL = 1.0
    STRAIN_REL_TOL = 1e-3
    STRAIN_ABS_TOL = 1e-8
    DEFAULT_REL_TOL = 1e-3
    DEFAULT_ABS_TOL = 1e-6
    
    if 'stress' in col.lower():
        return STRESS_REL_TOL, STRESS_ABS_TOL
    elif 'strain' in col.lower():
        return STRAIN_REL_TOL, STRAIN_ABS_TOL
    else:
        return DEFAULT_REL_TOL, DEFAULT_ABS_TOL

def check_value_tolerance(test_val, gold_val, rel_tol, abs_tol):
    """Check if values are within tolerance.
    
    Returns: (passed, error_message)
    """
    abs_diff = abs(test_val - gold_val)
    
    # For large values: check relative tolerance
    if abs(gold_val) > abs_tol:
        rel_diff = abs_diff / abs(gold_val)
        if rel_diff > rel_tol:
            return False, f"rel err: {rel_diff:.2e} (>{rel_tol:.2e})"
    # For small values: check absolute tolerance
    else:
        if abs_diff > abs_tol:
            return False, f"abs err: {abs_diff:.2e} (>{abs_tol:.2e})"
    
    return True, None

def compare_max_values(test_data, gold_data, cols, verbose):
    """Compare only maximum absolute values for each variable."""
    
    if verbose:
        print("Comparing maximum values only...")
    
    for col in cols:
        if col == 'time':
            continue
        
        # Get max absolute values
        test_vals = [abs(float(row[col])) for row in test_data]
        gold_vals = [abs(float(row[col])) for row in gold_data]
        
        test_max = max(test_vals)
        gold_max = max(gold_vals)
        
        # Get tolerance
        rel_tol, abs_tol = get_tolerance(col)
        
        # Check tolerance
        passed, error_msg = check_value_tolerance(test_max, gold_max, rel_tol, abs_tol)
        
        if not passed:
            if verbose:
                print(f"FAIL: {col}")
                print(f"  Test max: {test_max:.6e}")
                print(f"  Gold max: {gold_max:.6e}")
                print(f"  Error: {error_msg}")
            return False
        
        if verbose:
            abs_diff = abs(test_max - gold_max)
            if abs(gold_max) > abs_tol:
                rel_diff = abs_diff / abs(gold_max)
                print(f"PASS: {col}: {test_max:.6e} vs {gold_max:.6e} (rel: {rel_diff:.2e})")
            else:
                print(f"PASS: {col}: {test_max:.6e} vs {gold_max:.6e} (abs: {abs_diff:.2e})")
    
    return True

def compare_all_timesteps(test_data, gold_data, cols, verbose):
    """Compare every timestep (original behavior)."""
    
    # Check same number of rows
    if len(test_data) != len(gold_data):
        if verbose:
            print(f"Row count mismatch: {len(test_data)} vs {len(gold_data)}")
        return False
    
    if verbose:
        print(f"Comparing all {len(test_data)} timesteps...")
    
    max_error = 0
    worst_var = ""
    worst_timestep = 0
    
    for col in cols:
        if col == 'time':
            continue
        
        rel_tol, abs_tol = get_tolerance(col)
        
        for i, (test_row, gold_row) in enumerate(zip(test_data, gold_data)):
            test_val = float(test_row[col])
            gold_val = float(gold_row[col])
            
            passed, error_msg = check_value_tolerance(test_val, gold_val, rel_tol, abs_tol)
            
            if not passed:
                if verbose:
                    print(f"FAIL at timestep {i}, {col}:")
                    print(f"  Test: {test_val:.6e}")
                    print(f"  Gold: {gold_val:.6e}")
                    print(f"  {error_msg}")
                return False
            
            # Track worst error for reporting
            abs_diff = abs(test_val - gold_val)
            if abs(gold_val) > abs_tol:
                rel_diff = abs_diff / abs(gold_val)
                if rel_diff > max_error:
                    max_error = rel_diff
                    worst_var = col
                    worst_timestep = i
    
    if verbose and max_error > 0:
        print(f"Max error: {max_error:.2e} in {worst_var} at timestep {worst_timestep}")
    
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare MOOSE CSV output with gold file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare max values only (default)
  %(prog)s test_out.csv gold/test_out.csv
  
  # Compare all timesteps
  %(prog)s test_out.csv gold/test_out.csv --all-timesteps
  
  # Verbose output
  %(prog)s test_out.csv gold/test_out.csv -v
        """
    )
    parser.add_argument('test_file', type=Path, help='Test output CSV file')
    parser.add_argument('gold_file', type=Path, help='Gold reference CSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed comparison results')
    parser.add_argument('--all-timesteps', action='store_true',
                       help='Compare all timesteps instead of just max values')
    args = parser.parse_args()
    
    max_values_only = not args.all_timesteps
    
    if compare_csv(args.test_file, args.gold_file, args.verbose, max_values_only):
        if args.verbose:
            print("\nComparison PASSED")
        sys.exit(0)
    else:
        if args.verbose:
            print("\nComparison FAILED")
        sys.exit(1)

