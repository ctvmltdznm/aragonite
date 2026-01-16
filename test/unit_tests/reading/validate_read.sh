#!/bin/bash
# Quick check if exodus test worked

CSV="tiny_aragonite_out.csv"

echo "=== EXODUS TEST VALIDATION ==="

# Check 1: File exists
if [ -f "$CSV" ]; then
    echo "✓ Output CSV exists"
else
    echo "✗ No CSV output - simulation failed?"
    exit 1
fi

# Check 2: Euler angles loaded
euler_range=$(awk -F',' 'NR==2 {print $8 - $7}' "$CSV")
if (( $(echo "$euler_range > 10" | bc -l) )); then
    echo "✓ Euler angles varying (range = $euler_range°)"
else
    echo "✗ Euler angles not varying (all same?)"
fi

# Check 3: Stress variation
stress_range=$(awk -F',' 'NR==11 {print $11}' "$CSV")  # time=10
if (( $(echo "$stress_range > 100" | bc -l) )); then
    echo "✓ Stress varying per element (range = $stress_range MPa)"
else
    echo "✗ No stress variation (rotation not working?)"
fi

# Check 4: Plasticity occurred
max_plastic=$(awk -F',' 'END {print $16}' "$CSV")
if (( $(echo "$max_plastic > 0.001" | bc -l) )); then
    echo "✓ Plasticity occurred (max = $max_plastic)"
else
    echo "⚠ No plasticity yet (stopped early?)"
fi

echo "=== DONE ==="
