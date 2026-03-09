#!/bin/bash
# ============================================================================
# Plot Validation Results
# ============================================================================
# Generates comparison plots for all validation tests
#
# Usage:
#   ./plot_validation.sh              # Plot all categories
#   ./plot_validation.sh -c ortho     # Plot specific category
#   ./plot_validation.sh --compare    # Compare with gold files
#
# Output: test/validation/figures/*.png
# ============================================================================

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
TEST_ROOT="$REPO_ROOT/test/validation"
FIGURE_DIR="$TEST_ROOT/figures"

# Create figures directory
mkdir -p "$FIGURE_DIR"

# Colors
BLUE='\033[0;34m'
NC='\033[0m'

# Default settings
CATEGORY_FILTER=""
COMPARE_GOLD=false

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate validation plots for visual verification.

Options:
    -h, --help              Show this help
    -c, --category CAT      Plot only category CAT
    --compare               Overlay gold file data
    --clean                 Remove all figures

Examples:
    $0                      # Plot all tests
    $0 -c ortho             # Plot orthotropic tests only
    $0 --compare            # Compare with gold files

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help) show_usage; exit 0 ;;
        -c|--category) CATEGORY_FILTER="$2"; shift 2 ;;
        --compare) COMPARE_GOLD=true; shift ;;
        --clean) rm -rf "$FIGURE_DIR"/*; echo "Cleaned $FIGURE_DIR"; exit 0 ;;
        *) echo "Unknown: $1"; show_usage; exit 1 ;;
    esac
done

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

plot_category() {
    local category="$1"
    local category_dir="$TEST_ROOT/$category"
    
    if [ ! -d "$category_dir" ]; then
        return 0
    fi
    
    echo -e "${BLUE}Plotting: $category${NC}"
    
    # Call Python plotting script for this category
    python3 "$SCRIPT_DIR/plot_validation_category.py" \
        --category "$category" \
        --test-root "$TEST_ROOT" \
        --output-dir "$FIGURE_DIR" \
        $([ "$COMPARE_GOLD" = true ] && echo "--compare")
}

# ============================================================================
# MAIN
# ============================================================================

echo "========================================="
echo "Validation Plotting"
echo "========================================="
echo "Output: $FIGURE_DIR"
echo ""

# Check Python available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found"
    exit 1
fi

# Check matplotlib available
if ! python3 -c "import matplotlib" 2>/dev/null; then
    echo "ERROR: matplotlib not installed"
    echo "Install with: pip install matplotlib --break-system-packages"
    exit 1
fi

# Plot categories
CATEGORIES=("ortho" "asymmetry" "postyield" "grains" "czm")

for category in "${CATEGORIES[@]}"; do
    if [ -n "$CATEGORY_FILTER" ] && [ "$category" != "$CATEGORY_FILTER" ]; then
        continue
    fi
    plot_category "$category"
done

echo ""
echo "✓ Figures generated in: $FIGURE_DIR"
echo ""
ls -lh "$FIGURE_DIR"

