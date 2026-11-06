#!/bin/bash
# ============================================================================
# Validation Test Suite with Auto-Detection
# ============================================================================

# set -e

# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
TEST_ROOT="$REPO_ROOT/test/validation"

# ---- executable resolution --------------------------------------------------
# The app name comes from the Makefile, NOT from the directory name: this repo
# is "aragonite" on GitHub and "tuc-fe-moose" on GitLab, so basename() is not a
# reliable source. APP=/path/to/app-opt overrides everything.
#
# Several *-opt in the repo root is an ERROR, never a silent pick. The old code
# did `find ... | head -n 1`, which returns readdir order; with both
# aragonite-opt (Oct 23) and tuc_fe_edu-opt (Oct 8) present it ran the stale
# binary and reported only "Differs from gold" for four CZM tests. Symlinks are
# excluded by -type f, so a convenience alias is safe.
APP_NAME=$(awk -F':=' '/^[[:space:]]*APPLICATION_NAME/{gsub(/[[:space:]]/,"",$2); print $2}' \
             "$REPO_ROOT/Makefile" 2>/dev/null)
EXECUTABLE=${APP:-}
if [ -z "$EXECUTABLE" ]; then
    if [ -n "$APP_NAME" ] && [ -x "$REPO_ROOT/${APP_NAME}-opt" ]; then
        EXECUTABLE="$REPO_ROOT/${APP_NAME}-opt"
    else
        mapfile -t CAND < <(find "$REPO_ROOT" -maxdepth 1 -name '*-opt' -type f | sort)
        case ${#CAND[@]} in
          0) echo "ERROR: no *-opt executable in $REPO_ROOT." >&2
             echo "       Expected '${APP_NAME:-<APPLICATION_NAME>}-opt'. Build first: make -j" >&2
             exit 1 ;;
          1) EXECUTABLE="${CAND[0]}" ;;
          *) echo "ERROR: ${#CAND[@]} executables in $REPO_ROOT -- refusing to guess:" >&2
             ls -lt "${CAND[@]}" >&2
             echo "       Makefile says APPLICATION_NAME = ${APP_NAME:-<unset>}." >&2
             echo "       Delete the stale one(s), or set APP=<path>." >&2
             exit 1 ;;
        esac
    fi
fi
[ -x "$EXECUTABLE" ] || { echo "ERROR: $EXECUTABLE is not executable." >&2; exit 1; }

# ---- provenance, printed on EVERY run --------------------------------------
# A suite that does not record which binary and which source it ran cannot be
# compared across machines, and that is how a two-week-old library produced
# four "Differs from gold" failures that looked like a physics regression.
echo "executable : $EXECUTABLE"
echo "built      : $(date -r "$EXECUTABLE" '+%F %T')"
echo "app_name   : ${APP_NAME:-<unset>}"
echo "git        : $(git -C "$REPO_ROOT" describe --always --dirty 2>/dev/null || echo 'not a git repo')"
STALE=$(find "$REPO_ROOT/src" "$REPO_ROOT/include" \( -name '*.C' -o -name '*.h' \) \
          -newer "$EXECUTABLE" 2>/dev/null | head -5)
if [ -n "$STALE" ]; then
    echo "WARNING: source is NEWER than the executable -- rebuild before trusting results:"
    printf '         %s\n' $STALE
fi

COMPARE_SCRIPT="$SCRIPT_DIR/compare_with_gold.py"

# Test categories
declare -A TEST_CATEGORIES=(
    ["ortho"]="Orthotropic plasticity single-element tests"
    ["asymmetry"]="Tension-compression asymmetry tests"
    ["postyield"]="Post-yield behavior modes"
    ["grains"]="Multi-grain bicrystal tests"
    ["czm"]="Cohesive zone model tests"
    ["fabric"]="Fabric-scaling elasticity + plasticity tests"
    ["rve"]="Representative volume element tests"
)

# Default settings
PARALLEL_PROCS=1
RUN_COMPARISON=true
GENERATE_GOLD=false
CATEGORY_FILTER=""
VERBOSE=false
ALL_TIMESTEPS=false

# Results tracking
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0
SKIPPED_TESTS=0

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ============================================================================
# FUNCTIONS
# ============================================================================

print_header() {
    echo ""
    echo "========================================="
    echo "$1"
    echo "========================================="
}

print_category() {
    echo ""
    echo -e "${BLUE}Category: $1${NC}"
    echo -e "${BLUE}Description: $2${NC}"
    echo ""
}

print_test() {
    echo -n "  Testing: $1 ... "
}

print_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
}

print_skip() {
    echo -e "${YELLOW}[SKIP]${NC} $1"
}

print_summary() {
    echo ""
    print_header "Test Summary"
    echo "Total:   $TOTAL_TESTS"
    echo -e "Passed:  ${GREEN}$PASSED_TESTS${NC}"
    echo -e "Failed:  ${RED}$FAILED_TESTS${NC}"
    echo -e "Skipped: ${YELLOW}$SKIPPED_TESTS${NC}"
    echo ""
    
    if [ $FAILED_TESTS -gt 0 ]; then
        echo -e "${RED}Some tests failed.${NC}"
        return 1
    else
        echo -e "${GREEN}All tests passed!${NC}"
        return 0
    fi
}

check_prerequisites() {
    if [ ! -f "$EXECUTABLE" ]; then
        echo -e "${RED}ERROR: Executable not found${NC}"
        echo "Searched: $REPO_ROOT/*-opt"
        echo "Run 'make -j8' first."
        exit 1
    fi
    
    if [ ! -d "$TEST_ROOT" ]; then
        echo -e "${RED}ERROR: Test directory not found: $TEST_ROOT${NC}"
        exit 1
    fi
    
    if [ "$RUN_COMPARISON" = true ] && [ ! -f "$COMPARE_SCRIPT" ]; then
        echo -e "${YELLOW}WARNING: Comparison script not found${NC}"
        RUN_COMPARISON=false
    fi
    
    if [ "$RUN_COMPARISON" = true ] && ! command -v python3 &> /dev/null; then
        echo -e "${YELLOW}WARNING: python3 not found${NC}"
        RUN_COMPARISON=false
    fi
}

run_single_test() {
    local test_dir="$1"
    local test_name=$(basename "$test_dir")
    local input_file="$test_dir/${test_name}.i"
    local output_csv="${test_name}_out.csv"
    local gold_csv="$test_dir/gold/${output_csv}"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    if [ ! -f "$input_file" ]; then
        print_skip "Input file not found"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        return 0
    fi
    
    print_test "$test_name"
    
    cd "$test_dir"
    
    # Per-test rank count: a line "# validation_mpi_procs = N" in the input.
    # -p on the command line overrides it for every test. Without either, serial.
    # (A single global -p would also push the 2-element tests onto N ranks.)
    local procs=$PARALLEL_PROCS
    if [ "$PARALLEL_PROCS" -eq 1 ]; then
        local want
        want=$(grep -m1 -E '^#[[:space:]]*validation_mpi_procs[[:space:]]*=[[:space:]]*[0-9]+' "$input_file" \
                 | grep -oE '[0-9]+$')
        [ -n "$want" ] && procs=$want
    fi
    [ "$procs" -gt 1 ] && echo -n "[np=$procs] "

    # Run simulation (suppress output, ignore segfaults)
    if [ $procs -gt 1 ]; then
        mpirun -n $procs "$EXECUTABLE" -i "$input_file" > /dev/null 2>&1 || true
    else
        "$EXECUTABLE" -i "$input_file" > /dev/null 2>&1 || true
    fi
    
    if [ ! -f "$output_csv" ]; then
        print_fail "No output generated"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Generate gold file if requested
    if [ "$GENERATE_GOLD" = true ]; then
        mkdir -p "gold"
        cp "$output_csv" "$gold_csv"
        print_pass "Gold file generated"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    fi
    
    # Compare with gold
    if [ "$RUN_COMPARISON" = true ] && [ -f "$gold_csv" ]; then
        COMPARE_CMD="python3 \"$COMPARE_SCRIPT\" \"$output_csv\" \"$gold_csv\""
        
        if [ "$ALL_TIMESTEPS" = true ]; then
            COMPARE_CMD="$COMPARE_CMD --all-timesteps"
        fi
        
        if [ "$VERBOSE" = true ]; then
            COMPARE_CMD="$COMPARE_CMD -v"
        fi
        
        if eval $COMPARE_CMD > /dev/null 2>&1; then
            print_pass "Matches gold file"
            PASSED_TESTS=$((PASSED_TESTS + 1))
            return 0
        else
            print_fail "Differs from gold"
            if [ "$VERBOSE" = true ]; then
                eval $COMPARE_CMD || true
            fi
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    else
        print_pass "Simulation completed"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    fi
}

run_category_tests() {
    local category="$1"
    local category_dir="$TEST_ROOT/$category"
    
    if [ ! -d "$category_dir" ]; then
        return 0
    fi
    
    print_category "$category" "${TEST_CATEGORIES[$category]}"
    
    for test_dir in "$category_dir"/*; do
        # A test is a directory <name>/ containing <name>.i. Anything else
        # (__pycache__/, scratch dirs, check scripts' helpers) is not a test
        # and is ignored rather than counted as SKIPPED.
        if [ -d "$test_dir" ] && [ -f "$test_dir/$(basename "$test_dir").i" ]; then
            run_single_test "$test_dir"
        fi
    done
}

run_all_tests() {
    for category in "${!TEST_CATEGORIES[@]}"; do
        if [ -n "$CATEGORY_FILTER" ] && [ "$category" != "$CATEGORY_FILTER" ]; then
            continue
        fi
        run_category_tests "$category"
    done
}

show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    -h, --help              Show help
    -c, --category CAT      Run only category CAT
    -p, --parallel N        Use N MPI processes
    --no-compare            Skip gold comparison
    --generate-gold         Generate gold files
    --all-timesteps         Compare all timesteps (default: max values only)
    -v, --verbose           Show details
    --list                  List tests

Comparison modes:
    Default:        Compares max absolute values only
    --all-timesteps: Compares every timestep (stricter)

Examples:
    $0                          # Run all (max values)
    $0 -c ortho                 # Ortho tests only
    $0 -p 4                     # Parallel
    $0 --all-timesteps          # Compare all timesteps
    $0 --generate-gold -c ortho # Generate gold files

EOF
}

list_tests() {
    print_header "Available Tests"
    for category in "${!TEST_CATEGORIES[@]}"; do
        echo ""
        echo -e "${BLUE}$category:${NC} ${TEST_CATEGORIES[$category]}"
        category_dir="$TEST_ROOT/$category"
        if [ -d "$category_dir" ]; then
            for test_dir in "$category_dir"/*; do
                if [ -d "$test_dir" ]; then
                    test_name=$(basename "$test_dir")
                    if [ -f "$test_dir/${test_name}.i" ]; then
                        echo "  ✓ $test_name"
                    fi
                fi
            done
        fi
    done
    echo ""
}

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help) show_usage; exit 0 ;;
        -c|--category) CATEGORY_FILTER="$2"; shift 2 ;;
        -p|--parallel) PARALLEL_PROCS="$2"; shift 2 ;;
        --no-compare) RUN_COMPARISON=false; shift ;;
        --generate-gold) GENERATE_GOLD=true; shift ;;
        --all-timesteps) ALL_TIMESTEPS=true; shift ;;
        -v|--verbose) VERBOSE=true; shift ;;
        --list) list_tests; exit 0 ;;
        *) echo "Unknown: $1"; show_usage; exit 1 ;;
    esac
done

# ============================================================================
# MAIN
# ============================================================================

print_header "Validation Suite"

echo "Configuration:"
echo "  Repository:    $REPO_ROOT"
echo "  Executable:    $EXECUTABLE"
echo "  Test Root:     $TEST_ROOT"
echo "  MPI Processes: $PARALLEL_PROCS"
echo "  Gold Compare:  $RUN_COMPARISON"
if [ "$RUN_COMPARISON" = true ]; then
    if [ "$ALL_TIMESTEPS" = true ]; then
        echo "  Compare Mode:  All timesteps"
    else
        echo "  Compare Mode:  Max values only"
    fi
fi
echo "  Generate Gold: $GENERATE_GOLD"
[ -n "$CATEGORY_FILTER" ] && echo "  Category:      $CATEGORY_FILTER"

check_prerequisites
run_all_tests
print_summary
exit $?
