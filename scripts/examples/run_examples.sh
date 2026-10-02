#!/usr/bin/env bash
# ============================================================================
# run_examples.sh -- run the example inputs sequentially
#
#   ./run_examples.sh                 every example (57 runs)
#   ./run_examples.sh bone            only the bone examples
#   ./run_examples.sh coral           only the coral examples
#   ./run_examples.sh compact_ti      substring filter on the path
#   ./run_examples.sh shear           every shear case, bone and coral
#
#   APP=/path/to/aragonite-opt ./run_examples.sh
#   GOLD=1 ./run_examples.sh          also copy each CSV into gold/ afterwards
#   NOANALYSE=1 ./run_examples.sh     skip the analysis step
#
# Each run writes, next to its input:
#     <case>_out.csv    the postprocessor output MOOSE produces
#     <case>.log        stdout and stderr
# and the exit code is reported inline. At the end analyse_examples.py compares
# every CSV against the values predicted from the presets.
#
# Sequential on purpose: these are single- and two-element runs of a few
# seconds each, and a serial log is easier to read than interleaved output.
# ============================================================================
set -u

APP=${APP:-~/projects/aragonite/aragonite-opt}
ANALYSE=${ANALYSE:-./analyse_examples.py}
FILTER=${1:-}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1

if [ ! -x "$(eval echo "$APP")" ] && ! command -v "$APP" >/dev/null 2>&1; then
  echo "ERROR: MOOSE executable not found at '$APP'."
  echo "       Set it explicitly, e.g.  APP=~/projects/aragonite/aragonite-opt $0"
  exit 1
fi

mapfile -t INPUTS < <(find bone coral -name '*.i' 2>/dev/null | sort)
if [ -n "$FILTER" ]; then
  mapfile -t INPUTS < <(printf '%s\n' "${INPUTS[@]}" | grep -- "$FILTER")
fi

n=${#INPUTS[@]}
if [ "$n" -eq 0 ]; then
  echo "no inputs matched '${FILTER}'"
  exit 1
fi

echo "============================================================"
echo "running $n example(s)   APP=$APP"
[ -n "$FILTER" ] && echo "filter: $FILTER"
echo "============================================================"

fail=0
i=0
start_all=$SECONDS
for inp in "${INPUTS[@]}"; do
  i=$((i + 1))
  dir=$(dirname "$inp")
  base=$(basename "$inp" .i)
  printf '[%2d/%2d] %-52s ' "$i" "$n" "$inp"
  t0=$SECONDS
  ( cd "$dir" && eval "$APP" -i "$base.i" ) > "$dir/$base.log" 2>&1
  rc=$?
  dt=$((SECONDS - t0))

  if [ $rc -ne 0 ]; then
    fail=$((fail + 1))
    echo "FAILED (exit $rc, ${dt}s)"
    # the first error line is almost always the useful one
    grep -m1 -A3 '\*\*\* ERROR' "$dir/$base.log" | sed 's/^/          /'
  elif [ ! -f "$dir/${base}_out.csv" ]; then
    fail=$((fail + 1))
    echo "NO CSV (exit 0, ${dt}s) -- did [Outputs] csv = true survive?"
  else
    rows=$(($(wc -l < "$dir/${base}_out.csv") - 1))
    echo "ok (${dt}s, $rows rows)"
  fi

  # The startup self-checks abort on failure, so reaching this point means they
  # passed. Surface them anyway the first time, so it is visible they ran.
  if [ "$i" -eq 1 ]; then
    grep -m1 'FD check' "$dir/$base.log" | sed 's/^/          /'
    grep -m1 'Mandel round-trip' "$dir/$base.log" | sed 's/^/          /'
  fi

  if [ "${GOLD:-0}" = "1" ] && [ -f "$dir/${base}_out.csv" ]; then
    mkdir -p "$dir/gold"
    cp "$dir/${base}_out.csv" "$dir/gold/${base}_out.csv"
  fi
done

echo "============================================================"
echo "$((n - fail))/$n ran, $fail failed, $((SECONDS - start_all))s total"
[ "${GOLD:-0}" = "1" ] && echo "gold files written to <dir>/gold/"
echo "============================================================"

if [ "${NOANALYSE:-0}" != "1" ] && [ -f "$ANALYSE" ]; then
  echo
  python3 "$ANALYSE" ${FILTER:+--filter "$FILTER"}
fi

exit $([ $fail -eq 0 ] && echo 0 || echo 1)
