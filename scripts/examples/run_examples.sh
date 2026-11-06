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
#   APP=/path/to/app-opt ./run_examples.sh   override the executable
#                                     (default: <repo>/<APPLICATION_NAME>-opt,
#                                      i.e. aragonite-opt on GitHub,
#                                      tuc-fe-moose-opt on GitLab)
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

ANALYSE=${ANALYSE:-./analyse_examples.py}
FILTER=${1:-}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
cd "$HERE" || exit 1

# ---- executable resolution (same rules as scripts/run_validation.sh) --------
# The app name comes from the Makefile's APPLICATION_NAME, not from a hardcoded
# path: the repo builds aragonite-opt on GitHub and tuc-fe-moose-opt on GitLab.
# APP=<path> overrides; several *-opt in the repo root is an error, never a
# silent pick (a stale second binary once produced bogus gold mismatches).
APP_NAME=$(awk -F':=' '/^[[:space:]]*APPLICATION_NAME/{gsub(/[[:space:]]/,"",$2); print $2}' \
             "$REPO_ROOT/Makefile" 2>/dev/null)
APP=${APP:-}
APP=${APP/#\~/$HOME}                 # allow APP=~/... without eval
if [ -z "$APP" ]; then
  if [ -n "$APP_NAME" ] && [ -x "$REPO_ROOT/${APP_NAME}-opt" ]; then
    APP="$REPO_ROOT/${APP_NAME}-opt"
  else
    mapfile -t CAND < <(find "$REPO_ROOT" -maxdepth 1 -name '*-opt' -type f | sort)
    case ${#CAND[@]} in
      0) echo "ERROR: no *-opt executable in $REPO_ROOT." >&2
         echo "       Expected '${APP_NAME:-<APPLICATION_NAME>}-opt'. Build first: make -j2" >&2
         exit 1 ;;
      1) APP="${CAND[0]}" ;;
      *) echo "ERROR: ${#CAND[@]} executables in $REPO_ROOT -- refusing to guess:" >&2
         ls -lt "${CAND[@]}" >&2
         echo "       Makefile says APPLICATION_NAME = ${APP_NAME:-<unset>}." >&2
         echo "       Delete the stale one(s), or set APP=<path>." >&2
         exit 1 ;;
    esac
  fi
fi
if [ ! -x "$APP" ]; then
  echo "ERROR: '$APP' is not an executable file." >&2
  exit 1
fi

# provenance, so a log can always be traced to a binary and a source state
echo "executable : $APP"
echo "built      : $(date -r "$APP" '+%F %T')"
echo "git        : $(git -C "$REPO_ROOT" describe --always --dirty 2>/dev/null || echo 'not a git repo')"
STALE=$(find "$REPO_ROOT/src" "$REPO_ROOT/include" \( -name '*.C' -o -name '*.h' \) \
          -newer "$APP" 2>/dev/null | head -5)
if [ -n "$STALE" ]; then
  echo "WARNING: source is NEWER than the executable -- rebuild before trusting results:"
  printf '         %s\n' $STALE
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
  ( cd "$dir" && "$APP" -i "$base.i" ) > "$dir/$base.log" 2>&1
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
