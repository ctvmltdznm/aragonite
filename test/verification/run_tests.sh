#!/usr/bin/env bash
# ============================================================================
# run_tests.sh - tangent-work test driver
#
#   ./run_tests.sh fd            FD unit check of yield gradient/Hessian (1 step)
#   ./run_tests.sh A B C         single-element matrix (A regression, B rotated, C twin)
#   ./run_tests.sh Bt Bs Bf      B variants: tight local tol / SMALL strain / FDJ
#   ./run_tests.sh D             bone cube, production settings, 3 steps
#   NOREF=1 ./run_tests.sh D     ... without the expensive elastic reference
#   ./run_tests.sh Ddt           bone step-size study (dt = 0.025)
#   ./run_tests.sh Dtol          bone tolerance control (nl_abs_tol = 1e-8)
#   ./run_tests.sh V             every viscosity mode on the single element
#   MODE=exponential ./run_tests.sh Vr   rate-dependence + dt-independence check
#   ./run_tests.sh P             forces the Primal CPPA fallback (tests the 7x7 path)
#
# Every matrix test is run twice with the SAME binary:
#   *_ref   tangent_operator = elastic   (reference: tangent cannot change the answer)
#   *_full  code default (FULL)          (under test)
# then analyze_run.py prints solver health for both and compares the CSVs.
#
# Results land in ./results/<ID>_{ref,full}.{log,csv} and ./results/<ID>.summary
# ============================================================================
set -u

# ---- edit these -------------------------------------------------------------
APP=${APP:-~/projects/aragonite/aragonite-opt}                    # MOOSE executable
MPI_SMALL=${MPI_SMALL:-""}                     # single-element runs: serial is fine
MPI_BONE=${MPI_BONE:-"mpiexec -n 32"}           # bone cube
LH=${LH:-linear_hardening.i}
TWIN=${TWIN:-twin.i}
BONE=${BONE:-cube_tension_22.i}
BASELINE_LH=${BASELINE_LH:-baseline_lh.csv}    # = the linear_hardening_out.csv from Step 2
ANALYZE=${ANALYZE:-./analyse_run.py}
TOL=${TOL:-1e-5}                               # relative, per CSV column

# Production solver settings, applied to every bone run. These replace the
# values that were tuned around a broken tangent and are now all too loose.
# ksp_rtol lives inside petsc_options_value, so it is overridden as a vector.
# (array, because the petsc option string contains spaces)
PROD=(
  Materials/plasticity/absolute_tolerance=1e-8
  Executioner/line_search=bt
  Executioner/nl_abs_tol=1e-6
  "Executioner/petsc_options_value=gamg 2 gmres 300 200 1e-6"
  Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true
)
# -----------------------------------------------------------------------------

mkdir -p results
QUIET="Outputs/print_linear_residuals=false"
# FDJ=1 ./run_tests.sh ...  -> every run uses MOOSE's corotational
# finite-deformation Jacobian (needed for quadratic convergence with FINITE strain)
[ "${FDJ:-0}" = "1" ] && QUIET="$QUIET Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true"

run () {   # run <tag> <mpi> <input> <extra overrides...>
  local tag=$1 mpi=$2 input=$3; shift 3
  echo ">>> $tag"
  $mpi $APP -i "$input" "$@" $QUIET > "results/$tag.log" 2>&1
  echo "    exit $?  ($(wc -l < results/$tag.log) log lines)"
}

pair () { # pair <ID> <mpi> <input> <csv-override-key> <extra overrides...>
  # optional: set TMAX=<t> before calling to limit the CSV comparison
  local id=$1 mpi=$2 input=$3 csvkey=$4; shift 4
  if [ "${NOREF:-0}" = "1" ]; then
    run "${id}_full" "$mpi" "$input" "$csvkey=results/${id}_full" "$@"
    python3 "$ANALYZE" "results/${id}_full.log" | tee "results/${id}.summary"
    return
  fi
  run "${id}_ref"  "$mpi" "$input" Materials/stress/tangent_operator=elastic \
      "$csvkey=results/${id}_ref"  "$@"
  run "${id}_full" "$mpi" "$input" \
      "$csvkey=results/${id}_full" "$@"
  python3 "$ANALYZE" "results/${id}_full.log" --ref-log "results/${id}_ref.log" \
      --ref "results/${id}_ref.csv" --test "results/${id}_full.csv" --tol "$TOL" \
      ${TMAX:+--tmax $TMAX} \
      | tee "results/${id}.summary"
}

for T in "$@"; do
  echo
  echo "=================================== $T"
  case $T in

  fd)   # needs Patch F compiled in; one step is enough, the check fires at t_step 1
    run fd "$MPI_SMALL" "$LH" Executioner/num_steps=1 Outputs/file_base=results/fd
    grep -m1 "FD check" results/fd.log || echo "    no FD line - is Patch F compiled in?"
    grep -m1 -i "error" results/fd.log
    ;;

  A)    # sentinel: axis-aligned, no material-frame shear
    pair A "$MPI_SMALL" "$LH" Outputs/file_base
    if [ -f "$BASELINE_LH" ]; then
      echo "--- A_full vs Step-2 baseline (must be unchanged: no shear in this test)"
      python3 "$ANALYZE" results/A_full.log --ref "$BASELINE_LH" --test results/A_full.csv \
          --tol "$TOL" | sed -n '/CSV/,$p' | tee -a results/A.summary
    fi
    ;;

  B)    # key test: rotation + material-frame shear, no softening, no viscosity
    pair B "$MPI_SMALL" "$LH" Outputs/file_base \
      ICs/euler_phi1_0/value=30 ICs/euler_Phi_0/value=45 ICs/euler_phi2_0/value=60 \
      Materials/plasticity/viscosity_mode=rate_independent
    ;;

  V)    # item 2: every viscosity mode on the single element (A geometry).
        # M_FOR_MODE overrides the exponent where the default is unusable.
    for mode in linear exponential logarithmic polynomial powerlaw; do
      case $mode in
        powerlaw) MV=0.2 ;;   # visc ~ x^(1/m); 1/m = 1000 underflows to 0
        *)        MV=0.001 ;;
      esac
      echo "--- viscosity_mode = $mode  (m = $MV)"
      run "V_${mode}" "$MPI_SMALL" "$LH" Outputs/file_base=results/V_${mode} \
        Materials/plasticity/viscosity_mode=$mode Materials/plasticity/m=$MV \
        Materials/plasticity/absolute_tolerance=1e-8 \
        Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true
      python3 "$ANALYZE" "results/V_${mode}.log" | tee "results/V_${mode}.summary"
    done
    ;;

  Vr)   # item 2b: does the viscous term behave like a rate effect?
        # Same final strain (0.12) in every run, reached at different rates,
        # plus a dt refinement at fixed rate. MODE=<name> M=<value> to vary.
        #   - dt refinement at fixed rate MUST give the same stress
        #     (the term is eta*kappa_dot, so it may not depend on dt)
        #   - stress MUST increase monotonically with rate
    VMODE=${MODE:-linear}; VM=${M:-0.001}
    COMMON="Materials/plasticity/viscosity_mode=$VMODE Materials/plasticity/m=$VM
            Materials/plasticity/absolute_tolerance=1e-8
            Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true"
    #          tag          rate   function      end_time  dt
    for spec in "rate1     1      0.001*t       120       0.1" \
                "rate1_dt2 1      0.001*t       120       0.05" \
                "rate10    10     0.01*t        12        0.01" \
                "rate01    0.1    0.0001*t      1200      1.0" \
                "rateindep -      0.001*t       120       0.1"; do
      set -- $spec; tag=$1; fn=$3; et=$4; d=$5
      extra="Materials/plasticity/viscosity_mode=$VMODE"
      [ "$tag" = "rateindep" ] && extra="Materials/plasticity/viscosity_mode=rate_independent"
      run "Vr_${VMODE}_${tag}" "$MPI_SMALL" "$LH" \
        Outputs/file_base=results/Vr_${VMODE}_${tag} \
        "BCs/pull_x/function=$fn" Executioner/end_time=$et Executioner/dt=$d \
        Materials/plasticity/m=$VM Materials/plasticity/absolute_tolerance=1e-8 \
        Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true $extra
    done
    python3 - "$VMODE" <<'PY' | tee "results/Vr_${VMODE}.summary"
import csv, sys
mode = sys.argv[1]
print(f"--- viscosity_mode = {mode}: final stress_xx at strain_xx = 0.12")
base = None
for tag, label in [("rateindep", "rate-independent  "),
                   ("rate01",    "rate x 0.1        "),
                   ("rate1",     "rate x 1          "),
                   ("rate1_dt2", "rate x 1, dt/2    "),
                   ("rate10",    "rate x 10         ")]:
    try:
        r = list(csv.DictReader(open(f"results/Vr_{mode}_{tag}.csv")))[-1]
    except Exception as e:
        print(f"    {label} MISSING ({e})"); continue
    s, e_xx = float(r["stress_xx"]), float(r["strain_xx"])
    if tag == "rateindep":
        base = s
    over = "" if base is None else f"   overstress {100*(s-base)/base:+7.3f}%"
    print(f"    {label} stress_xx = {s:12.6f}   strain_xx = {e_xx:.6f}{over}")
print("    CHECK 1  'rate x 1' and 'rate x 1, dt/2' must agree to ~1e-6 relative")
print("    CHECK 2  overstress must increase monotonically with rate")
PY
    ;;

  P)    # forces the Primal CPPA fallback by starving the stage-2 Newton, so the
        # 7x7 path (normally never entered) is actually exercised. Runs the
        # axis-aligned and the rotated case; both must match the normal runs.
    for cfg in "P_axis   -                                                    -" \
               "P_rot    ICs/euler_phi1_0/value=30 ICs/euler_Phi_0/value=45 ICs/euler_phi2_0/value=60"; do
      set -- $cfg; tag=$1; shift
      rot=""; [ "$1" != "-" ] && rot="$*"
      run "$tag" "$MPI_SMALL" "$LH" Outputs/file_base=results/$tag \
        Materials/plasticity/max_iterations_newton=1 \
        Materials/plasticity/absolute_tolerance=1e-8 \
        Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true $rot
      ref=results/A_full.csv; [ "$tag" = "P_rot" ] && ref=results/B_full.csv
      echo "--- $tag (Primal CPPA forced) vs $ref"
      python3 "$ANALYZE" "results/$tag.log" --ref "$ref" --test "results/$tag.csv" \
        --tol 1e-6 | tee "results/$tag.summary"
    done
    ;;

  Bs)   # B with SMALL strain: removes MOOSE's approximate finite-strain Jacobian
    pair Bs "$MPI_SMALL" "$LH" Outputs/file_base \
      ICs/euler_phi1_0/value=30 ICs/euler_Phi_0/value=45 ICs/euler_phi2_0/value=60 \
      Materials/plasticity/viscosity_mode=rate_independent \
      Physics/SolidMechanics/QuasiStatic/all/strain=SMALL
    ;;

  Bf)   # B with FINITE strain + MOOSE's corotational finite-deformation Jacobian
    pair Bf "$MPI_SMALL" "$LH" Outputs/file_base \
      ICs/euler_phi1_0/value=30 ICs/euler_Phi_0/value=45 ICs/euler_phi2_0/value=60 \
      Materials/plasticity/viscosity_mode=rate_independent \
      Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true
    ;;

  C)    # rotation + shear + softening + viscosity
    # Softening makes the grain split a bifurcation: the seed that selects
    # the localising grain, and WHEN it localises, is amplified from round-off
    # and depends on the Newton path. After ~t=90 neither grain nor global
    # columns are a valid tangent check. Verdict: every column up to t=90.
    TMAX=90 pair C "$MPI_SMALL" "$TWIN" Outputs/file_base
    echo "--- C global columns, full run (INFORMATIONAL: post-bifurcation, seed-dependent)"
    python3 "$ANALYZE" results/C_full.log --ref results/C_ref.csv --test results/C_full.csv \
        --tol "$TOL" --cols global | sed -n '/CSV/,$p' | tee -a results/C.summary
    ;;

  Bt)   # diagnostic: B with a 1e6x tighter local return-map tolerance (FULL only)
    run Bt_full "$MPI_SMALL" "$LH" Outputs/file_base=results/Bt_full \
      ICs/euler_phi1_0/value=30 ICs/euler_Phi_0/value=45 ICs/euler_phi2_0/value=60 \
      Materials/plasticity/viscosity_mode=rate_independent \
      Materials/plasticity/absolute_tolerance=1e-12 \
      Materials/plasticity/max_iterations_newton=200
    python3 "$ANALYZE" results/Bt_full.log --ref-log results/B_full.log \
        --ref results/B_ref.csv --test results/Bt_full.csv --tol "$TOL" | tee results/Bt.summary
    ;;

  D)    # bone cube, production settings, first 3 steps
        # NOREF=1 skips the (expensive) elastic reference
    pair D "$MPI_BONE" "$BONE" Outputs/csv/file_base \
      Outputs/exo_out/file_base=results/D_exo Executioner/end_time=0.3 "${PROD[@]}"
    ;;

  Ddt)  # step-size study: same total strain as D step 1, in 4 increments.
        # Tests whether the line-search damping is caused by the step size.
    run Ddt_full "$MPI_BONE" "$BONE" Outputs/csv/file_base=results/Ddt_full \
      Outputs/exo_out/file_base=results/Ddt_exo \
      Executioner/dt=0.025 Executioner/end_time=0.1 "${PROD[@]}"
    python3 "$ANALYZE" results/Ddt_full.log --ref-log results/D_full.log \
      | tee results/Ddt.summary
    ;;

  Dtol) # tolerance control: same tangent, 100x tighter nonlinear tolerance.
        # If D vs Dtol differs as much as D_full vs D_ref, the difference is
        # solver tolerance, not the tangent.
    run Dtol_full "$MPI_BONE" "$BONE" Outputs/csv/file_base=results/Dtol_full \
      Outputs/exo_out/file_base=results/Dtol_exo \
      Executioner/end_time=0.3 "${PROD[@]}" Executioner/nl_abs_tol=1e-8
    python3 "$ANALYZE" results/Dtol_full.log --ref-log results/D_full.log \
      --ref results/D_full.csv --test results/Dtol_full.csv --tol "$TOL" \
      | tee results/Dtol.summary
    ;;

  *) echo "unknown test '$T' (use: fd A B Bt Bs Bf V Vr P C D Ddt Dtol)";;
  esac
done
