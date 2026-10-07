# Examples: material-model flags

Self-contained MOOSE inputs, one per test, each runnable with no arguments.

    aragonite-opt -i examples/bone/trabecular_iso/tension_xx.i
    aragonite-opt -i examples/coral/coral_czm_mode_I_opening.i

Reference CSVs are captured by the runner:
 
    GOLD=1 ./run_examples.sh
 
## Layout
 
    examples/
      run_examples.sh          run them all, sequentially
      analyse_examples.py      compare results against the preset predictions
      selftest_analyser.py     check the analyser itself, no MOOSE needed
      generate_examples.py     rebuild every input below
      expected_peaks.py        predicted yield points, standalone
      verify_presets.py        preset formulas, shared by the scripts above
      bone/
        trabecular_iso/           tension_{xx,yy,zz}.i
        trabecular_ti/            compression_{xx,yy,zz}.i
        trabecular_fabric_ti/     shear_{xy,xz,yz}.i
        trabecular_fabric_ortho/
        compact_iso/
        compact_ti/
      coral/
        coral_czm_mode_I_opening.i
        coral_czm_mode_II_shear.i
        coral_czm_mixed_mode.i
 
54 bone inputs (6 flags x 9 load cases) and 3 coral inputs.
 
## One flag, both responses
 
Each input sets the model once, in `[GlobalParams]`:
 
    [GlobalParams]
      material_model = trabecular_fabric_ortho
      density_rho    = 0.2
      fabric_m1 = 0.85
      fabric_m2 = 0.95
      fabric_m3 = 1.20
    []
 
`ComputeFabricElasticityTensor` and `OrthotropicPlasticityStressUpdate` both
read `material_model`, so the elastic and plastic responses cannot drift apart.
Setting `elastic_model` or `plastic_model` inside one block overrides it there
and prints a warning; mixing is still supported, it is just no longer the
accident it used to be.
 
## The six bone flags
 
IDs are the UMAT `PROPS(1)` values from `UMAT_QUADRIC_PRIMAL_Major.f`.
 
| ID | flag | elastic | yield surface | post-yield | rho used here |
|---|---|---|---|---|---|
| 0 | `trabecular_iso` | isotropic | isotropic | `simple_softening` | 0.20 |
| 1 | `trabecular_ti` | transversely isotropic | TI | `simple_softening` | 0.20 |
| 2 | `trabecular_fabric_ti` | fabric, TI | fabric, TI | `simple_softening` | 0.20 |
| 3 | `trabecular_fabric_ortho` | fabric, orthotropic | fabric, orthotropic | `simple_softening` | 0.20 |
| 4 | `compact_iso` | isotropic | isotropic | `exp_hardening` | 0.90 |
| 5 | `compact_ti` | transversely isotropic | TI | `exp_hardening` | 0.90 |
 
The trabecular presets run at BV/TV = 0.20 and the compact ones at 0.90.
A compact preset at BV/TV = 0.20 is not a meaningful material. The TI and
fabric-TI presets use `main_direction = 3`, so `zz` is axial and `xx`/`yy` are
transverse.
 
## The nine load cases
 
**`tension_ii` and `compression_ii`** put symmetry BCs on the three minus faces
and drive one plus face. The other two plus faces stay traction free, so the
state is **uniaxial stress** and the peak maps straight onto the directional
strength.
 
**`shear_ij`** uses affine (Taylor) BCs on every face: `u_i = gamma * x_j * t`,
with the other two components held at zero everywhere. That prescribes the whole
deformation gradient, so the field is **homogeneous** and every quadrature point
sees the same state. Homogeneity is what makes the element-average CSV columns
checkable at all: `phi(sigma) = r(kappa)` holds pointwise, and with a uniform
field it also holds for the averages, to about 1e-5.
 
**The shear state is not pure shear, and cannot be made so on one element.** The
quadric surface carries a linear term, so it is pressure sensitive, and the
associated flow for a pure shear stress has normal components (`N_11`, `N_22`,
`N_33` run from 0.22 to 0.62 of the flow direction depending on the preset).
Affine BCs forbid that plastic dilatation, so a normal reaction stress develops,
up to about -40 to -44 MPa for he compact presets and about -1 MPa for the trabecular ones at 2% shear.
The analyser therefore reports the **state purity** at the peak and treats the
peak comparison as informational, rather than asserting a relation that does not
hold.
 
Freeing the lateral directions to allow the dilatation does not work, and this
was tried: prescribing only the driven component affinely and holding the other
two on one face each. Those constraints do not only hold the normal strains,
they also fix `d(u_y)/d(x)` and the whole `u_z` field, so without them the
element escapes into other shear modes. A `shear_xy` run came back carrying an
`xz` stress 18% of the driven component, the field stopped being homogeneous,
and the consistency check collapsed from 1e-5 to 0.4 -- not because the model
was wrong, but because `phi(average sigma) != r(average kappa)` once the field
varies within the element. Homogeneity is worth more here than mode purity.
Getting both would need traction-free normal faces plus multi-point constraints
tying `d(u_y)/d(x)`, which is right for a real RVE and far too much machinery
for a single-element example.
 
All cases reach 2% nominal strain (or engineering shear) at t = 1, in 100 steps.
 
## Yields land at 0.7 x strength, on purpose
 
The UMAT's strengths are **ultimate** strengths and `RADK` starts at
`RDY = 0.7`, so every bone preset carries `initial_yield_ratio = 0.7` and first
yield occurs at 0.7 x strength. This is the single most common surprise when
reading the output. See `MATERIAL_MODEL_THEORY.md` section 4.
 
## What to expect
 
Computed from the presets by `expected_peaks.py`, independently of MOOSE.
"margin" is how far past first yield each run goes by t = 1.
 
```
model                    case              yield stress  yield strain  margin @2%
----------------------------------------------------------------------------------
trabecular_iso           tension_xx               2.851       0.00460        4.3x
trabecular_iso           compression_xx          -4.142      -0.00669        3.0x
trabecular_iso           tension_yy               2.851       0.00460        4.3x
trabecular_iso           compression_yy          -4.142      -0.00669        3.0x
trabecular_iso           tension_zz               2.851       0.00460        4.3x
trabecular_iso           compression_zz          -4.142      -0.00669        3.0x
trabecular_iso           shear_xy                 2.191       0.00882        2.3x
trabecular_iso           shear_xz                 2.191       0.00882        2.3x
trabecular_iso           shear_yz                 2.191       0.00882        2.3x
trabecular_ti            tension_xx               2.335       0.00490        4.1x
trabecular_ti            compression_xx          -3.139      -0.00659        3.0x
trabecular_ti            tension_yy               2.335       0.00490        4.1x
trabecular_ti            compression_yy          -3.139      -0.00659        3.0x
trabecular_ti            tension_zz               4.551       0.00336        6.0x
trabecular_ti            compression_zz          -7.646      -0.00565        3.5x
trabecular_ti            shear_xy                 1.561       0.00868        2.3x
trabecular_ti            shear_xz                 2.380       0.00877        2.3x
trabecular_ti            shear_yz                 2.380       0.00877        2.3x
trabecular_fabric_ti     tension_xx               2.456       0.00420        4.8x
trabecular_fabric_ti     compression_xx          -3.678      -0.00629        3.2x
trabecular_fabric_ti     tension_yy               2.456       0.00420        4.8x
trabecular_fabric_ti     compression_yy          -3.678      -0.00629        3.2x
trabecular_fabric_ti     tension_zz               4.493       0.00408        4.9x
trabecular_fabric_ti     compression_zz          -6.730      -0.00611        3.3x
trabecular_fabric_ti     shear_xy                 1.887       0.00793        2.5x
trabecular_fabric_ti     shear_xz                 2.108       0.00781        2.6x
trabecular_fabric_ti     shear_yz                 2.108       0.00781        2.6x
trabecular_fabric_ortho  tension_xx               2.164       0.00420        4.8x
trabecular_fabric_ortho  compression_xx          -3.241      -0.00629        3.2x
trabecular_fabric_ortho  tension_yy               2.733       0.00415        4.8x
trabecular_fabric_ortho  compression_yy          -4.094      -0.00622        3.2x
trabecular_fabric_ortho  tension_zz               4.464       0.00406        4.9x
trabecular_fabric_ortho  compression_zz          -6.687      -0.00608        3.3x
trabecular_fabric_ortho  shear_xy                 1.543       0.00788        2.5x
trabecular_fabric_ortho  shear_xz                 1.972       0.00779        2.6x
trabecular_fabric_ortho  shear_yz                 2.217       0.00774        2.6x
compact_iso              tension_xx              84.769       0.00521        3.8x
compact_iso              compression_xx        -137.200      -0.00843        2.4x
compact_iso              tension_yy              84.769       0.00521        3.8x
compact_iso              compression_yy        -137.200      -0.00843        2.4x
compact_iso              tension_zz              84.769       0.00521        3.8x
compact_iso              compression_zz        -137.200      -0.00843        2.4x
compact_iso              shear_xy                60.704       0.01002        2.0x
compact_iso              shear_xz                60.704       0.01002        2.0x
compact_iso              shear_yz                60.704       0.01002        2.0x
compact_ti               tension_xx              33.137       0.00261        7.7x
compact_ti               compression_xx        -117.965      -0.00929        2.2x
compact_ti               tension_yy              33.137       0.00261        7.7x
compact_ti               compression_yy        -117.965      -0.00929        2.2x
compact_ti               tension_zz             103.383       0.00499        4.0x
compact_ti               compression_zz        -157.067      -0.00759        2.6x
compact_ti               shear_xy                36.451       0.00839        2.4x
compact_ti               shear_xz                48.380       0.00873        2.3x
compact_ti               shear_yz                48.380       0.00873        2.3x
----------------------------------------------------------------------------------
smallest margin over all 54 cases: 2.0x  OK, every case yields well before t=1
```
 
## The coral examples
 
Two HEX8 grains with different crystal orientations, bonded by a
`HomogenizedExponentialCZM` interface. The smallest complete version of the
aragonite RVE setup. Three files, identical except for the boundary conditions,
so they diff cleanly against each other.
 
| file | loading | predicted `T_peak` | predicted `delta_0_eff` |
|---|---|---|---|
| `coral_czm_mode_I_opening.i` | pull along x (the interface normal) | 626.0 MPa | 1.91e-4 µm |
| `coral_czm_mode_II_shear.i` | slide along y, x held | 374.0 MPa | 2.17e-4 µm |
| `coral_czm_mixed_mode.i` | 45 deg, equal opening and sliding | 515.6 MPa |2.03e-4 µm |
 
The interface normal is x, so jump component 0 is opening and components 1 and 2
are sliding. The model mixes the modes **twice over**:
 
- **Peak traction**, by an elliptic interaction on the direction cosines of the
  jump vector: `T_peak = sqrt((sigma_n*rn)^2 + (tau_s*rs)^2 + (tau_t*rt)^2)`.
- **Characteristic opening**, by the Wang 2025 mixed-mode form, running from
  `delta_0_normal` at pure opening to `delta_0_tangent` at pure sliding.
The table above is what those two formulas give for the shipped parameters, and
it is what each run should be checked against. The mixed case is the informative
one: its peak must land between the two pure modes, and `normal_jump` and
`tangent_jump` must stay equal to each other, which is the check that the 45 deg
path held.
 
All three write the same postprocessor set, so mode purity is visible directly:
in mode I the tangential columns stay at zero, in mode II the normal ones do, and
in the mixed case neither does.
 
**Boundary conditions differ by mode, deliberately.** Mode I keeps symmetry
conditions so the bulk can contract freely. Modes II and mixed prescribe all
three components on both end faces, because any unconstrained normal motion
would let the interface open and contaminate the mode. The bulk deforms by about
1e-6 µm before the interface reaches its peak, three orders below the applied
displacement, so in every case the jump is essentially the applied displacement.
 
**`delta_0` and the length unit.** Lengths are in micrometres. `delta_0_normal`
= 1.91e-4 and `delta_0_tangent` = 2.17e-4 are the MD values, 0.19 nm and
0.22 nm, used unchanged. They are separations beyond the equilibrium interface
(the organic layer itself is about 1 nm thick); the zero-thickness cohesive
element represents only that extra separation. The values are provisional.

Do not enter these openings in millimetres. The class smooths the opening
with a fixed width of 1e-6 in mesh units, which must stay far below `delta_0`;
in millimetres `delta_0` would be 1.9e-7 and the interface would start out
damaged. Keep `delta_0` at or above 1e-4 in mesh units.

The element edge here is 1e-4 µm. That is not a grain size: these are
single-interface tests, and the bulk elements are only there to carry the
load to the interface. What governs the numerics is the displacement
increment per step relative to `delta_0`, which these inputs keep small. The
peak traction is physical; the post-peak branch depends on `eta`, which is a
regularisation choice.
 
`zeta12/13/23` have no default and must be set: zero is not an acceptable
placeholder and no calibrated coral value exists yet. The 0.27 in the examples is
provisional, carried over from the bone value.
 
Hexahedral and tetrahedral meshes both work with the cohesive interface. One
failure has been seen: tetrahedra made by MOOSE's hex-to-tet conversion left
the interface tractions at the regularisation floor. Mesh with tetrahedra
directly instead of converting.
 
## Running them
 
    ./run_examples.sh                 # all 57, sequentially
    ./run_examples.sh bone            # only the bone examples
    ./run_examples.sh coral           # only the coral examples
    ./run_examples.sh compact_ti      # substring filter on the path
    ./run_examples.sh shear           # every shear case, bone and coral
 
    APP=/path/to/aragonite-opt ./run_examples.sh
    GOLD=1 ./run_examples.sh          # also copy each CSV into <dir>/gold/
    NOANALYSE=1 ./run_examples.sh     # skip the analysis step
 
Each run writes `<case>_out.csv` and `<case>.log` next to its input, and the
exit code is reported inline; a failure prints the first MOOSE error line
immediately rather than making you open the log. Sequential on purpose: these
are single- and two-element runs of a few seconds each, and a serial log reads
better than interleaved output. When the runs finish, `analyse_examples.py`
runs automatically.
 
## Checking the results
 
    python3 analyse_examples.py                 # everything
    python3 analyse_examples.py --filter coral
    python3 analyse_examples.py --gold          # also diff against gold/
 
The analyser compares each CSV against values **predicted from the presets**,
not against a previous run, so it is meaningful the very first time and does not
depend on a gold file existing. Three checks:
 
**1. Elastic slope.** The pre-yield slope of stress against strain, against
E_i = 1/(C^-1)_ii for the normal cases and G_ij for the shear cases. Uses no
plasticity at all, so it isolates the elasticity side.
 
**2. Yield-surface consistency — the main test.** At every row where the
material is plastic, the stress must sit on the current yield surface:
 
    phi(sigma) = sqrt(s:F:s) + f_lin.s  ==  r(kappa)
 
Both sides come out of the CSV: the six stress components and `plastic_strain`,
which is kappa. F, f_lin and r(.) come from the preset formulas. Nothing here
depends on how far the run went, on the time step, or on a stored reference,
which is what makes it the strongest of the three. A violation means the return
map converged off the surface, or the preset resolved to the wrong constants.
 
**3. First yield and peak, reported rather than asserted.** First yield is where
kappa leaves zero, linearly interpolated between the bracketing rows, against
r0 x strength. The raw first plastic row overshoots by up to one elastic
increment, about 4% of the yield stress at these settings, so even the
interpolated value is a discretisation estimate. The peak is compared with
r(kappa) evaluated at the row where the peak occurs, because the peak is
wherever r happens to be largest within the strain range, not a material
constant.
 
Note that **peak is not strength** for the bone presets. UMAT strengths are
ultimate strengths, r(0) = 0.7, and these runs stop at 2% strain, which leaves
r near 0.93 to 0.94 rather than 1.
 
For the coral cases it checks the peak traction against the elliptic-interaction
prediction, mode purity (the off-mode traction as a fraction of the loaded one),
that damage never decreases, and for the mixed case that the two jumps stay
equal, which is the check that the 45 deg path held.
 
**Tolerances.** `--tol-elastic 2e-2`, `--tol-yield 2e-3`, `--tol-peak 3e-2`.
The yield tolerance was tightened from 1e-2 after the first clean baseline:
the worst observed `|phi/r - 1|` is 4.5e-4. The floor on that check is set by
FINITE strain: the CSV reports Cauchy stress, and at 2% strain that differs
from the measure the return map uses in the fourth digit.
 
### Is the analyser itself right?
 
    python3 selftest_analyser.py
 
It builds synthetic CSVs that satisfy the model exactly by construction, so a
correct analyser must report round-off, then perturbs one case by 5% and
confirms it is caught. Current result: worst `|phi/r - 1|` 4.4e-16, worst
elastic-slope error 3.6e-16, injected error caught. This matters because the
analyser hard-codes the stress 6-vector order and the slot mapping, and getting
either wrong would produce confident nonsense rather than an obvious crash.
 
## Regenerating
 
    python3 generate_examples.py --out examples
    python3 expected_peaks.py            # needs verify_presets.py alongside
 
Regenerate after any change to the flags or the presets rather than editing the
54 inputs by hand.
 
## Self-checks
 
Every example has `debug_checks = true` on the plasticity block. At startup it
finite-difference tests the yield gradient, the Hessian and the softening
derivative at three stress states with shear, and round-trips the elasticity
tensor through the Mandel helpers. All of them abort the run on failure, so a
clean start is itself a result. Look for the `FD check ... OK` and
`Mandel round-trip check ... OK` lines. Turn it off for production runs.
 
