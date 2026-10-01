# Examples: material-model flags

Self-contained MOOSE inputs, one per test, each runnable with no arguments.

    aragonite-opt -i examples/bone/trabecular_iso/tension_xx.i
    aragonite-opt -i examples/coral/coral_czm_two_element.i

Gold files:

    aragonite-opt -i examples/bone/trabecular_iso/tension_xx.i --generate-gold

## Layout

    examples/
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

**`shear_ij`** uses affine (Taylor) BCs on every face: `u_i = gamma * x_j * t`.
In the material frame both the orthotropic stiffness and the quadric yield
surface are block diagonal, so an affine simple shear produces a **pure shear
stress** and the peak maps onto the corresponding tau. At finite gamma there is
a second-order normal stress; at the 2% used here it is about 1e-4 of the shear
stress.

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
| `coral_czm_mode_I_opening.i` | pull along x (the interface normal) | 626.0 MPa | 1.91e-4 mm |
| `coral_czm_mode_II_shear.i` | slide along y, x held | 374.0 MPa | 2.17e-4 mm |
| `coral_czm_mixed_mode.i` | 45 deg, equal opening and sliding | 515.6 MPa | 2.03e-4 mm |

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
1e-6 mm before the interface reaches its peak, three orders below the applied
displacement, so in every case the jump is essentially the applied displacement.

**`delta_0` and the mesh scale.** This is the parameter that confuses people.
`delta_0` is a characteristic opening, and what governs the numerics is its
ratio to the element size. The MD values are about 0.19 nm, three orders below
any mesh that resolves a coral grain; used directly they put the interface in
softening at the first increment and Newton stalls. Production runs regularise
`delta_0` up to roughly the element size. Here the element edge is 1e-4 mm
(0.1 um) and `delta_0` is near 2e-4 mm, a ratio near 2. The peak traction, which
sets failure initiation, is physical. The fracture energy is not.

`zeta12/13/23` have no default and must be set: zero is not an acceptable
placeholder and no calibrated coral value exists yet. The 0.27 in the examples is
provisional, carried over from the bone value.

HEX8 is required for the cohesive interface. MOOSE hex-to-tet splitting leaves
interface tractions frozen at the regularisation floor.

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
