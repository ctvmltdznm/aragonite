# Orthotropic Plasticity Model for Mineralised Tissues

A multiscale finite element framework for modeling coral aragonite biomechanics and trabecular bone, implemented in the [MOOSE](https://mooseframework.inl.gov/) framework.

[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL_v2.1-blue.svg)](https://spdx.org/licenses/LGPL-2.1-only.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21258386.svg)](https://doi.org/10.5281/zenodo.21258386)

## Overview

This framework provides elastic-viscoplastic constitutive models for mineralised tissues (coral aragonite, trabecular and compact bone) and a cohesive zone model for grain-boundary interfaces parameterised from molecular dynamics (MD). Key features:

- **Material-model flag** (`material_model`) that selects elasticity, yield surface and post-yield law together:
  - six bone models identical to the material sets of the reference Abaqus UMAT (Schwiedrzik & Zysset): `trabecular_iso`, `trabecular_ti`, `trabecular_fabric_ti`, `trabecular_fabric_ortho`, `compact_iso`, `compact_ti`
  - `coral`: orthotropic aragonite single crystals with explicit constants
  - `legacy` (default): every constant set by hand; older inputs run unchanged
- **Orthotropic quadric yield surface** (Schwiedrzik et al. 2013) with tension-compression asymmetry and associated flow
- **Six post-yield laws**: perfect plasticity, linear and exponential hardening, and three softening laws
- **Continuous Perzyna viscosity** (five viscosity functions) or rate-independent
- **Two-stage return mapping** (Newton, then primal closest-point projection) with a **consistent tangent**: quadratic convergence of the global Newton iteration, verified up to 2.4 million degrees of freedom
- **Density-fabric scaling** for bone (Zysset-Curnier elasticity, fabric-based strength)
- **Homogenized cohesive zone model (CZM)** for grain boundaries with MD-derived parameters and Gauss-Hermite homogenization of contact quality (Wang et al. 2025)
- **Per-element Euler angles** for polycrystalline representative volume elements (RVEs)

The documentation covers the equations, the code and its verification. Comparison with experiments is not part of this repository.

## Physical Systems

### Coral Aragonite

- **Dense aragonite needles:** orthotropic single crystals / 'dry' grains
- **Interfaces:** protein-water layers between needles and between grains
- **MD data:** interfaces are 22–49× weaker than bulk aragonite (Kvashin et al. 2026)

### Trabecular and Compact Bone

- **Input:** bone volume fraction and, for the fabric models, fabric eigenvalues from μCT
- **Constants:** tissue properties, exponents and interaction parameters come from the selected preset and can be overridden individually

## Installation

### Prerequisites

- C++17 compiler (GCC 9+, Clang 10+)
- [conda](https://docs.conda.io) (Miniconda) — quick setup below
- Python 3.6+ with matplotlib (for validation plotting)

Quick conda setup:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow prompts, answer 'yes' to initialization
source ~/.bashrc
```

### Build

This is a MOOSE application: it builds against a local clone of the **MOOSE source tree**, which must sit in the *same parent folder* as this repository - the app's Makefile looks for `../moose/framework/build.mk`. So clone both repos side by side:

```bash
cd ~   # or any parent directory

# 1. MOOSE source tree (provides framework/build.mk)
git clone https://github.com/idaholab/moose.git

# 2. This application, as a sibling of moose/
git clone git clone https://github.com/ctvmltdznm/aragonite.git

# 3. MOOSE build environment (conda)
conda config --add channels https://conda.software.inl.gov/public
conda create -n moose moose-dev -c conda-forge
conda activate moose

# 4. Build the app
cd tuc-fe-moose
make -j8
```

Resulting layout:
```
~/
├── moose/          # MOOSE framework source
└── tuc-fe-moose/   # this application
```

If MOOSE is cloned elsewhere, set `MOOSE_DIR=/path/to/moose` before `make`.

**MOOSE version.** The app builds against current MOOSE and, in principle, tracks
the latest framework. For an exact, reproducible build it was tested with:

> MOOSE (`moose-dev`) build **2026-03-28**, commit `73af6c9`.

To reproduce that exact build, pin the MOOSE checkout before step 4:
```bash
cd ~/moose && git checkout 73af6c9 && cd -
```

**Alternative conda setup (`environment.yml`).** The repo also ships an
`environment.yml` declaring the same dependencies:
```bash
conda env create -f environment.yml
conda activate moose
```
If this errors with `database is locked` (a conda repodata-cache issue on some
conda versions), use the `conda create` command in step 3 instead.

## Quick Start

After building, run the three test sets. Together they take a few minutes on one core (the RVE tests of the gold-file suite about ten minutes more).

```bash
# 1. Gold-file suite: 25 tests, detects unintended changes
cd scripts
chmod +x *.sh *.py              # first time only
./run_validation.sh

# 2. Example suite: 57 inputs, every material model and three interface
#    modes, each compared with a value predicted from the parameters
cd examples
./run_examples.sh               # ends with "VERDICT: PASS"

# 3. Offline preset check against the UMAT formulas (no MOOSE needed)
python3 verify_presets.py
```

### Gold-File Suite (`scripts/run_validation.sh`)

| Category | Tests | Purpose |
|----------|-------|---------|
| **ortho** | 6 | Elastic constants and yield stresses in six loading directions |
| **asymmetry** | 2 | Tension-compression asymmetry |
| **postyield** | 5 | Post-yield evolution modes |
| **grains** | 4 | Multi-grain orientation effects |
| **czm** | 4 | Interface failure (Mode I, II, mixed, complete separation) |
| **fabric** | 2 | Fabric-based vs equivalent explicit tensor |
| **rve** | 2 | Polycrystalline RVE without and with cohesive interfaces |

```bash
./run_validation.sh -c ortho    # one category
./run_validation.sh -p 8        # 8 MPI ranks
./run_validation.sh --list      # list tests
./plot_validation.sh            # figures in test/validation/figures/
```

A test passes if it shows `[PASS] Matches gold file`. Segfaults during cleanup are harmless. Gold files detect change, not error; after modifying the material classes use the tests below first, then regenerate gold files with `--generate-gold`.

The `rve` tests use short runs (elastic regime only). For a full deformation run, increase `end_time` in `test/validation/rve/grain_interfaces/grain_interfaces.i`.

### Example Suite (`scripts/examples/`)

```bash
./run_examples.sh bone          # the 54 bone inputs (6 models x 9 load cases)
./run_examples.sh coral         # the 3 cohesive-interface inputs
./run_examples.sh compact_ti    # substring filter
```

The inputs are also the maintained starting points for your own simulations: one small input per model and load case.

### Tangent and Solver Verification (`test/verification/run_tests.sh`)

Developer tests. Every test is run with the consistent tangent and with an elastic-tangent reference of the same binary; the results must agree and the iteration counts show the quality of the tangent.

```bash
cd test/verification
./run_tests.sh fd               # derivatives against finite differences
./run_tests.sh A B Bs Bf C      # single element, rotated, two grains
./run_tests.sh V P M            # viscosity modes, primal stage, presets
NOREF=1 ./run_tests.sh D        # bone cube, 2.4M DOFs (32 MPI ranks)
```

## Usage Examples

### Example 1: Bone Model from a Preset

```
[GlobalParams]
  displacements  = 'disp_x disp_y disp_z'
  material_model = trabecular_fabric_ortho
  density_rho    = 0.2          # BV/TV
  fabric_m1      = 0.85         # sum of the three = 3
  fabric_m2      = 0.95
  fabric_m3      = 1.20
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    absolute_tolerance    = 1e-8
    max_iterations_newton = 50
  []
[]
```

Constants, post-yield law and viscosity come from the preset. In the bone models the preset strengths are ultimate strengths; first yield occurs at 0.7 times that value.

### Example 2: Coral Grains (Aragonite)

```
[GlobalParams]
  displacements  = 'disp_x disp_y disp_z'
  material_model = coral
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor     # C_ijkl defaults to aragonite
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate # strengths default to aragonite
    zeta12 = 0.27     # required; provisional value
    zeta13 = 0.27
    zeta23 = 0.27
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
    absolute_tolerance    = 1e-8
    max_iterations_newton = 50
  []
[]
```

Euler angles follow the MOOSE convention, which is the inverse of the Bunge convention of EBSD software: enter Bunge angles (φ₁, Φ, φ₂) as (180° − φ₂, Φ, 180° − φ₁).

### Example 3: Explicit Constants (Legacy Mode)

```
[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    sigma_xx_tension = 4980    # MPa
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    zeta12 = <set>             # initialised to 0 if omitted; set all three
    zeta13 = <set>
    zeta23 = <set>
    postyield_mode    = exp_softening
    residual_strength = 0.7
    kmax   = 0.001
    kslope = 30
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
    absolute_tolerance    = 1e-8
    max_iterations_newton = 50
  []
[]
```

Porous material: add `density_rho` with `elasticity_density_exponent` (elasticity) and `yield_density_exponent` (plasticity).

### Example 4: Grain Boundary Interfaces (CZM)

```
[Materials]
  [czm_interface]
    type = HomogenizedExponentialCZM
    boundary = 'grain_boundaries'

    # MD-derived (Kvashin et al. 2026); provisional parameter set
    normal_strength  = 626.0       # MPa
    shear_strength_s = 374.0       # MPa
    shear_strength_t = 374.0       # MPa
    delta_0_normal   = 1.91e-4     # um (0.191 nm)
    delta_0_tangent  = 2.17e-4     # um (0.217 nm)

    mu  = 0.92                     # loading exponent
    eta = 0.27                     # softening exponent (regularisation)

    quality_std_dev         = 0.10   # within a quadrature point
    spatial_quality_std_dev = 0.15   # between quadrature points
    spatial_random_seed     = 1234

    damage_viscosity        = 2.5    # simulation time units; keep dt <= this
  []
[]
```

**Use micrometres as the mesh length unit with this class.** It contains a regularisation length of 1e-6 mesh units; with the openings above entered in millimetres the law is not resolved.

The homogenization is deterministic (fixed quadrature, no random sampling) and independent of the element area.

### Solver Settings

Measured on a 2.4 million DOF bone specimen; use as the starting point:

```
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    use_finite_deform_jacobian = true   # required with FINITE strain
  []
[]

[Executioner]
  type = Transient
  solve_type  = NEWTON
  line_search = bt                      # required
  nl_rel_tol  = 1e-8
  nl_abs_tol  = 1e-6
  petsc_options_iname = '-pc_type -pc_gamg_agg_nsmooths -ksp_type -ksp_gmres_restart -ksp_max_it -ksp_rtol'
  petsc_options_value = 'gamg 2 gmres 300 200 1e-4'
  dt = 0.1
[]
```

For single-element problems `-pc_type lu` is sufficient.

## Documentation

Complete technical documentation: [`doc/documentation.pdf`](doc/documentation.pdf)

- **Chapter 1:** Introduction, material models, application domains
- **Chapter 2:** Mathematical formulation (kinematics, elasticity, yield surface, post-yield laws, viscosity, return mapping, consistent tangent, orientation)
- **Chapter 3:** Cohesive zone model
- **Chapters 4 & 5:** MOOSE installation and building
- **Chapter 6:** Material models, preset constants, parameter reference, input templates, solver settings
- **Chapter 7:** Verification (automated tests, earlier studies, example inputs)
- **Chapter 8:** Troubleshooting

## Material Parameters

### Aragonite (Dense Crystal, `coral` Defaults)

| Property | Value |
|----------|-------|
| E₁, E₂, E₃ | 140.4, 70.3, 63.4 GPa |
| G₁₂, G₁₃, G₂₃ | 46.6, 31.1, 42.1 GPa |
| σ₁, σ₂, σ₃ (tensile strength) | 4980, 4100, 5340 MPa |
| τ₁₂, τ₁₃, τ₂₃ (shear strength) | 4510, 5080, 5060 MPa |

### Protein-Aragonite Interface (CZM)

| Property | Protein | Water (001/100) |
|----------|---------|-----------------|
| Normal strength σₙ | 626 MPa | 742 MPa |
| Shear strength τₛ | 374 MPa | 354 MPa |
| Characteristic opening δ₀,ₙ | 0.19 nm | 0.224 nm |
| Characteristic opening δ₀,ₜ | 0.22 nm | 0.182 nm |
| Loading exponent μ | 0.916 | 1.224 |

Water-010 face excluded (non-equivalent electrostatics). δ₀ is the separation beyond the equilibrium interface thickness (about 1 nm). The set is provisional.

### Bone

The constants of the six bone models are tabulated in Chapter 6 of the documentation and defined in `include/materials/MaterialModelPresets.h`.

**References:** [Schwiedrzik et al. (2013)](https://doi.org/10.1007/s10237-013-0472-5), [Zysset & Curnier (1995)](https://doi.org/10.1016/0167-6636(95)00018-6), [Wolfram et al. (2012)](https://doi.org/10.1016/j.jmbbm.2012.07.005), [Rincón-Kohli & Zysset (2009)](https://doi.org/10.1007/s10237-008-0128-z), [Wang et al. (2025)](https://doi.org/10.1016/j.conbuildmat.2025.142454), [Kvashin et al. (2026)](https://doi.org/10.1016/j.jmbbm.2026.107403).

## Implementation Details

### Material Classes

| Class | Purpose |
|-------|---------|
| `OrthotropicPlasticityStressUpdate` | Quadric plasticity, return mapping, consistent tangent |
| `ComputeFabricElasticityTensor` | Density- and fabric-based elasticity; model flags |
| `ComputeElasticityTensorCoupled` | Explicit orthotropic elasticity with per-element Euler angles |
| `HomogenizedExponentialCZM` | Grain-boundary failure with Gauss-Hermite homogenization |

Header-only: `MaterialModelPresets.h` (preset constants), `MandelHelpers.h` (tensor ↔ six-component conversion).

### Verification Status

Verified with the current code:

- ✅ Derivatives of the yield function and post-yield laws against finite differences
- ✅ Preset constants of all models against the UMAT formulas (64 flag combinations)
- ✅ Elastic modulus, yield onset and peak stress of every bone model (single element)
- ✅ Consistent tangent against an elastic-tangent reference, with rotation, softening and viscosity
- ✅ Both stages of the return mapping against each other
- ✅ Quadratic global convergence on a bone cube with 2.4 million DOFs
- ✅ Viscosity: rate dependence and step-size independence
- ✅ CZM peak traction, mode mix and damage on a single interface

Not verified: bulk damage coupling; a bone preset in a multi-element run; the interface model under non-proportional loading. Not attempted here: comparison with experiments or with a run of the reference UMAT.

## Project Structure

```
tuc-fe-moose/
├── include/
│   ├── base/                # Application base
│   └── materials/           # Material headers (.h)
├── src/
│   ├── base/                # Application implementation
│   ├── main.C               # Main entry point
│   └── materials/           # Material implementation (.C)
├── test/
│   ├── validation/          # Gold-file suite (run_validation.sh)
│   │   ├── ortho/           # 6 tests
│   │   ├── asymmetry/       # 2 tests
│   │   ├── postyield/       # 5 tests
│   │   ├── grains/          # 4 tests
│   │   ├── czm/             # 4 tests
│   │   ├── fabric/          # 2 tests
│   │   ├── rve/             # 2 tests
│   │   └── figures/         # Generated plots
│   └── verification/        # Tangent and solver tests
│       ├── run_tests.sh
│       ├── analyse_run.py
│       └── *.i              # linear_hardening, twin, bone cube, model flags
├── scripts/
│   ├── run_validation.sh
│   ├── plot_validation.sh
│   ├── plot_validation_category.py
│   ├── compare_with_gold.py
│   ├── compute_fabric_reference.py
│   └── examples/            # Example suite
│       ├── run_examples.sh
│       ├── analyse_examples.py
│       ├── generate_examples.py
│       ├── verify_presets.py
│       ├── bone/            # 6 models x 9 load cases
│       └── coral/           # 3 interface load cases
├── doc/
│   └── documentation.pdf
├── environment.yml
├── CITATION.cff
├── codemeta.json
├── Makefile
├── LICENSE
└── README.md
```

## Citation

If you use this software in your research, please cite it via its Zenodo record
(metadata in [`CITATION.cff`](CITATION.cff)):

> Kvashin, N., & Wolfram, U. (2026). *TUC-FE-MOOSE: a MOOSE application for hierarchical biomineral mechanics* (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.21258387

```bibtex
@software{kvashin_tuc_fe_moose_2026,
  author    = {Kvashin, Nikolai and Wolfram, Uwe},
  title     = {TUC-FE-MOOSE: a MOOSE application for hierarchical biomineral mechanics},
  year      = {2026},
  version   = {1.0.1},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21258387},
  url       = {https://doi.org/10.5281/zenodo.21258387}
}
```

Please also cite the MD study the interface parameters derive from:
[Kvashin et al. (2026)](https://doi.org/10.1016/j.jmbbm.2026.107403).

## License

This project is licensed under the GNU Lesser General Public License v2.1 (LGPL-2.1-only) — see the [LICENSE](LICENSE) file for details.

---

**For detailed usage, see the [documentation](doc/documentation.pdf).**
