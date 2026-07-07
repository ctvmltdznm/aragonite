# Orthotropic Plasticity Model for Mineralised Tissues

A multiscale finite element framework for modeling coral aragonite biomechanics and trabecular bone, implemented in the [MOOSE](https://mooseframework.inl.gov/) framework.

[![License: LGPL](https://img.shields.io/badge/License-LGPL-blue.svg)](https://opensource.org/licenses/LGPL-2.1)
<!-- Add after Zenodo upload: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX) -->

Fork to create a new MOOSE-based application.

For more information see: [https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app](https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app)

## Overview

This framework bridges molecular dynamics (MD) simulations to continuum finite element analysis (FEA) for modelling porous biominerals and biological composites. Key features:

- **Three modeling paradigms:**
  - Explicit orthotropic (dense crystalline materials: aragonite, minerals)
  - Porosity-scaled orthotropic (porous materials with known full-density properties)
  - Fabric tensor-based (trabecular bone, foam-like structures with μCT-derived anisotropy)

- **Orthotropic plasticity** with quadric yield surface (Schwiedrzik formulation)
- **Tension-compression asymmetry** via linear term in yield function
- **Five post-yield evolution modes**: perfect plasticity, exponential/linear hardening, exponential/piecewise softening
- **Density-fabric scaling** for trabecular bone (Zysset-Curnier elasticity, Schwiedrzik plasticity)
- **Homogenized cohesive zone model (CZM)** for grain boundary interfaces with MD-derived parameters and Gauss-Hermite probabilistic homogenization (Wang et el. 2025)
- **Euler angle-based orientation** for polycrystalline representative volume elements (RVEs)

## Physical Systems

### Primary: Coral Aragonite Biomechanics

Hierarchical coral biomineral structure:
- **Dense aragonite needles:** Orthotropic single crystals / 'dry' grains
- **Porous sclerodermites:** Aragonite aggregates with varying porosity
- **Interfaces:** Protein-water layers between needles (CZ1) and between grains (CZ2)
- **MD data:** Interfaces are 22–49× weaker than bulk aragonite (Kvashin et al. 2026)

### Secondary: Trabecular Bone Mechanics

Skeletal tissue modeling:
- **Base tissue properties:** E₀ = 10–20 GPa, ν₀ = 0.3
- **Fabric tensor:** From μCT analysis (structural anisotropy)
- **Apparent properties:** Density-fabric scaled constants at continuum level
- **Applications:** Vertebral bodies, femoral heads, bone quality assessment

## Installation

### Prerequisites

- [MOOSE Framework](https://mooseframework.inl.gov/getting_started/installation/index.html) (conda install recommended)
- C++17 compiler (GCC 9+, Clang 10+)
- Python 3.6+ with matplotlib (for validation plotting)

### Build

```bash
# Activate MOOSE conda environment
conda activate moose

# Clone this repository
git clone https://github.com/ctvmltdznm/aragonite.git
cd aragonite

# Build
make -j8

# Verify
./aragonite-opt --version
```

## Quick Start

### Single-Element Validation Test

```bash
# Run X-direction tension test
cd test/validation/ortho/x_tension
../../../../aragonite-opt -i x_tension.i

# View results
paraview x_tension_out.e
```

### Complete Validation Suite

**25 tests across 7 categories:**

| Category | Tests | Purpose |
|----------|-------|---------|
| **ortho** | 6 | Elastic constants (E₁, E₂, E₃, G₁₂, G₁₃, G₂₃) and yield stresses |
| **asymmetry** | 2 | Tension-compression asymmetry (σ_c ≠ σ_t) |
| **postyield** | 5 | Five r(κ) evolution modes |
| **grains** | 4 | Multi-grain orientation effects |
| **czm** | 4 | Interface failure (Mode I, II, mixed, cycling) |
| **fabric** | 2 | `fabric_plasticity` (Zysset-Curnier) vs `explicit_plasticity` (equivalent `symmetric9`): verifies both give identical elastic slope and yield stress |
| **rve** | 2 | Polycrystalline RVE with GH cohesive interfaces (`grain_interfaces`); run + gold-compared |

```bash
cd scripts

# Make scripts executable (first time only)
chmod +x *.sh *.py

# Run all 25 validation tests
./run_validation.sh

# Run specific category
./run_validation.sh -c ortho
./run_validation.sh -c czm
./run_validation.sh -c fabric

# Parallel execution
./run_validation.sh -p 8

# List available tests
./run_validation.sh --list
```

**Test passes if:**
- Shows `[PASS] Output matches gold file`
- Segfaults during cleanup are harmless

**Generate validation plots:**
```bash
cd scripts
./plot_validation.sh           # Generate all plots
./plot_validation.sh -c ortho  # Plot specific category
```

**Output:** Figures saved to `test/validation/figures/`
- `ortho.png` — All 6 orthotropic directions on one panel
- `asymmetry.png` — Tension-compression asymmetry + biaxial loading
- `postyield.png` — All 5 post-yield evolution modes
- `grains.png` — Multi-grain orientation effects
- `czm.png` — All 4 CZM modes (Mode I/II/mixed/complete separation)
- `fabric.png` — Fabric tensor vs explicit tensor overlay (elastic slope + yield)

## Usage Examples

### Example 1: Explicit Orthotropic (Dense Aragonite)

For dense crystalline materials with measured orthotropic constants:

```
[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    fill_method = symmetric9
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []

  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    use_fabric_scaling = false

    # Yield strengths (MPa)
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060

    # Post-yield behavior
    postyield_mode = 'exp_softening'
    residual_strength = 0.7
    kmin = 0.002
    kslope = 100.0
  []
[]
```

### Example 2: Porosity-Scaled Orthotropic (Porous Coral)

For porous materials with known full-density constants:

```
[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    fill_method = symmetric9
    density_rho = 0.85
    elasticity_density_exponent = 2.0
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []

  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    use_fabric_scaling = false
    sigma_xx_tension = 4980
    # ... (all yield strengths)
    density_rho = 0.85
    yield_density_exponent = 1.8
  []
[]
```

### Example 3: Fabric-Based (Trabecular Bone)

For trabecular structures with fabric tensor from μCT:

```
[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    E_0 = 10000.0
    nu_0 = 0.3
    fabric_m1 = 1.1        # Fabric eigenvalues (sum = 3)
    fabric_m2 = 1.0
    fabric_m3 = 0.9
    density_rho = 0.2      # BV/TV from μCT
    exponent_k = 2.0
    exponent_l = 1.0
  []

  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    use_fabric_scaling = true

    # Base tissue yield strengths [MPa] (Wolfram et al. 2012)
    sigma_0_tension = 74.589
    sigma_0_compression = 111.724
    tau_0 = 47.33
    zeta_0 = 0.22

    # Fabric tensor (must match elasticity)
    fabric_m1 = 1.1
    fabric_m2 = 1.0
    fabric_m3 = 0.9
    density_rho = 0.2
    exponent_p = 1.686
    exponent_q = 1.02

    postyield_mode = 'exp_softening'
    residual_strength = 0.7
  []
[]
```

### Example 4: Grain Boundary Interfaces (CZM)

Protein-aragonite interfaces with Gauss-Hermite probabilistic homogenization:

```
[Materials]
  [czm_interface]
    type = HomogenizedExponentialCZM
    boundary = 'grain_boundaries'

    # MD-derived strengths (Kvashin et al. 2026)
    normal_strength  = 626.0       # MPa
    shear_strength_s = 374.0       # MPa
    shear_strength_t = 374.0       # MPa

    # Mode-specific critical openings
    delta_0_normal  = 1.91e-7      # mm (191 pm)
    delta_0_tangent = 2.17e-7      # mm (217 pm)

    # Shape exponents (eta < mu required)
    mu  = 0.92
    eta = 0.27

    # Gauss-Hermite homogenization
    quality_std_dev         = 0.10   # Within-QP quality factor CV
    spatial_quality_std_dev = 0.15   # Across-QP spatial variation
    spatial_random_seed     = 1234   # Deterministic; no random_seed needed

    # Regularization
    damage_viscosity        = 2.5    # s
    damage_activation_width = 0.02
  []
[]
```

**Key properties of the GH model:**
- Deterministic — identical results across runs and MPI partition counts
- Area-independent — mesh independence by construction
- 2000× faster than previous Monte Carlo implementation

### Example 5: Polycrystalline RVE

See `test/validation/rve/` for complete input files. The `grain_interfaces` test runs a polycrystalline mesh with GH cohesive interfaces active and compares CSV output against gold files.

### NOTE on RVE Tests (Validation vs Production)

The `rve/dry` and `rve/grain_interfaces` tests use **short runs** for validation:

```bash
# Validation (default)
./run_validation.sh -c rve      # ~10 min, elastic regime only

# Production
cd test/validation/rve/grain_interfaces
vim grain_interfaces.i          # increase end_time for full deformation
../../../aragonite-opt -i grain_interfaces.i
```

## Documentation

Complete technical documentation: [`doc/aragonite_plasticity_documentation.pdf`](doc/aragonite_plasticity_documentation.pdf)

**Contents:**
- **Chapter 1:** Introduction and physical system
- **Chapter 2:** Mathematical formulation (yield surface, fabric tensor, return mapping)
- **Chapter 3:** Cohesive zone model (GH homogenization, damage evolution, quality factor distributions)
- **Chapters 4 & 5:** MOOSE installation and building
- **Chapter 6:** Complete parameter reference
- **Chapter 7:** Validation workflow and examples
- **Chapter 8:** Troubleshooting

## Material Parameters

### Aragonite (Dense Crystal)

| Property | Value |
|----------|-------|
| E₁, E₂, E₃ | 140.4, 70.3, 63.4 GPa |
| G₁₂, G₁₃, G₂₃ | 46.6, 31.1, 42.1 GPa |
| σ₁, σ₂, σ₃ (tension yield) | 4980, 4100, 5340 MPa |
| τ₁₂, τ₁₃, τ₂₃ (shear yield) | 4510, 5080, 5060 MPa |

### Protein-Aragonite Interface (CZM)

| Property | Protein | Water (001/100) |
|----------|---------|-----------------|
| Normal strength σₙ | 626 MPa | 742 MPa |
| Shear strength τₛ | 374 MPa | 354 MPa |
| Critical opening δ₀,ₙ | 191 pm | 224 pm |
| Critical opening δ₀,ₜ | 217 pm | 182 pm |
| Loading exponent μ | 0.916 | 1.224 |

Water-010 face excluded (non-equivalent electrostatics).

### Trabecular Bone (Fabric-Based)

**Base tissue properties:** E₀ = 10–20 GPa, ν₀ = 0.3, σ₀⁺ = 32–75 MPa, σ₀⁻ = 48–112 MPa

**Scaling exponents:** k = 2.0, l = 1.0 (elasticity); p = 1.28–1.69, q = 0.5–1.5 (plasticity)

**References:** [Wolfram et al. (2012)](https://doi.org/10.1016/j.jmbbm.2012.07.005), [Rincón-Kohli & Zysset (2009)](https://doi.org/10.1007/s10237-008-0128-z), [Zysset & Curnier (1995)](https://doi.org/10.1016/0167-6636(95)00018-6), [Schwiedrzik et al. (2013)](https://doi.org/10.1007/s10237-013-0472-5), [Wang et el. (2025)](https://doi.org/10.1016/j.conbuildmat.2025.142454), [Kvashin et al. (2026)](https://doi.org/10.1016/j.jmbbm.2026.107403).

## Implementation Details

### Material Classes

| Class | Purpose | Model |
|-------|---------|-------|
| `ComputeElasticityTensorCoupled` | Explicit orthotropic elasticity | Euler angle rotation |
| `ComputeFabricElasticityTensor` | Fabric-based elasticity | Zysset-Curnier (1995) |
| `OrthotropicPlasticityStressUpdate` | Anisotropic plasticity | Quadric yield surface |
| `HomogenizedExponentialCZM` | Grain boundary failure | GH probabilistic homogenization |

### Validation Status

- ✅ All 6 orthotropic elastic moduli validated against analytical solutions
- ✅ Tension-compression asymmetry verified with biaxial tests
- ✅ Five post-yield evolution modes tested independently
- ✅ Multi-grain load redistribution demonstrated
- ✅ CZM Mode I/II/mixed validated; mesh independence confirmed (512× element range)
- ✅ Fabric tensor (Zysset-Curnier) vs equivalent explicit symmetric9 tensor: stress-strain curves coincide at elastic slope and yield
- ✅ Complete RVE framework (grain boundaries + interfaces; two-run protocol)

## Project Structure

```
aragonite/
├── include/
│   ├── base/                # Application base
│   └── materials/           # Material headers (.h)
├── src/
│   ├── base/                # Application implementation
│   ├── main.C               # Main entry point
│   └── materials/           # Material implementation (.C)
├── test/
│   └── validation/
│       ├── ortho/           # 6 tests
│       ├── asymmetry/       # 2 tests
│       ├── postyield/       # 5 tests
│       ├── grains/          # 3 tests
│       ├── czm/             # 4 tests
│       ├── fabric/          # 2 tests (fabric_plasticity, explicit_plasticity)
│       ├── rve/             # 2 test  (dry, grain_interfaces)
│       └── figures/         # Generated plots
├── scripts/
│   ├── run_validation.sh
│   ├── plot_validation.sh
│   ├── plot_validation_category.py
│   ├── compare_with_gold.py
│   └── compute_fabric_reference.py
├── doc/
│   └── aragonite_plasticity_documentation.pdf
├── Makefile
├── LICENSE
└── README.md
```

## Citation

If you use this software in your research, please cite:

<!--
```bibtex
@article{yourname2026,
  title={Multiscale Finite Element Framework for Coral Aragonite Biomechanics},
  author={Your Name},
  journal={Computer Methods in Applied Mechanics and Engineering},
  year={2026},
  note={In preparation}
}
```
-->
<!-- After Zenodo upload, add DOI citation -->

## License

This project is licensed under the GNU Library General Public License v2.1 — see the [LICENSE](LICENSE) file for details.

---

**For detailed usage, see [documentation](doc/aragonite_plasticity_documentation.pdf).**
