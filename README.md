# Orthotropic Plasticity Model for Aragonite Crystals

A multiscale finite element framework for modeling coral aragonite biomechanics and trabecular bone, implemented in the [MOOSE](https://mooseframework.inl.gov/) framework.

[![License: LGPL](https://img.shields.io/badge/License-LGPL-blue.svg)](https://opensource.org/licenses/LGPL-2.1)
<!-- Add after Zenodo upload: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX) -->

Fork "aragonite" to create a new MOOSE-based application.

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
- **Homogenized cohesive zone model (CZM)** for grain boundary interfaces with MD-derived parameters
- **Euler angle-based orientation** for polycrystalline representative volume elements (RVEs)

## Physical Systems

### Primary: Coral Aragonite Biomechanics

Hierarchical coral biomineral structure:
- **Dense aragonite needles:** Orthotropic single crystals / 'dry' grains
- **Porous sclerodermites:** Aragonite aggregates with varying porosity
- **Interfaces:** Protein-water layers between needles (CZ1) and between grains (CZ2)
- **MD data:** Interfaces are 22–49× weaker than bulk aragonite

### Secondary: Trabecular Bone Mechanics
 
Skeletal tissue modeling:
- **Base tissue properties:** E₀ = 10–20 GPa, ν₀ = 0.3
- **Fabric tensor:** From μCT analysis (structural anisotropy)
- **Apparent properties:** Density-fabric scaled constants at continuum level
- **Applications:** Vertebral bodies, femoral heads, bone quality assessment
 
## Installation

### Prerequisites

- [MOOSE Framework](https://mooseframework.inl.gov/getting_started/installation/index.html) (development version)
- C++17 compiler (GCC 9+, Clang 10+)
- Python 3.6+ with matplotlib (for validation plotting)

### Build

```bash
# Clone MOOSE (if not already installed)
git clone https://github.com/idaholab/moose.git
cd moose
git checkout master
./scripts/update_and_rebuild_libmesh.sh
cd ..

# Clone this repository
git clone https://github.com/ctvmltdznm/aragonite.git
cd aragonite

# Set MOOSE_DIR
export MOOSE_DIR=/path/to/moose

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

**21 tests across 6 categories:**

| Category | Tests | Purpose |
|----------|-------|---------|
| **Orthotropic** | 6 | Elastic constants (E₁, E₂, E₃, G₁₂, G₁₃, G₂₃) and yield stresses |
| **Asymmetry** | 2 | Tension-compression asymmetry (σ_c ≠ σ_t) |
| **Post-yield** | 5 | Five r(κ) evolution modes |
| **Grains** | 3-4 | Multi-grain orientation effects |
| **CZM** | 4 | Interface failure (Mode I, II, mixed, cycling) |
| **RVE** | 2 | Complete multiscale framework |

```bash
cd scripts

# After cloning, make the validation scripts executable:
chmod +x *.sh *.py

# Quick validation, specific category
cd scripts
./run_validation.sh -c ortho

# Run all 21 validation tests
./run_validation.sh

# Parallel execution
./run_validation.sh -p 4

# List available tests
./run_validation.sh --list
```

**Test passes if:**
- Shows `[PASS] Output matches gold file`
- Segfaults during cleanup are harmless

**Generate validation plots:**
```bash
cd scripts
./plot_validation.sh          # Generate all plots
./plot_validation.sh -c ortho # Plot specific category
```

**Output:** Figures saved to `test/validation/figures/`
- `ortho.png` - All 6 orthotropic directions on one panel
- `asymmetry.png` - Tension-compression asymmetry + biaxial loading
- `postyield.png` - All 5 post-yield evolution modes on one panel
- `grains.png` - Multi-grain orientation effects (solid=Grain 1, dashed=Grain 2)
- `czm.png` - All 4 CZM modes (Mode I/II/mixed/complete separation)

## Usage Examples
### Example 1: Explicit Orthotropic (Dense Aragonite)
 
For dense crystalline materials with measured orthotropic constants:
 
```bash
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
    use_fabric_scaling = false  # Explicit mode
    
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
 
```bash
[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    fill_method = symmetric9
    density_rho = 0.85                   # 85% density
    elasticity_density_exponent = 2.0    # Gibson-Ashby scaling
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
 
```bash
[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    E_0 = 10000.0          # Base tissue modulus [MPa]
    nu_0 = 0.3             # Base Poisson ratio
    fabric_m1 = 1.1        # Fabric eigenvalues (sum = 3)
    fabric_m2 = 1.0
    fabric_m3 = 0.9
    density_rho = 0.2      # BV/TV from μCT
    exponent_k = 2.0       # Elasticity density exponent
    exponent_l = 1.0       # Fabric exponent
  []
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    use_fabric_scaling = true  # Fabric mode
    
    # Base tissue yield strengths [MPa]
    sigma_0_tension = 74.589
    sigma_0_compression = 111.724
    tau_0 = 47.33
    zeta_0 = 0.22
    
    # Fabric tensor (must match elasticity)
    fabric_m1 = 1.1
    fabric_m2 = 1.0
    fabric_m3 = 0.9
    
    # Density and exponents
    density_rho = 0.2
    exponent_p = 1.686     # Yield density exponent
    exponent_q = 1.02      # Yield fabric exponent
    
    # Post-yield behavior
    postyield_mode = 'exp_softening'
    residual_strength = 0.7
  []
[]
```

### Example 4: Polycrystalline RVE

```bash
[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'tiny_aragonite_blocks.e'  # Exodus mesh with orientations
  []
[]

# Create boundaries using coordinate expressions
# ADJUST the coordinate thresholds based on your mesh size!
# OPtionally add cohesive zone interfaces
# See [`test/validation/rve'] for details

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
  []
[]

# Add AUX variables and kernels

[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    # Reads Euler angles from exodus file
  []
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
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
...

# Run simulation
aragonite-opt -i input_file_name.i
```

### NOTE on RVE Tests (Validation vs Production)

The `rve/dry` and `rve/grain_interfaces` tests use **SHORT RUNS** for validation:

**Validation mode (default):**
```bash
./run_validation.sh         # Includes RVE tests (1-2 min each)
```
- Loads to ~5-20% of yield, then stops
- Validates: exodus mesh reading, Euler angle assignment, CZM setup
- Runtime: ~10 minutes per test
- Sufficient for framework verification

**Production mode (for research):**
```bash
cd test/validation/rve/grain_interfaces

# Edit input file
vim grain_interfaces.i

# Change this line:
end_time = 1.0    # Change this

# To this:
end_time = 100.0  # Or higher for more deformation

# Run full simulation
../../../aragonite-opt -i grain_interfaces.i

# Output: grain_interfaces_out.e (for ParaView)
```

**Production gives you:**
- Complete plastic deformation in all grains
- Interface debonding and failure (if CZM enabled)
- Realistic polycrystal stress-strain curves

**Note:** Validation gold files reflect SHORT RUN (elastic only). Production results will differ (plastic regime).

## Documentation

Complete technical documentation: [`doc/aragonite_plasticity_documentation.pdf`](doc/aragonite_plasticity_documentation.pdf)

**Contents:**
- **Chapter 1:** Introduction and physical system
- **Chapter 2:** Mathematical formulation (yield surface, fabric tensor theory, return mapping)
- **Chapter 3:** Cohesive zone model formulation
- **Chapters 4 & 5:** MOOSE installation and building your application
- **Chapter 6:** Complete parameter reference (explicit + fabric modes)
- **Chapter 7:** Validation workflow and examples
- **Chapter 8:** Troubleshooting common issues
- **Chapter 9:** Summary

## Material Parameters

### Aragonite (Dense Crystal)

**Aragonite elastic constants:**
- E₁ = 140.4 GPa, E₂ = 70.3 GPa, E₃ = 113.3 GPa
- G₁₂ = 46.6 GPa, G₁₃ = 46.9 GPa, G₂₃ = 42.1 GPa
- ν₁₂ = 0.31, ν₁₃ = 0.30, ν₂₃ = 0.39

**Yield strengths:**
- Tension: σ₁ = 4980 MPa, σ₂ = 4100 MPa, σ₃ = 5340 MPa
- Compression: σ₁ = 3900 MPa, σ₂ = 3100 MPa, σ₃ = 4200 MPa
- Shear: τ₁₂ = 4510 MPa, τ₁₃ = 5080 MPa, τ₂₃ = 5060 MPa

### Trabecular Bone (Fabric-Based)

**Base tissue properties:**
- E₀ = 10000–20000 MPa
- ν₀ = 0.3
- σ₀⁺ = 32–75 MPa (Rincón-Kohli vs Wolfram)
- σ₀⁻ = 48–112 MPa

**Scaling exponents:**
- Elasticity: k = 2.0, l = 1.0
- Plasticity: p = 1.28–1.69, q = 0.5–1.5

**References:** [Wolfram et al. (2012)](https://doi.org/10.1016/j.jmbbm.2012.07.005), [Rincón-Kohli & Zysset (2009)](https://doi.org/10.1007/s10237-008-0128-z), [Zysset & Curnier (1995)](https://doi.org/10.1016/0167-6636(95)00018-6), [Schwiedrzik et al. (2013)](https://doi.org/10.1007/s10237-013-0472-5).

## Implementation Details

### Key Materials

| Material | Purpose | Model |
|----------|---------|-------|
| `ComputeElasticityTensorCoupled` | Explicit orthotropic elasticity | Euler angle rotation |
| `ComputeFabricElasticityTensor` | Fabric-based elasticity | Zysset-Curnier (1995) |
| `OrthotropicPlasticityStressUpdate` | Anisotropic plasticity | Quadric yield surface |
| `HomogenizedExponentialCZM` | Dual-exponential CZM | MD-derived parameters |

### Validation Features

- ✅ All 6 orthotropic elastic moduli validated against analytical solutions
- ✅ Tension-compression asymmetry verified with biaxial tests
- ✅ Five post-yield evolution modes tested independently
- ✅ Multi-grain load redistribution demonstrated (elastic unloading in fixed grain)
- ✅ CZM Mode I/II/mixed validated against MD traction-separation curves
- ✅ Complete RVE framework with grain boundaries and interfaces

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
│   └── validation/          # Validation test suite
│       ├── ortho/           # 6 tests
│       ├── asymmetry/       # 2 tests
│       ├── postyield/       # 5 tests
│       ├── grains/          # 4 tests
│       ├── czm/             # 4 tests
│       ├── rve/             # 2 tests
│       └── figures/         # Generated plots
├── scripts/
│   ├── run_validation.sh    # Run test suite
│   ├── plot_validation.sh   # Generate plots
│   ├── plot_validation_category.py
│   └── compare_with_gold.py
├── doc/
│   └── aragonite_plasticity_documentation.pdf
├── Makefile                 # Build system
├── LICENSE                  # LGPL license
└── README.md                # This file
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

This project is licensed under the GNU Library General Public License v2.1 - see the [LICENSE](LICENSE) file for details.


---

**For detailed usage, see [documentation](doc/aragonite_plasticity_documentation.pdf).**

