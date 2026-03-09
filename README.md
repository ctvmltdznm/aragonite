# Orthotropic Plasticity Model for Aragonite Crystals

A multiscale finite element framework for modeling coral aragonite biomechanics, implemented in the [MOOSE](https://mooseframework.inl.gov/) framework.

<!-- Add after Zenodo upload: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX) -->

Fork "aragonite" to create a new MOOSE-based application.

For more information see: [https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app](https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app)

## Overview

This framework bridges molecular dynamics (MD) simulations to continuum finite element analysis (FEA) for modeдling coral aragonite biomechanics. Key features:

- **Orthotropic plasticity** with quadric yield surface (Schwiedrzik formulation)
- **Tension-compression asymmetry** via linear term in yield function
- **Five post-yield evolution modes**: perfect plasticity, exponential/linear hardening, exponential/piecewise softening
- **Homogenized cohesive zone model (CZM)** for grain boundary interfaces with MD-derived parameters
- **Euler angle-based orientation** for polycrystalline representative volume elements (RVEs)

## Physical System

Hierarchical coral biomineral structure:
- **Sclerodermites** (grains) containing aragonite needle crystals
- **Interfaces**: Protein-water layers between needles (CZ1) and between grains (CZ2)
- **MD data**: Interfaces are 22–49× weaker than bulk aragonite

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

```bash
cd scripts

# Run all 21 validation tests
./run_validation.sh

# Generate comparison plots
./plot_validation.sh

# View plots
cd ../test/validation/figures
ls *.png
```

## Usage

### Example: Polycrystalline RVE

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
aragonite-opt -i polycrystal_10grain.i
```

## Validation Test Suite

**21 tests across 6 categories:**

| Category | Tests | Purpose |
|----------|-------|---------|
| **Orthotropic** | 6 | Elastic constants (E₁, E₂, E₃, G₁₂, G₁₃, G₂₃) and yield stresses |
| **Asymmetry** | 2 | Tension-compression asymmetry (σ_c ≠ σ_t) |
| **Post-yield** | 5 | Five r(κ) evolution modes |
| **Grains** | 3-4 | Multi-grain orientation effects |
| **CZM** | 4 | Interface failure (Mode I, II, mixed, cycling) |
| **RVE** | 2 | Complete multiscale framework |

**Run validation:**
```bash
# Quick validation, specific category
cd scripts
./run_validation.sh -c ortho

# All tests
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



### RVE Tests (Validation vs Production)

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
- Publication-quality results

**Note:** Validation gold files reflect SHORT RUN (elastic only). 
Production results will differ (plastic regime).

## Documentation

Complete technical documentation: [`doc/aragonite_plasticity_documentation.pdf`](doc/aragonite_plasticity_documentation.pdf)

**Contents:**
- Chapter 1: Introduction and physical system
- Chapter 2: Mathematical formulation (yield surface, return mapping)
- Chapter 3: CZM formulation
- Chapters 4&5: MOOSE installation and building
- Chapter 6: Input file parameter reference
- Chapter 7: Troubleshooting most common issues
- Chapter 9: Summary

## Material Parameters

**Aragonite elastic constants (from MD):**
- E₁ = 140.4 GPa, E₂ = 70.3 GPa, E₃ = 113.3 GPa
- G₁₂ = 46.6 GPa, G₁₃ = 46.9 GPa, G₂₃ = 42.1 GPa
- ν₁₂ = 0.31, ν₁₃ = 0.30, ν₂₃ = 0.39

**Yield strengths:**
- Tension: σ₁ = 4980 MPa, σ₂ = 4100 MPa, σ₃ = 5340 MPa
- Compression: σ₁ = 3900 MPa, σ₂ = 3100 MPa, σ₃ = 4200 MPa
- Shear: τ₁₂ = 4510 MPa, τ₁₃ = 5080 MPa, τ₂₃ = 5060 MPa

## Project Structure

```
aragonite/
├── include/materials/    # Header files (.h)
├── src/materials/        # Implementation (.C)
├── test/validation/      # Validation test suite
├── scripts/              # Helper scripts
├── examples/             # Example input files
├── doc/                  # LaTeX documentation
└── Makefile              # Build system
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

This project is licensed under the GNU Library Public License - see the [LICENSE](LICENSE) file for details.

---

**For detailed usage, see [documentation](doc/aragonite_plasticity_documentation.pdf).**

