// ComputeFabricElasticityTensor.h
// Fabric-based anisotropic elasticity following Zysset & Curnier (1995)
// Computes orthotropic stiffness from isotropic base properties + fabric tensor

#pragma once

#include "ComputeElasticityTensorBase.h"

/**
 * ComputeFabricElasticityTensor computes an orthotropic elasticity tensor
 * from base isotropic properties scaled by density and fabric tensor.
 * 
 * Based on: Zysset PK, Curnier A (1995) "An alternative model for anisotropic 
 * elasticity based on fabric tensors." Mech Mater 21:243-250.
 * 
 * The compliance tensor is transformed as:
 *   S = (1/ρᵏ) M⁻ˡ S₀ M⁻ˡ
 * 
 * Where S₀ is isotropic compliance and M is the fabric tensor (diagonal in
 * principal axes with eigenvalues m₁, m₂, m₃).
 * 
 * This yields directional properties:
 *   Eᵢ = E₀ × ρᵏ × mᵢ^(2l)
 *   Gᵢⱼ = G₀ × ρᵏ × mᵢˡ × mⱼˡ  
 *   νᵢⱼ = ν₀ × (mᵢ/mⱼ)ˡ
 */
class ComputeFabricElasticityTensor : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();
  ComputeFabricElasticityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;
  
  /// Build 6x6 stiffness matrix from orthotropic constants
  void buildOrthotropicStiffness(Real E1, Real E2, Real E3,
                                  Real G12, Real G13, Real G23,
                                  Real nu12, Real nu13, Real nu23);
  
  // ==========================================================================
  // BASE ISOTROPIC PROPERTIES
  // ==========================================================================
  const Real _E_0;          ///< Base Young's modulus [MPa]
  const Real _nu_0;         ///< Base Poisson's ratio [-]
  Real _G_0;                ///< Base shear modulus [MPa] (computed or input)
  
  // ==========================================================================
  // FABRIC AND DENSITY PARAMETERS
  // ==========================================================================
  const Real _density_rho;  ///< Relative density (BV/TV), 0-1
  const Real _fabric_m1;    ///< Fabric eigenvalue in direction 1
  const Real _fabric_m2;    ///< Fabric eigenvalue in direction 2
  const Real _fabric_m3;    ///< Fabric eigenvalue in direction 3
  
  // ==========================================================================
  // SCALING EXPONENTS
  // ==========================================================================
  const Real _exponent_k;   ///< Density exponent (typically 1.8-2.0)
  const Real _exponent_l;   ///< Fabric exponent (typically 1.0)
  
  // ==========================================================================
  // COMPUTED ORTHOTROPIC CONSTANTS (for output/debugging)
  // ==========================================================================
  Real _E1, _E2, _E3;           ///< Directional Young's moduli
  Real _G12, _G13, _G23;        ///< Shear moduli
  Real _nu12, _nu13, _nu23;     ///< Poisson's ratios
  Real _nu21, _nu31, _nu32;     ///< Conjugate Poisson's ratios
};
