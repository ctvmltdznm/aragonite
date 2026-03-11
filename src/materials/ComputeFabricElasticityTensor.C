// ComputeFabricElasticityTensor.C
// Fabric-based anisotropic elasticity following Zysset & Curnier (1995)

#include "ComputeFabricElasticityTensor.h"

registerMooseObject("aragoniteApp", ComputeFabricElasticityTensor);

InputParameters
ComputeFabricElasticityTensor::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  
  params.addClassDescription(
    "Computes orthotropic elasticity tensor from base isotropic properties "
    "scaled by density and fabric tensor (Zysset & Curnier 1995). "
    "Use with OrthotropicPlasticityStressUpdate in fabric mode for "
    "consistent trabecular bone modeling.");
  
  // ==========================================================================
  // BASE ISOTROPIC PROPERTIES
  // ==========================================================================
  params.addRequiredParam<Real>("E_0", 
    "Base (isotropic reference) Young's modulus E₀ [MPa]");
  params.addRequiredParam<Real>("nu_0", 
    "Base Poisson's ratio ν₀ [-]. Typical: 0.3 for bone tissue.");
  params.addParam<Real>("G_0", 
    "Base shear modulus G₀ [MPa]. Default: E₀/(2(1+ν₀))");
  
  // ==========================================================================
  // FABRIC TENSOR EIGENVALUES
  // ==========================================================================
  params.addParam<Real>("fabric_m1", 1.0, 
    "Fabric tensor eigenvalue m₁ (direction 1). Normalized: m₁+m₂+m₃=3");
  params.addParam<Real>("fabric_m2", 1.0, 
    "Fabric tensor eigenvalue m₂ (direction 2)");
  params.addParam<Real>("fabric_m3", 1.0, 
    "Fabric tensor eigenvalue m₃ (direction 3)");
  
  // ==========================================================================
  // DENSITY
  // ==========================================================================
  params.addParam<Real>("density_rho", 1.0, 
    "Relative density ρ (BV/TV for bone). Range 0-1. "
    "Trabecular: 0.05-0.30, Cortical: 0.85-0.95, Dense: 1.0");
  
  // ==========================================================================
  // SCALING EXPONENTS
  // ==========================================================================
  params.addParam<Real>("exponent_k", 2.0, 
    "Density exponent k for elasticity. Typical: 1.8-2.0 for trabecular bone. "
    "Higher k = stronger density dependence.");
  params.addParam<Real>("exponent_l", 1.0, 
    "Fabric exponent l. Typical: 1.0. "
    "Higher l = stronger anisotropy from fabric.");
  
  // Cortical bone correction
  params.addParam<Real>("delta_cortical", 1.0,
    "Cortical bone density correction factor (UMAT DELTA). "
    "Only affects ρ > 0.5. Set to 1.0 to disable.");
  
  return params;
}

ComputeFabricElasticityTensor::ComputeFabricElasticityTensor(
    const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _E_0(getParam<Real>("E_0")),
    _nu_0(getParam<Real>("nu_0")),
    _density_rho(getParam<Real>("density_rho")),
    _fabric_m1(getParam<Real>("fabric_m1")),
    _fabric_m2(getParam<Real>("fabric_m2")),
    _fabric_m3(getParam<Real>("fabric_m3")),
    _exponent_k(getParam<Real>("exponent_k")),
    _exponent_l(getParam<Real>("exponent_l")),
    _delta_cortical(getParam<Real>("delta_cortical"))
{
  // Compute or read base shear modulus
  if (isParamValid("G_0"))
    _G_0 = getParam<Real>("G_0");
  else
    _G_0 = _E_0 / (2.0 * (1.0 + _nu_0));
  
  // Validate inputs
  if (_E_0 <= 0.0)
    mooseError("E_0 must be positive");
  if (_nu_0 < -1.0 || _nu_0 > 0.5)
    mooseError("nu_0 must be in range (-1, 0.5) for stability");
  if (_density_rho <= 0.0 || _density_rho > 1.0)
    mooseError("density_rho must be in range (0, 1]");
  if (_fabric_m1 <= 0.0 || _fabric_m2 <= 0.0 || _fabric_m3 <= 0.0)
    mooseError("Fabric eigenvalues must be positive");
  
  // Check fabric tensor normalization
  Real fabric_sum = _fabric_m1 + _fabric_m2 + _fabric_m3;
  if (std::abs(fabric_sum - 3.0) > 0.01)
    mooseWarning("Fabric eigenvalues should sum to 3.0, got ", fabric_sum, 
                 ". Consider normalizing.");
  
  // ==========================================================================
  // COMPUTE ORTHOTROPIC CONSTANTS (Zysset-Curnier 1995)
  // ==========================================================================
  
  // Density scaling with cortical correction (TSFU function from UMAT)
  Real rho_k = computeTSFU(_density_rho, _exponent_k, _delta_cortical);
  
  // Fabric powers
  Real m1_l = std::pow(_fabric_m1, _exponent_l);
  Real m2_l = std::pow(_fabric_m2, _exponent_l);
  Real m3_l = std::pow(_fabric_m3, _exponent_l);
  Real m1_2l = std::pow(_fabric_m1, 2.0 * _exponent_l);
  Real m2_2l = std::pow(_fabric_m2, 2.0 * _exponent_l);
  Real m3_2l = std::pow(_fabric_m3, 2.0 * _exponent_l);
  
  // Young's moduli: Eᵢ = E₀ × ρᵏ × mᵢ^(2l)
  _E1 = _E_0 * rho_k * m1_2l;
  _E2 = _E_0 * rho_k * m2_2l;
  _E3 = _E_0 * rho_k * m3_2l;
  
  // Shear moduli: Gᵢⱼ = G₀ × ρᵏ × mᵢˡ × mⱼˡ
  _G12 = _G_0 * rho_k * m1_l * m2_l;
  _G13 = _G_0 * rho_k * m1_l * m3_l;
  _G23 = _G_0 * rho_k * m2_l * m3_l;
  
  // Poisson's ratios: νᵢⱼ = ν₀ × (mᵢ/mⱼ)ˡ
  // These satisfy the symmetry requirement: νᵢⱼ/Eᵢ = νⱼᵢ/Eⱼ
  _nu12 = _nu_0 * std::pow(_fabric_m1 / _fabric_m2, _exponent_l);
  _nu13 = _nu_0 * std::pow(_fabric_m1 / _fabric_m3, _exponent_l);
  _nu23 = _nu_0 * std::pow(_fabric_m2 / _fabric_m3, _exponent_l);
  _nu21 = _nu_0 * std::pow(_fabric_m2 / _fabric_m1, _exponent_l);
  _nu31 = _nu_0 * std::pow(_fabric_m3 / _fabric_m1, _exponent_l);
  _nu32 = _nu_0 * std::pow(_fabric_m3 / _fabric_m2, _exponent_l);
  
  // Check positive definiteness condition
  // For orthotropic materials: 1 - ν₁₂ν₂₁ - ν₂₃ν₃₂ - ν₃₁ν₁₃ - 2ν₂₁ν₃₂ν₁₃ > 0
  Real delta = 1.0 - _nu12*_nu21 - _nu23*_nu32 - _nu31*_nu13 
               - 2.0*_nu21*_nu32*_nu13;
  if (delta <= 0.0)
    mooseError("Computed orthotropic constants violate positive definiteness. "
               "delta = ", delta, ". Try reducing nu_0 or adjusting fabric values.");
  
  // Report computed values
  Moose::out << "\n=== FABRIC-BASED ELASTICITY (Zysset-Curnier 1995) ===" << std::endl;
  Moose::out << "Input: E₀=" << _E_0 << " ν₀=" << _nu_0 
             << " G₀=" << _G_0 << " MPa" << std::endl;
  Moose::out << "Fabric: m₁=" << _fabric_m1 << " m₂=" << _fabric_m2 
             << " m₃=" << _fabric_m3 << " (sum=" << fabric_sum << ")" << std::endl;
  Moose::out << "Density: ρ=" << _density_rho << " k=" << _exponent_k 
             << " l=" << _exponent_l << " δ=" << _delta_cortical << std::endl;
  Moose::out << "TSFU(ρ,k,δ)=" << rho_k << " (effective density scaling)" << std::endl;
  Moose::out << "Computed orthotropic constants:" << std::endl;
  Moose::out << "  E₁=" << _E1 << " E₂=" << _E2 << " E₃=" << _E3 << " MPa" << std::endl;
  Moose::out << "  G₁₂=" << _G12 << " G₁₃=" << _G13 << " G₂₃=" << _G23 << " MPa" << std::endl;
  Moose::out << "  ν₁₂=" << _nu12 << " ν₁₃=" << _nu13 << " ν₂₃=" << _nu23 << std::endl;
  Moose::out << "  ν₂₁=" << _nu21 << " ν₃₁=" << _nu31 << " ν₃₂=" << _nu32 << std::endl;
  Moose::out << "  Δ=" << delta << " (positive definiteness check)" << std::endl;
  Moose::out << "======================================================\n" << std::endl;
}

void
ComputeFabricElasticityTensor::computeQpElasticityTensor()
{
  // Build orthotropic stiffness tensor from pre-computed constants
  buildOrthotropicStiffness(_E1, _E2, _E3, _G12, _G13, _G23, 
                            _nu12, _nu13, _nu23);
}

void
ComputeFabricElasticityTensor::buildOrthotropicStiffness(
    Real E1, Real E2, Real E3,
    Real G12, Real G13, Real G23,
    Real nu12, Real nu13, Real nu23)
{
  // Compute conjugate Poisson's ratios from symmetry: νᵢⱼ/Eᵢ = νⱼᵢ/Eⱼ
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;
  
  // Positive definiteness denominator
  Real delta = 1.0 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2.0*nu21*nu32*nu13;
  
  // Build stiffness matrix components
  // C = S⁻¹ where S is compliance matrix
  //
  // For orthotropic material in Voigt notation:
  // [σ₁₁]   [C₁₁ C₁₂ C₁₃  0   0   0 ] [ε₁₁]
  // [σ₂₂]   [C₁₂ C₂₂ C₂₃  0   0   0 ] [ε₂₂]
  // [σ₃₃] = [C₁₃ C₂₃ C₃₃  0   0   0 ] [ε₃₃]
  // [σ₂₃]   [ 0   0   0  C₄₄  0   0 ] [2ε₂₃]
  // [σ₁₃]   [ 0   0   0   0  C₅₅  0 ] [2ε₁₃]
  // [σ₁₂]   [ 0   0   0   0   0  C₆₆] [2ε₁₂]
  
  Real C11 = E1 * (1.0 - nu23*nu32) / delta;
  Real C22 = E2 * (1.0 - nu13*nu31) / delta;
  Real C33 = E3 * (1.0 - nu12*nu21) / delta;
  
  Real C12 = E1 * (nu21 + nu31*nu23) / delta;  // = E2*(nu12 + nu32*nu13)/delta
  Real C13 = E1 * (nu31 + nu21*nu32) / delta;  // = E3*(nu13 + nu12*nu23)/delta
  Real C23 = E2 * (nu32 + nu12*nu31) / delta;  // = E3*(nu23 + nu21*nu13)/delta
  
  Real C44 = G23;  // yz shear
  Real C55 = G13;  // xz shear
  Real C66 = G12;  // xy shear
  
  // Fill the elasticity tensor
  _elasticity_tensor[_qp].zero();
  
  // Use fillFromInputVector with symmetric9 ordering
  // symmetric9: C11, C12, C13, C22, C23, C33, C44, C55, C66
  std::vector<Real> Cvals = {C11, C12, C13, C22, C23, C33, C44, C55, C66};
  _elasticity_tensor[_qp].fillFromInputVector(Cvals, RankFourTensor::symmetric9);
}

Real
ComputeFabricElasticityTensor::computeTSFU(Real rho, Real exponent, Real delta) const
{
  // UMAT lines 2036-2045: TSFU function for density scaling with cortical correction
  // For ρ ≤ 0.5: simple power law
  // For ρ > 0.5: adds cortical bone correction term
  //
  // TSFU(ρ, p, δ) = ρᵖ                           if ρ ≤ 0.5
  //               = ρᵖ + (δ-1)((ρ-0.5)/0.5)ᵖ     if ρ > 0.5
  //
  // Set δ = 1.0 to disable cortical correction (reduces to simple power law)
  
  Real rho_p = std::pow(rho, exponent);
  
  if (rho > 0.5 && std::abs(delta - 1.0) > 1e-12)
  {
    // Cortical correction for dense bone
    Real correction = (delta - 1.0) * std::pow((rho - 0.5) / 0.5, exponent);
    return rho_p + correction;
  }
  
  return rho_p;
}