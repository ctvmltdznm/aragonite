// HomogenizedExponentialCZM.h
// PROOF OF CONCEPT: Area-based homogenization with identical contacts
// Purpose: Test that 1 μm² interface = 10,000 × (10 nm²) MD contacts

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

class HomogenizedExponentialCZM : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  HomogenizedExponentialCZM(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeInterfaceTractionAndDerivatives() override;
  
  /// Compute peak strength for mixed-mode loading (Mohr-Coulomb)
  Real computePeakStrength(const Real delta_n, const Real delta_s, const Real delta_t);
  
  /// Compute peak strength for a specific contact (with distributed strengths)
  Real computePeakStrengthForContact(const Real delta_n, 
                                     const Real delta_s, 
                                     const Real delta_t,
                                     const Real sigma_n, 
                                     const Real tau_s, 
                                     const Real tau_t);

  /// Compute damage variable for a single contact
  Real computeDamage(Real delta_eff);
  
  /// Compute damage for a specific contact (delta_0 variation)
  Real computeDamageForContact(Real delta_eff, Real delta_0_contact, Real phi);
  
  /// Compute derivative of damage
  Real computeDamageDeriv(Real delta_eff, Real delta_0_use, Real phi);
  Real computeDamageDerivWrtDelta0(Real delta_eff, Real delta_0_mixed, Real eta);

  void computeDamageGradient(Real delta_eff, Real delta_0_contact,
                            Real delta_n_pos, Real delta_s, Real delta_t,
                            Real phi, Real dD_grad[3]);

  /// Mixed mode separation - following Wang25
  Real computeMixedModeDelta0(Real phi, Real delta_0_n, Real delta_0_t) const;
  Real getEffectiveDelta0() const;
  
  /// Solve for failure displacement from traction threshold
  Real solveForFailureDisplacement(Real traction_ratio, Real exponent);
  
  /// Compute interface area at quadrature point
  Real computeInterfaceArea();

  // ============================================================================
  // MD-DERIVED PARAMETERS (per contact)
  // ============================================================================
  
  /// Peak strengths from MD (MPa)
  const Real _normal_strength;
  const Real _shear_strength_s;
  const Real _shear_strength_t;
  
  /// Characteristic length from MD (mm) - this is the MD scale!
  const Real _delta_0;
  const Real _delta_0_normal;   // Mode I (normal) peak separation
  const Real _delta_0_tangent;  // Mode II (tangential) peak separation

  /// Failure displacement
  Real _delta_c;
  
  /// Exponential coefficients
  const Real _mu;
  const Real _eta;
  
  /// Traction ratio for failure threshold
  const Real _failure_traction_ratio;
  
  /// Viscoplastic time scale (seconds)
  const Real _damage_viscosity;
  
  /// Normal gap tolerance
  const Real _normal_gap_tol;
  
  // ============================================================================
  // HOMOGENIZATION PARAMETERS
  // ============================================================================
  
  /// Area of single MD contact (mm²) - typically (10 nm)² = 1e-10 mm²
  const Real _md_contact_area;
  
  /// Maximum number of contacts to simulate (for computational efficiency)
  const unsigned int _max_contacts;
  
  // ============================================================================
  // STORAGE FOR HOMOGENIZATION
  // ============================================================================
  
  /// Number of MD contacts at this quadrature point
  MaterialProperty<Real> & _n_contacts;
  const MaterialProperty<Real> & _n_contacts_old;
  
  /// Damage state for each contact (vector per QP)
  MaterialProperty<std::vector<Real>> & _contact_damage;
  const MaterialProperty<std::vector<Real>> & _contact_damage_old;
  
  // ============================================================================
  // DISTRIBUTION PARAMETERS
  // ============================================================================

  std::vector<Real> _contact_strength_n;  // Normal strength for each contact
  std::vector<Real> _contact_strength_s;  // Shear-s strength for each contact
  std::vector<Real> _contact_strength_t;  // Shear-t strength for each contact
  std::vector<Real> _contact_delta0;     // delta_0 for each contact

  // Add parameters for distribution
  Real _strength_std_dev;  // Standard deviation (as fraction of mean)
  Real _initial_damage_max;      // Maximum initial damage (0 to 1)
  Real _delta0_std_dev;          // NEW: Standard deviation for delta_0

  // Seed for ditributions
  const unsigned int _random_seed;

  // ============================================================================
  // OUTPUT PROPERTIES
  // ============================================================================
  
  /// Average damage across all contacts (for AuxVariable)
  MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;
  
  /// Effective opening (for AuxVariable)
  MaterialProperty<Real> & _delta_eff;
  const MaterialProperty<Real> & _delta_eff_old;
  
  /// Interface area at this QP (for output/debugging)
  MaterialProperty<Real> & _interface_area;
  
  /// Number of contacts at this QP (exposed as material property for output)
  MaterialProperty<Real> & _n_contacts_prop;
};
