// CORRECTED Orthotropic Plasticity - UMAT Style (Unnormalized)
// CRITICAL FIX: f = sqrt(σ:F:σ) + f_lin·σ - r(κ)*σ_0 (NOT sqrt of sum!)

#pragma once
#include "StressUpdateBase.h"

class OrthotropicPlasticityStressUpdate : public StressUpdateBase
{
public:
  static InputParameters validParams();
  OrthotropicPlasticityStressUpdate(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;
  virtual bool requiresIsotropicTensor() override { return false; }

protected:
  virtual void updateState(RankTwoTensor & strain_increment,
                          RankTwoTensor & inelastic_strain_increment,
                          const RankTwoTensor & rotation_increment,
                          RankTwoTensor & stress_new,
                          const RankTwoTensor & stress_old,
                          const RankFourTensor & elasticity_tensor,
                          const RankTwoTensor & elastic_strain_old,
                          bool compute_full_tangent_operator,
                          RankFourTensor & tangent_operator) override;

  // Core functions
  void computeRotationMatrix(Real phi1, Real Phi, Real phi2,
                            RankTwoTensor & R, RankTwoTensor & R_inv) const;
  RankTwoTensor rotateToMaterial(const RankTwoTensor & t, const RankTwoTensor & R) const;
  RankTwoTensor rotateToGlobal(const RankTwoTensor & t, const RankTwoTensor & R) const;

  RankFourTensor rotateElasticityTensor(const RankFourTensor & C, 
                                     const RankTwoTensor & R) const;

  // CORRECTED yield function (UMAT line 892)
  Real computeYieldFunction(const RankTwoTensor & stress, Real kappa,
                           Real delta_kappa, Real dt) const;
  RankTwoTensor computeYieldGradient(const RankTwoTensor & stress) const;
  RankFourTensor computeYieldHessian(const RankTwoTensor & stress) const;
  
  // Dyadic product: C_ijkl = A_ij * B_kl (UMAT: VECDYAD)
  RankFourTensor dyadicProduct(const RankTwoTensor & A, const RankTwoTensor & B) const;

  Real computeSofteningFactor(Real kappa) const;
  Real computeSofteningDerivative(Real kappa) const;
  Real computeDamage(Real kappa) const;
  Real computeDamageDerivative(Real kappa) const;

  // Viscosity computation
  Real computeViscosity(Real dkappa, Real HI, Real dt) const;
  void computeViscosityDerivatives(Real dkappa, Real HI, 
                                   const RankTwoTensor & dHI_ds, Real dHI_dk, Real dt,
                                   RankTwoTensor & dvisc_ds, Real & dvisc_dk) const;

  // Two-stage return mapping
  bool performNewtonRaphson(const RankTwoTensor & stress_trial,
                           const RankFourTensor & C_inv,
                           Real kappa_old, Real damage_old, Real dt,
                           RankTwoTensor & stress_new,
                           Real & delta_kappa, Real & kappa_new, Real & damage_new,
                           unsigned int & iterations);

  bool performPrimalCPP(const RankTwoTensor & stress_trial,
                       const RankFourTensor & C_inv,
                       Real kappa_old, Real damage_old, Real dt,
                       RankTwoTensor & stress_new,
                       Real & delta_kappa, Real & kappa_new, Real & damage_new,
                       unsigned int & iterations);

  // Parameters - Yield (tension)
  const Real _sigma_xx_tension, _sigma_yy_tension, _sigma_zz_tension;
  const Real _tau_xy_max, _tau_xz_max, _tau_yz_max;

  // Parameters - Yield (compression, defaults to tension)
  const Real _sigma_xx_compression, _sigma_yy_compression, _sigma_zz_compression;

  // Euler angles
  const VariableValue & _euler_angle_1, & _euler_angle_2, & _euler_angle_3;
  //const Real _euler_angle_1, _euler_angle_2, _euler_angle_3;

  // Quadric yield surface
  std::vector<std::vector<Real>> _F_matrix;  // 6×6
  std::vector<Real> _f_lin_vector;           // 6×1 linear term

  /// Post-yield behavior mode
  enum class PostYieldMode
  {
    PERFECT_PLASTICITY,
    EXP_HARDENING,
    LINEAR_HARDENING,
    EXP_SOFTENING,
    PIECEWISE_SOFTENING
  };
  
  PostYieldMode _postyield_mode;
  
  /// Post-yield parameters
  Real _residual_strength;  // Residual strength ratio (0-1)
  Real _kslope;            // Hardening/softening rate
  Real _kmax;              // Start of softening transition
  Real _kmin;              // End of softening transition

  /// Viscosity mode
  enum class ViscosityMode
  {
    RATE_INDEPENDENT,
    LINEAR,
    EXPONENTIAL,
    LOGARITHMIC,
    POLYNOMIAL,
    POWERLAW
  };
  
  ViscosityMode _viscosity_mode;

  const Real _eta;
  const Real _m;  // Viscosity exponent (for exponential, logarithmic, etc.)

  // Damage (optional, D=0 default)
  const bool _use_damage;
  const Real _damage_critical, _damage_rate;

  // Numerical
  const Real _absolute_tolerance;  // Dimensionless (yield function is normalized!)
  const unsigned int _max_iterations_newton, _max_iterations_primal;
  const bool _use_primal_cpp;
  const Real _line_search_beta, _line_search_eta;

  // State variables
  MaterialProperty<Real> & _equivalent_plastic_strain;
  const MaterialProperty<Real> & _equivalent_plastic_strain_old;
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;

  // Diagnostics
  MaterialProperty<Real> & _return_mapping_stage;
  MaterialProperty<Real> & _return_mapping_iterations;
};
