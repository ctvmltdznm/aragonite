//* This file is part of the MOOSE framework
//* https://www.mooseframework.org

#pragma once

#include "StressUpdateBase.h"

class OrthotropicPlasticityStressUpdate : public StressUpdateBase
{
public:
  static InputParameters validParams();
  OrthotropicPlasticityStressUpdate(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void updateState(RankTwoTensor & strain_increment,
                          RankTwoTensor & inelastic_strain_increment,
                          const RankTwoTensor & rotation_increment,
                          RankTwoTensor & stress_new,
                          const RankTwoTensor & stress_old,
                          const RankFourTensor & elasticity_tensor,
                          const RankTwoTensor & elastic_strain_old,
                          bool compute_full_tangent_operator = false,
                          RankFourTensor & tangent_operator = _identityTensor) override;

  virtual TangentCalculationMethod getTangentCalculationMethod() override
  {
    return TangentCalculationMethod::ELASTIC;
  }

  virtual bool requiresIsotropicTensor() override
  {
    return false;
  }

  Real computeYieldFunction(const RankTwoTensor & stress);
  RankTwoTensor computeFlowDirection(const RankTwoTensor & stress);
  Real getCurrentYieldStress(unsigned int direction);

  // Yield stresses in material coordinate system (MPa)
  const Real _sigma_xx_max;
  const Real _sigma_yy_max;
  const Real _sigma_zz_max;
  const Real _tau_xy_max;
  const Real _tau_xz_max;
  const Real _tau_yz_max;

  // Softening parameters
  const Real _hardening_modulus;
  const Real _max_plastic_strain;
  const Real _absolute_tolerance;
  const unsigned int _max_iterations;

  // State variables
  MaterialProperty<Real> & _equivalent_plastic_strain;
  const MaterialProperty<Real> & _equivalent_plastic_strain_old;
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  unsigned int _active_direction;
  static RankFourTensor _identityTensor;
};