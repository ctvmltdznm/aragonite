// COMPLETE UMAT IMPLEMENTATION FOR MOOSE
// Orthotropic plasticity with damage, viscosity, and Primal CPPA
// Based on UMAT_QUADRIC_PRIMAL_Major_LINEAR_HARDENING.f
// For production use in research

#include "OrthotropicPlasticityStressUpdate.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "VoigtHelpers.h"

using namespace libMesh;

// ============================================================================
// MISSING MOOSE FRAMEWORK INTEGRATION
// Add these to the TOP of your .C file (after #includes, before other functions)
// ============================================================================

registerMooseObject("aragoniteApp", OrthotropicPlasticityStressUpdate);

InputParameters
OrthotropicPlasticityStressUpdate::validParams()
{
  InputParameters params = StressUpdateBase::validParams();
  
  params.addClassDescription("Orthotropic plasticity with UMAT-style implementation");
  
  // Yield strengths - Tension
  params.addRequiredParam<Real>("sigma_xx_tension", "Yield strength in xx direction (tension) [MPa]");
  params.addRequiredParam<Real>("sigma_yy_tension", "Yield strength in yy direction (tension) [MPa]");
  params.addRequiredParam<Real>("sigma_zz_tension", "Yield strength in zz direction (tension) [MPa]");
  params.addRequiredParam<Real>("tau_xy_max", "Maximum shear stress in xy plane [MPa]");
  params.addRequiredParam<Real>("tau_xz_max", "Maximum shear stress in xz plane [MPa]");
  params.addRequiredParam<Real>("tau_yz_max", "Maximum shear stress in yz plane [MPa]");
  
  // Yield strengths - Compression (optional, defaults to tension)
  params.addParam<Real>("sigma_xx_compression", "Yield strength in xx direction (compression) [MPa]");
  params.addParam<Real>("sigma_yy_compression", "Yield strength in yy direction (compression) [MPa]");
  params.addParam<Real>("sigma_zz_compression", "Yield strength in zz direction (compression) [MPa]");
  
  // Euler angles for material orientation
  params.addRequiredCoupledVar("euler_angle_1", "First Euler angle (phi1) in degrees");
  params.addRequiredCoupledVar("euler_angle_2", "Second Euler angle (Phi) in degrees");
  params.addRequiredCoupledVar("euler_angle_3", "Third Euler angle (phi2) in degrees");
  //params.addRequiredParam<Real>("euler_angle_1", "First Euler angle (phi1) in degrees");
  //params.addRequiredParam<Real>("euler_angle_2", "Second Euler angle (Phi) in degrees");
  //params.addRequiredParam<Real>("euler_angle_3", "Third Euler angle (phi2) in degrees");

  // params.addParam<Real>("kappa_start", 0.002, 
  //                        "Plastic strain at which softening begins");
  //params.addParam<Real>("softening_range", 0.015, 
  //                        "Plastic strain range over which linear softening occurs");
  
  // Viscoplasticity
  //params.addParam<bool>("use_viscoplasticity", false, "Enable viscoplastic regularization");
  params.addParam<Real>("eta", 1e-3, "Inverse viscosity parameter (MPa*s)^-1");
  
  // Damage
  params.addParam<bool>("use_damage", false, "Enable damage evolution");
  params.addParam<Real>("damage_critical", 0.05, "Critical plastic strain for damage onset");
  params.addParam<Real>("damage_rate", 0.1, "Damage evolution rate");
  
  // Numerical parameters
  params.addParam<Real>("absolute_tolerance", 1e-6, "Absolute tolerance for yield function");
  params.addParam<unsigned int>("max_iterations_newton", 10, "Maximum iterations for Newton-Raphson");
  params.addParam<unsigned int>("max_iterations_primal", 1000, "Maximum iterations for Primal CPPA");
  params.addParam<bool>("use_primal_cpp", true, "Enable Primal CPPA as fallback");
  params.addParam<Real>("line_search_beta", 1e-4, "Line search beta parameter (Armero 2002)");
  params.addParam<Real>("line_search_eta", 0.25, "Line search eta parameter");

  // Post-yield mode selection
  MooseEnum postyield_mode("perfect exp_hardening linear_hardening "
                          "simple_softening exp_softening piecewise_softening",
                          "exp_softening");
  params.addParam<MooseEnum>("postyield_mode", postyield_mode,
                            "Post-yield behavior mode");
  
  // Post-yield parameters (mode-specific)
  params.addParam<Real>("residual_strength", 0.7, 
    "Residual strength ratio (0-1). For softening: final strength. "
    "For hardening: maximum additional strength.");
  params.addParam<Real>("kslope", 10.0, "Hardening/softening rate");
  params.addParam<Real>("kmax", 0.001, "Start of softening transition");
  params.addParam<Real>("kmin", 0.015, "End of softening transition");

  // Viscocity
  MooseEnum viscosity_mode("rate_independent linear exponential "
                          "logarithmic polynomial powerlaw",
                          "linear");
  params.addParam<MooseEnum>("viscosity_mode", viscosity_mode,
                            "Viscoplastic regularization mode");
  params.addParam<Real>("m", 0.001,
                    "Viscosity exponent for non-linear viscosity models");

  return params;
}

OrthotropicPlasticityStressUpdate::OrthotropicPlasticityStressUpdate(
    const InputParameters & parameters)
  : StressUpdateBase(parameters),
    // Yield strengths - tension
    _sigma_xx_tension(getParam<Real>("sigma_xx_tension")),
    _sigma_yy_tension(getParam<Real>("sigma_yy_tension")),
    _sigma_zz_tension(getParam<Real>("sigma_zz_tension")),
    _tau_xy_max(getParam<Real>("tau_xy_max")),
    _tau_xz_max(getParam<Real>("tau_xz_max")),
    _tau_yz_max(getParam<Real>("tau_yz_max")),
    
    // Yield strengths - compression (default to tension if not provided)
    _sigma_xx_compression(isParamValid("sigma_xx_compression") ? 
                          getParam<Real>("sigma_xx_compression") : _sigma_xx_tension),
    _sigma_yy_compression(isParamValid("sigma_yy_compression") ? 
                          getParam<Real>("sigma_yy_compression") : _sigma_yy_tension),
    _sigma_zz_compression(isParamValid("sigma_zz_compression") ? 
                          getParam<Real>("sigma_zz_compression") : _sigma_zz_tension),
    
    // Euler angles
    _euler_angle_1(coupledValue("euler_angle_1")),
    _euler_angle_2(coupledValue("euler_angle_2")),
    _euler_angle_3(coupledValue("euler_angle_3")),
    //_euler_angle_1(getParam<Real>("euler_angle_1")),
    //_euler_angle_2(getParam<Real>("euler_angle_2")),
    //_euler_angle_3(getParam<Real>("euler_angle_3")),
    
    // Softening
    //_softening_exponent(getParam<Real>("softening_exponent")),
    //_residual_strength_ratio(getParam<Real>("residual_strength_ratio")),
    //_max_plastic_strain(getParam<Real>("max_plastic_strain")),

    // Softening configuration
    //_softening_type(getParam<std::string>("softening_type")),
    //_kappa_start(getParam<Real>("kappa_start")),
    //_softening_range(getParam<Real>("softening_range")),

    // Plastic flow
    _residual_strength(getParam<Real>("residual_strength")),
    _kslope(getParam<Real>("kslope")),
    _kmax(getParam<Real>("kmax")),
    _kmin(getParam<Real>("kmin")),
    
    // Viscoplasticity
    //_use_viscoplasticity(getParam<bool>("use_viscoplasticity")),
    _eta(getParam<Real>("eta")),
    // viscosity exponent
    _m(getParam<Real>("m")),
    
    // Damage
    _use_damage(getParam<bool>("use_damage")),
    _damage_critical(getParam<Real>("damage_critical")),
    _damage_rate(getParam<Real>("damage_rate")),
    
    // Numerical
    _absolute_tolerance(getParam<Real>("absolute_tolerance")),
    _max_iterations_newton(getParam<unsigned int>("max_iterations_newton")),
    _max_iterations_primal(getParam<unsigned int>("max_iterations_primal")),
    _use_primal_cpp(getParam<bool>("use_primal_cpp")),
    _line_search_beta(getParam<Real>("line_search_beta")),
    _line_search_eta(getParam<Real>("line_search_eta")),
    
    // State variables
    _equivalent_plastic_strain(declareProperty<Real>("effective_plastic_strain")),
    _equivalent_plastic_strain_old(getMaterialPropertyOld<Real>("effective_plastic_strain")),
    _plastic_strain(declareProperty<RankTwoTensor>("plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("plastic_strain")),
    _damage(declareProperty<Real>("damage_variable")),
    _damage_old(getMaterialPropertyOld<Real>("damage_variable")),
    
    // Diagnostics
    _return_mapping_stage(declareProperty<Real>("return_mapping_stage")),
    _return_mapping_iterations(declareProperty<Real>("return_mapping_iterations"))
{
  // Parse post-yield mode
  MooseEnum mode = getParam<MooseEnum>("postyield_mode");
  if (mode == "perfect")
    _postyield_mode = PostYieldMode::PERFECT_PLASTICITY;
  else if (mode == "exp_hardening")
    _postyield_mode = PostYieldMode::EXP_HARDENING;
  else if (mode == "linear_hardening")
    _postyield_mode = PostYieldMode::LINEAR_HARDENING;
  else if (mode == "exp_softening")
    _postyield_mode = PostYieldMode::EXP_SOFTENING;
  else if (mode == "piecewise_softening")
    _postyield_mode = PostYieldMode::PIECEWISE_SOFTENING;

  // Viscocity
  MooseEnum visc_mode = getParam<MooseEnum>("viscosity_mode");
  if (visc_mode == "rate_independent")
    _viscosity_mode = ViscosityMode::RATE_INDEPENDENT;
  else if (visc_mode == "linear")
    _viscosity_mode = ViscosityMode::LINEAR;
  else if (visc_mode == "exponential")
    _viscosity_mode = ViscosityMode::EXPONENTIAL;
  else if (visc_mode == "logarithmic")
    _viscosity_mode = ViscosityMode::LOGARITHMIC;
  else if (visc_mode == "polynomial")
    _viscosity_mode = ViscosityMode::POLYNOMIAL;
  else if (visc_mode == "powerlaw")
    _viscosity_mode = ViscosityMode::POWERLAW;

  // Initialize F matrix and f_lin vector (6x6 and 6x1)
  _F_matrix.resize(6, std::vector<Real>(6, 0.0));
  _f_lin_vector.resize(6, 0.0);
  
  // Build quadric yield surface from yield strengths
  // F matrix for orthotropic yield (Voigt notation)
  // Based on Hill-type criterion with tension/compression asymmetry
  
  // Compute F matrix components
  Real sigma_t_xx = _sigma_xx_tension;
  Real sigma_t_yy = _sigma_yy_tension;
  Real sigma_t_zz = _sigma_zz_tension;
  Real sigma_c_xx = _sigma_xx_compression;
  Real sigma_c_yy = _sigma_yy_compression;
  Real sigma_c_zz = _sigma_zz_compression;
  
  // Average strengths for quadric part
  Real sigma_xx_avg = 0.5 * (sigma_t_xx + sigma_c_xx);
  Real sigma_yy_avg = 0.5 * (sigma_t_yy + sigma_c_yy);
  Real sigma_zz_avg = 0.5 * (sigma_t_zz + sigma_c_zz);
  
  // Normalize to get dimensionless F matrix
  Real norm = 1.0;  // Keep dimensional for now
  
  // Diagonal terms (normal stresses)
  _F_matrix[0][0] = norm / (sigma_xx_avg * sigma_xx_avg);  // xx
  _F_matrix[1][1] = norm / (sigma_yy_avg * sigma_yy_avg);  // yy
  _F_matrix[2][2] = norm / (sigma_zz_avg * sigma_zz_avg);  // zz
  
  // Shear terms
  _F_matrix[3][3] = norm / (_tau_yz_max * _tau_yz_max);  // yz, not sure about 4.0
  _F_matrix[4][4] = norm / (_tau_xz_max * _tau_xz_max);  // xz
  _F_matrix[5][5] = norm / (_tau_xy_max * _tau_xy_max);  // xy
  
  //debug
  /*Moose::out << "*** norm = " << norm << "\n";
  Moose::out << "*** tau_xy_max = " << _tau_xy_max << "\n";
  Moose::out << "*** F_matrix[5][5] = " << _F_matrix[5][5] << "\n";
  Moose::out << "*** 1/sqrt(F_66) = " << 1.0/std::sqrt(_F_matrix[5][5]) << "\n";
  Moose::out << "*** CHECKING MATRIX PRODUCTS:\n";
  Real test_stress = 2427.0;
  Real test_quadric = test_stress * test_stress * _F_matrix[5][5];
  Moose::out << "*** If stress=2427: quadric=" << test_quadric << ", phi=" << std::sqrt(test_quadric) << "\n";
  */
  // Linear term for tension/compression asymmetry
  _f_lin_vector[0] = (sigma_c_xx - sigma_t_xx) / (sigma_c_xx * sigma_t_xx);  // xx
  _f_lin_vector[1] = (sigma_c_yy - sigma_t_yy) / (sigma_c_yy * sigma_t_yy);  // yy
  _f_lin_vector[2] = (sigma_c_zz - sigma_t_zz) / (sigma_c_zz * sigma_t_zz);  // zz
  // Shear components are zero (symmetric in shear)
}

void
OrthotropicPlasticityStressUpdate::initQpStatefulProperties()
{
  _equivalent_plastic_strain[_qp] = 0.0;
  _plastic_strain[_qp].zero();
  _damage[_qp] = 0.0;
  _return_mapping_stage[_qp] = -1.0;
  _return_mapping_iterations[_qp] = 0.0;
}

void
OrthotropicPlasticityStressUpdate::propagateQpStatefulProperties()
{
  _equivalent_plastic_strain[_qp] = _equivalent_plastic_strain_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _damage[_qp] = _damage_old[_qp];
}

// ============================================================================
// ROTATION FUNCTIONS
// ============================================================================

void
OrthotropicPlasticityStressUpdate::computeRotationMatrix(
    Real phi1, Real Phi, Real phi2,
    RankTwoTensor & R, RankTwoTensor & R_inv) const
{
  // Convert degrees to radians
  Real phi1_rad = phi1 * M_PI / 180.0;
  Real Phi_rad = Phi * M_PI / 180.0;
  Real phi2_rad = phi2 * M_PI / 180.0;
  
  // Compute rotation matrix using ZXZ Euler angles (Bunge convention)
  Real c1 = std::cos(phi1_rad);
  Real s1 = std::sin(phi1_rad);
  Real c = std::cos(Phi_rad);
  Real s = std::sin(Phi_rad);
  Real c2 = std::cos(phi2_rad);
  Real s2 = std::sin(phi2_rad);
  
  // R = Rz(phi1) * Rx(Phi) * Rz(phi2)
  R.zero();
  R(0, 0) = c1*c2 - s1*s2*c;
  R(0, 1) = -c1*s2 - s1*c2*c;
  R(0, 2) = s1*s;
  R(1, 0) = s1*c2 + c1*s2*c;
  R(1, 1) = -s1*s2 + c1*c2*c;
  R(1, 2) = -c1*s;
  R(2, 0) = s2*s;
  R(2, 1) = c2*s;
  R(2, 2) = c;
  
  // R_inv = R^T (orthogonal matrix)
  R_inv = R.transpose();
}

RankTwoTensor
OrthotropicPlasticityStressUpdate::rotateToMaterial(
    const RankTwoTensor & t, const RankTwoTensor & R) const
{
  // t_material = R^T * t_global * R
  return R.transpose() * t * R;
}

RankTwoTensor
OrthotropicPlasticityStressUpdate::rotateToGlobal(
    const RankTwoTensor & t, const RankTwoTensor & R) const
{
  // t_global = R * t_material * R^T
  return R * t * R.transpose();
}

RankFourTensor
OrthotropicPlasticityStressUpdate::rotateElasticityTensor(
    const RankFourTensor & C, const RankTwoTensor & R) const
{
  // Rotate 4th-order elasticity tensor: C'_ijkl = R_im R_jn C_mnop R_ko R_lp
  RankFourTensor C_rotated;
  C_rotated.zero();
  
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int n = 0; n < 3; n++)
              for (unsigned int o = 0; o < 3; o++)
                for (unsigned int p = 0; p < 3; p++)
                  C_rotated(i, j, k, l) += R(i, m) * R(j, n) * C(m, n, o, p) * R(k, o) * R(l, p);
  
  return C_rotated;
}

// ============================================================================
// MAIN UPDATE STATE FUNCTION
// ============================================================================

void
OrthotropicPlasticityStressUpdate::updateState(
    RankTwoTensor & strain_increment,
    RankTwoTensor & inelastic_strain_increment,
    const RankTwoTensor & rotation_increment,
    RankTwoTensor & stress_new,
    const RankTwoTensor & stress_old,
    const RankFourTensor & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old,
    bool compute_full_tangent_operator,
    RankFourTensor & tangent_operator)
{
  // Get Euler angles at current quadrature point
  Real phi1 = _euler_angle_1[_qp];
  Real Phi = _euler_angle_2[_qp];
  Real phi2 = _euler_angle_3[_qp];
  
  //Real phi1 = -_euler_angle_1;
  //Real Phi = -_euler_angle_2;
  //Real phi2 = -_euler_angle_3;
  
  // Compute rotation matrices
  RankTwoTensor R, R_inv;
  computeRotationMatrix(phi1, Phi, phi2, R, R_inv);

  // Rotaion in plastic regime is consistent with MOOSE elastic rotations
  std::swap(R, R_inv);

  // Rotate stress to material coordinates
  RankTwoTensor stress_old_material = rotateToMaterial(stress_old, R);

  // CRITICAL: Rotate elasticity tensor to material coordinates!
  RankFourTensor elasticity_tensor_material = rotateElasticityTensor(elasticity_tensor, R_inv);
  //RankFourTensor elasticity_tensor_material = rotateElasticityTensor(elasticity_tensor, R);

  // CRITICAL: Rotate strain increment to material coordinates too!
  RankTwoTensor strain_increment_material = rotateToMaterial(strain_increment, R);

  // Compute trial stress in material coordinates (now consistent!)
  RankTwoTensor stress_trial_material = stress_old_material + 
                                        elasticity_tensor_material * strain_increment_material;

  // Get compliance tensor in material coordinates
  RankFourTensor C_inv_material = elasticity_tensor_material.invSymm();
  //

  // Instead of commented block above - NO ROTATION - work in global coordinates (same as material for angles=0)
  //RankTwoTensor stress_trial = stress_old + elasticity_tensor * strain_increment;
  //RankFourTensor C_inv = elasticity_tensor.invSymm();

  // Get old state
  Real kappa_old = _equivalent_plastic_strain_old[_qp];
  Real damage_old = _damage_old[_qp];
  Real dt = _t - _t_old;
  
  // Return mapping
  RankTwoTensor stress_return;
  Real delta_kappa = 0.0;
  Real kappa_new = kappa_old;
  Real damage_new = damage_old;
  unsigned int iterations = 0;
  
  bool success = false;
  
  // Try Newton-Raphson first
  success = performNewtonRaphson(stress_trial_material, C_inv_material,
                                 kappa_old, damage_old, dt,
                                 stress_return, delta_kappa, 
                                 kappa_new, damage_new, iterations);
  
  // If Newton fails and Primal is enabled, try Primal CPPA
  if (!success && _use_primal_cpp) {
    mooseWarning("Newton failed at qp=", _qp, ", switching to Primal CPPA");
    iterations = 0;
    success = performPrimalCPP(stress_trial_material, C_inv_material,
                               kappa_old, damage_old, dt,
                               stress_return, delta_kappa,
                               kappa_new, damage_new, iterations);
  }
  
  if (!success) {
    mooseError("Both Newton and Primal CPPA failed at qp=", _qp);
  }
  
  // Rotate stress back to global coordinates - for constistency with angles in elastic
  stress_new = rotateToGlobal(stress_return, R);

  // without rotation, stay in material frame
  // stress_new = stress_return;
  
  // Compute inelastic strain increment in material coordinates
  RankTwoTensor plastic_strain_increment = delta_kappa * computeYieldGradient(stress_return);
  
  // Rotate plastic strain to global coordinates
  inelastic_strain_increment = rotateToGlobal(plastic_strain_increment, R);
  // without rotation
  //inelastic_strain_increment = plastic_strain_increment;
  
  // Update state variables
  _equivalent_plastic_strain[_qp] = kappa_new;
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + plastic_strain_increment;
  _damage[_qp] = damage_new;
  _return_mapping_iterations[_qp] = iterations;
  
  // Tangent operator (elastic for now - could add consistent tangent)
  if (compute_full_tangent_operator)
    tangent_operator = elasticity_tensor;
}

// Initiatilization pass

// Alternative: Manual LU decomposition if libMesh version unavailable
void manualInvert7x7(const std::vector<std::vector<Real>> & A,
                     std::vector<std::vector<Real>> & Ainv)
{
  const int n = 7;
  std::vector<std::vector<Real>> mat = A;  // Copy
  std::vector<std::vector<Real>> inv(n, std::vector<Real>(n, 0.0));
  std::vector<int> indx(n);
  
  // Initialize inverse to identity
  for (int i = 0; i < n; i++)
    inv[i][i] = 1.0;
  
  // LU decomposition
  for (int i = 0; i < n; i++) {
    // Find pivot
    int imax = i;
    Real amax = std::abs(mat[i][i]);
    for (int k = i + 1; k < n; k++) {
      if (std::abs(mat[k][i]) > amax) {
        amax = std::abs(mat[k][i]);
        imax = k;
      }
    }
    
    if (amax < 1e-14)
      mooseError("Singular matrix in 7x7 inversion");
    
    // Swap rows
    if (imax != i) {
      std::swap(mat[i], mat[imax]);
      std::swap(inv[i], inv[imax]);
    }
    
    // Eliminate column
    for (int k = i + 1; k < n; k++) {
      Real factor = mat[k][i] / mat[i][i];
      for (int j = i; j < n; j++)
        mat[k][j] -= factor * mat[i][j];
      for (int j = 0; j < n; j++)
        inv[k][j] -= factor * inv[i][j];
    }
  }
  
  // Back substitution
  for (int i = n - 1; i >= 0; i--) {
    for (int j = 0; j < n; j++) {
      inv[i][j] /= mat[i][i];
      for (int k = 0; k < i; k++)
        inv[k][j] -= mat[k][i] * inv[i][j];
    }
  }
  
  Ainv = inv;
}

// ============================================================================
// DYADIC PRODUCT: C_ijkl = A_ij * B_kl (UMAT: VECDYAD)
// ============================================================================
RankFourTensor
OrthotropicPlasticityStressUpdate::dyadicProduct(const RankTwoTensor & A,
                                                  const RankTwoTensor & B) const
{
  RankFourTensor C;
  C.zero();
  
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          C(i, j, k, l) = A(i, j) * B(k, l);
  
  return C;
}

// ============================================================================
// YIELD HESSIAN: ∂²f/∂σ∂σ (UMAT: DDSY, lines 1508-1509, 1548-1549)
// ============================================================================
RankFourTensor
OrthotropicPlasticityStressUpdate::computeYieldHessian(const RankTwoTensor & stress) const
{
  // Voigt notation
  std::vector<Real> s(6);
  s[0] = stress(0,0); s[1] = stress(1,1); s[2] = stress(2,2);
  s[3] = stress(1,2); s[4] = stress(0,2); s[5] = stress(0,1);


  // Compute Fσ (FFS in UMAT)
  std::vector<Real> Fs(6, 0.0);
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 6; j++)
      Fs[i] += _F_matrix[i][j] * s[j];
  
  // Compute σ:F:σ (SFFS in UMAT)
  Real SFFS = 0.0;
  for (unsigned int i = 0; i < 6; i++)
    SFFS += s[i] * Fs[i];
  
  if (SFFS < 1e-16)
    SFFS = 1e-16;  // Regularize
  
  // UMAT formula: DDSY = -1/(SFFS)^1.5 * VECDYAD(FFS,FFS) + 1/sqrt(SFFS) * FFFF
  Real inv_sffs_sqrt = 1.0 / std::sqrt(SFFS);
  Real inv_sffs_3_2 = -1.0 / std::pow(SFFS, 1.5);
  
  // Build Fσ as RankTwoTensor
  RankTwoTensor Fs_tensor;
  Fs_tensor.zero();
  Fs_tensor(0,0) = Fs[0]; Fs_tensor(1,1) = Fs[1]; Fs_tensor(2,2) = Fs[2];
  Fs_tensor(1,2) = Fs_tensor(2,1) = Fs[3];
  Fs_tensor(0,2) = Fs_tensor(2,0) = Fs[4];
  Fs_tensor(0,1) = Fs_tensor(1,0) = Fs[5];
  
  // First term: -1/(SFFS)^1.5 * (Fσ ⊗ Fσ)
  RankFourTensor dyadic_term = dyadicProduct(Fs_tensor, Fs_tensor);
  dyadic_term *= inv_sffs_3_2;
  
  // Second term: 1/sqrt(SFFS) * F
  RankFourTensor F_tensor;
  F_tensor.zero();
  
  // Voigt to tensor conversion
  auto voigt_to_tensor = [](int v, int & i, int & j) {
    const int map[6][2] = {{0,0}, {1,1}, {2,2}, {1,2}, {0,2}, {0,1}};
    i = map[v][0]; j = map[v][1];
  };
  
  for (unsigned int v1 = 0; v1 < 6; v1++) {
    for (unsigned int v2 = 0; v2 < 6; v2++) {
      int i, j, k, l;
      voigt_to_tensor(v1, i, j);
      voigt_to_tensor(v2, k, l);
      
      Real val = _F_matrix[v1][v2];
      F_tensor(i, j, k, l) = val;
      
      // Enforce tensor symmetry
      if (i != j) F_tensor(j, i, k, l) = val;
      if (k != l) F_tensor(i, j, l, k) = val;
      if (i != j && k != l) F_tensor(j, i, l, k) = val;
    }
  }
  
  F_tensor *= inv_sffs_sqrt;
  
  return dyadic_term + F_tensor;
}

// ============================================================================
// FULL NEWTON-RAPHSON (UMAT lines 1452-1522)
// ============================================================================
bool
OrthotropicPlasticityStressUpdate::performNewtonRaphson(
    const RankTwoTensor & stress_trial,
    const RankFourTensor & C_inv,
    Real kappa_old,
    Real damage_old,
    Real dt,
    RankTwoTensor & stress_new,
    Real & delta_kappa,
    Real & kappa_new,
    Real & damage_new,
    unsigned int & iterations)
{
  const Real TOL = _absolute_tolerance;
  
  // Initialize
  stress_new = stress_trial;
  delta_kappa = 0.0;
  kappa_new = kappa_old;
  damage_new = damage_old;
  iterations = 0;
  
  // Check if elastic
  Real f_trial = computeYieldFunction(stress_trial, kappa_old, 0.0, dt);
  if (f_trial <= TOL) {
    _return_mapping_stage[_qp] = -1;
    return true;
  }
  
  _return_mapping_stage[_qp] = 1;  // Plastic
  
  // Newton iteration variables
  RankTwoTensor stress_i = stress_trial;
  Real dkappa_i = 0.0;
  Real NORMRR = 1e10;
  Real ABSY = 1e10;
  
  while ((NORMRR > TOL || ABSY > TOL) && iterations < _max_iterations_newton) {
    iterations++;
    
    // Update state
    Real kappa_i = kappa_old + dkappa_i;
    Real r_i = computeSofteningFactor(kappa_i);
    Real dr_i = computeSofteningDerivative(kappa_i);  // Note: returns positive dr/dκ
    Real D_i = computeDamage(kappa_i);
    Real dD_i = computeDamageDerivative(kappa_i);
    
    // Compute gradient, Hessian (UMAT lines 1504-1512)
    RankTwoTensor DSY = computeYieldGradient(stress_i);
    RankFourTensor DDSY = computeYieldHessian(stress_i);
    
    Real HI = DSY.L2norm();
    if (HI < 1e-14) {
      mooseWarning("Newton: Singular yield gradient");
      return false;
    }
    
    RankTwoTensor NP = DSY / HI;
    
    // Compute yield function and residual
    Real f = computeYieldFunction(stress_i, kappa_i, dkappa_i, dt);
    
    RankTwoTensor trial_elastic = C_inv * stress_trial;
    RankTwoTensor RR = -(C_inv * (stress_i - stress_trial)) / (1.0 - D_i);
    
    if (_use_damage && D_i > 1e-6 && D_i < 0.99) {
      Real D_0 = damage_old;
      RR -= (D_i - D_0) / (1.0 - D_i) * trial_elastic;
    }
    
    RR -= dkappa_i * NP;
    
    NORMRR = RR.L2norm();
    ABSY = std::abs(f);
    
    if (NORMRR <= TOL && ABSY <= TOL) {
      stress_new = stress_i;
      delta_kappa = dkappa_i;
      kappa_new = kappa_i;
      damage_new = D_i;
      _return_mapping_stage[_qp] = 2;
      return true;
    }
    
    // ===== COMPUTE JACOBIAN (UMAT lines 1463-1472) =====
    
    // DHDS = 1/HI * DDSY * DSY
    RankTwoTensor DHDS = (DDSY * DSY) / HI;
    
    // DYDS = DSY + viscous correction
    RankTwoTensor DYDS = DSY;
    Real DYDK = dr_i;  // UMAT: DRAD = -dr/dκ, but my dr_i = +dr/dκ, so use -dr_i
    
    if (_viscosity_mode != ViscosityMode::RATE_INDEPENDENT && dt > 1e-16 && HI > 1e-14) {
      RankTwoTensor dvisc_ds;
      Real dvisc_dk;
      computeViscosityDerivatives(dkappa_i, HI, DHDS, 0.0, dt, dvisc_ds, dvisc_dk);
      DYDS += dvisc_ds;
      DYDK += dvisc_dk;
    }
    
    // DNPDS = (DDSY*HI - dyadic(DSY,DHDS)) / HI^2
    RankFourTensor DNPDS = (DDSY * HI - dyadicProduct(DSY, DHDS)) / (HI * HI);
    
    // DRRDS = -C/(1-D) - dkappa*DNPDS
    RankFourTensor DRRDS = -C_inv / (1.0 - D_i) - DNPDS * dkappa_i;
    
    // DRRDK = -dD/(1-D)^2 * [C*(σ-σ_tr) + (1-D0)*ε_tr] - NP
    RankTwoTensor DRRDK = -NP;
    if (_use_damage && D_i < 0.99 && dD_i > 1e-14) {
      RankTwoTensor damage_contrib = C_inv * (stress_i - stress_trial);
      Real D_0 = damage_old;
      if (D_0 < 0.99)
        damage_contrib += (1.0 - D_0) * trial_elastic;
      DRRDK -= (dD_i / ((1.0 - D_i) * (1.0 - D_i))) * damage_contrib;
    }
    
    // Invert DRRDS → SSSA (UMAT line 1474)
    RankFourTensor SSSA = (-DRRDS).invSymm();
    
    // Compute Newton update (UMAT lines 1476-1480)
    Real DYDS_norm = DYDS.L2norm();
    if (DYDS_norm < 1e-14)
      DYDS_norm = 1e-14;
    
    Real numerator = -(f / DYDS_norm + (NP * (SSSA * RR)).trace());
    RankTwoTensor SSSA_DRRDK = SSSA * DRRDK;
    Real denominator = (NP * SSSA_DRRDK).trace() + DYDK / DYDS_norm;
    
    if (std::abs(denominator) < 1e-14) {
      mooseWarning("Newton: Singular Jacobian");
      return false;
    }
    
    Real DDK1 = numerator / denominator;
    RankTwoTensor DSS1 = SSSA * (RR + DRRDK * DDK1);
    
    // Update
    stress_i += DSS1;
    dkappa_i += DDK1;
    
    // Enforce dκ >= 0
    if (dkappa_i < 0.0) {
      stress_i = stress_trial;
      dkappa_i = 0.0;
    }
  }
  
  // Failed to converge
  mooseWarning("Newton failed: iter=", iterations, ", NORMRR=", NORMRR, ", ABSY=", ABSY);
  return false;
}

// ============================================================================
// PRIMAL CPPA WITH LINE SEARCH (UMAT lines 1525-1856)
// ============================================================================
// All Voigt notation conversions done properly

bool
OrthotropicPlasticityStressUpdate::performPrimalCPP(
    const RankTwoTensor & stress_trial,
    const RankFourTensor & C_inv,
    Real kappa_old,
    Real damage_old,
    Real dt,
    RankTwoTensor & stress_new,
    Real & delta_kappa,
    Real & kappa_new,
    Real & damage_new,
    unsigned int & iterations)
{
  const Real TOL = _absolute_tolerance;
  const Real BETA = _line_search_beta;
  const Real ETAL = _line_search_eta;
  const unsigned int MAXITER = _max_iterations_primal;
  
  // Initialize
  RankTwoTensor stress_i = stress_trial;
  Real dkappa_i = 0.0;
  Real kappa_i = kappa_old;
  iterations = 0;
  
  // Check yield
  Real f = computeYieldFunction(stress_i, kappa_i, dkappa_i, dt);
  if (f <= TOL) {
    stress_new = stress_i;
    delta_kappa = 0.0;
    kappa_new = kappa_old;
    damage_new = damage_old;
    _return_mapping_stage[_qp] = -1;
    return true;
  }
  
  _return_mapping_stage[_qp] = 3;  // Primal CPPA
  
  // Initial state
  Real D_i = computeDamage(kappa_i);
  Real dD_i = computeDamageDerivative(kappa_i);
  Real r_i = computeSofteningFactor(kappa_i);
  Real dr_i = computeSofteningDerivative(kappa_i);
  
  RankTwoTensor DSY = computeYieldGradient(stress_i);
  Real HI = DSY.L2norm();
  if (HI < 1e-14) HI = 1e-14;
  RankTwoTensor NP = DSY / HI;
  
  RankTwoTensor trial_elastic = C_inv * stress_trial;
  RankTwoTensor RR = -(C_inv * (stress_i - stress_trial)) / (1.0 - D_i);
  if (_use_damage && D_i > 1e-6 && D_i < 0.99) {
    RR -= (D_i - damage_old) / (1.0 - D_i) * trial_elastic;
  }
  RR -= dkappa_i * NP;
  
  Real NORMRR = RR.L2norm();
  Real ABSY = std::abs(f);
  
  // Store initial residual RRI (7x1)
  std::vector<Real> RRI(7);
  tensorToVoigt(RR, RRI);  // First 6 components
  RRI[6] = f;  // 7th component
  
  // ===== PRIMAL CPPA LOOP =====
  while ((NORMRR > TOL || ABSY > TOL) && iterations < MAXITER) {
    iterations++;
    
    // Store current state for line search
    RankTwoTensor stress_old_iter = stress_i;
    Real dkappa_old_iter = dkappa_i;
    
    // Update state variables
    kappa_i = kappa_old + dkappa_i;
    r_i = computeSofteningFactor(kappa_i);
    dr_i = computeSofteningDerivative(kappa_i);
    D_i = computeDamage(kappa_i);
    dD_i = computeDamageDerivative(kappa_i);
    
    // Compute gradient and Hessian
    DSY = computeYieldGradient(stress_i);
    RankFourTensor DDSY = computeYieldHessian(stress_i);
    
    HI = DSY.L2norm();
    if (HI < 1e-14) {
      mooseWarning("Primal: Singular gradient at qp=", _qp);
      return false;
    }
    NP = DSY / HI;
    
    // ===== COMPUTE JACOBIAN TERMS =====
    
    // DHDS = 1/HI * DDSY * DSY
    RankTwoTensor DHDS = contractRankFourTwo(DDSY, DSY) / HI;
    
    // DYDS = DSY + viscous correction
    RankTwoTensor DYDS = DSY;
    Real DYDK = dr_i;
    if (_viscosity_mode != ViscosityMode::RATE_INDEPENDENT && dt > 1e-16 && HI > 1e-14) {
      RankTwoTensor dvisc_ds;
      Real dvisc_dk;
      computeViscosityDerivatives(dkappa_i, HI, DHDS, 0.0, dt, dvisc_ds, dvisc_dk);
      DYDS += dvisc_ds;
      DYDK += dvisc_dk;
    }
//      DYDS += -(_eta / dt) / (HI * HI) * DHDS;
//      DYDK += -(_eta / dt) / HI;
//    }
    
    // DNPDS = (DDSY*HI - dyadic(DSY,DHDS)) / HI^2
    RankFourTensor DNPDS = (DDSY * HI - dyadicProduct(DSY, DHDS)) / (HI * HI);
    
    // DRRDS = -C/(1-D) - dkappa*DNPDS
    RankFourTensor DRRDS = -C_inv / (1.0 - D_i) - DNPDS * dkappa_i;
    
    // DRRDK = damage terms - NP
    RankTwoTensor DRRDK = -NP;
    if (_use_damage && D_i < 0.99 && dD_i > 1e-14) {
      RankTwoTensor damage_contrib = C_inv * (stress_i - stress_trial);
      if (damage_old < 0.99)
        damage_contrib += (1.0 - damage_old) * trial_elastic;
      DRRDK -= (dD_i / ((1.0 - D_i) * (1.0 - D_i))) * damage_contrib;
    }
    
    // ===== BUILD 7x7 JACOBIAN SYSTEM (PROPER VOIGT) =====
    
    // Convert to Voigt notation
    std::vector<std::vector<Real>> DRRDS_voigt;
    rankFourToVoigt(DRRDS, DRRDS_voigt);  // 6x6
    
    std::vector<Real> DRRDK_voigt, DYDS_voigt;
    tensorToVoigt(DRRDK, DRRDK_voigt);  // 6x1
    tensorToVoigt(DYDS, DYDS_voigt);    // 6x1
    
    std::vector<Real> RR_voigt;
    tensorToVoigt(RR, RR_voigt);  // 6x1
    
    // Build 7x7 Jacobian: JJJ = [DRRDS  DRRDK]
    //                           [DYDS   DYDK ]
    std::vector<std::vector<Real>> JJJ(7, std::vector<Real>(7, 0.0));
    
    // J11 = DRRDS (6x6)
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
        JJJ[i][j] = DRRDS_voigt[i][j];
    
    // J12 = DRRDK (6x1)
    for (int i = 0; i < 6; i++)
      JJJ[i][6] = DRRDK_voigt[i];
    
    // J21 = DYDS (1x6)
    for (int j = 0; j < 6; j++)
      JJJ[6][j] = DYDS_voigt[j];
    
    // J22 = DYDK (scalar)
    JJJ[6][6] = DYDK;
    
    // Build 7x1 residual: RRR = [RR]
    //                           [f ]
    std::vector<Real> RRR(7);
    for (int i = 0; i < 6; i++)
      RRR[i] = RR_voigt[i];
    RRR[6] = f;
    
    // ===== INVERT 7x7 JACOBIAN =====
    std::vector<std::vector<Real>> INVJ;
    manualInvert7x7(JJJ, INVJ);
    
    // ===== DETERMINE SEARCH DIRECTION (PRIMAL vs NEWTON) =====
    
    // Check if at elastic trial (UMAT line 1612-1615)
    Real AUX = 0.0;
    if (std::abs(dkappa_i) < 1e-14) {
      AUX = (NP * (stress_i - stress_trial)).trace();
    }
    
    int FLAG = 0;
    std::vector<Real> DDD(7);
    
    if (AUX <= 0.0) {
      // NEWTON DIRECTION (UMAT line 1619-1623)
      FLAG = 1;
      
      // DDD = -INVJ * RRR
      matVecMult7(INVJ, RRR, DDD);
      for (int i = 0; i < 7; i++)
        DDD[i] = -DDD[i];
    }
    else {
      // PRIMAL DIRECTION (UMAT line 1625-1637)
      FLAG = 0;
      
      // DDJ = INVJ * INVJ^T
      std::vector<std::vector<Real>> DDJ(7, std::vector<Real>(7, 0.0));
      for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
          for (int k = 0; k < 7; k++)
            DDJ[i][j] += INVJ[i][k] * INVJ[j][k];  // Note: INVJ[j][k] is transpose
      
      // Zero out kappa row/column (enforce constraint)
      for (int i = 0; i < 6; i++) {
        DDJ[6][i] = 0.0;
        DDJ[i][6] = 0.0;
      }
      
      // DDD = -DDJ * JJJ^T * RRR
      std::vector<Real> temp1(7), temp2(7);
      
      // temp1 = JJJ^T * RRR
      for (int i = 0; i < 7; i++) {
        temp1[i] = 0.0;
        for (int j = 0; j < 7; j++)
          temp1[i] += JJJ[j][i] * RRR[j];
      }
      
      // temp2 = DDJ * temp1
      matVecMult7(DDJ, temp1, temp2);
      
      // DDD = -temp2
      for (int i = 0; i < 7; i++)
        DDD[i] = -temp2[i];
    }
    
    // ===== EXTRACT STRESS AND KAPPA UPDATES =====
    
    RankTwoTensor DSS1;
    voigtToTensor(DDD, DSS1);  // First 6 components → stress tensor
    Real DDK1 = DDD[6];         // 7th component → kappa
    
    // ===== TENTATIVE UPDATE =====
    
    stress_i = stress_old_iter + DSS1;
    dkappa_i = dkappa_old_iter + DDK1;
    
    // Store state vectors for line search
    std::vector<Real> XX1(7), XXI(7);
    tensorToVoigt(stress_i, XX1);
    XX1[6] = dkappa_i;
    tensorToVoigt(stress_old_iter, XXI);
    XXI[6] = dkappa_old_iter;
    
    // Enforce constraint: dκ >= 0
    if (dkappa_i < 0.0)
      dkappa_i = 0.0;
    
    kappa_i = kappa_old + dkappa_i;
    
    // ===== EVALUATE AT NEW POINT =====
    
    D_i = computeDamage(kappa_i);
    r_i = computeSofteningFactor(kappa_i);
    DSY = computeYieldGradient(stress_i);
    HI = DSY.L2norm();
    if (HI < 1e-14) HI = 1e-14;
    NP = DSY / HI;
    
    f = computeYieldFunction(stress_i, kappa_i, dkappa_i, dt);
    
    RR = -(C_inv * (stress_i - stress_trial)) / (1.0 - D_i);
    if (_use_damage && D_i > 1e-6 && D_i < 0.99)
      RR -= (D_i - damage_old) / (1.0 - D_i) * trial_elastic;
    RR -= dkappa_i * NP;
    
    // Update residual vector
    tensorToVoigt(RR, RR_voigt);
    for (int i = 0; i < 6; i++)
      RRR[i] = RR_voigt[i];
    RRR[6] = f;
    
    // ===== MERIT FUNCTION (UMAT line 1682, 1702) =====
    
    Real MKI = 0.5 * (RR * RR).trace() + 0.5 * f * f;
    
    // Compute merit derivative (for line search)
    Real DMK = 0.0;
    if (FLAG == 1) {
      // Newton direction: DMK = -2*MKI
      DMK = -2.0 * MKI;
    }
    else {
      // Primal direction: DMK = RRR^T * JJJ * DDD
      std::vector<Real> temp(7);
      matVecMult7(JJJ, DDD, temp);
      DMK = dotProduct7(RRR, temp);
    }
    
    Real MK1 = 0.5 * (RR * RR).trace() + 0.5 * f * f;
    
    // ===== LINE SEARCH (ARMERO 2002, UMAT lines 1704-1773) =====
    
    Real ALPHA = 1.0;
    Real UPLIM = 0.0;
    
    // Compute acceptance threshold (UMAT line 1707-1711)
    if (FLAG == 1 && dkappa_i >= 0.0) {
      // Newton direction
      UPLIM = (1.0 - 2.0 * BETA * ALPHA) * MKI;
    }
    else {
      // Primal direction: UPLIM = MKI + BETA * RRR^T * JJJ * (XX1 - XXI)
      std::vector<Real> diff(7), temp(7);
      for (int i = 0; i < 7; i++)
        diff[i] = XX1[i] - XXI[i];
      
      matVecMult7(JJJ, diff, temp);
      Real derivative = dotProduct7(RRR, temp);
      UPLIM = MKI + BETA * derivative;
    }
    
    // Perform line search if merit increased
    if (MK1 > UPLIM) {
      unsigned int ITERL = 0;
      
      while (MK1 > UPLIM && ALPHA >= BETA && ITERL < 20) {
        ITERL++;
        
        // Compute new step length (quadratic interpolation, UMAT line 1723-1730)
        Real ALPHA1 = ETAL * ALPHA;  // Linear backtrack
        Real ALPHA2 = -ALPHA * ALPHA * DMK / (2.0 * (MK1 - MKI - ALPHA * DMK));  // Quadratic
        
        ALPHA = (ALPHA1 >= ALPHA2) ? ALPHA1 : ALPHA2;
        
        // Update with reduced step
        dkappa_i = dkappa_old_iter + ALPHA * DDK1;
        if (dkappa_i < 0.0)
          dkappa_i = 0.0;
        
        kappa_i = kappa_old + dkappa_i;
        stress_i = stress_old_iter + DSS1 * ALPHA;
        
        // Update state vector
        tensorToVoigt(stress_i, XX1);
        XX1[6] = dkappa_i;
        
        // Evaluate at new point
        D_i = computeDamage(kappa_i);
        r_i = computeSofteningFactor(kappa_i);
        DSY = computeYieldGradient(stress_i);
        HI = DSY.L2norm();
        if (HI < 1e-14) HI = 1e-14;
        NP = DSY / HI;
        
        f = computeYieldFunction(stress_i, kappa_i, dkappa_i, dt);
        RR = -(C_inv * (stress_i - stress_trial)) / (1.0 - D_i);
        if (_use_damage && D_i > 1e-6 && D_i < 0.99)
          RR -= (D_i - damage_old) / (1.0 - D_i) * trial_elastic;
        RR -= dkappa_i * NP;
        
        // Update merit function
        MK1 = 0.5 * (RR * RR).trace() + 0.5 * f * f;
        
        // Recompute acceptance threshold
        if (FLAG == 1 && dkappa_i >= 0.0) {
          UPLIM = (1.0 - 2.0 * BETA * ALPHA) * MKI;
        }
        else {
          std::vector<Real> diff(7), temp(7);
          for (int i = 0; i < 7; i++)
            diff[i] = XX1[i] - XXI[i];
          matVecMult7(JJJ, diff, temp);
          Real derivative = dotProduct7(RRI, temp);  // Use original RRI
          UPLIM = MKI + BETA * derivative;
        }
      }
      
      // After line search, recompute full state (UMAT lines 1775-1799)
      D_i = computeDamage(kappa_i);
      dD_i = computeDamageDerivative(kappa_i);
      r_i = computeSofteningFactor(kappa_i);
      dr_i = computeSofteningDerivative(kappa_i);
      
      DSY = computeYieldGradient(stress_i);
      DDSY = computeYieldHessian(stress_i);
      HI = DSY.L2norm();
      if (HI < 1e-14) HI = 1e-14;
      NP = DSY / HI;
      
      f = computeYieldFunction(stress_i, kappa_i, dkappa_i, dt);
      RR = -(C_inv * (stress_i - stress_trial)) / (1.0 - D_i);
      if (_use_damage && D_i > 1e-6 && D_i < 0.99)
        RR -= (D_i - damage_old) / (1.0 - D_i) * trial_elastic;
      RR -= dkappa_i * NP;
      
      // Update RRI for next iteration
      tensorToVoigt(RR, RRI);
      RRI[6] = f;
    }
    
    // ===== CHECK CONVERGENCE =====
    
    NORMRR = RR.L2norm();
    ABSY = std::abs(f);
    
    if (NORMRR <= TOL && ABSY <= TOL) {
      stress_new = stress_i;
      delta_kappa = dkappa_i;
      kappa_new = kappa_i;
      damage_new = D_i;
      _return_mapping_stage[_qp] = 4;  // Primal converged
      return true;
    }
  }
  
  // Failed to converge
  mooseWarning("Primal CPPA failed after ", iterations, " iterations, NORMRR=", NORMRR, ", ABSY=", ABSY);
  return false;
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// ============================================================================
// CONFIGURABLE SOFTENING - BOTH EXPONENTIAL AND PIECEWISE LINEAR
// EXPONENTIAL SOFTENING FACTOR: r(κ) = r_res + (1 - r_res) * exp(-β*κ)
// ============================================================================

Real
OrthotropicPlasticityStressUpdate::computeSofteningFactor(Real kappa) const
{
  Real r = 1.0;
  
  switch (_postyield_mode)
  {
    case PostYieldMode::PERFECT_PLASTICITY:
      r = 1.0;
      break;
    
    case PostYieldMode::EXP_HARDENING:
      // Start at 1.0, harden to (1 + hardening_amount)
      // r = _residual_strength + (1.0 - _residual_strength) * (1.0 - std::exp(-_kslope * kappa));
      r = 1.0 + _residual_strength * (1.0 - std::exp(-_kslope * kappa));
      break;
    
    case PostYieldMode::LINEAR_HARDENING:
      // Start at 1.0, harden linearly
      // kslope controls hardening rate
      r = 1.0 + _kslope * kappa;
      break;

    case PostYieldMode::EXP_SOFTENING:
      //r = _rdy + (1.0 - _rdy) * std::exp(-_kslope * kappa);
      if (kappa < _kmax)
        r = 1.0;  // Perfect plasticity before kmax
      else
      {
        Real kappa_shifted = kappa - _kmax;
        r = _residual_strength + (1.0 - _residual_strength) * std::exp(-_kslope * kappa_shifted);
      }
      break;

    case PostYieldMode::PIECEWISE_SOFTENING:
      if (kappa < _kmax)
      {
        r = 1.0;  // Perfect plasticity before softening
      }
      else if (kappa < _kmin)
      {
        // Simple LINEAR drop from 1.0 to gmin
        Real progress = (kappa - _kmax) / (_kmin - _kmax);
        r = 1.0 - (1.0 - _residual_strength) * progress;
      }
      else
      {
        r = _residual_strength;  // Residual plateau
      }
      break;  // No scaling!
  }
  
  return r;
}

// ============================================================================
// computeSofteningDerivative() - NEW VERSION WITH BOTH TYPES
// EXPONENTIAL SOFTENING DERIVATIVE: dr/dκ = -(1 - r_res) * β * exp(-β*κ)
// ============================================================================

Real
OrthotropicPlasticityStressUpdate::computeSofteningDerivative(Real kappa) const
{
  Real dr_dk = 0.0;
  
  switch (_postyield_mode)
  {
    case PostYieldMode::PERFECT_PLASTICITY:
      dr_dk = 0.0;
      break;
    
    case PostYieldMode::EXP_HARDENING:
      dr_dk = _residual_strength * _kslope * std::exp(-_kslope * kappa);
      break;
    
    case PostYieldMode::LINEAR_HARDENING:
      dr_dk = _kslope;
      break;
    
    case PostYieldMode::EXP_SOFTENING:
      if (kappa < _kmax)
        dr_dk = 0.0;
      else
      {
        Real kappa_shifted = kappa - _kmax;
        dr_dk = -(1.0 - _residual_strength) * _kslope * std::exp(-_kslope * kappa_shifted);
      }
      break;
    
    case PostYieldMode::PIECEWISE_SOFTENING:
      if (kappa < _kmax)
      {
        dr_dk = 0.0;  // No softening yet
      }
      else if (kappa < _kmin)
      {
        // Constant negative slope during linear drop
        dr_dk = -(1.0 - _residual_strength) / (_kmin - _kmax);
      }
      else
      {
        dr_dk = 0.0;  // Plateau - no more softening
      }
      break;
  }
  
  return dr_dk;
}

// ============================================================================
// VISCOSITY FUNCTION
// ============================================================================
Real
OrthotropicPlasticityStressUpdate::computeViscosity(
    Real dkappa, Real HI, Real dt) const
{
  // Rate-independent case
  if (_viscosity_mode == ViscosityMode::RATE_INDEPENDENT || dt < 1e-16 || HI < 1e-14)
    return 0.0;
  
  Real eta_dt = _eta / dt;
  Real visc = 0.0;
  
  switch (_viscosity_mode)
  {
    case ViscosityMode::RATE_INDEPENDENT:
      visc = 0.0;
      break;
    
    case ViscosityMode::LINEAR:
      // VISCFL = 1: -η/Δt × Δκ/HI
      visc = -eta_dt * dkappa / HI;
      break;
    
    case ViscosityMode::EXPONENTIAL:
      // VISCFL = 2: -log(1 + η/Δt × Δκ/HI)/m
      visc = -std::log(1.0 + eta_dt * dkappa / HI) / _m;
      break;
    
    case ViscosityMode::LOGARITHMIC:
      // VISCFL = 3: -(exp(η/Δt × Δκ/HI) - 1)/m
      visc = -(std::exp(eta_dt * dkappa / HI) - 1.0) / _m;
      break;
    
    case ViscosityMode::POLYNOMIAL:
      // VISCFL = 4: 0.5×m - sqrt((m²)/4 + η/Δt × Δκ/HI)
      visc = 0.5 * _m - std::sqrt(0.25 * _m * _m + eta_dt * dkappa / HI);
      break;
    
    case ViscosityMode::POWERLAW:
      // VISCFL = 5: -sign(Δκ) × (η/Δt × |Δκ|/HI)^(1/m)
      {
        Real sign = (dkappa > 0) ? 1.0 : -1.0;
        visc = -sign * std::pow(eta_dt * std::abs(dkappa) / HI, 1.0 / _m);
      }
      break;
  }
  
  return visc;
}

// ============================================================================
// VISCOSITY DERIVATIVES
// ============================================================================
void
OrthotropicPlasticityStressUpdate::computeViscosityDerivatives(
    Real dkappa, Real HI, const RankTwoTensor & dHI_ds, Real dHI_dk, Real dt,
    RankTwoTensor & dvisc_ds, Real & dvisc_dk) const
{
  dvisc_ds.zero();
  dvisc_dk = 0.0;
  
  // Rate-independent case
  if (_viscosity_mode == ViscosityMode::RATE_INDEPENDENT || dt < 1e-16 || HI < 1e-14)
    return;
  
  Real eta_dt = _eta / dt;
  
  switch (_viscosity_mode)
  {
    case ViscosityMode::RATE_INDEPENDENT:
      // Already zero
      break;
    
    case ViscosityMode::LINEAR:
      // ∂visc/∂σ = (1/HI²) × dHI/dσ × η/Δt × Δκ
      dvisc_ds = (1.0 / (HI * HI)) * dHI_ds * eta_dt * dkappa;
      // ∂visc/∂κ = -η/Δt / HI
      dvisc_dk = -eta_dt / HI;
      break;
    
    case ViscosityMode::EXPONENTIAL:
      {
        Real factor = eta_dt * dkappa / (_m * HI * HI) / (1.0 + eta_dt * dkappa / HI);
        dvisc_ds = factor * dHI_ds;
        dvisc_dk = -(eta_dt / _m) / (HI * (1.0 + eta_dt * dkappa / HI));
      }
      break;
    
    case ViscosityMode::LOGARITHMIC:
      {
        Real exp_term = std::exp(eta_dt * dkappa / HI);
        dvisc_ds = (eta_dt * dkappa / (_m * HI * HI)) * dHI_ds * exp_term;
        dvisc_dk = -(eta_dt / _m) / HI * exp_term;
      }
      break;
    
    case ViscosityMode::POLYNOMIAL:
      {
        Real denom = std::sqrt(0.25 * _m * _m + eta_dt * dkappa / HI);
        dvisc_ds = 0.5 * dHI_ds * eta_dt * std::abs(dkappa) / (HI * HI) / denom;
        dvisc_dk = -0.5 * eta_dt / HI / denom;
      }
      break;
    
    case ViscosityMode::POWERLAW:
      {
        Real sign = (dkappa > 0) ? 1.0 : -1.0;
        Real abs_dk = std::abs(dkappa);
        Real base = eta_dt * abs_dk / HI;
        Real power = std::pow(base, 1.0 / _m);
        
        dvisc_ds = -(1.0 / (_m * HI)) * dHI_ds * sign * power;
        dvisc_dk = -(1.0 / _m) * sign / HI * power;
      }
      break;
  }
}


// ============================================================================
// DAMAGE: D(κ) - Linear damage evolution
// ============================================================================
Real
OrthotropicPlasticityStressUpdate::computeDamage(Real kappa) const
{
  if (!_use_damage)
    return 0.0;
  
  // Linear damage evolution (UMAT style)
  // D = 0                           if κ < κ_crit
  // D = D_rate * (κ - κ_crit)       if κ ≥ κ_crit
  
  if (kappa < _damage_critical)
    return 0.0;
  
  Real D = _damage_rate * (kappa - _damage_critical);
  
  // Clamp to [0, 1) - damage cannot exceed 1
  if (D > 0.99)
    D = 0.99;
  
  return D;
}

// ============================================================================
// DAMAGE DERIVATIVE: dD/dκ
// ============================================================================
Real
OrthotropicPlasticityStressUpdate::computeDamageDerivative(Real kappa) const
{
  if (!_use_damage)
    return 0.0;
  
  // Derivative of linear damage
  // dD/dκ = 0         if κ < κ_crit
  // dD/dκ = D_rate    if κ ≥ κ_crit (constant)
  
  if (kappa < _damage_critical)
    return 0.0;
  
  // Check if we're at saturation
  Real D = _damage_rate * (kappa - _damage_critical);
  if (D > 0.99)
    return 0.0;  // No more damage growth
  
  return _damage_rate;
}

// ============================================================================
// YIELD FUNCTION: f = √(σ:F:σ) + f_lin·σ - r(κ)*σ₀ + viscosity
// ============================================================================
Real
OrthotropicPlasticityStressUpdate::computeYieldFunction(
    const RankTwoTensor & stress,
    Real kappa,
    Real delta_kappa,
    Real dt) const
{
  // Convert to Voigt notation
  std::vector<Real> s(6);
  s[0] = stress(0,0); s[1] = stress(1,1); s[2] = stress(2,2);
  s[3] = stress(1,2); s[4] = stress(0,2); s[5] = stress(0,1);
  


  // Quadric part: √(σ:F:σ)
  std::vector<Real> Fs(6, 0.0);
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 6; j++)
      Fs[i] += _F_matrix[i][j] * s[j];
  
  Real quadric = 0.0;
  for (unsigned int i = 0; i < 6; i++)
    quadric += s[i] * Fs[i];
  
  if (quadric < 0.0)
    quadric = 0.0;  // Regularize for tension-compression asymmetry
  
  Real phi_quadric = std::sqrt(quadric);
  
  // Linear part: f_lin·σ
  Real phi_linear = 0.0;
  for (unsigned int i = 0; i < 6; i++)
    phi_linear += _f_lin_vector[i] * s[i];
  
  // Total phi
  Real phi = phi_quadric + phi_linear;
  
  // Softening factor
  Real r = computeSofteningFactor(kappa);
  
  // Yield function: f = φ - r
  Real f = phi - r;
  
  // Add viscoplastic regularization (UMAT VISCY)
  RankTwoTensor grad = computeYieldGradient(stress);
  Real HI = grad.L2norm();
  Real visc = computeViscosity(delta_kappa, HI, dt);
  f += visc;

  // debug
   /*if (delta_kappa > 0 && _qp == 0) {
     Moose::out << "*** YIELD: stress(0,1) = " << stress(0,1) 
                << ", phi = " << phi << ", r = " << r << "\n";
   }*/

  // Debug output
  /*if (_t_step % 25 == 0 && _qp == 0) {  // Print every 25 steps, first qp
    mooseWarning("=== YIELD FUNCTION DEBUG ===\n",
                 "stress(0,1) = ", stress(0,1), " MPa\n",
                 "s[5] = ", s[5], " MPa\n",
                 "F_66 = ", _F_matrix[5][5], "\n",
                 "quadric term = ", s[5] * s[5] * _F_matrix[5][5], "\n",
                 "phi = ", std::sqrt(s[5] * s[5] * _F_matrix[5][5]), "\n",
                 "kappa = ", kappa);
  }*/

  // Comprehensive diagnostics
  if (_qp == 0 && (_t_step % 10 == 0 || delta_kappa > 0)) {
    Moose::out << "t=" << _t << " step=" << _t_step 
               << " | σ_xx=" << s[0] << " σ_yy=" << s[1] << " σ_zz=" << s[2]
               << " | σ_yz=" << s[3] << " σ_xz=" << s[4] << " σ_xy=" << s[5]
               << " | phi=" << phi << " r=" << r << " Δκ=" << delta_kappa << "\n";
  }

  return f;
}

// ============================================================================
// YIELD GRADIENT: ∂f/∂σ = (1/√(σ:F:σ)) * F·σ + f_lin
// ============================================================================
RankTwoTensor
OrthotropicPlasticityStressUpdate::computeYieldGradient(const RankTwoTensor & stress) const
{
  // Convert to Voigt
  std::vector<Real> s(6);
  s[0] = stress(0,0); s[1] = stress(1,1); s[2] = stress(2,2);
  s[3] = stress(1,2); s[4] = stress(0,2); s[5] = stress(0,1);
  
  // Compute F·σ
  std::vector<Real> Fs(6, 0.0);
  for (unsigned int i = 0; i < 6; i++)
    for (unsigned int j = 0; j < 6; j++)
      Fs[i] += _F_matrix[i][j] * s[j];
  
  // Compute σ:F:σ
  Real SFFS = 0.0;
  for (unsigned int i = 0; i < 6; i++)
    SFFS += s[i] * Fs[i];
  
  if (SFFS < 1e-16)
    SFFS = 1e-16;  // Regularize
  
  Real inv_sqrt = 1.0 / std::sqrt(SFFS);
  
  // ∂f/∂σ = (1/√(σ:F:σ)) * F·σ + f_lin
  std::vector<Real> grad(6);
  for (unsigned int i = 0; i < 6; i++)
    grad[i] = inv_sqrt * Fs[i] + _f_lin_vector[i];
  
  // Convert back to tensor
  RankTwoTensor grad_tensor;
  grad_tensor.zero();
  grad_tensor(0,0) = grad[0];
  grad_tensor(1,1) = grad[1];
  grad_tensor(2,2) = grad[2];
  grad_tensor(1,2) = grad_tensor(2,1) = grad[3];
  grad_tensor(0,2) = grad_tensor(2,0) = grad[4];
  grad_tensor(0,1) = grad_tensor(1,0) = grad[5];
  
  return grad_tensor;
}