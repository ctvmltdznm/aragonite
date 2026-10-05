// COMPLETE UMAT IMPLEMENTATION FOR MOOSE
// Orthotropic plasticity with damage, viscosity, and Primal CPPA
// Based on UMAT_QUADRIC_PRIMAL_Major_LINEAR_HARDENING.f
// For production use in research

#include "OrthotropicPlasticityStressUpdate.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "MooseException.h"
#include "MandelHelpers.h"

using namespace libMesh;

// ============================================================================
// MISSING MOOSE FRAMEWORK INTEGRATION
// Add these to the TOP of your .C file (after #includes, before other functions)
// ============================================================================

registerMooseObject("aragoniteApp", OrthotropicPlasticityStressUpdate);

using MaterialModelPresets::Model;

InputParameters
OrthotropicPlasticityStressUpdate::validParams()
{
  InputParameters params = StressUpdateBase::validParams();
  
  params.addClassDescription("Orthotropic plasticity with UMAT-style implementation. "
    "Supports both explicit directional strengths and fabric-based orthotropy "
    "(Schwiedrzik et al. 2013).");
  
  // ==========================================================================
  // MATERIAL MODEL FLAG (UMAT PROPS(1) + coral + legacy)
  // ==========================================================================
  params.addParam<MooseEnum>("material_model", MaterialModelPresets::modelEnum(),
                             MaterialModelPresets::sharedModelDoc());
  params.addParam<MooseEnum>("plastic_model", MaterialModelPresets::modelEnum(),
                             "Plastic model for THIS block only, overriding material_model. "
                             "Leave unset to follow material_model. " +
                             MaterialModelPresets::modelDoc() +
                             " Bone presets also set the UMAT post-yield and viscosity "
                             "defaults (initial_yield_ratio = 0.7, simple_softening or "
                             "exp_hardening, linear viscosity with eta = 1e-4).");
  params.addRangeCheckedParam<unsigned int>(
      "main_direction", 1, "main_direction>=1 & main_direction<=3",
      "Material axis of transverse isotropy (1,2,3), UMAT PROPS(6). Used by "
      "trabecular_ti, compact_ti and trabecular_fabric_ti.");

  // TI axial strengths (preset defaults for trabecular_ti / compact_ti)
  params.addParam<Real>("sigma_a_tension", "Axial tensile strength (TI models) [MPa]");
  params.addParam<Real>("sigma_a_compression", "Axial compressive strength (TI models) [MPa]");
  params.addParam<Real>("tau_a",
      "Shear strength of the planes containing the axis (TI models) [MPa]");
  params.addParam<Real>("zeta_a",
      "Axial interaction parameter (TI models). Referenced to the AXIAL F as in the "
      "UMAT; converted internally to the lower-index convention of _F_matrix.");

  params.addParam<bool>("debug_checks", false,
      "Run the in-code self-checks at startup: finite-difference check of the yield "
      "gradient, Hessian and softening derivative at stress states with shear, and a "
      "Mandel round-trip check of the elasticity tensor on the first stress update. "
      "Cheap (a few hundred yield evaluations), but prints; off in production.");

  // ==========================================================================
  // YIELD INPUT MODE (legacy only)
  // ==========================================================================
  params.addParam<bool>("use_fabric_scaling", false,
    "Enable fabric-based orthotropy (Schwiedrzik et al. 2013). When true, "
    "directional strengths are computed from base strengths, fabric eigenvalues, "
    "and density. When false (default), use explicit directional strengths.");
  
  // ==========================================================================
  // EXPLICIT MODE: Direct yield strengths (original interface)
  // ==========================================================================
  // Yield strengths - Tension
  params.addParam<Real>("sigma_xx_tension", "Yield strength in xx direction (tension) [MPa]");
  params.addParam<Real>("sigma_yy_tension", "Yield strength in yy direction (tension) [MPa]");
  params.addParam<Real>("sigma_zz_tension", "Yield strength in zz direction (tension) [MPa]");
  params.addParam<Real>("tau_xy_max", "Maximum shear stress in xy plane [MPa]");
  params.addParam<Real>("tau_xz_max", "Maximum shear stress in xz plane [MPa]");
  params.addParam<Real>("tau_yz_max", "Maximum shear stress in yz plane [MPa]");
  
  // Yield strengths - Compression (optional, defaults to tension)
  params.addParam<Real>("sigma_xx_compression", "Yield strength in xx direction (compression) [MPa]");
  params.addParam<Real>("sigma_yy_compression", "Yield strength in yy direction (compression) [MPa]");
  params.addParam<Real>("sigma_zz_compression", "Yield strength in zz direction (compression) [MPa]");
  
  // Yield surface coupling for explicit mode
  params.addParam<Real>("zeta12", 0.0, "X-Y coupling parameter (12 interaction)");
  params.addParam<Real>("zeta13", 0.0, "X-Z coupling parameter (13 interaction)");
  params.addParam<Real>("zeta23", 0.0, "Y-Z coupling parameter (23 interaction)");

  // Density scaling for explicit mode
  params.addParam<Real>("yield_density_exponent", 0.0,
    "Density exponent for yield scaling in EXPLICIT mode. "
    "σ_scaled = σ_input × ρ^p. Set to 0 to disable. "
    "Typical values: 1.5-2.0 for trabecular bone. "
    "Note: In fabric mode, use 'exponent_p' parameter instead.");
  
  // ==========================================================================
  // FABRIC MODE: Base strengths + fabric tensor (Schwiedrzik et al. 2013)
  // ==========================================================================
  // Base yield strengths (isotropic reference)
  params.addParam<Real>("sigma_0_tension", 
    "Base tensile yield strength σ₀⁺ [MPa] (fabric mode)");
  params.addParam<Real>("sigma_0_compression", 
    "Base compressive yield strength σ₀⁻ [MPa] (fabric mode). Defaults to sigma_0_tension.");
  params.addParam<Real>("tau_0", 
    "Base shear yield strength τ₀ [MPa] (fabric mode)");
  params.addParam<Real>("zeta_0", 0.22,
    "Base interaction parameter ζ₀ (fabric mode). Typical: 0.22 for bone.");
  
  // Fabric tensor eigenvalues (normalized so m1+m2+m3=3)
  params.addParam<Real>("fabric_m1", 1.0, 
    "Fabric tensor eigenvalue m₁ (fabric mode). Normalized: m₁+m₂+m₃=3");
  params.addParam<Real>("fabric_m2", 1.0, 
    "Fabric tensor eigenvalue m₂ (fabric mode)");
  params.addParam<Real>("fabric_m3", 1.0, 
    "Fabric tensor eigenvalue m₃ (fabric mode)");
  
  // Density and exponents
  params.addParam<Real>("density_rho", 1.0, 
    "Relative density ρ (BV/TV for bone, 0-1). Set to 1.0 for full density.");
  params.addParam<Real>("exponent_p", 1.6, 
    "Density exponent p (fabric mode). Typical: 1.6 for trabecular bone.");
  params.addParam<Real>("exponent_q", 1.0, 
    "Fabric exponent q (fabric mode). Typical: 1.0 for trabecular bone.");
  
  // Cortical bone correction (UMAT DELTA parameter for ρ > 0.5)
  params.addParam<Real>("delta_cortical", 1.0,
    "Cortical bone density correction factor (UMAT DELTA). "
    "Only affects ρ > 0.5. Set to 1.0 to disable.");

  // ==========================================================================
  // EULER ANGLES (both modes)
  // ==========================================================================
  params.addRequiredCoupledVar("euler_angle_1", "First Euler angle (phi1) in degrees");
  params.addRequiredCoupledVar("euler_angle_2", "Second Euler angle (Phi) in degrees");
  params.addRequiredCoupledVar("euler_angle_3", "Third Euler angle (phi2) in degrees");

  // params.addParam<Real>("kappa_start", 0.002, 
  //                        "Plastic strain at which softening begins");
  //params.addParam<Real>("softening_range", 0.015, 
  //                        "Plastic strain range over which linear softening occurs");
  
  // Viscoplasticity
  //params.addParam<bool>("use_viscoplasticity", false, "Enable viscoplastic regularization");
  params.addParam<Real>("eta", 1e-3, "Viscosity parameter [s/MPa]");
  
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
  params.addParam<Real>("kwidth", 8.0,
      "simple_softening (UMAT PYFL=2, KWIDTH): width of the post-yield peak. "
      "The peak sits at kappa = kmax and r decays back to initial_yield_ratio.");
  params.addParam<Real>("initial_yield_ratio", 1.0,
      "r(0), the UMAT's RDY. With r(0) < 1 the input strengths are ULTIMATE "
      "strengths and yield onsets at initial_yield_ratio * strength. Bone presets "
      "use 0.7. Only exp_hardening and simple_softening read it; the other "
      "post-yield modes start at r(0) = 1.");

  // Viscocity
  MooseEnum viscosity_mode("rate_independent linear exponential "
                          "logarithmic polynomial powerlaw",
                          "linear");
  params.addParam<MooseEnum>("viscosity_mode", viscosity_mode,
                            "Viscoplastic regularization mode");
  params.addParam<Real>("m", 0.001,
      "Viscosity shape parameter. Meaning depends on viscosity_mode: divisor "
      "for exponential and logarithmic, curvature for polynomial, RECIPROCAL "
      "EXPONENT for powerlaw (visc ~ x^(1/m)). The default 0.001 is sensible "
      "only for the divisor modes; powerlaw needs m >~ 1 or it underflows to "
      "rate-independent behaviour.");

  return params;
}

OrthotropicPlasticityStressUpdate::OrthotropicPlasticityStressUpdate(
    const InputParameters & parameters)
  : StressUpdateBase(parameters),
    // Yield input mode
    _use_fabric_scaling(getParam<bool>("use_fabric_scaling")),
    _model(resolveModel()),
    _main_direction(getParam<unsigned int>("main_direction")),
    _debug_checks(getParam<bool>("debug_checks")),
    
    // Initialize yield strengths to zero (will be set in constructor body)
    _sigma_xx_tension(0.0), _sigma_yy_tension(0.0), _sigma_zz_tension(0.0),
    _tau_xy_max(0.0), _tau_xz_max(0.0), _tau_yz_max(0.0),
    _sigma_xx_compression(0.0), _sigma_yy_compression(0.0), _sigma_zz_compression(0.0),
    _zeta12(0.0), _zeta13(0.0), _zeta23(0.0),
    
    // Fabric parameters (read if fabric mode, defaults otherwise)
    _sigma_0_tension(isParamValid("sigma_0_tension") ? getParam<Real>("sigma_0_tension") : 0.0),
    _sigma_0_compression(isParamValid("sigma_0_compression") ? 
                         getParam<Real>("sigma_0_compression") : _sigma_0_tension),
    _tau_0(isParamValid("tau_0") ? getParam<Real>("tau_0") : 0.0),
    _zeta_0(getParam<Real>("zeta_0")),
    _sigma_a_tension(0.0), _sigma_a_compression(0.0), _tau_a(0.0), _zeta_a(0.0),
    _fabric_m1(getParam<Real>("fabric_m1")),
    _fabric_m2(getParam<Real>("fabric_m2")),
    _fabric_m3(getParam<Real>("fabric_m3")),
    _density_rho(getParam<Real>("density_rho")),
    _exponent_p(getParam<Real>("exponent_p")),
    _exponent_q(getParam<Real>("exponent_q")),
    _delta_cortical(getParam<Real>("delta_cortical")),
    _yield_density_exponent(getParam<Real>("yield_density_exponent")),
    
    // Euler angles
    _euler_angle_1(coupledValue("euler_angle_1")),
    _euler_angle_2(coupledValue("euler_angle_2")),
    _euler_angle_3(coupledValue("euler_angle_3")),

    // Softening configuration
    _residual_strength(getParam<Real>("residual_strength")),
    _kslope(getParam<Real>("kslope")),
    _kmax(getParam<Real>("kmax")),
    _kmin(getParam<Real>("kmin")),
    _kwidth(getParam<Real>("kwidth")),
    _initial_yield_ratio(getParam<Real>("initial_yield_ratio")),
    
    // Viscoplasticity
    _eta(getParam<Real>("eta")),
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
    _yield_function(declareProperty<Real>("yield_function")),
    _return_mapping_stage(declareProperty<Real>("return_mapping_stage")),
    _return_mapping_iterations(declareProperty<Real>("return_mapping_iterations"))//,
{
  // ==========================================================================
  // COMPUTE EFFECTIVE YIELD PARAMETERS
  // ==========================================================================
  
  if (_model == Model::LEGACY)
    warnIgnored({"sigma_a_tension", "sigma_a_compression", "tau_a", "zeta_a", "main_direction"});

  if (_model != Model::LEGACY)
  {
    if (_use_fabric_scaling)
      paramError("use_fabric_scaling",
                 "use_fabric_scaling is legacy-only. With plastic_model set, the model "
                 "itself decides whether fabric scaling applies.");
    if (_model == Model::CORAL)
      applyCoralPreset();
    else
      applyBonePreset();
  }
  else if (_use_fabric_scaling)
  {
    // ========================================================================
    // FABRIC-BASED MODE (Schwiedrzik et al. 2013, Section 2.3)
    // ========================================================================
    
    // Validate fabric parameters
    if (!isParamValid("sigma_0_tension"))
      mooseError("Fabric mode requires 'sigma_0_tension' parameter");
    if (!isParamValid("tau_0"))
      mooseError("Fabric mode requires 'tau_0' parameter");
    
    // Check fabric tensor normalization (m1 + m2 + m3 should = 3)
    Real fabric_sum = _fabric_m1 + _fabric_m2 + _fabric_m3;
    if (std::abs(fabric_sum - 3.0) > 0.01)
      mooseWarning("Fabric eigenvalues should sum to 3.0, got ", fabric_sum, 
                   ". Consider normalizing.");
    
    // Compute density scaling: TSFU(ρ, p, δ)
    Real rho_p = computeTSFU(_density_rho, _exponent_p, _delta_cortical);
    
    // Fabric eigenvalue powers
    Real m1_2q = std::pow(_fabric_m1, 2.0 * _exponent_q);
    Real m2_2q = std::pow(_fabric_m2, 2.0 * _exponent_q);
    Real m3_2q = std::pow(_fabric_m3, 2.0 * _exponent_q);
    Real m1_q = std::pow(_fabric_m1, _exponent_q);
    Real m2_q = std::pow(_fabric_m2, _exponent_q);
    Real m3_q = std::pow(_fabric_m3, _exponent_q);
    
    // Directional yield strengths (Eq. 43 in paper)
    // σ⁺ᵢᵢ = σ⁺₀ · ρᵖ · mᵢ^(2q)
    _sigma_xx_tension = _sigma_0_tension * rho_p * m1_2q;
    _sigma_yy_tension = _sigma_0_tension * rho_p * m2_2q;
    _sigma_zz_tension = _sigma_0_tension * rho_p * m3_2q;
    
    _sigma_xx_compression = _sigma_0_compression * rho_p * m1_2q;
    _sigma_yy_compression = _sigma_0_compression * rho_p * m2_2q;
    _sigma_zz_compression = _sigma_0_compression * rho_p * m3_2q;
    
    // Shear yield strengths: τᵢⱼ = τ₀ · ρᵖ · mᵢ^q · mⱼ^q
    _tau_xy_max = _tau_0 * rho_p * m1_q * m2_q;  // 12 plane
    _tau_xz_max = _tau_0 * rho_p * m1_q * m3_q;  // 13 plane  
    _tau_yz_max = _tau_0 * rho_p * m2_q * m3_q;  // 23 plane
    
    // Interaction parameters scaled by fabric (Eq. 44)
    // ζᵢⱼ = ζ₀ · (mᵢ/mⱼ)^(2q)
    _zeta12 = _zeta_0 * std::pow(_fabric_m1 / _fabric_m2, 2.0 * _exponent_q);
    _zeta13 = _zeta_0 * std::pow(_fabric_m1 / _fabric_m3, 2.0 * _exponent_q);
    _zeta23 = _zeta_0 * std::pow(_fabric_m2 / _fabric_m3, 2.0 * _exponent_q);
    
    // Report computed values
    Moose::out << "\n=== FABRIC-BASED YIELD PARAMETERS (Schwiedrzik et al. 2013) ===\n";
    Moose::out << "Input: σ₀⁺=" << _sigma_0_tension << " σ₀⁻=" << _sigma_0_compression 
               << " τ₀=" << _tau_0 << " ζ₀=" << _zeta_0 << "\n";
    Moose::out << "Fabric: m₁=" << _fabric_m1 << " m₂=" << _fabric_m2 
               << " m₃=" << _fabric_m3 << " (sum=" << fabric_sum << ")\n";
    Moose::out << "Density: ρ=" << _density_rho << " p=" << _exponent_p 
               << " q=" << _exponent_q << " δ=" << _delta_cortical << "\n";
    Moose::out << "Computed strengths:\n";
    Moose::out << "  σ_xx: +" << _sigma_xx_tension << " / -" << _sigma_xx_compression << "\n";
    Moose::out << "  σ_yy: +" << _sigma_yy_tension << " / -" << _sigma_yy_compression << "\n";
    Moose::out << "  σ_zz: +" << _sigma_zz_tension << " / -" << _sigma_zz_compression << "\n";
    Moose::out << "  τ_xy=" << _tau_xy_max << " τ_xz=" << _tau_xz_max 
               << " τ_yz=" << _tau_yz_max << "\n";
    Moose::out << "  ζ₁₂=" << _zeta12 << " ζ₁₃=" << _zeta13 << " ζ₂₃=" << _zeta23 << "\n";
    Moose::out << "===============================================================\n\n";
  }
  else
  {
    // ========================================================================
    // EXPLICIT MODE (original interface - backward compatible)
    // ========================================================================
    // legacy only: reached through the else-branch of the model dispatch above
    // Validate that explicit parameters are provided
    if (!isParamValid("sigma_xx_tension"))
      mooseError("Explicit mode requires 'sigma_xx_tension' parameter. "
                 "Set 'use_fabric_scaling = true' for fabric-based mode.");
    if (!isParamValid("sigma_yy_tension"))
      mooseError("Explicit mode requires 'sigma_yy_tension' parameter");
    if (!isParamValid("sigma_zz_tension"))
      mooseError("Explicit mode requires 'sigma_zz_tension' parameter");
    if (!isParamValid("tau_xy_max"))
      mooseError("Explicit mode requires 'tau_xy_max' parameter");
    if (!isParamValid("tau_xz_max"))
      mooseError("Explicit mode requires 'tau_xz_max' parameter");
    if (!isParamValid("tau_yz_max"))
      mooseError("Explicit mode requires 'tau_yz_max' parameter");

    // Read explicit yield strengths
    _sigma_xx_tension = getParam<Real>("sigma_xx_tension");
    _sigma_yy_tension = getParam<Real>("sigma_yy_tension");
    _sigma_zz_tension = getParam<Real>("sigma_zz_tension");
    _tau_xy_max = getParam<Real>("tau_xy_max");
    _tau_xz_max = getParam<Real>("tau_xz_max");
    _tau_yz_max = getParam<Real>("tau_yz_max");

    // Compression strengths (default to tension if not provided)
    _sigma_xx_compression = isParamValid("sigma_xx_compression") ? 
                            getParam<Real>("sigma_xx_compression") : _sigma_xx_tension;
    _sigma_yy_compression = isParamValid("sigma_yy_compression") ? 
                            getParam<Real>("sigma_yy_compression") : _sigma_yy_tension;
    _sigma_zz_compression = isParamValid("sigma_zz_compression") ? 
                            getParam<Real>("sigma_zz_compression") : _sigma_zz_tension;

    // Yield surface coupling (explicit values)
    _zeta12 = getParam<Real>("zeta12");
    _zeta13 = getParam<Real>("zeta13");
    _zeta23 = getParam<Real>("zeta23");

    // ========================================================================
    // OPTIONAL DENSITY SCALING FOR EXPLICIT MODE
    // ========================================================================
    if (_yield_density_exponent > 0.0 && _density_rho < 1.0)
    {
      Real rho_p = std::pow(_density_rho, _yield_density_exponent);
      
      // Scale all yield strengths by ρ^p
      _sigma_xx_tension *= rho_p;
      _sigma_yy_tension *= rho_p;
      _sigma_zz_tension *= rho_p;
      _sigma_xx_compression *= rho_p;
      _sigma_yy_compression *= rho_p;
      _sigma_zz_compression *= rho_p;
      _tau_xy_max *= rho_p;
      _tau_xz_max *= rho_p;
      _tau_yz_max *= rho_p;
  
      Moose::out << "\n=== DENSITY-SCALED EXPLICIT YIELD PARAMETERS ===" << std::endl;
      Moose::out << "Relative density: ρ = " << _density_rho << std::endl;
      Moose::out << "Yield exponent: p = " << _yield_density_exponent << std::endl;
      Moose::out << "Scale factor: ρ^p = " << rho_p << std::endl;
      Moose::out << "Scaled yield strengths:" << std::endl;
      Moose::out << "  σ_xx: +" << _sigma_xx_tension << " / -" << _sigma_xx_compression << " MPa" << std::endl;
      Moose::out << "  σ_yy: +" << _sigma_yy_tension << " / -" << _sigma_yy_compression << " MPa" << std::endl;
      Moose::out << "  σ_zz: +" << _sigma_zz_tension << " / -" << _sigma_zz_compression << " MPa" << std::endl;
      Moose::out << "  τ_xy=" << _tau_xy_max << " τ_xz=" << _tau_xz_max 
                 << " τ_yz=" << _tau_yz_max << " MPa" << std::endl;
      Moose::out << "=================================================\n" << std::endl;
    }
  }
  // Post-yield mode, viscosity mode and their scalars. For a preset model the
  // defaults come from MaterialModelPresets.h unless the user set them.
  resolvePostYieldAndViscosity();

  // Initialize F matrix and f_lin vector (6x6 and 6x1)
  _F_matrix.resize(6, std::vector<Real>(6, 0.0));
  _f_lin_vector.resize(6, 0.0);
  
  // Build quadric yield surface from yield strengths
  // PLAIN COMPONENT 6-vector, NOT Mandel: _F_matrix acts on [s_xx s_yy s_zz
  // s_yz s_xz s_xy] with F[5][5] = 1/tau_xy^2, so SFFS is a plain 6-sum and
  // comes out correct. The conversion of the RESULT to a tensor is where the
  // sqrt(2)/2 factors appear (see below). Do not replace with the Mandel
  // helpers.
  // Based on Hill-type criterion with tension/compression asymmetry
  
  // Compute F matrix components
  Real sigma_t_xx = _sigma_xx_tension;
  Real sigma_t_yy = _sigma_yy_tension;
  Real sigma_t_zz = _sigma_zz_tension;
  Real sigma_c_xx = _sigma_xx_compression;
  Real sigma_c_yy = _sigma_yy_compression;
  Real sigma_c_zz = _sigma_zz_compression;
  
  // F matrix diagonal terms from yield condition first principles
  // For normal stresses: F_ii = [(σ_t + σ_c)/(2 σ_t σ_c)]²
  // For shear stresses: F_ij = 1/τ²
  // Derivation: From φ = √(σ:F:σ) + f·σ = 1 at yield
  
  Real F_xx_sqrt = (sigma_t_xx + sigma_c_xx) / (2.0 * sigma_t_xx * sigma_c_xx);
  Real F_yy_sqrt = (sigma_t_yy + sigma_c_yy) / (2.0 * sigma_t_yy * sigma_c_yy);
  Real F_zz_sqrt = (sigma_t_zz + sigma_c_zz) / (2.0 * sigma_t_zz * sigma_c_zz);
  
  // Normal stress terms (squared)
  _F_matrix[0][0] = F_xx_sqrt * F_xx_sqrt;
  _F_matrix[1][1] = F_yy_sqrt * F_yy_sqrt;
  _F_matrix[2][2] = F_zz_sqrt * F_zz_sqrt;
  
  // Shear terms
  _F_matrix[3][3] = 1.0 / (_tau_yz_max * _tau_yz_max);
  _F_matrix[4][4] = 1.0 / (_tau_xz_max * _tau_xz_max);
  _F_matrix[5][5] = 1.0 / (_tau_xy_max * _tau_xy_max);
  
  // Off-diagonal coupling: Schwiedrzik, Wolfram & Zysset (2013), Eq. 47 (general orthotropy).
  // F_ij = -zeta_ij * F_ii, reference index = lower-numbered index of the pair (i < j).
  // NOTE: this is NOT symmetric in i,j as written -- F_12 uses F_11, F_13 uses F_11,
  // F_23 uses F_22. This matches Eq. 55/56 convexity bounds exactly; it is only
  // numerically equivalent to a sqrt(F_ii*F_jj) form for the fabric-power-law special
  // case (Eq. 43-44), which is why fabric mode below is unaffected by this change.
  if (std::abs(_zeta12) > 1e-12 || std::abs(_zeta13) > 1e-12 || std::abs(_zeta23) > 1e-12) {

    // X-Y coupling (1-2), referenced to F_11
    _F_matrix[0][1] = -_zeta12 * _F_matrix[0][0];
    _F_matrix[1][0] = _F_matrix[0][1];  // Symmetric

    // X-Z coupling (1-3), referenced to F_11
    _F_matrix[0][2] = -_zeta13 * _F_matrix[0][0];
    _F_matrix[2][0] = _F_matrix[0][2];  // Symmetric

    // Y-Z coupling (2-3), referenced to F_22
    _F_matrix[1][2] = -_zeta23 * _F_matrix[1][1];
    _F_matrix[2][1] = _F_matrix[1][2];  // Symmetric
  }

  checkYieldSurfaceConvexity();

  // Linear term for tension/compression asymmetry
  _f_lin_vector[0] = (sigma_c_xx - sigma_t_xx) / (2.0 * sigma_c_xx * sigma_t_xx);  // xx
  _f_lin_vector[1] = (sigma_c_yy - sigma_t_yy) / (2.0 * sigma_c_yy * sigma_t_yy);  // yy
  _f_lin_vector[2] = (sigma_c_zz - sigma_t_zz) / (2.0 * sigma_c_zz * sigma_t_zz);  // zz
  // Shear components are zero (symmetric in shear)
  
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

  // Self-checks last: they evaluate the finished yield surface (_F_matrix AND
  // _f_lin_vector) and the resolved post-yield law.
  if (_debug_checks)
    runDebugChecks();
}

OrthotropicPlasticityStressUpdate::~OrthotropicPlasticityStressUpdate()
{
  // Deliberately empty. An earlier version printed the fallback totals here and
  // the line never reached the log: MOOSE destroys the material objects after
  // the output system is gone. Verified on the P runs -- five gated
  // "switching to Primal CPPA" warnings, no totals line. The report lives in
  // timestepSetup() instead.
}

void
OrthotropicPlasticityStressUpdate::timestepSetup()
{
  // Cumulative fallback totals, printed only when they have changed since the
  // last report. The per-occurrence warnings are gated to the first five, so
  // counting warning lines undercounts any run that falls back more than that.
  //
  // timestepSetup() runs at the START of a step, so the numbers cover
  // everything through the previous step; whatever happens in the final step is
  // not reported. That is acceptable for a diagnostic, and it is the only hook
  // a MaterialBase has that is guaranteed to reach the log.
  if (_primal_fallbacks != _reported_primal || _tangent_fallbacks != _reported_tangent)
  {
    _reported_primal = _primal_fallbacks;
    _reported_tangent = _tangent_fallbacks;
    Moose::out << "OrthotropicPlasticityStressUpdate totals: Newton->Primal "
               << _primal_fallbacks << ", C_ep fallbacks " << _tangent_fallbacks
               << " (cumulative, through step " << (_t_step > 0 ? _t_step - 1 : 0) << ")\n";
  }
}


void
OrthotropicPlasticityStressUpdate::checkYieldSurfaceConvexity() const
{
  // sqrt(s:F:s) is convex iff F is positive semi-definite. F is block diagonal
  // (normal 3x3 block; shear block diagonal and positive by construction), so
  // the test reduces to the normal block
  //   [ F11          -zeta12 F11   -zeta13 F11 ]
  //   [ -zeta12 F11   F22          -zeta23 F22 ]
  //   [ -zeta13 F11  -zeta23 F22    F33        ]
  // PSD <=> every principal minor >= 0 (Sylvester, semi-definite form).
  //
  // This replaces the earlier bounds |zeta_ij| <= F_jj/F_ii (Eq. 55) plus the
  // Eq. 56 cubic. Those are not the PSD conditions of THIS matrix: the exact
  // 2x2 bound is |zeta_ij| <= sqrt(F_jj/F_ii), and the cubic differs from the
  // determinant of the block in its squares and in the sign of the cross term.
  // They are also not invariant under relabelling the axes, which is exactly
  // what main_direction does: they rejected the convex UMAT compact_ti surface
  // for main_direction = 2 and 3. Verified against eigenvalues for every
  // preset and 20000 random zeta triples (verify_presets.py).
  const Real F11 = _F_matrix[0][0], F22 = _F_matrix[1][1], F33 = _F_matrix[2][2];
  const Real F12 = _F_matrix[0][1], F13 = _F_matrix[0][2], F23 = _F_matrix[1][2];
  const Real scale = std::max({F11, F22, F33});
  const Real tol = 1e-10;

  if (F11 <= 0.0 || F22 <= 0.0 || F33 <= 0.0)
    mooseError("OrthotropicPlasticityStressUpdate: non-positive normal quadric term "
               "(F11, F22, F33 = ", F11, ", ", F22, ", ", F33, ")");

  const Real m12 = F11 * F22 - F12 * F12;
  const Real m13 = F11 * F33 - F13 * F13;
  const Real m23 = F22 * F33 - F23 * F23;
  const Real det = F11 * m23 - F12 * (F12 * F33 - F23 * F13) + F13 * (F12 * F23 - F22 * F13);

  if (m12 < -tol * scale * scale)
    mooseError("OrthotropicPlasticityStressUpdate: zeta12 = ", _zeta12,
               " gives a non-convex yield surface; requires |zeta12| <= sqrt(F22/F11) = ",
               std::sqrt(F22 / F11));
  if (m13 < -tol * scale * scale)
    mooseError("OrthotropicPlasticityStressUpdate: zeta13 = ", _zeta13,
               " gives a non-convex yield surface; requires |zeta13| <= sqrt(F33/F11) = ",
               std::sqrt(F33 / F11));
  if (m23 < -tol * scale * scale)
    mooseError("OrthotropicPlasticityStressUpdate: zeta23 = ", _zeta23,
               " gives a non-convex yield surface; requires |zeta23| <= sqrt(F33/F22) = ",
               std::sqrt(F33 / F22));
  if (det < -tol * scale * scale * scale)
    mooseError("OrthotropicPlasticityStressUpdate: zeta12/13/23 = ", _zeta12, ", ", _zeta13,
               ", ", _zeta23, " give a non-convex yield surface (determinant of the normal "
               "block = ", det, ")");
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

// Forward declaration — defined below after manualInvert7x7
static void manualInvert6x6(const std::vector<std::vector<Real>> & A,
                             std::vector<std::vector<Real>> & Ainv);


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
  //RankTwoTensor stress_trial_material = rotateToMaterial(stress_new, R);

  // Self-check (debug_checks): the Mandel helpers must reproduce invSymm().
  // Needs an elasticity tensor, so it cannot run in the constructor.
  if (_debug_checks && !_mandel_checked)
  {
    const Real err = mandelRoundTripError(elasticity_tensor_material);
    Moose::out << "Mandel round-trip check: max rel err " << err
               << (err < 1e-10 ? "  OK" : "  FAILED") << "\n";
    if (err >= 1e-10)
      mooseError("Mandel round-trip check failed (", err,
                 "): rankFourToMandel/invert/mandelToRankFour does not reproduce invSymm()");
    _mandel_checked = true;
  }

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
    if (_primal_fallbacks++ < 5)
      mooseWarning("Newton failed at qp=", _qp, " (elem ", _current_elem->id(),
                   "), switching to Primal CPPA");
    iterations = 0;
    success = performPrimalCPP(stress_trial_material, C_inv_material,
                               kappa_old, damage_old, dt,
                               stress_return, delta_kappa,
                               kappa_new, damage_new, iterations);
  }
  
  if (!success) {
    // MooseException instead of mooseError: MOOSE catches this, marks the solve
    // failed and cuts the time step, rather than aborting the job. This is the
    // equivalent of the UMAT's PNEWDT = 0.5 (UMAT line 1855).
    throw MooseException("OrthotropicPlasticityStressUpdate: return mapping failed at qp=",
                         _qp, " (elem ", _current_elem->id(),
                         "), |stress_trial| = ", stress_trial_material.L2norm(),
                         ", kappa_old = ", kappa_old, ", dt = ", dt);
  }

  // Rotate stress back to global coordinates - for constistency with angles in elastic
  stress_new = rotateToGlobal(stress_return, R);

  // without rotation, stay in material frame
  // stress_new = stress_return;
  
  // Compute inelastic strain increment in material coordinates
  //RankTwoTensor plastic_strain_increment = delta_kappa * computeYieldGradient(stress_return);

  // Normalize at final stress
  RankTwoTensor DSY_final = computeYieldGradient(stress_return);
  Real HI_final = DSY_final.L2norm();
  RankTwoTensor NP_final = DSY_final / HI_final;
  RankTwoTensor plastic_strain_increment = delta_kappa * NP_final;  // Normalized!
  
  // Rotate plastic strain to global coordinates
  inelastic_strain_increment = rotateToGlobal(plastic_strain_increment, R);
  // without rotation
  //inelastic_strain_increment = plastic_strain_increment;
  
  // Update state variables
  _equivalent_plastic_strain[_qp] = kappa_new;
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + plastic_strain_increment;
  _damage[_qp] = damage_new;
  _return_mapping_iterations[_qp] = iterations;

  // ── Consistent elastoplastic tangent (UMAT lines 1840-1848) ────────────────
  // Computed in the same call in which it is requested (see
  // ComputeMultipleInelasticStressBase::computeAdmissibleState). All 6x6/6x1
  // algebra is MANDEL, so matvec, inversion, dot and norm are the tensor ops.
  if (compute_full_tangent_operator)
  {
    if (delta_kappa > 0.0)
    {
      const Real dr_new = computeSofteningDerivative(kappa_new);
      const Real D_new  = damage_new;

      const RankFourTensor DDSY_final = computeYieldHessian(stress_return);
      const RankTwoTensor  DHDS_final = (DDSY_final * DSY_final) / HI_final;

      // Consistency row, built exactly as in the return map: yield gradient
      // plus viscous corrections. On a perfect-plasticity plateau (dr = 0)
      // the viscous term is the only resistance to plastic flow; without it
      // the tangent is singular along the flow direction.
      RankTwoTensor DYDS_final = DSY_final;
      Real          DYDK_final = -dr_new;
      if (_viscosity_mode != ViscosityMode::RATE_INDEPENDENT && dt > 1e-16 && HI_final > 1e-14)
      {
        RankTwoTensor dvisc_ds;
        Real dvisc_dk;
        computeViscosityDerivatives(delta_kappa, HI_final, DHDS_final, 0.0, dt,
                                    dvisc_ds, dvisc_dk);
        DYDS_final += dvisc_ds;
        DYDK_final += dvisc_dk;
      }

      const RankFourTensor DNPDS_final =
          (DDSY_final * HI_final - dyadicProduct(DSY_final, DHDS_final)) / (HI_final * HI_final);
      const RankFourTensor DRRDS_final =
          -C_inv_material / (1.0 - D_new) - DNPDS_final * delta_kappa;

      // Damage contributions. D = D(kappa_old + delta_kappa), so both the
      // 1/(1-D) factor and the (D-D0)/(1-D) term depend on delta_kappa:
      //   dRR/d(dk)   -= dD/(1-D)^2 * [ C^-1(sigma - sigma_tr) + (1-D0)*C^-1 sigma_tr ]
      //   dRR/dsigma_tr = C^-1 * (1 - D + D0)/(1 - D)   ->  prefactor g on C_ep
      // Gated exactly as in the return map so the tangent matches the
      // residual it differentiates.
      // UNTESTED: use_damage = false in every test run so far. Verify with a
      // finite-difference check of C_ep at a state with D > 0 before trusting
      // it (perturb strain_increment, re-run the return map, compare columns).
      RankTwoTensor DRRDK_final = -NP_final;
      Real g_damage = 1.0;
      if (_use_damage && D_new > 1e-6 && D_new < 0.99)
      {
        const Real dD_new = computeDamageDerivative(kappa_new);
        const Real omd = 1.0 - D_new;
        const RankTwoTensor A = C_inv_material * (stress_return - stress_trial_material);
        const RankTwoTensor B = C_inv_material * stress_trial_material;
        DRRDK_final -= (dD_new / (omd * omd)) * (A + (1.0 - damage_old) * B);
        g_damage = (1.0 - D_new + damage_old) / omd;
      }

      std::vector<std::vector<Real>> neg_DRRDS_m;
      std::vector<std::vector<Real>> SSSA_m(6, std::vector<Real>(6, 0.0));
      rankFourToMandel(-DRRDS_final, neg_DRRDS_m);
      const bool inverted = invertMatrix6x6(neg_DRRDS_m, SSSA_m);

      std::vector<Real> DRRDK_m(6), DYDS_m(6);
      tensorToMandel(DRRDK_final, DRRDK_m);
      tensorToMandel(DYDS_final,  DYDS_m);

      // dsigma = SSSA.deps + SSSA.DRRDK dDk,   DYDS:dsigma + DYDK dDk = 0
      // => C_ep = SSSA - (SSSA.DRRDK) x (SSSA^T.DYDS) / (DYDS.SSSA.DRRDK + DYDK)
      // SSSA is not symmetric (DNPDS is a one-sided projection), hence the
      // transpose. Uses DYDS directly rather than UMAT 1836's NP and
      // DYDK/|DYDS|, which is exact only when DYDS is parallel to NP (i.e.
      // rate-independent); with viscosity it is not.
      std::vector<Real> SSSA_dRRdk(6, 0.0), SSSAt_DYDS(6, 0.0);
      matVecMult6(SSSA_m, DRRDK_m, SSSA_dRRdk);
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          SSSAt_DYDS[b] += SSSA_m[a][b] * DYDS_m[a];

      const Real denom_a = dotProduct6(DYDS_m, SSSA_dRRdk);
      const Real denom   = denom_a + DYDK_final;

      // No symmetrisation: for rate-independent associated flow the exact
      // C_ep is symmetric by itself; with viscosity it genuinely is not.
      std::vector<std::vector<Real>> C_el_m;
      rankFourToMandel(elasticity_tensor_material, C_el_m);
      Real c_el_max = 0.0;
      for (int a = 0; a < 6; a++)
        for (int b = 0; b < 6; b++)
          c_el_max = std::max(c_el_max, std::abs(C_el_m[a][b]));

      const char * reason = nullptr;
      if (!inverted)
        reason = "singular -DRRDS";
      else if (!std::isfinite(denom) ||
               std::abs(denom) <= 1e-12 * (std::abs(denom_a) + std::abs(DYDK_final)))
        reason = "vanishing denominator";

      std::vector<std::vector<Real>> C_ep_m = SSSA_m;
      if (!reason)
        for (int a = 0; a < 6; a++)
          for (int b = 0; b < 6; b++)
            C_ep_m[a][b] -= SSSA_dRRdk[a] * SSSAt_DYDS[b] / denom;

      // dRR/dsigma_trial = C^-1 * (1 - D + D0)/(1 - D), so the whole tangent
      // carries that factor. g_damage is 1 unless use_damage is on.
      if (g_damage != 1.0)
        for (int a = 0; a < 6; a++)
          for (int b = 0; b < 6; b++)
            C_ep_m[a][b] *= g_damage;

      for (int a = 0; a < 6 && !reason; a++)
        for (int b = 0; b < 6 && !reason; b++)
        {
          if (!std::isfinite(C_ep_m[a][b]))
            reason = "non-finite entry";
          else if (std::abs(C_ep_m[a][b]) > 10.0 * c_el_max)
            reason = "entry exceeds 10x elastic stiffness";
        }

      if (!reason)
      {
        RankFourTensor C_ep_material;
        mandelToRankFour(C_ep_m, C_ep_material);
        tangent_operator = rotateElasticityTensor(C_ep_material, R);
      }
      else
      {
        tangent_operator = elasticity_tensor;
        if (_tangent_fallbacks++ < 5)
          mooseWarning("C_ep fallback to elastic (", reason, ") at elem ",
                       _current_elem->id(), " qp ", _qp,
                       ", delta_kappa = ", delta_kappa, ", denom = ", denom);
      }
    }
    else
    {
      // Elastic step. Must never be skipped: the caller's tangent storage is a
      // plain member, not per-qp, so an unwritten tangent reuses the last qp's.
      tangent_operator = elasticity_tensor;
    }
  }
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
      throw MooseException("OrthotropicPlasticityStressUpdate: singular 7x7 Jacobian in "
                           "return mapping");

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

// General 6×6 LU inversion — no symmetry assumed (UMAT: CALL MIGS(-DRRDS,6,SSSA))
// invSymm() would be WRONG for SSSA because DRRDS lacks major symmetry.
static void manualInvert6x6(const std::vector<std::vector<Real>> & A,
                             std::vector<std::vector<Real>> & Ainv)
{
  const int n = 6;
  std::vector<std::vector<Real>> mat = A;
  std::vector<std::vector<Real>> inv(n, std::vector<Real>(n, 0.0));
  for (int i = 0; i < n; i++) inv[i][i] = 1.0;
  for (int i = 0; i < n; i++) {
    int imax = i; Real amax = std::abs(mat[i][i]);
    for (int k = i+1; k < n; k++)
      if (std::abs(mat[k][i]) > amax) { amax = std::abs(mat[k][i]); imax = k; }
    if (amax < 1e-14) mooseError("manualInvert6x6: singular matrix in tangent");
    if (imax != i) { std::swap(mat[i], mat[imax]); std::swap(inv[i], inv[imax]); }
    for (int k = i+1; k < n; k++) {
      Real f = mat[k][i] / mat[i][i];
      for (int j = i; j < n; j++) mat[k][j] -= f * mat[i][j];
      for (int j = 0; j < n; j++) inv[k][j] -= f * inv[i][j];
    }
  }
  for (int i = n-1; i >= 0; i--) {
    for (int j = 0; j < n; j++) {
      inv[i][j] /= mat[i][i];
      for (int k = 0; k < i; k++) inv[k][j] -= mat[k][i] * inv[i][j];
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
  // PLAIN COMPONENT 6-vector
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
  // Same convention as computeYieldGradient: off-diagonal slots get half.
  Fs_tensor(1,2) = Fs_tensor(2,1) = 0.5 * Fs[3];
  Fs_tensor(0,2) = Fs_tensor(2,0) = 0.5 * Fs[4];
  Fs_tensor(0,1) = Fs_tensor(1,0) = 0.5 * Fs[5];
  
  // First term: -1/(SFFS)^1.5 * (Fσ ⊗ Fσ)
  RankFourTensor dyadic_term = dyadicProduct(Fs_tensor, Fs_tensor);
  dyadic_term *= inv_sffs_3_2;
  
  // Second term: 1/sqrt(SFFS) * F
  RankFourTensor F_tensor;
  F_tensor.zero();
  
  // 6-vector
  auto voigt_to_tensor = [](int v, int & i, int & j) {
    const int map[6][2] = {{0,0}, {1,1}, {2,2}, {1,2}, {0,2}, {0,1}};
    i = map[v][0]; j = map[v][1];
  };
  
  for (unsigned int v1 = 0; v1 < 6; v1++) {
    for (unsigned int v2 = 0; v2 < 6; v2++) {
      int i, j, k, l;
      voigt_to_tensor(v1, i, j);
      voigt_to_tensor(v2, k, l);
      
      // _F_matrix acts on plain 6-vector components. The equivalent
      // RankFourTensor satisfies sigma:F:sigma = s.F.s only if each shear
      // index divides by 2.
      const Real w1 = (v1 < 3) ? 1.0 : 2.0;
      const Real w2 = (v2 < 3) ? 1.0 : 2.0;
      Real val = _F_matrix[v1][v2] / (w1 * w2);
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
    _yield_function[_qp] = f_trial;
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
      _yield_function[_qp] = f;
      _return_mapping_stage[_qp] = 2;
      return true;
    }
    
    // ===== COMPUTE JACOBIAN (UMAT lines 1463-1472) =====
    
    // DHDS = 1/HI * DDSY * DSY
    RankTwoTensor DHDS = (DDSY * DSY) / HI;
    
    // DYDS = DSY + viscous correction
    RankTwoTensor DYDS = DSY;
    Real DYDK = -dr_i;  // UMAT: DRAD = -dr/dκ, but my dr_i = +dr/dκ, so use -dr_i
    //Real DYDK = dr_i; // OLD--wrong apparently

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
    
    // Condensation of
    //   DRRDS:dsigma + DRRDK*ddk = -RR
    //   DYDS :dsigma + DYDK *ddk = -f
    // with SSSA = (-DRRDS)^-1, giving dsigma = SSSA*(RR + DRRDK*ddk) and
    //   ddk = -(f + DYDS:(SSSA*RR)) / (DYDS:(SSSA*DRRDK) + DYDK)
    // The previous form divided f and DYDK by ||DYDS|| and contracted with NP,
    // which is exact only when DYDS is parallel to NP (rate-independent).
    // With any viscosity mode DYDS = DSY + dvisc_ds is not.
    const RankTwoTensor SSSA_RR    = SSSA * RR;
    const RankTwoTensor SSSA_DRRDK = SSSA * DRRDK;

    const Real denom_a   = DYDS.doubleContraction(SSSA_DRRDK);
    const Real numerator = -(f + DYDS.doubleContraction(SSSA_RR));
    const Real denominator = denom_a + DYDK;

    // Relative guard: removing the 1/||DYDS|| scaling changed the magnitude
    // of this quantity, so an absolute threshold is no longer meaningful.
    if (std::abs(denominator) <= 1e-14 * (std::abs(denom_a) + std::abs(DYDK))) {
      mooseWarning("Newton: singular Jacobian (vanishing consistency denominator)");
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
  if (_primal_fallbacks < 5)
    mooseWarning("Newton failed: iter=", iterations, ", NORMRR=", NORMRR, ", ABSY=", ABSY);
  return false;
}

// ============================================================================
// PRIMAL CPPA WITH LINE SEARCH (UMAT lines 1525-1856)
// ============================================================================

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
  _yield_function[_qp] = f;

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
  tensorToMandel(RR, RRI);  // First 6 components
  RRI[6] = f;  // 7th component
  
  // ===== PRIMAL CPPA LOOP =====
  while ((NORMRR > TOL || ABSY > TOL) && iterations < MAXITER) {
    iterations++;
    
    // Store current state for line search
    RankTwoTensor stress_old_iter = stress_i;
    Real dkappa_old_iter = dkappa_i;

    // Merit at the CURRENT iterate, before the step is applied. NORMRR and
    // ABSY are up to date here (set before the loop and at the end of each
    // iteration). Previously MKI and MK1 were both evaluated AFTER the update,
    // so MK1 == MKI identically: the Armijo test fired every iteration and the
    // step was damped every time. UMAT computes these at lines 1682 and 1702,
    // i.e. on either side of the update.
    const Real MKI = 0.5 * NORMRR * NORMRR + 0.5 * ABSY * ABSY;
    
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
    // Real DYDK = dr_i; // OLD--wrong apparently
    Real DYDK = -dr_i;
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
    
    // ===== BUILD 7x7 JACOBIAN SYSTEM =====
    
    // Convert
    std::vector<std::vector<Real>> DRRDS_voigt;
    rankFourToMandel(DRRDS, DRRDS_voigt);  // 6x6
    
    std::vector<Real> DRRDK_voigt, DYDS_voigt;
    tensorToMandel(DRRDK, DRRDK_voigt);  // 6x1
    tensorToMandel(DYDS, DYDS_voigt);    // 6x1
    
    std::vector<Real> RR_voigt;
    tensorToMandel(RR, RR_voigt);  // 6x1
    
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
    
    // Directional derivative of the merit along DDD, evaluated with the RRR
    // and JJJ of the CURRENT iterate — both are overwritten below.
    Real DMK = 0.0;
    if (FLAG == 1) {
      // Newton direction: DMK = RRR^T*J*(-J^-1*RRR) = -||RRR||^2 = -2*MKI
      DMK = -2.0 * MKI;
    }
    else {
      std::vector<Real> temp(7);
      matVecMult7(JJJ, DDD, temp);
      DMK = dotProduct7(RRR, temp);
    }

    // ===== EXTRACT STRESS AND KAPPA UPDATES =====
    
    RankTwoTensor DSS1;
    mandelToTensor(DDD, DSS1);  // First 6 components → stress tensor
    Real DDK1 = DDD[6];         // 7th component → kappa
    
    // ===== TENTATIVE UPDATE =====
    
    stress_i = stress_old_iter + DSS1;
    dkappa_i = dkappa_old_iter + DDK1;
    
    // Store state vectors for line search
    std::vector<Real> XX1(7), XXI(7);
    tensorToMandel(stress_i, XX1);
    XX1[6] = dkappa_i;
    tensorToMandel(stress_old_iter, XXI);
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
    tensorToMandel(RR, RR_voigt);
    for (int i = 0; i < 6; i++)
      RRR[i] = RR_voigt[i];
    RRR[6] = f;
    
    // ===== MERIT AT THE NEW POINT (UMAT line 1702) =====

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
      // XX1 - XXI = ALPHA * DDD, so RRR^T*JJJ*(XX1-XXI) = ALPHA * DMK,
      // with DMK already evaluated at the current iterate.
      UPLIM = MKI + BETA * ALPHA * DMK;
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
        tensorToMandel(stress_i, XX1);
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
          UPLIM = MKI + BETA * ALPHA * DMK;
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
      tensorToMandel(RR, RRI);
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
      _yield_function[_qp] = f;
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
      // r = r0 + (1 - r0 + h) (1 - exp(-kslope kappa)),  h = _residual_strength
      //   r0 = 1 : legacy form  1 + h (1 - exp(-kslope kappa))
      //   h  = 0 : UMAT PYFL=1  RDY + (1 - RDY)(1 - exp(-KSLOPE kappa))
      // Derivative below must stay the exact derivative of this line.
      r = _initial_yield_ratio +
          (1.0 - _initial_yield_ratio + _residual_strength) * (1.0 - std::exp(-_kslope * kappa));
      break;

    case PostYieldMode::SIMPLE_SOFTENING:
      // UMAT PYFL=2 (RADK, UMAT 2075-2078):
      //   r = RDY + (1-RDY) [ exp(-(k-kmax)^2 / (kwidth kmax^2))
      //                       - exp(-1/kwidth - kslope k) ]
      // r(0) = RDY exactly (the two exponentials cancel at kappa = 0), rises to
      // ~1 at kappa = kmax, then decays back to RDY.
      r = _initial_yield_ratio +
          (1.0 - _initial_yield_ratio) *
              (std::exp(-(kappa - _kmax) * (kappa - _kmax) / (_kwidth * _kmax * _kmax)) -
               std::exp(-1.0 / _kwidth - _kslope * kappa));
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
      dr_dk = (1.0 - _initial_yield_ratio + _residual_strength) * _kslope *
              std::exp(-_kslope * kappa);
      break;

    case PostYieldMode::SIMPLE_SOFTENING:
      // exact derivative of the SIMPLE_SOFTENING branch above; matches UMAT
      // DRADK (UMAT 2123-2127). Checked by FD with debug_checks = true.
      dr_dk = (1.0 - _initial_yield_ratio) *
              (-2.0 * (kappa - _kmax) / (_kwidth * _kmax * _kmax) *
                   std::exp(-(kappa - _kmax) * (kappa - _kmax) / (_kwidth * _kmax * _kmax)) +
               _kslope * std::exp(-1.0 / _kwidth - _kslope * kappa));
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
      // VISCFL = 2: -(exp(η/Δt × Δκ/HI) - 1)/m
      visc = -(std::exp(eta_dt * dkappa / HI) - 1.0) / _m;
      break;
    
    case ViscosityMode::LOGARITHMIC:
      // VISCFL = 3: -log(1 + η/Δt × Δκ/HI)/m
      visc = -std::log(1.0 + eta_dt * dkappa / HI) / _m;
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
        Real exp_term = std::exp(eta_dt * dkappa / HI);
        dvisc_ds = (eta_dt * dkappa / (_m * HI * HI)) * dHI_ds * exp_term;
        dvisc_dk = -(eta_dt / _m) / HI * exp_term;
      }
      break;
    
    case ViscosityMode::LOGARITHMIC:
      {
        Real factor = eta_dt * dkappa / (_m * HI * HI) / (1.0 + eta_dt * dkappa / HI);
        dvisc_ds = factor * dHI_ds;
        dvisc_dk = -(eta_dt / _m) / (HI * (1.0 + eta_dt * dkappa / HI));
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

        // d(visc)/dHI = +sign*(1/m)*power/HI   (UMAT VDYDS, VISCFL=5)
        dvisc_ds = (1.0 / (_m * HI)) * dHI_ds * sign * power;
        // d(visc)/d(dkappa) = -(1/m)*power/|dkappa|, independent of sign.
        // The UMAT writes this as (eta/dt)^(1/m) * (|dk|/HI)^(1/m - 1) / HI,
        // which is the same expression. Diverges as dkappa -> 0 when m > 1
        // (exponent 1/m < 1), so guard the division.
        dvisc_dk = (abs_dk > 1e-30) ? -(1.0 / _m) * power / abs_dk : 0.0;
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
  // 6-vector
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
/*  if (_qp == 0 && (_t_step % 100 == 0 || delta_kappa > 0)) {
    Moose::out << "t=" << _t << " step=" << _t_step 
               << " | σ_xx=" << s[0] << " σ_yy=" << s[1] << " σ_zz=" << s[2]
               << " | σ_yz=" << s[3] << " σ_xz=" << s[4] << " σ_xy=" << s[5]
               << " | phi=" << phi << " r=" << r << " Δκ=" << delta_kappa << "\n";
  }
*/
  return f;
}

// ============================================================================
// YIELD GRADIENT: ∂f/∂σ = (1/√(σ:F:σ)) * F·σ + f_lin
// ============================================================================
RankTwoTensor
OrthotropicPlasticityStressUpdate::computeYieldGradient(const RankTwoTensor & stress) const
{
  // PLAIN COMPONENT 6-vector, NOT Mandel: s[5] = sigma_12 and _F_matrix[5][5]
  // = 1/tau_xy^2, so s.F.s equals sigma:F:sigma with no sqrt(2) factors. The
  // conversion back to a tensor at the end is where the halving appears.
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
  // grad[p] = dY/ds_p with s_p the plain component sigma_ij. For symmetric
  // tensors dY = G:dsigma requires 2*G_ij = dY/ds_p for i != j, so each
  // off-diagonal slot receives half. Matches the UMAT, whose Mandel DSY is
  // the Mandel image of this tensor.
  grad_tensor(1,2) = grad_tensor(2,1) = 0.5 * grad[3];
  grad_tensor(0,2) = grad_tensor(2,0) = 0.5 * grad[4];
  grad_tensor(0,1) = grad_tensor(1,0) = 0.5 * grad[5];
  
  return grad_tensor;
}

// ============================================================================
// TSFU: Tissue Function for Density Scaling (UMAT lines 2036-2045)
// ============================================================================
// Computes ρ^exp with optional cortical bone correction for ρ > 0.5
// TSFU(ρ, exp, δ) = ρ^exp                                      if ρ ≤ 0.5
//                 = ρ^exp + (δ-1) * ((ρ-0.5)/0.5)^exp          if ρ > 0.5
// Setting δ = 1.0 disables the cortical correction.
// ============================================================================
Real
OrthotropicPlasticityStressUpdate::computeTSFU(Real rho, Real exponent, Real delta) const
{
  if (rho <= 0.0)
    return 0.0;
  
  if (rho <= 0.5)
  {
    // Standard power law scaling
    return std::pow(rho, exponent);
  }
  else
  {
    // Cortical bone correction (UMAT line 2042)
    // For dense bone (ρ > 0.5), add correction term
    Real base = std::pow(rho, exponent);
    Real correction = (delta - 1.0) * std::pow((rho - 0.5) / 0.5, exponent);
    return base + correction;
  }
}

// ============================================================================
// MATERIAL-MODEL PRESETS
// ============================================================================
// Every preset only fills the effective-strength slots. The quadric itself is
// assembled by the single block at the end of the constructor, in PLAIN
// COMPONENT form (F[5][5] = 1/tau_xy^2 acting on sigma_12). A preset must
// never write _F_matrix in Mandel form: the yield surface would still be
// right, but the gradient, and therefore the plastic flow direction, would be
// wrong by a factor of 2 in shear, and every uniaxial test would still pass.

MaterialModelPresets::Model
OrthotropicPlasticityStressUpdate::resolveModel() const
{
  // Precedence: the block-specific flag, then the shared one, then legacy.
  // MOOSE applies [GlobalParams] through InputParameters::applyParameter, which
  // only fills a parameter the block did not set and copies the "set by user"
  // flag across, so isParamSetByUser() is true for a GlobalParams value and a
  // per-block value always wins. That is exactly the precedence we want.
  const bool spec_set = isParamSetByUser("plastic_model");
  const bool shared_set = isParamSetByUser("material_model");
  const MooseEnum spec = getParam<MooseEnum>("plastic_model");
  const MooseEnum shared = getParam<MooseEnum>("material_model");

  if (spec_set)
  {
    if (shared_set && static_cast<int>(spec) != static_cast<int>(shared))
      mooseWarning("plastic_model = ", spec, " overrides material_model = ", shared,
                   " in this block, so the elastic and plastic responses use different "
                   "models. That is supported, but it is rarely intended: drop "
                   "plastic_model to follow material_model.");
    return static_cast<MaterialModelPresets::Model>(static_cast<int>(spec));
  }
  if (shared_set)
    return static_cast<MaterialModelPresets::Model>(static_cast<int>(shared));
  return MaterialModelPresets::Model::LEGACY;
}

Real
OrthotropicPlasticityStressUpdate::resolve(const std::string & name, Real preset) const
{
  return isParamSetByUser(name) ? getParam<Real>(name) : preset;
}

void
OrthotropicPlasticityStressUpdate::warnIgnored(const std::vector<std::string> & names) const
{
  for (const auto & n : names)
    if (isParamSetByUser(n))
      mooseWarning("Parameter '", n, "' is ignored for the plastic model ",
                   isParamSetByUser("plastic_model") ? getParam<MooseEnum>("plastic_model")
                                                     : getParam<MooseEnum>("material_model"));
}

void
OrthotropicPlasticityStressUpdate::applyBonePreset()
{
  using namespace MaterialModelPresets;
  const PlasticPreset pre = plasticPreset(_model);
  const MooseEnum model_name = isParamSetByUser("plastic_model")
                                   ? getParam<MooseEnum>("plastic_model")
                                   : getParam<MooseEnum>("material_model");

  // explicit-mode inputs never apply to a bone preset
  warnIgnored({"sigma_xx_tension", "sigma_yy_tension", "sigma_zz_tension",
               "sigma_xx_compression", "sigma_yy_compression", "sigma_zz_compression",
               "tau_xy_max", "tau_xz_max", "tau_yz_max", "zeta12", "zeta13", "zeta23",
               "yield_density_exponent"});
  if (isIso(_model))
    warnIgnored({"sigma_a_tension", "sigma_a_compression", "tau_a", "zeta_a", "tau_0",
                 "fabric_m1", "fabric_m2", "fabric_m3", "exponent_q", "main_direction"});
  else if (isTI(_model))
    warnIgnored({"tau_0", "fabric_m1", "fabric_m2", "fabric_m3", "exponent_q"});
  else if (_model == Model::TRABECULAR_FABRIC_ORTHO)
    warnIgnored({"sigma_a_tension", "sigma_a_compression", "tau_a", "zeta_a", "main_direction"});
  else // fabric TI
    warnIgnored({"sigma_a_tension", "sigma_a_compression", "tau_a", "zeta_a"});

  if (!isParamSetByUser("density_rho"))
    paramError("density_rho", "Bone plastic models require density_rho (BV/TV)");
  if (isFabric(_model) && !(isParamSetByUser("fabric_m1") && isParamSetByUser("fabric_m2") &&
                            isParamSetByUser("fabric_m3")))
    mooseError("plastic_model = ", model_name, " requires fabric_m1, fabric_m2, fabric_m3");

  _sigma_0_tension = resolve("sigma_0_tension", pre.s0p);
  _sigma_0_compression = resolve("sigma_0_compression", pre.s0n);
  _zeta_0 = resolve("zeta_0", pre.zeta0);
  _exponent_p = resolve("exponent_p", pre.p);
  _delta_cortical = resolve("delta_cortical", pre.delta);
  if (isFabric(_model))
  {
    _tau_0 = resolve("tau_0", pre.tau0);
    _exponent_q = resolve("exponent_q", pre.q);
  }
  if (isTI(_model))
  {
    _sigma_a_tension = resolve("sigma_a_tension", pre.sap);
    _sigma_a_compression = resolve("sigma_a_compression", pre.san);
    _tau_a = resolve("tau_a", pre.taua);
    _zeta_a = resolve("zeta_a", pre.zetaa);
  }

  // density scaling, UMAT TSFU(RHO,PP,DELTA) applied to every strength
  const Real t = computeTSFU(_density_rho, _exponent_p, _delta_cortical);
  auto S = [](Real sp, Real sn) { return (sp + sn) / (2.0 * sp * sn); };
  // shear strength of an ISOTROPIC plane: not an input in the UMAT, it follows
  // from the isotropy constraint (UMAT 413): TAUD0 = sqrt(0.5/S0^2/(1+ZETA0)).
  const Real tau_iso =
      1.0 / (S(_sigma_0_tension, _sigma_0_compression) * std::sqrt(2.0 * (1.0 + _zeta_0)));
  const unsigned int a = _main_direction - 1;
  // plane order as used by _F_matrix: 0 -> 12 (xy), 1 -> 13 (xz), 2 -> 23 (yz).
  // NOTE the UMAT's slot order is 4 = xy, 5 = xz, 6 = yz, i.e. slots 4 and 6
  // are interchanged relative to MOOSE (UMAT_TO_MOOSE_MAPPING.md section 1).
  // The presets are written in terms of PLANES, not slot numbers, so that
  // swap cannot leak in here.
  const unsigned int pl[3][2] = {{0, 1}, {0, 2}, {1, 2}};

  Real sT[3], sC[3], tau[3], zeta[3];

  if (isIso(_model))
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      sT[i] = _sigma_0_tension * t;
      sC[i] = _sigma_0_compression * t;
      tau[i] = tau_iso * t;
      zeta[i] = _zeta_0;
    }
  }
  else if (isTI(_model))
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      sT[i] = (i == a ? _sigma_a_tension : _sigma_0_tension) * t;
      sC[i] = (i == a ? _sigma_a_compression : _sigma_0_compression) * t;
    }
    for (unsigned int p = 0; p < 3; ++p)
    {
      const unsigned int i = pl[p][0], j = pl[p][1];
      const bool axial_plane = (i == a || j == a);
      tau[p] = (axial_plane ? _tau_a : tau_iso) * t;
      // The UMAT writes F_aj = -zeta_a * F_aa, i.e. referenced to the AXIAL
      // entry. _F_matrix uses F_ij = -zeta_ij * F_ii with i < j, so
      //   zeta_ij = zeta_a * F_aa / F_ii,
      // which is the identity only when i == a (main_direction = 1). Without
      // this conversion main_direction = 2 or 3 gives the wrong surface.
      if (axial_plane)
      {
        const Real Faa = std::pow(S(sT[a], sC[a]), 2);
        const Real Fii = std::pow(S(sT[i], sC[i]), 2);
        zeta[p] = _zeta_a * Faa / Fii;
      }
      else
        zeta[p] = _zeta_0;
    }
  }
  else // fabric models (Schwiedrzik et al. 2013, Eq. 43/44)
  {
    Real m[3] = {_fabric_m1, _fabric_m2, _fabric_m3};
    const Real fabric_sum = m[0] + m[1] + m[2];
    if (std::abs(fabric_sum - 3.0) > 0.01)
      mooseWarning("Fabric eigenvalues should sum to 3.0, got ", fabric_sum, ".");

    const bool ti = (_model == Model::TRABECULAR_FABRIC_TI);
    if (ti)
    {
      // UMAT averages the two transverse eigenvalues (UMAT 669-670)
      const unsigned int b = (a + 1) % 3, c = (a + 2) % 3;
      m[b] = m[c] = 0.5 * (m[b] + m[c]);
      _fabric_m1 = m[0]; _fabric_m2 = m[1]; _fabric_m3 = m[2];
    }
    const Real q = _exponent_q;
    for (unsigned int i = 0; i < 3; ++i)
    {
      sT[i] = _sigma_0_tension * t * std::pow(m[i], 2.0 * q);
      sC[i] = _sigma_0_compression * t * std::pow(m[i], 2.0 * q);
    }
    for (unsigned int p = 0; p < 3; ++p)
    {
      const unsigned int i = pl[p][0], j = pl[p][1];
      // fabric TI: the transverse plane must stay isotropic, so it uses the
      // derived tau_iso instead of tau_0 (UMAT 723-725)
      const bool transverse = ti && i != a && j != a;
      tau[p] = (transverse ? tau_iso : _tau_0) * t * std::pow(m[i] * m[j], q);
      zeta[p] = _zeta_0 * std::pow(m[i] / m[j], 2.0 * q);
    }
  }

  _sigma_xx_tension = sT[0]; _sigma_yy_tension = sT[1]; _sigma_zz_tension = sT[2];
  _sigma_xx_compression = sC[0]; _sigma_yy_compression = sC[1]; _sigma_zz_compression = sC[2];
  _tau_xy_max = tau[0]; _tau_xz_max = tau[1]; _tau_yz_max = tau[2];
  _zeta12 = zeta[0]; _zeta13 = zeta[1]; _zeta23 = zeta[2];

  Moose::out << "\n=== YIELD SURFACE: model = " << model_name << " ===\n"
             << "rho=" << _density_rho << " p=" << _exponent_p << " delta=" << _delta_cortical
             << " TSFU=" << t;
  if (usesMainDirection(_model))
    Moose::out << " main_direction=" << _main_direction;
  if (isFabric(_model))
    Moose::out << "\nm (used)=" << _fabric_m1 << " " << _fabric_m2 << " " << _fabric_m3
               << " q=" << _exponent_q << " tau_0=" << _tau_0;
  Moose::out << "\nsigma_0+=" << _sigma_0_tension << " sigma_0-=" << _sigma_0_compression
             << " zeta_0=" << _zeta_0;
  if (isTI(_model))
    Moose::out << " sigma_a+=" << _sigma_a_tension << " sigma_a-=" << _sigma_a_compression
               << " tau_a=" << _tau_a << " zeta_a=" << _zeta_a;
  Moose::out << "\nEffective strengths (material frame):\n"
             << "  sigma_xx: +" << _sigma_xx_tension << " / -" << _sigma_xx_compression << "\n"
             << "  sigma_yy: +" << _sigma_yy_tension << " / -" << _sigma_yy_compression << "\n"
             << "  sigma_zz: +" << _sigma_zz_tension << " / -" << _sigma_zz_compression << "\n"
             << "  tau_xy=" << _tau_xy_max << " tau_xz=" << _tau_xz_max
             << " tau_yz=" << _tau_yz_max << "\n"
             << "  zeta12=" << _zeta12 << " zeta13=" << _zeta13 << " zeta23=" << _zeta23
             << "  (F_ij = -zeta_ij F_ii, i<j)\n";
}

void
OrthotropicPlasticityStressUpdate::applyCoralPreset()
{
  const MaterialModelPresets::CoralStrengths d;

  warnIgnored({"density_rho", "fabric_m1", "fabric_m2", "fabric_m3", "exponent_p", "exponent_q",
               "delta_cortical", "sigma_0_tension", "sigma_0_compression", "tau_0", "zeta_0",
               "sigma_a_tension", "sigma_a_compression", "tau_a", "zeta_a", "main_direction",
               "yield_density_exponent"});

  _sigma_xx_tension = resolve("sigma_xx_tension", d.sxx);
  _sigma_yy_tension = resolve("sigma_yy_tension", d.syy);
  _sigma_zz_tension = resolve("sigma_zz_tension", d.szz);
  _tau_xy_max = resolve("tau_xy_max", d.txy);
  _tau_xz_max = resolve("tau_xz_max", d.txz);
  _tau_yz_max = resolve("tau_yz_max", d.tyz);
  _sigma_xx_compression = resolve("sigma_xx_compression", _sigma_xx_tension);
  _sigma_yy_compression = resolve("sigma_yy_compression", _sigma_yy_tension);
  _sigma_zz_compression = resolve("sigma_zz_compression", _sigma_zz_tension);

  // No default: zeta = 0 is not acceptable and no calibrated coral value exists.
  for (const std::string n : {"zeta12", "zeta13", "zeta23"})
    if (!isParamSetByUser(n))
      paramError(n, "the coral plastic model requires zeta12, zeta13 and zeta23 to be set "
                    "explicitly: there is no calibrated coral value yet and zeta = 0 is "
                    "not an acceptable placeholder.");
  _zeta12 = getParam<Real>("zeta12");
  _zeta13 = getParam<Real>("zeta13");
  _zeta23 = getParam<Real>("zeta23");
  if (_zeta12 == 0.0 || _zeta13 == 0.0 || _zeta23 == 0.0)
    mooseWarning("coral plastic model with a zero zeta_ij: the normal stresses are "
                 "uncoupled in that pair.");

  Moose::out << "\n=== YIELD SURFACE: model = coral ===\n"
             << "  sigma_xx: +" << _sigma_xx_tension << " / -" << _sigma_xx_compression << "\n"
             << "  sigma_yy: +" << _sigma_yy_tension << " / -" << _sigma_yy_compression << "\n"
             << "  sigma_zz: +" << _sigma_zz_tension << " / -" << _sigma_zz_compression << "\n"
             << "  tau_xy=" << _tau_xy_max << " tau_xz=" << _tau_xz_max
             << " tau_yz=" << _tau_yz_max << "\n"
             << "  zeta12=" << _zeta12 << " zeta13=" << _zeta13 << " zeta23=" << _zeta23 << "\n"
             << "Note: coral RVEs normally also carry HomogenizedExponentialCZM on the "
                "grain/needle interfaces (see MATERIAL_MODELS.md).\n";
}

void
OrthotropicPlasticityStressUpdate::resolvePostYieldAndViscosity()
{
  MooseEnum py = getParam<MooseEnum>("postyield_mode");
  MooseEnum visc = getParam<MooseEnum>("viscosity_mode");

  if (_model != Model::LEGACY)
  {
    const MaterialModelPresets::PlasticPreset pre = MaterialModelPresets::plasticPreset(_model);

    // The preset scalars belong to the preset's post-yield LAW: kmax is a peak
    // position for simple_softening and a softening onset for exp_softening,
    // so they are only applied while that law is in use. Overriding
    // postyield_mode falls back to the plain parameter defaults. Same for
    // viscosity_mode and eta/m.
    if (!isParamSetByUser("postyield_mode"))
    {
      py = pre.postyield;
      _initial_yield_ratio = resolve("initial_yield_ratio", pre.rdy);
      _residual_strength = resolve("residual_strength", pre.residual);
      _kslope = resolve("kslope", pre.kslope);
      _kmax = resolve("kmax", pre.kmax);
      if (pre.kmin > 0.0)
        _kmin = resolve("kmin", pre.kmin);
      _kwidth = resolve("kwidth", pre.kwidth);
    }
    if (!isParamSetByUser("viscosity_mode"))
    {
      visc = pre.viscosity;
      _eta = resolve("eta", pre.eta);
      // m is UNUSED by linear (every bone preset and coral use linear); it is
      // carried only so an override of viscosity_mode alone lands on the UMAT
      // value. Only powerlaw (reciprocal exponent) and the divisor modes read it.
      _m = resolve("m", pre.m);
    }
  }

  if (py == "perfect")
    _postyield_mode = PostYieldMode::PERFECT_PLASTICITY;
  else if (py == "exp_hardening")
    _postyield_mode = PostYieldMode::EXP_HARDENING;
  else if (py == "linear_hardening")
    _postyield_mode = PostYieldMode::LINEAR_HARDENING;
  else if (py == "simple_softening")
    _postyield_mode = PostYieldMode::SIMPLE_SOFTENING;
  else if (py == "exp_softening")
    _postyield_mode = PostYieldMode::EXP_SOFTENING;
  else if (py == "piecewise_softening")
    _postyield_mode = PostYieldMode::PIECEWISE_SOFTENING;
  else
    mooseError("Unhandled postyield_mode '", py, "'");

  if (visc == "rate_independent")
    _viscosity_mode = ViscosityMode::RATE_INDEPENDENT;
  else if (visc == "linear")
    _viscosity_mode = ViscosityMode::LINEAR;
  else if (visc == "exponential")
    _viscosity_mode = ViscosityMode::EXPONENTIAL;
  else if (visc == "logarithmic")
    _viscosity_mode = ViscosityMode::LOGARITHMIC;
  else if (visc == "polynomial")
    _viscosity_mode = ViscosityMode::POLYNOMIAL;
  else if (visc == "powerlaw")
    _viscosity_mode = ViscosityMode::POWERLAW;
  else
    mooseError("Unhandled viscosity_mode '", visc, "'");

  if (_initial_yield_ratio <= 0.0 || _initial_yield_ratio > 1.0)
    paramError("initial_yield_ratio", "must be in (0, 1]");
  const bool uses_r0 = _postyield_mode == PostYieldMode::EXP_HARDENING ||
                       _postyield_mode == PostYieldMode::SIMPLE_SOFTENING;
  if (!uses_r0 && _initial_yield_ratio != 1.0)
    mooseWarning("initial_yield_ratio = ", _initial_yield_ratio,
                 " is ignored by postyield_mode = ", py, " (that law starts at r(0) = 1)");
  if (_postyield_mode == PostYieldMode::SIMPLE_SOFTENING && (_kmax <= 0.0 || _kwidth <= 0.0))
    mooseError("simple_softening requires kmax > 0 and kwidth > 0");
  if (_viscosity_mode == ViscosityMode::POWERLAW && _m < 0.5)
    mooseWarning("viscosity_mode = powerlaw with m = ", _m, ": visc ~ x^(1/m) underflows "
                 "to rate-independent behaviour for small m (see UMAT_TO_MOOSE_MAPPING.md 6.5)");

  Moose::out << "Post-yield: " << py << " r0=" << _initial_yield_ratio
             << " kslope=" << _kslope << " kmax=" << _kmax << " kmin=" << _kmin
             << " kwidth=" << _kwidth << " residual_strength=" << _residual_strength
             << "\nViscosity: " << visc << " eta=" << _eta << " m=" << _m
             << (_viscosity_mode == ViscosityMode::LINEAR ? "  (m unused by linear)" : "")
             << "\n==============================================================\n\n";
}

// ============================================================================
// SELF-CHECKS (debug_checks = true)
// ============================================================================
// Runs once per material object at construction, so every preset is checked
// with its own _F_matrix, at stress states WITH shear. A pure-normal state
// cannot see the shear-slot factor-of-2 class of bug (see
// UMAT_TO_MOOSE_MAPPING.md section 6.3), which is why the states below always
// carry shear.

void
OrthotropicPlasticityStressUpdate::fdCheckAt(const RankTwoTensor & stress,
                                             Real & grad_err, Real & hess_err) const
{
  // symmetric perturbation E^(ij) = (e_i x e_j + e_j x e_i)/2, so that
  // dphi/dh along E equals G:E = G_ij, and dG/dh along E equals H:E.
  const int IJ[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
  const Real scale = stress.L2norm();
  const Real h = 1e-6 * (scale > 0.0 ? scale : 1.0);

  const RankTwoTensor G = computeYieldGradient(stress);
  const RankFourTensor H = computeYieldHessian(stress);

  Real gmax = 0.0, hmax = 0.0, gdiff = 0.0, hdiff = 0.0;
  for (int p = 0; p < 6; p++)
  {
    const int i = IJ[p][0], j = IJ[p][1];
    RankTwoTensor E;
    E.zero();
    E(i, j) = E(j, i) = (i == j) ? 1.0 : 0.5;

    RankTwoTensor sp = stress + E * h, sm = stress - E * h;

    // kappa = 0, delta_kappa = 0, dt = 0: r(0) is a constant that cancels in
    // the difference, and computeViscosity returns 0 for dt < 1e-16, so this
    // differentiates exactly phi = sqrt(s:F:s) + f_lin.s.
    const Real fp = computeYieldFunction(sp, 0.0, 0.0, 0.0);
    const Real fm = computeYieldFunction(sm, 0.0, 0.0, 0.0);
    const Real g_fd = (fp - fm) / (2.0 * h);
    gdiff = std::max(gdiff, std::abs(g_fd - G(i, j)));
    gmax = std::max(gmax, std::abs(G(i, j)));

    const RankTwoTensor Gp = computeYieldGradient(sp);
    const RankTwoTensor Gm = computeYieldGradient(sm);
    const RankTwoTensor HE = contractRankFourTwo(H, E);
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
      {
        const Real h_fd = (Gp(k, l) - Gm(k, l)) / (2.0 * h);
        hdiff = std::max(hdiff, std::abs(h_fd - HE(k, l)));
        hmax = std::max(hmax, std::abs(HE(k, l)));
      }
  }
  grad_err = (gmax > 0.0) ? gdiff / gmax : gdiff;
  hess_err = (hmax > 0.0) ? hdiff / hmax : hdiff;
}

void
OrthotropicPlasticityStressUpdate::runDebugChecks() const
{
  // Stress states scaled to the yield surface of THIS preset, all with shear.
  const Real s0 = 0.5 / std::sqrt(_F_matrix[0][0]);
  const Real t12 = 0.5 / std::sqrt(_F_matrix[5][5]);
  const Real t13 = 0.5 / std::sqrt(_F_matrix[4][4]);
  const Real t23 = 0.5 / std::sqrt(_F_matrix[3][3]);

  std::vector<RankTwoTensor> states(3);
  for (auto & S : states)
    S.zero();
  // mixed normal + full shear
  states[0](0, 0) = 0.6 * s0;  states[0](1, 1) = -0.4 * s0; states[0](2, 2) = 0.25 * s0;
  states[0](0, 1) = states[0](1, 0) = 0.35 * t12;
  states[0](0, 2) = states[0](2, 0) = 0.20 * t13;
  states[0](1, 2) = states[0](2, 1) = 0.30 * t23;
  // shear dominated
  states[1](0, 0) = 0.05 * s0;
  states[1](0, 1) = states[1](1, 0) = 0.70 * t12;
  states[1](1, 2) = states[1](2, 1) = 0.45 * t23;
  // compression + one shear
  states[2](0, 0) = -0.5 * s0; states[2](1, 1) = -0.3 * s0; states[2](2, 2) = -0.2 * s0;
  states[2](0, 2) = states[2](2, 0) = 0.40 * t13;

  Real gmax = 0.0, hmax = 0.0;
  for (const auto & S : states)
  {
    Real ge = 0.0, he = 0.0;
    fdCheckAt(S, ge, he);
    gmax = std::max(gmax, ge);
    hmax = std::max(hmax, he);
  }

  // Softening derivative. r(kappa) is only PIECEWISE smooth: exp_softening has
  // a derivative jump at kmax, and piecewise_softening at kmax and kmin. A
  // central difference straddling such a kink returns the mean of the two
  // one-sided derivatives, so at kappa = kmax it reports exactly half the true
  // value -- a relative error of 0.5, which is what the first version of this
  // check reported for coral (exp_softening, residual 0.7, kslope 30: analytic
  // -9, central difference -4.5). The derivative code was correct.
  //
  // So compare the two ONE-SIDED differences first. Where they agree the point
  // is smooth and the analytic derivative must match there; where they disagree
  // the point is a kink, which is reported and skipped rather than failed.
  Real rmax = 0.0;
  unsigned int n_kinks = 0;
  Real first_kink = 0.0;
  const Real k_ref = (_kmax > 0.0) ? _kmax : 0.01;
  for (const Real k : {0.2 * k_ref, 0.9 * k_ref, 1.0 * k_ref, 1.5 * k_ref, 4.0 * k_ref})
  {
    const Real h = 1e-6 * k_ref;
    const Real r0 = computeSofteningFactor(k);
    const Real fwd = (computeSofteningFactor(k + h) - r0) / h;
    const Real bwd = (r0 - computeSofteningFactor(k - h)) / h;
    const Real an = computeSofteningDerivative(k);
    const Real sc = std::max({std::abs(an), std::abs(fwd), std::abs(bwd), 1.0});

    if (std::abs(fwd - bwd) > 1e-4 * sc)
    {
      // derivative discontinuity: no finite difference is meaningful here
      if (!n_kinks)
        first_kink = k;
      n_kinks++;
      continue;
    }
    rmax = std::max(rmax, std::abs(0.5 * (fwd + bwd) - an) / sc);
  }

  // r ITSELF must be continuous even where its derivative is not. A kink is a
  // convergence nuisance; a jump in r moves the yield surface discontinuously
  // and is a genuine defect. Check across both switch points a law can carry.
  Real r_jump = 0.0;
  for (const Real ks : {_kmax, _kmin})
    if (ks > 0.0)
    {
      const Real h = 1e-6 * std::max(k_ref, ks);
      r_jump = std::max(r_jump, std::abs(computeSofteningFactor(ks + h) -
                                         computeSofteningFactor(ks - h)));
    }

  const bool ok = gmax < 1e-5 && hmax < 1e-5 && rmax < 1e-5 && r_jump < 1e-4;
  const MooseEnum fd_model = isParamSetByUser("plastic_model")
                                 ? getParam<MooseEnum>("plastic_model")
                                 : getParam<MooseEnum>("material_model");
  Moose::out << "FD check (" << fd_model
             << ", 3 states with shear): yield gradient " << gmax << ", Hessian " << hmax
             << ", softening derivative " << rmax << ", r continuity " << r_jump
             << (ok ? "  OK" : "  FAILED") << "\n";
  if (n_kinks)
    Moose::out << "  note: postyield_mode = " << getParam<MooseEnum>("postyield_mode")
               << " has a derivative kink at kappa = " << first_kink << " ("
               << n_kinks << " sample point(s) skipped). r stays continuous there, so this "
                  "is valid, but the Jacobian jumps as a quadrature point crosses it.\n";

  if (gmax >= 1e-5 || hmax >= 1e-5)
    mooseError("FD check failed: analytic yield gradient/Hessian disagree with finite "
               "differences (", gmax, ", ", hmax, "). A shear-slot conversion factor is "
               "the usual cause.");
  if (rmax >= 1e-5)
    mooseError("FD check failed: computeSofteningDerivative is not the derivative of "
               "computeSofteningFactor for postyield_mode = ",
               getParam<MooseEnum>("postyield_mode"), " (", rmax,
               "), at a point where the one-sided differences agree, so this is not a kink.");
  if (r_jump >= 1e-4)
    mooseError("FD check failed: computeSofteningFactor is discontinuous at kmax or kmin "
               "for postyield_mode = ", getParam<MooseEnum>("postyield_mode"),
               " (jump ", r_jump, "). The yield surface would move discontinuously.");
}
