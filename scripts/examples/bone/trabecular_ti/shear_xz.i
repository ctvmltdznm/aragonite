# ============================================================================
# trabecular_ti -- shear_xz
#
# UMAT flag 1. Transversely isotropic trabecular bone (TIRGWTB).
# main_direction = 3, so zz is the axial direction and xx/yy are
# transverse. Exercises the zeta_a -> zeta_ij reference-index
# conversion (see MATERIAL_MODEL_THEORY.md section 2.3).
#
# LOAD CASE: simple SHEAR in the xz plane, engineering gamma = 2.0% at t = 1.
#            Peak stress_xz = tau_xz_max x r(0).
#
# Single HEX8 element, 1 x 1 x 1 mm. Run as is, no arguments:
#     aragonite-opt -i shear_xz.i
# Gold file:
#     aragonite-opt -i shear_xz.i --generate-gold
#
# The elastic and plastic responses are BOTH driven by one flag,
# material_model, set once in [GlobalParams]. Setting elastic_model or
# plastic_model in a single block would override it there and warn.
# ============================================================================

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmax = 1
    ymax = 1
    zmax = 1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

  # ---- the material model, and the microstructure it needs -----------------
  # One place. ComputeFabricElasticityTensor and
  # OrthotropicPlasticityStressUpdate both read these, so the elastic and
  # plastic responses cannot drift apart.
  material_model = trabecular_ti
  density_rho    = 0.2
  main_direction = 3
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    # REQUIRED with FINITE strain and solve_type = NEWTON. Without it MOOSE
    # drops the rotation-increment terms and convergence stays linear however
    # exact the material tangent is.
    use_finite_deform_jacobian = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz
                       strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'
  []
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    # model and microstructure come from [GlobalParams]
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    # model and microstructure come from [GlobalParams]
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    absolute_tolerance = 1e-8
    max_iterations_newton = 50
    # Startup self-checks: finite-difference test of the yield gradient, Hessian
    # and softening derivative at stress states with shear, and a Mandel
    # round-trip of the elasticity tensor. Cheap; turn off for production.
    debug_checks = true
  []
[]

[Functions]
  [shear_fn]
    type = ParsedFunction
    expression = '0.02 * z * t'
  []
[]

[BCs]
  # Affine (Taylor) boundary conditions on EVERY face:
  #     u_x = gamma * z * t,   the other two components zero.
  # The deformation is homogeneous simple shear with engineering shear
  # gamma = 0.02 at t = 1. Both the orthotropic stiffness and the quadric
  # yield surface are block diagonal in the material frame, so the stress
  # stays pure shear and the peak is tau_xz x r(0).
  [u_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right bottom top back front'
    function = shear_fn
  []
  [u_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left right bottom top back front'
    value = 0
  []
  [u_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'left right bottom top back front'
    value = 0
  []
[]

[Postprocessors]
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [stress_xz]
    type = ElementAverageValue
    variable = stress_xz
  []
  [stress_yz]
    type = ElementAverageValue
    variable = stress_yz
  []
  [strain_xz]
    type = ElementAverageValue
    variable = strain_xz
  []
  [plastic_strain]
    type = ElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = bt
  dt = 0.01
  end_time = 1.0
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-8
  nl_max_its = 50
[]

[Outputs]
  csv = true
[]
