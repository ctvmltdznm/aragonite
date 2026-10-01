# ============================================================================
# compact_iso -- compression_yy
#
# UMAT flag 4. Isotropic compact (cortical) bone (ISORBCT, JJS 2013).
# Post-yield is exp_hardening (UMAT PYFL=1), not softening.
#
# LOAD CASE: uniaxial COMPRESSION along y, 2.0% nominal strain at t = 1.
#            Peak |stress_yy| = sigma_yy_compression x r(0).
#
# Single HEX8 element, 1 x 1 x 1 mm. Run as is, no arguments:
#     aragonite-opt -i compression_yy.i
# Gold file:
#     aragonite-opt -i compression_yy.i --generate-gold
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
  material_model = compact_iso
  density_rho    = 0.9
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

[BCs]
  # Symmetry on the three minus faces: removes rigid-body motion and lets
  # the material contract freely on the opposite faces.
  [sym_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [sym_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [sym_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  # Driven face. The other two plus faces are traction free, so the state
  # is uniaxial stress and the peak is the directional strength x r(0).
  [drive]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = '-0.02*t'
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
  [strain_yy]
    type = ElementAverageValue
    variable = strain_yy
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
