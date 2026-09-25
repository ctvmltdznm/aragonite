# ============================================================
# Bone micro-FE — quadric yield surface (our OrthotropicPlasticityStressUpdate)
# Specimen  : cube
# L         : 5.6520 mm
# BV/TV     : 0.1248
# Load case : tension_22
#
# Elasticity: ComputeFabricElasticityTensor, m=[1,1,1], rho=1 → isotropic E=12700.0 MPa, nu=0.32
# Plasticity: OrthotropicPlasticityStressUpdate, quadric yield surface (Schwiedrzik 2013)
#   With m=[1,1,1], rho=1 → isotropic cast-iron: sigma_t=52.07 MPa, sigma_c=105.41 MPa
#   tau_0=37.04 MPa (Mohr-Coulomb cohesion), zeta_0=0.27 (Wolfram 2012)
#   linear_hardening: r = 1 + kslope*kappa, kslope = H/sigma_t = 635/52.07 = 12.2
#     → direct translation of Wolfram bilinear tissue (H = 0.05*E = 635 MPa)
#
# Timing: ~60 steps × ~5 min/step on 64 cores ≈ 5.0 h
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = FileMesh
  file = cube.exo
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz
                       strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz
                       vonmises_stress hydrostatic_stress'
  []
[]

[AuxVariables]
  [euler_angle_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle_2]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle_3]
    order = CONSTANT
    family = MONOMIAL
  []
  [plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [euler1_ic]
    type = ConstantIC
    variable = euler_angle_1
    value = 0.0
  []
  [euler2_ic]
    type = ConstantIC
    variable = euler_angle_2
    value = 0.0
  []
  [euler3_ic]
    type = ConstantIC
    variable = euler_angle_3
    value = 0.0
  []
[]

[AuxKernels]
  [plastic_strain_aux]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain
    execute_on = TIMESTEP_END
  []
[]

[Materials]
  # Isotropic tissue elasticity via fabric tensor with m=[1,1,1], rho=1.
  # With rho=1, m1=m2=m3=1: reduces exactly to isotropic E=12700.0 MPa, nu=0.32.
  [elasticity]
    type = ComputeFabricElasticityTensor
    E_0         = 12700.0
    nu_0        = 0.32
    fabric_m1   = 1.0
    fabric_m2   = 1.0
    fabric_m3   = 1.0
    density_rho = 1.0
    exponent_k  = 1.0
    exponent_l  = 1.0
    delta_cortical = 1.0
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
    max_iterations = 50
    absolute_tolerance = 1e-8
  []
  # Quadric yield surface (Schwiedrzik 2013), isotropic reduction: m=[1,1,1], rho=1.
  # sigma_i± = sigma_0±, tau_ij = tau_0 for all axes/planes.
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    use_fabric_scaling      = true
    sigma_0_tension         = 52.07
    sigma_0_compression     = 105.41
    tau_0                   = 37.0428761707295
    zeta_0                  = 0.27    # Wolfram 2012: zeta_0=0.27 for trabecular bone tissue
    fabric_m1               = 1.0
    fabric_m2               = 1.0
    fabric_m3               = 1.0
    density_rho             = 1.0
    exponent_p              = 1.0
    exponent_q              = 1.0
    delta_cortical          = 1.0
    postyield_mode          = linear_hardening
    kslope                  = 12.2
    viscosity_mode          = rate_independent
    euler_angle_1           = euler_angle_1
    euler_angle_2           = euler_angle_2
    euler_angle_3           = euler_angle_3
    absolute_tolerance      = 1e-4
    max_iterations_newton   = 100
    use_primal_cpp          = true
    max_iterations_primal   = 1000
  []
[]

[BCs]
  # --- Fixed face: yminus (all DOFs pinned) ---
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = nodeset_3_yminus
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = nodeset_3_yminus
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = nodeset_3_yminus
    value = 0
  []

  # --- Driven face: yplus, ramp u_y ---
  [drive_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = nodeset_4_yplus
    function = '0.113040 * t'
  []
  # Lateral DOFs on driven face are free (allows Poisson contraction)
[]

[Postprocessors]
  # ---------------------------------------------------------------
  # Volume-averaged (apparent) stresses: sigma_bar_ij = (1/V) ∫ sigma_ij dV
  # BV/TV correction applied in post-processing: sigma_apparent = sigma_avg * BV/TV
  # ---------------------------------------------------------------
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

  # ---------------------------------------------------------------
  # Volume-averaged (apparent) strains — stored for diagnostics only.
  # DO NOT divide by BV/TV in post-processing: ElementAverageValue over
  # bone-only elements gives tissue-average strain ≈ ε_macro/BV/TV due
  # to strain concentration. Dividing again by BV/TV would give ε_macro/BV/TV²
  # (inflated ~260× for BV/TV≈0.062), making the 0.2% offset criterion wrong.
  # The correct macroscopic apparent strain used for yield identification is
  #   ε_m = MACRO_STRAIN * time        (normal cases)
  #   ε_m = MACRO_STRAIN * time / √2   (shear cases, γ → tensor → Frobenius)
  # ---------------------------------------------------------------
  [strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  []
  [strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  []
  [strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  []
  [strain_xy]
    type = ElementAverageValue
    variable = strain_xy
  []
  [strain_xz]
    type = ElementAverageValue
    variable = strain_xz
  []
  [strain_yz]
    type = ElementAverageValue
    variable = strain_yz
  []

  # ---------------------------------------------------------------
  # Element-level von Mises extremes (used for linear-mode yield scaling)
  # ---------------------------------------------------------------
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
  # Volume-averaged plastic strain (onset indicator for nonlinear mode)
  [plastic_strain_avg]
    type = ElementAverageValue
    variable = plastic_strain
  []
[]


[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_gamg_agg_nsmooths -ksp_type -ksp_gmres_restart -ksp_max_it -ksp_rtol'
  petsc_options_value = 'gamg     2                     gmres     300                200        1e-4'

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-2
  line_search = l2
  nl_max_its = 100
  l_max_its  = 200

  dt       = 0.1
  end_time = 1.5
[]

[Outputs]
  [csv]
    type = CSV
    file_base = cube_tension_22
  []
  [exo_out]
    type = Exodus
    file_base = cube_tension_22
    time_step_interval = 5
  []
  print_linear_residuals = false
  perf_graph             = false
  [console]
    type = Console
    verbose            = false
  []
[]
