# Two-Grain Aragonite CZM Test - Tensile load-unload-reload test

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 1
    nz = 1
    xmin = 0
    xmax = 2.0
    ymin = 0
    ymax = 1.0
    zmin = 0
    zmax = 1.0
  []
  
  [block1]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    block_name = 'grain_1'
    bottom_left = '0 0 0'
    top_right = '1.0 1.0 1.0'
  []
  
  [block2]
    type = SubdomainBoundingBoxGenerator
    input = block1
    block_id = 2
    block_name = 'grain_2'
    bottom_left = '1.0 0 0'
    top_right = '2.0 1.0 1.0'
  []
  
  [break_mesh]
    type = BreakMeshByBlockGenerator
    input = block2
    split_interface = true
  []

  # Newer MOOSE refuses a displaced mesh together with the disjoint neighbour
  # pairs BreakMeshByBlockGenerator creates, and aborts on GlobalParams/
  # displacements below. strain = SMALL needs no displaced mesh, so this is
  # both correct and version-portable. Without it a clean MOOSE checkout fails
  # with "A displaced mesh cannot be used with a mesh that has disjoint
  # neighbor boundary pairs".
  use_displaced_mesh = false
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics]
  [QuasiStatic]
    [bulk]
      strain = SMALL
      incremental = true
      add_variables = true
      generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
    []
  []
    
  [CohesiveZone]
    [./czm]
      boundary = 'grain_1_grain_2'
      strain = SMALL
      generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
    [../]
  []
[]

[AuxVariables]
  [plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [grain_id]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_phi1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_Phi]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_phi2]
    order = CONSTANT
    family = MONOMIAL
  []
  [damage]
    order = CONSTANT
    family = MONOMIAL
  []
  [delta_eff]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [phi1_1]
    type = ConstantIC
    variable = euler_phi1
    value = 0
    block = grain_1
  []
  [Phi_1]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = grain_1
  []
  [phi2_1]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = grain_1
  []
  
  [phi1_2]
    type = ConstantIC
    variable = euler_phi1
    value = 0#-31.1
    block = grain_2
  []
  [Phi_2]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = grain_2
  []
  [phi2_2]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = grain_2
  []
  
  [grain1_id]
    type = ConstantIC
    variable = grain_id
    value = 1
    block = grain_1
  []
  [grain2_id]
    type = ConstantIC
    variable = grain_id
    value = 2
    block = grain_2
  []
[]

[AuxKernels]
  [plastic_strain]
    type = MaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = 'TIMESTEP_END'
  []
  # Interface properties - only exist on CZM boundary!
  [damage_aux]
    type = MaterialRealAux
    variable = damage
    property = damage
    boundary = 'grain_1_grain_2'  # Only on interface
    execute_on = 'TIMESTEP_END'
  []
  [delta_eff_aux]
    type = MaterialRealAux
    variable = delta_eff
    property = delta_eff
    boundary = 'grain_1_grain_2'  # Only on interface
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    postyield_mode = 'exp_softening'
    residual_strength = 0.7
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0
    viscosity_mode = 'linear'
    eta = 0.002
    m = 0.001
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000
    line_search_beta = 1e-4
    line_search_eta = 0.25
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
  []
  
  # HOMOGENIZED INTERFACE
  [czm_homogenized]
    type = HomogenizedExponentialCZM
    boundary = 'grain_1_grain_2'

    # ========================================
    # MD-DERIVED STRENGTHS
    # ========================================
    normal_strength = 700.0       # σc (MPa)
    shear_strength_s = 300.0      # τs (MPa)
    shear_strength_t = 500.0      # τt (MPa)
    
    # ========================================
    # CHARACTERISTIC SEPARATION FROM MD
    # ========================================
    delta_0_tangent = 1e-4
    delta_0_normal = 1.2e-4

    # ========================================
    # HOMOGENIZATION PARAMETERS
    # ========================================
    # REGRESSION TEST: the stochastic quality field is switched OFF so this
    # case is deterministic and mesh-numbering independent. f_spatial is drawn
    # from _current_elem->id(), so any remesh reshuffles it; and the GH spread
    # makes the result depend on five contact classes. The heterogeneous
    # behaviour is exercised by the assertion-based suite, not by a gold file.
    quality_std_dev = 0            # within-QP GH spread  (0 = single contact)
    spatial_quality_std_dev = 0    # across-QP spatial heterogeneity (0 = uniform)
    spatial_random_seed = 1234     # spatial quality field seed
    
    # ========================================
    # EXPONENTIAL COEFFICIENTS
    # ========================================
    mu = 1.0                      # Loading exponent
    eta = 0.15                   # Softening exponent

    # ========================================
    # VISCOPLASTIC REGULARIZATION
    # Helps with convergence by smoothing damage evolution
    # ========================================
    # damage_viscosity is in SIMULATION time units, and dtmax here is 0.01.
    # The old value 1.25 was therefore a 125-STEP lag, not a mild
    # regularisation: D fell ~0.04 behind D_target, and since the post-peak
    # traction is T_peak*r*(1-D) with T_peak*r unbounded, that lag became a
    # factor ~200 in traction (3727 MPa against a model ceiling of 900 MPa).
    # Rate-independent is the correct reference: the identity
    #   1 - D_law = e*r^(eta-1)*exp(-r^eta)
    # then makes T_peak*r*(1-D) reproduce the exponential envelope exactly.
    damage_viscosity = 0          # rate-independent

    # Pinned so a future change to the code default cannot move this gold.
    damage_activation_width = 0.02

    # With eta = 0.15 the softening tail is very fat: the traction is still 27%
    # of peak at r = 5000, so the solve for delta_c from
    # failure_traction_ratio = 0.01 needs r ~ 7.7e5 and hits the internal 5000
    # clamp. Setting delta_c explicitly keeps the previous effective value,
    # documents it, and silences the clamp warning. Complete failure in this
    # model is driven by D -> 1, not by reaching delta_c.
    delta_c = 0.6
  []
[]

[Functions]
  [load_unload_reload]
    type = PiecewiseLinear
    x = '0  2  4  12'         # Time
    y = '0  0.02  0  0.08'       # Displacement
    # Load → Unload → Reload
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []

  [pull]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = load_unload_reload
  []
[]

[Postprocessors]
  # Bulk behavior
  [stress_xx_global]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx_global]
    type = ElementAverageValue
    variable = strain_xx
  []
  
  # Interface behavior
  [max_normal_traction]
    type = ElementExtremeValue
    variable = normal_traction
  []
  [max_normal_jump]
    type = ElementExtremeValue
    variable = normal_jump
  []
  [avg_separation]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'grain_1_grain_2'
  []
  [avg_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'grain_1_grain_2'
  []

  # Damage tracking
  [avg_damage]
    type = SideAverageValue
    variable = damage
    boundary = 'grain_1_grain_2'
  []
  [max_damage]
    type = ElementExtremeValue
    variable = damage
  []
  [avg_delta_eff]
    type = SideAverageValue
    variable = delta_eff
    boundary = 'grain_1_grain_2'
  []
  
  [newton_its]
    type = NumNonlinearIterations
  []
[]

[Executioner]
  type = Transient
  dtmax = 0.01
  dtmin = 1e-10
  dt = 0.001
  end_time = 12
  
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_abs_tol = 1e-8
  nl_max_its = 200
  
  line_search = 'l2'
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 12
    iteration_window = 5
    growth_factor = 1.2
    cutback_factor = 0.1
  []
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
