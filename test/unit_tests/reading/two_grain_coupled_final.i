# Two-Grain Test: 0° vs 90° with Coupled Elasticity
# Both elastic AND plastic responses correct!

[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 5
    nz = 5
    xmin = 0
    xmax = 2.0
    ymin = 0
    ymax = 1.0
    zmin = 0
    zmax = 1.0
  []
  
  [split]
    type = SubdomainBoundingBoxGenerator
    input = generated
    block_id = 2
    block_name = 'grain_2'
    bottom_left = '1.0 0 0'
    top_right = '2.0 1.0 1.0'
  []
  
  [rename]
    type = RenameBlockGenerator
    input = split
    old_block = '0'
    new_block = 'grain_1'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
  []
[]

[AuxVariables]
  [plastic_strain_equiv]
    order = CONSTANT
    family = MONOMIAL
  []
  [grain_id]
    order = CONSTANT
    family = MONOMIAL
  []
  
  # Euler angles - per element!
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
[]

[AuxKernels]
  [plastic_strain_equiv]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_equiv
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  [grain_1_ic]
    type = ConstantIC
    variable = grain_id
    value = 1
    block = 'grain_1'
  []
  [grain_2_ic]
    type = ConstantIC
    variable = grain_id
    value = 2
    block = 'grain_2'
  []
  
  # Set Euler angles per grain
  [euler_phi1_grain1]
    type = ConstantIC
    variable = euler_phi1
    value = 0
    block = 'grain_1'
  []
  [euler_phi1_grain2]
    type = ConstantIC
    variable = euler_phi1
    value = 90
    block = 'grain_2'
  []
  
  [euler_Phi_grain1]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = 'grain_1'
  []
  [euler_Phi_grain2]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = 'grain_2'
  []
  
  [euler_phi2_grain1]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = 'grain_1'
  []
  [euler_phi2_grain2]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = 'grain_2'
  []
[]

[Materials]
  #========================================
  # GRAIN 1 & 2: Single material definition!
  # Reads per-element angles from AuxVariables
  #========================================
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    
    # Couple to Euler angle AuxVariables
    # Note: "coupled_" prefix for elasticity!
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
    
    # Couple to SAME AuxVariables
    # Note: NO "coupled_" prefix for plasticity!
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
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
  
  [pull_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = '0.001*t'
  []
[]

[Postprocessors]
  # ========== GLOBAL AVERAGES ==========
  [stress_xx_global]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx_global]
    type = ElementAverageValue
    variable = strain_xx
  []
  [vonmises_global]
    type = ElementAverageValue
    variable = vonmises_stress
  []
  [plastic_strain_global]
    type = ElementAverageValue
    variable = plastic_strain_equiv
  []
  
  # ========== GRAIN 1 AVERAGES ==========
  [stress_xx_grain1]
    type = ElementAverageValue
    variable = stress_xx
    block = 'grain_1'
  []
  [strain_xx_grain1]
    type = ElementAverageValue
    variable = strain_xx
    block = 'grain_1'
  []
  [vonmises_grain1]
    type = ElementAverageValue
    variable = vonmises_stress
    block = 'grain_1'
  []
  [plastic_strain_grain1]
    type = ElementAverageValue
    variable = plastic_strain_equiv
    block = 'grain_1'
  []
  
  # ========== GRAIN 2 AVERAGES ==========
  [stress_xx_grain2]
    type = ElementAverageValue
    variable = stress_xx
    block = 'grain_2'
  []
  [strain_xx_grain2]
    type = ElementAverageValue
    variable = strain_xx
    block = 'grain_2'
  []
  [vonmises_grain2]
    type = ElementAverageValue
    variable = vonmises_stress
    block = 'grain_2'
  []
  [plastic_strain_grain2]
    type = ElementAverageValue
    variable = plastic_strain_equiv
    block = 'grain_2'
  []
  
  # ========== ELASTIC MODULI (for validation) ==========
  [E_grain1]
    type = ParsedPostprocessor
    pp_names = 'stress_xx_grain1 strain_xx_grain1'
    expression = 'stress_xx_grain1 / strain_xx_grain1'
  []
  [E_grain2]
    type = ParsedPostprocessor
    pp_names = 'stress_xx_grain2 strain_xx_grain2'
    expression = 'stress_xx_grain2 / strain_xx_grain2'
  []
  [E_ratio]
    type = ParsedPostprocessor
    pp_names = 'E_grain1 E_grain2'
    expression = 'E_grain2 / E_grain1'
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 150
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
  line_search = 'basic'
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 1.2
    cutback_factor = 0.5
  []
[]

[Outputs]
  csv = true
  console = true
  print_linear_residuals = false
  
  [exodus]
    type = Exodus
    file_base = two_grain_coupled_out
    time_step_interval = 10
  []
[]
