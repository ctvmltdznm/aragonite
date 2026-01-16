# DIAGNOSTIC TEST: Two-Grain Stress Reporting
# 
# Purpose: Check if stress is being computed correctly per block
#          or if there's an issue with averaging/reporting
#
# Test: Use EXTREME difference (0° vs 90°) to make problem obvious

[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4   # Small mesh for debugging
    ny = 2
    nz = 2
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
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
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
  # Add element-wise stress tracking
  [elem_stress_xx]
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
  # Track stress at element level
  [elem_stress_xx_kernel]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = elem_stress_xx
    index_i = 0
    index_j = 0
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
[]

[Materials]
  #========================================
  # GRAIN 1: Baseline (0,0,0)
  #========================================
  [elasticity_grain1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 'grain_1'
  []
  
  [stress_grain1]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity_grain1'
    block = 'grain_1'
  []
  
  [plasticity_grain1]
    type = OrthotropicPlasticityStressUpdate
    
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    
    postyield_mode = 'perfect'  # Perfect plasticity for now
    residual_strength = 1.0
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0
    
    viscosity_mode = 'rate_independent'  # No viscosity for clarity
    eta = 0.002
    m = 0.001
    
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000
    
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    
    block = 'grain_1'
  []
  
  #========================================
  # GRAIN 2: 90° rotation (0,90,0) -> Y-axis
  #========================================
  [elasticity_grain2]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    euler_angle_1 = 90
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 'grain_2'
  []
  
  [stress_grain2]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity_grain2'
    block = 'grain_2'
  []
  
  [plasticity_grain2]
    type = OrthotropicPlasticityStressUpdate
    
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    
    postyield_mode = 'perfect'
    residual_strength = 1.0
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0
    
    viscosity_mode = 'rate_independent'
    eta = 0.002
    m = 0.001
    
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000
    
    euler_angle_1 = 90
    euler_angle_2 = 0
    euler_angle_3 = 0
    
    block = 'grain_2'
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
    function = '0.0005*t'  # Smaller for elastic only
  []
[]

[Postprocessors]
  # Global
  [stress_xx_avg]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx_avg]
    type = ElementAverageValue
    variable = strain_xx
  []
  
  # Grain 1 - using stress_xx variable
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
  
  # Grain 2 - using stress_xx variable
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
  
  # Alternative: using elem_stress_xx auxiliary variable
  [elem_stress_grain1]
    type = ElementAverageValue
    variable = elem_stress_xx
    block = 'grain_1'
  []
  [elem_stress_grain2]
    type = ElementAverageValue
    variable = elem_stress_xx
    block = 'grain_2'
  []
  
  # Element extremes to check variation
  [stress_xx_min]
    type = ElementExtremeValue
    variable = stress_xx
    value_type = min
  []
  [stress_xx_max]
    type = ElementExtremeValue
    variable = stress_xx
    value_type = max
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 20.0  # Keep elastic
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  console = true
  file_base = diagnostic_out
  
  [exodus]
    type = Exodus
    file_base = diagnostic_out
  []
[]
