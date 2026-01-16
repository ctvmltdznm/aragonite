# Unit Test: YZ-Shear
# Expected: Yield at 5060 MPa, soften to 3542 MPa

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 1
  ymax = 1
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_yz strain_xz'
  []
[]

[AuxVariables]
  [plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plastic_strain]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  [elasticity]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
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
    
    # ========== NEW ENUM-BASED POST-YIELD SYSTEM ==========
    # Options: perfect, exp_hardening, linear_hardening,
    #          exp_softening, piecewise_softening
    postyield_mode = 'perfect'
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
    max_iterations_primal = 1000     # Primal CPPA max iterations
    line_search_beta = 1e-4          # Armero line search parameter
    line_search_eta = 0.25           # Line search backtracking factor
    
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  []
[]

[BCs]
  # ============================================
  # YZ SHEAR BCs: Fix back (z=0), shear front (z=1) in y
  # ============================================
  
  # Back face (z=0): completely fixed
  [fix_back_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0
  []
  [fix_back_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0
  []
  [fix_back_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  
  # Front face (z=1): apply shear in y, fix z
  [shear_front_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'front'
    # γ_yz = u_y / h = 0.001*t, so σ_yz = G*γ = 42100*0.001*t = 42.1*t MPa
    # Yield at σ_yz = 5060: t = 5060/42.1 ≈ 120.2 s
    function = 'if(t<=240, 0.001*t, 0.001*(480-t))'  # Load to t=55, unload to t=10
  []
  [fix_front_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
    value = 0
  []
  # Note: disp_x on front is FREE
  
  # Left/right faces: constrain x for plane strain
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [fix_right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
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
  [strain_yz]
    type = ElementAverageValue
    variable = strain_yz
  []
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 480.0
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  file_base = yz_shear_out_cyclic_perfect
  print_linear_residuals = false
[]
