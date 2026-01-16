# Unit Test: Y-Compression
# Expected: Yield at -4100 MPa (or user-defined compression strength)

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
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
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
    
    # Yield stresses (tension)
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    
    # Compression strengths (if different from tension)
    # Uncomment and adjust if you have asymmetry:
    # sigma_xx_compression = 6000
    # sigma_yy_compression = 5000
    # sigma_zz_compression = 6500
    
    # ========== NEW ENUM-BASED POST-YIELD SYSTEM ==========
    # Options: perfect, exp_hardening, linear_hardening,
    #          exp_softening, piecewise_softening
    postyield_mode = 'piecewise_softening'
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
  [compress_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = '-0.001*t'
  []
[]

[Postprocessors]
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  []
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 90.0
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  file_base = y_compression_out_linear_piecewise_softening
  print_linear_residuals = false
[]
