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
    strain = FINITE
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
  [plastic_strain]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  [euler_phi1_0]
    type = ConstantIC
    variable = euler_phi1
    value = 0
    block = '0'
  []
  [euler_Phi_0]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = '0'
  []
  [euler_phi2_0]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = '0'
  []
[]

[Materials]
  [elasticity_tensor]
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
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
    postyield_mode = perfect
    absolute_tolerance = 1e-6
    viscosity_mode = 'linear'
    eta = 0.002
    m = 0.001
    max_iterations_newton = 10
    use_primal_cpp = true
    max_iterations_primal = 1000     # Primal CPPA max iterations
    line_search_beta = 1e-4          # Armero line search parameter
    line_search_eta = 0.25           # Line search backtracking factor
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
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 120.0
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
