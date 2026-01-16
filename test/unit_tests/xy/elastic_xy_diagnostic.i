# Simple Elastic XY Shear Test (NO PLASTICITY)
# Purpose: Verify elastic response and identify any scaling issues

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

[Materials]
  [elasticity]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  []
  
  # ONLY ELASTIC - NO PLASTICITY!
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[BCs]
  [fix_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [fix_y_bottom]
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
  
  # Pure shear
  [shear_top]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    function = '0.001*t'
  []
  [shear_right]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'right'
    function = '0.001*t'
  []
[]

[Postprocessors]
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [strain_xy]
    type = ElementAverageValue
    variable = strain_xy
  []
[]

[Executioner]
  type = Transient
  dt = 1.0
  end_time = 10.0
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  file_base = elastic_xy_diagnostic
  print_linear_residuals = false
[]
