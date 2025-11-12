//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "aragoniteTestApp.h"
#include "aragoniteApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
aragoniteTestApp::validParams()
{
  InputParameters params = aragoniteApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

aragoniteTestApp::aragoniteTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  aragoniteTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

aragoniteTestApp::~aragoniteTestApp() {}

void
aragoniteTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  aragoniteApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"aragoniteTestApp"});
    Registry::registerActionsTo(af, {"aragoniteTestApp"});
  }
}

void
aragoniteTestApp::registerApps()
{
  registerApp(aragoniteApp);
  registerApp(aragoniteTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
aragoniteTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  aragoniteTestApp::registerAll(f, af, s);
}
extern "C" void
aragoniteTestApp__registerApps()
{
  aragoniteTestApp::registerApps();
}
