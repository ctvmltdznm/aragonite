#include "aragoniteApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Custom materials
#include "OrthotropicPlasticityStressUpdate.h"

InputParameters
aragoniteApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

aragoniteApp::aragoniteApp(const InputParameters & parameters) : MooseApp(parameters)
{
  aragoniteApp::registerAll(_factory, _action_factory, _syntax);
}

aragoniteApp::~aragoniteApp() {}

void
aragoniteApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<aragoniteApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"aragoniteApp"});
  Registry::registerActionsTo(af, {"aragoniteApp"});

  /* register custom execute flags, action syntax, etc. here */
  // Register custom materials
  registerMooseObject("aragoniteApp", OrthotropicPlasticityStressUpdate);

}

void
aragoniteApp::registerApps()
{
  registerApp(aragoniteApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
aragoniteApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  aragoniteApp::registerAll(f, af, s);
}
extern "C" void
aragoniteApp__registerApps()
{
  aragoniteApp::registerApps();
}
