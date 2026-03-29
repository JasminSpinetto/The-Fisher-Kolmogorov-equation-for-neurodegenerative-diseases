#include "NDProblem.hpp"
#include "InitialConditions.hpp"
#include "FiberFields.hpp"
#include "NDThetaSolver.hpp"
#include "NDConfig.hpp"

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // This is now a command-line driven program with no hard-coded defaults
  NDConfig config;
  config.parse(argc, argv);

  // --- THIS IS THE KEY CHANGE ---
  // Create our new "smart" initial condition.
  // It uses the mesh path provided via the command line to center itself.
  // CenteredExponentialInitialCondition<1> initial_condition(config.mesh);
  // in main()

  CenteredExponentialInitialCondition<1> initial_condition(config.mesh, config.C_0);
  // ExponentialInitialCondition<1> initial_condition;

  // The fiber field for a 1D problem is simple
  RadialFiberField<1> fiber_field;

  // Create the problem definition
  NDProblem<1> problem(config.mesh, config.alpha, config.d_ext, config.d_axn, initial_condition, fiber_field);
  problem.export_problem(std::string(config.output_dir) + config.output_filename + ".problem");

  // Create the numerical solver
  NDBackwardEulerSolver<1> solver(problem, config.deltat, config.T, config.degree, config.output_dir, config.output_filename);
  solver.setup();
  solver.solve();

  return EXIT_SUCCESS;
}