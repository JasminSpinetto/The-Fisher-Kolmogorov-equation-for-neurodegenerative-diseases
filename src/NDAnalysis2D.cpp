#include "NDProblem.hpp"
#include "NDConfig.hpp"
#include "InitialConditions.hpp"
#include "FiberFields.hpp"
#include "NDThetaSolver.hpp"
#include "SeedingRegions.hpp"


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  NDConfig config;
  config.parse(argc, argv);

  // --- Both factories now correctly return smart pointers ---
  auto initial_condition = SeedingRegion<2>::create(config.seeding_region_type, config.C_0);

  Point<2> sagittal_origin(70.0, 73.0);
  Point<2> default_semi_axes(35.0, 15.0);
  auto fiber_field = FiberFieldFactory<2>::create(
    config.fiber_field_type,
    sagittal_origin,
    default_semi_axes
  );

  // --- Both smart pointers are now dereferenced with '*' ---
  NDProblem<2> problem(config.mesh, config.alpha, config.d_ext, config.d_axn, *initial_condition, *fiber_field, config.gray_matter_distance_threshold);
  
  problem.export_problem(config.output_dir + config.output_filename + ".problem");

  NDBackwardEulerSolver<2> solver(
    problem, 
    config.deltat, 
    config.T, 
    config.degree, 
    config.output_dir, 
    config.output_filename,
    config.preconditioner_type
  );
  
  solver.setup();
  solver.solve();

  return EXIT_SUCCESS;
}

