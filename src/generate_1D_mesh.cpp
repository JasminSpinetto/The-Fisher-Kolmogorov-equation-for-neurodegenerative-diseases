#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <iostream>
#include <vector> // Make sure to include vector

using namespace dealii;

int main()
{
  const int dim = 1; // We are creating a 1D mesh

  // The domain is a line from P1=(-1) to P2=(+1)
  const Point<dim> p1(-1.0);
  const Point<dim> p2(1.0);

  // We want 200 elements along the line
  const unsigned int n_elements = 200;

  // This vector tells the generator how many times to subdivide in each dimension.
  // For 1D, it's just one number.
  const std::vector<unsigned int> subdivisions = {n_elements};

  // Create a Triangulation object to hold the mesh
  Triangulation<dim> line_tria;

  // --- THE FIX ---
  // Use subdivided_hyper_rectangle to create the line with the desired number of elements.
  GridGenerator::subdivided_hyper_rectangle(line_tria, subdivisions, p1, p2);

  // Define the output filename
  std::string filename = "../meshes/mesh-1D-centered.msh";

  // Create an output stream and write the mesh in Gmsh (.msh) format
  std::ofstream out(filename);
  GridOut grid_out;
  grid_out.write_msh(line_tria, out);

  std::cout << "Successfully generated mesh file: " << filename << std::endl;

  return 0;
}