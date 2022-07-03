
#include "sobolev_gradient_mpi.hpp"

int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run();
  }
  return EXIT_SUCCESS;
}