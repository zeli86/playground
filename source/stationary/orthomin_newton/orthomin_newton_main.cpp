
#include "orthomin_newton.hpp"

int main(int argc, char* argv[])
{
  deallog.depth_console(0);

  BreedSolver_1::MySolver<1> solver("params.xml");
  solver.run();
  return EXIT_SUCCESS;
}