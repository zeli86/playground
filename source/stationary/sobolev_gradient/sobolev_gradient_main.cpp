
#include "sobolev_gradient.hpp"

int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  solver::stationary::CSobolevGradient<1> solver("params_one.xml");
  solver.run();
  return EXIT_SUCCESS;
}