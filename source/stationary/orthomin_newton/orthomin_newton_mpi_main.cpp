
#include "orthomin_newton_mpi.hpp"

int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  MyParameterHandler params("params.xml");
  int                dim = 0;

  try
  {
    params.get("parameter.spatial_dimension", dim);
  }
  catch (mu::Parser::exception_type& e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    switch (dim)
    {
      case 2:
      {
        BreedSolver::MySolver<2> solver("params.xml");
        solver.run2b();
        break;
      }
      case 3:
      {
        BreedSolver::MySolver<3> solver("params.xml");
        solver.run2b();
        break;
      }
      default:
        std::cout << "You have found a new dimension!" << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
