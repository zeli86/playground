
#include <cassert>
#include <complex>
#include <cmath>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

namespace MyRealTools
{
  namespace MPI
  {
    namespace LA
    {
      using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;

    template <int dim>
    void AssembleSystem_tangent(const DoFHandler<dim>&,
                                const FE_Q<dim>&,
                                const AffineConstraints<double>&,
                                const LA::MPI::Vector&,
                                const Function<dim>& Potential,
                                const double mu,
                                const double gs,
                                LA::MPI::SparseMatrix& matrix,
                                LA::MPI::Vector& rhs);


    template <int dim>
    void compute_stepsize(MPI_Comm,
                          const DoFHandler<dim>&,
                          const FE_Q<dim>& fe,
                          const Function<dim>&,
                          const LA::MPI::Vector&,
                          const LA::MPI::Vector&,
                          const double,
                          const double,
                          double&);


    template <int dim>
    void orthonormalize(MPI_Comm mpi_communicator,
                        const DoFHandler<dim>&,
                        const FE_Q<dim>&,
                        const  AffineConstraints<double>&,
                        const LA::MPI::Vector&,
                        const LA::MPI::Vector&,
                        LA::MPI::Vector&);
  }
}
