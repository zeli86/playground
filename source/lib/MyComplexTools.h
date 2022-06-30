
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/numerics/matrix_tools.h>


namespace MyComplexTools
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  void Interpolate_R_to_C(const DoFHandler<dim>& dof_handler, const FE_Q<dim>&,
                          const Vector<double>& vec,
                          const DoFHandler<dim>& dof_handler_2,
                          const FESystem<dim>& _2,
                          const dealii::AffineConstraints<double>& constraints,
                          Vector<double>& ret);

} // end of namespace

namespace MyComplexTools
{
  namespace MPI
  {
    namespace LA
    {
      using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;


    template <int dim>
    std::complex<double> L2_dot_product(MPI_Comm,
                                        const DoFHandler<dim>&,
                                        const FESystem<dim>&, const LA::MPI::Vector&,
                                        const LA::MPI::Vector&);

    template <int dim>
    void Interpolate_R_to_C(MPI_Comm mpi_communicator,
                            const DoFHandler<dim>& dof_handler, const FE_Q<dim>&,
                            const LA::MPI::Vector& vec,
                            const DoFHandler<dim>& dof_handler_2,
                            const FESystem<dim>& _2,
                            const dealii::AffineConstraints<double>& constraints,
                            LA::MPI::Vector& ret);

  }
}
