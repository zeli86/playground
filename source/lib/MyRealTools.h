
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
  using namespace std;
  using namespace dealii;

  template<int dim>
  double Particle_Number(const DoFHandler<dim>&, const FE_Q<dim>&, const Vector<double>&);

  template <int dim>
  void compute_stepsize(const DoFHandler<dim>& dof_handler,
                        const FE_Q<dim>& fe,
                        const Function<dim>& Potential,
                        const Vector<double>& psi,
                        const Vector<double>& direction,
                        const double mu,
                        const double gs,
                        double& retval);
}


namespace MyRealTools
{
  namespace MPI
  {
    namespace LA
    {
      using namespace dealii::LinearAlgebraPETSc;
    }

    using namespace dealii;

    template<int dim>
    double Particle_Number(MPI_Comm, const DoFHandler<dim>&, const FE_Q<dim>&, const LA::MPI::Vector&);

    template<int dim>
    void Expectation_value_position(MPI_Comm, const DoFHandler<dim>&, const FE_Q<dim>&, const LA::MPI::Vector&, vector<double>& );

    template<int dim>
    void Expectation_value_width(MPI_Comm, const DoFHandler<dim>&, const FE_Q<dim>&, const LA::MPI::Vector&, const vector<double>&, vector<double>&);

    template <int dim>
    void AssembleSystem_Jacobian(const DoFHandler<dim>&, const FE_Q<dim>& fe, const  AffineConstraints<double>&,
                                 const LA::MPI::Vector&,
                                 const Function<dim>&,
                                 const double,
                                 const double,
                                 LA::MPI::SparseMatrix&);

    template <int dim>
    void AssembleRHS_L2gradient(const DoFHandler<dim>&,
                                const FE_Q<dim>& fe,
                                const AffineConstraints<double>&,
                                const LA::MPI::Vector&,
                                const Function<dim>&,
                                const double,
                                const double,
                                double& res,
                                LA::MPI::Vector&);

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


    template<int dim>
    void Compute_mu(MPI_Comm mpi_communicator, const DoFHandler<dim>&, const FE_Q<dim>&, const Function<dim>&, const LA::MPI::Vector&, double&);

    template<int dim>
    void Compute_E_lin(MPI_Comm mpi_communicator, const DoFHandler<dim>&, const FE_Q<dim>&,const  AffineConstraints<double>&, const Function<dim>&, const LA::MPI::Vector&, double&, double&, double&);


    template <int dim>
    void compute_stepsize(MPI_Comm,
                          const DoFHandler<dim>&,
                          const FE_Q<dim>& fe,
                          const Function<dim>& ,
                          const LA::MPI::Vector& ,
                          const LA::MPI::Vector& ,
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
