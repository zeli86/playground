
#include "IBase.hpp"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

namespace utils
{
  namespace complex_wavefunction
  {
    using namespace dealii;

    template<int iDim>
    using IComplexWavefunction = IBase<DoFHandler<iDim>, FESystem<iDim>, AffineConstraints<double>>;

    template <int dim>
    void assemble_mulvz(IComplexWavefunction<dim>*, const Vector<double>&, const std::complex<double>, SparseMatrix<double>&, Vector<double>&);

    template <int dim>
    void assemble_mulvz(IComplexWavefunction<dim>*, const LA::MPI::Vector&, const std::complex<double> z, LA::MPI::SparseMatrix&, LA::MPI::Vector&);

    template <int dim>
    void assemble_lin_step(IComplexWavefunction<dim>*, const Vector<double>&, const Function<dim>&, const double, SparseMatrix<double>&, Vector<double>&);

    template <int dim>
    void assemble_lin_step(IComplexWavefunction<dim>*, const Vector<double>&, const double, SparseMatrix<double>&, Vector<double>&);

    template <int dim>
    void assemble_lin_step(IComplexWavefunction<dim>*, const Vector<double>&, const double, SparseMatrix<double>&, Vector<double>&);

    template <int dim>
    void assemble_nl_step(IComplexWavefunction<dim>*, const Vector<double>&, const Function<dim>&, const double, const double, SparseMatrix<double>&, Vector<double>&);

    //template <int dim>
    //void Interpolate_R_to_C(IRealWavefunction<dim>*, IComplexWavefunction<dim>*, const Vector<double>&, Vector<double>&);
  }
}