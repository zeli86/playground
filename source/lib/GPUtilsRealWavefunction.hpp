
#include "IBase.hpp"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/affine_constraints.h"
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

#include <vector>
#include <tuple>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

namespace utils
{
  namespace real_wavefunction
  {
    using namespace dealii;

    template<int iDim>
    using IBaseRealWavefunction = IBase<DoFHandler<iDim>, FE_Q<iDim>, AffineConstraints<double>>;

    template<int iDim>
    std::tuple<double, double, double>
    GP(IBaseRealWavefunction<iDim>*, const std::vector<double>&,  const Function<iDim>&);

    template<int iDim>
    std::tuple<double, double, double>
    GP(IBaseRealWavefunction<iDim>*, const LA::MPI::Vector&,  const Function<iDim>&, MPI_Comm);

    template<int iDim>
    double
    mu(IBaseRealWavefunction<iDim>*, const std::vector<double>&,  const Function<iDim>&);

    template<int iDim>
    double
    mu(IBaseRealWavefunction<iDim>*, const LA::MPI::Vector&,  const Function<iDim>&, MPI_Comm);

    template<int iDim>
    SparseMatrix<double>
    assemble_jacobian(IBaseRealWavefunction<iDim>*, std::vector<double>&, const Function<iDim>&, const double mu, const double gs);

    template<int iDim>
    void
    assemble_jacobian(IBaseRealWavefunction<iDim>*, const LA::MPI::Vector&, const Function<iDim>&, LA::MPI::SparseMatrix&, const double mu, const double gs);

    template<int iDim>
    double
    assemble_L2gradient(const Vector<double>&, const Function<iDim>&, const double, const double, Vector<double>&);

    template<int iDim>
    double
    assemble_L2gradient(const LA::MPI::Vector&, const Function<iDim>&, const double, const double, LA::MPI::Vector&);
  }
}