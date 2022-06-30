
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

    template<int iDim>
    double particle_number(IComplexWavefunction<iDim>*, const Vector<double>&);

    template<int iDim>
    double particle_number(IComplexWavefunction<iDim>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    Point<iDim> expectation_value_position(IComplexWavefunction<iDim>*, const Vector<double>&);

    template <int iDim>
    Point<iDim> expectation_value_position(IComplexWavefunction<iDim>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    Point<iDim> expectation_value_momentum(IComplexWavefunction<iDim>*, const Vector<double>&);

    template <int iDim>
    Point<iDim> expectation_value_momentum(IComplexWavefunction<iDim>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    Point<iDim> expectation_value_width(IComplexWavefunction<iDim>*, const Vector<double>&, const Point<iDim>&);

    template <int iDim>
    Point<iDim> expectation_value_width(IComplexWavefunction<iDim>*, const LA::MPI::Vector& vec, const Point<iDim>&, MPI_Comm);
  }
}