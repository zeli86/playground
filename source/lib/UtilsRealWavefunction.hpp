
#include "IBase.hpp"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/affine_constraints.h"
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <vector>
#include <array>

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
    using IRealWavefunction = IBase<DoFHandler<iDim>, FE_Q<iDim>, AffineConstraints<double>>;

    template<int iDim>
    double particle_number(IRealWavefunction<iDim>*, const Vector<double>&);

    template<int iDim>
    double particle_number(IRealWavefunction<iDim>*, const LA::MPI::Vector&, MPI_Comm);

    template<int iDim>
    std::array<double, iDim> expectation_value_position(IRealWavefunction<iDim>*, const Vector<double>&);

    template<int iDim>
    std::array<double, iDim> expectation_value_position(IRealWavefunction<iDim>*, const LA::MPI::Vector&, MPI_Comm);

    template<int iDim>
    std::array<double, iDim> expectation_value_width(IRealWavefunction<iDim>*, const Vector<double>&, const Point<iDim>&);

    template<int iDim>
    std::array<double, iDim> expectation_value_width(IRealWavefunction<iDim>*, const LA::MPI::Vector&, const Point<iDim>&, MPI_Comm);
  }
}