
#include "IBase.hpp"
#include "deal.II/fe/fe_q.h"

namespace utils
{
  using namespace dealii;
  namespace LA
  {
    using namespace dealii::LinearAlgebraPETSc;
  }

  template<int dim>
  double particle_number(const DoFHandler<dim>& dof_handler, const FE_Q<dim>& fe, const Vector<double>& vec);

  template<int dim>
  double particle_number(const DoFHandler<dim>&, const FE_Q<dim>&, const LA::MPI::Vector&, MPI_Comm);

  //template <int dim>
  //double Particle_Number(MPI_Comm, IBase<DoFHandler<dim>, FE_Q<dim>, const Vector<double>&, Fe, tConstraints>* );
}