
#pragma once

// #include "deal.II/base/index_set.h"
// #include "deal.II/distributed/tria.h"
// #include "deal.II/dofs/dof_handler.h"
// #include "deal.II/grid/tria.h"
// #include "deal.II/lac/generic_linear_algebra.h"
// #include "deal.II/lac/sparse_matrix.h"
// #include "deal.II/lac/vector.h"

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
//#include <deal.II/numerics/error_estimator.h>
//#include <deal.II/numerics/derivative_approximation.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, MAXITER, SINGULAR };

template<int iDim>
using tMpiTriangulation = dealii::parallel::distributed::Triangulation<iDim>;
using tMpiVector = dealii::LinearAlgebraPETSc::MPI::Vector;
using tMpiMatrix = dealii::LinearAlgebraPETSc::MPI::SparseMatrix;

template<int iDim>
using tTriangulation = dealii::Triangulation<iDim>;
using tVector = dealii::Vector<double>;
using tMatrix = dealii::SparseMatrix<double>;

using tConstraints = dealii::AffineConstraints<double>;

template<int iDim>
using tDoFHandler = dealii::DoFHandler<iDim, iDim>;

using tIndexSet = dealii::IndexSet;

template<typename tDoFHandler, typename tFe, typename tConstraints>
class IBase
{
public:
  virtual tDoFHandler& get_dof_handler() = 0;
  virtual tFe& get_fe() = 0;
  virtual tConstraints& get_constraints() = 0;
};