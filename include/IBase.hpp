
#pragma once

#include "deal.II/base/index_set.h"
#include "deal.II/distributed/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/grid/tria.h"
#include "deal.II/lac/generic_linear_algebra.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/lac/vector.h"

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