//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
//
// This file is part of atus-pro testing.
//
// atus-pro testing is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// atus-pro testing is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
//

/** Želimir Marojević
 * Soll mal eines Tages das Megnetfeld eines Stromdurchfloßenen Leiters berechnen
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

#include "mpi.h"

namespace realtime_propagation
{
enum Status { SUCCESS, FAILED };

using namespace std;
using namespace dealii;

double Heaviside( double p ) { return 0.5*(1.0+p/fabs(p)); };
double rect( double x0, double L, double x ) { return Heaviside(L*L-(x-x0)*(x-x0)); };

class RightHandSide :  public AutoDerivativeFunction<3>
{
public:
  RightHandSide ();
  virtual double value (const Point<3> &p, const unsigned int component=0 ) const;
  virtual void vector_value (const Point<3> &p, Vector<double> &values) const;
protected:
  double A;
  double B;
  
};

RightHandSide::RightHandSide () : AutoDerivativeFunction (0.001, 3)
{
  A=0.1;
  B=0.1;
}

inline void RightHandSide::vector_value (const Point<3> &p, Vector<double> &values)  const 
{
  values(0) = A*exp(-(p(1)*p(1)+p(2)*p(2))/B);
  values(1) = 0;
  values(2) = 0;
}

inline double RightHandSide::value (const Point<3> &p, const unsigned int component ) const
{
  if( component==1 || component==2 ) return 0;  
return A*exp(-(p(1)*p(1)+p(2)*p(2))/B);
}

class MySolver
{
public:
  MySolver();
  ~MySolver();

  void run ();

protected:
  void make_grid();
  void setup_system();
  void assemble_B();
  void assemble_E();
  void assemble_0();

  double dotprod(const Tensor<1,3> &A, const Tensor<1,3> &B) const;
  double dotprod(const Tensor<1,3> &A, const Vector<double> &B) const;

  void solve();
  void output_results ( string );
  void setup_boundary_ids();

  MPI_Comm mpi_communicator;
  parallel::distributed::Triangulation<3> triangulation;
  FE_Nedelec<3> fe;
  DoFHandler<3> dof_handler;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  ConstraintMatrix constraints;

  LA::MPI::SparseMatrix system_matrix;
  LA::MPI::Vector system_rhs;
  LA::MPI::Vector sol;
  LA::MPI::Vector m_B; 
  LA::MPI::Vector m_E; 
  LA::MPI::Vector m_J; // Vektorpotential
  LA::MPI::Vector m_rotJ; // Vektorpotential
  Vector<double> m_error_per_cell;

  ConditionalOStream pcout;

  bool m_root;
  double Lh;
};

MySolver::MySolver ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator, typename Triangulation<3>::MeshSmoothing(Triangulation<3>::smoothing_on_refinement|Triangulation<3>::smoothing_on_coarsening)),
  fe (0),
  dof_handler (triangulation),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  Lh = 2;
  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
}

MySolver::~MySolver ()
{
  dof_handler.clear ();
}

void MySolver::make_grid ()
{
  Point<3,double> pt1( -Lh, -Lh, -Lh );
  Point<3,double> pt2( Lh, Lh, Lh );

  GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
  triangulation.refine_global(5);
  
  setup_boundary_ids();
}

void MySolver::setup_system()
{
  dof_handler.distribute_dofs (fe);

  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

  sol.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  m_B.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  m_E.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  m_J.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  m_rotJ.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  m_error_per_cell.reinit(triangulation.n_active_cells());

  int myrank;
  MPI_Comm_rank( mpi_communicator, &myrank );
  cout << "(" << myrank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;

  system_rhs = 0;

  vector<bool> mask (dof_handler.get_fe().n_components(), true );

  RightHandSide right_hand_side;
  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  VectorTools::project_boundary_values_curl_conforming(dof_handler, 0, ZeroFunction<3>(fe.n_components()), 0, constraints);
  VectorTools::project_boundary_values_curl_conforming(dof_handler, 0, right_hand_side, 1, constraints);
  constraints.close ();

  CompressedSimpleSparsityPattern csp (locally_relevant_dofs);

  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

  system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);
}

void MySolver::setup_boundary_ids()
{
  parallel::distributed::Triangulation<3>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  for(; cell!=endc; ++cell)
  {
    for (unsigned int f=0; f < GeometryInfo<3>::faces_per_cell; ++f)
    {
      const Point<3> face_center = cell->face(f)->center();
      if( cell->face(f)->at_boundary() && !(fabs(face_center[0])==Lh) )    
      {
	cell->face(f)->set_all_boundary_indicators(1);
      }
    }
  }
}
  

double MySolver::dotprod(const Tensor<1,3> &A, const Tensor<1,3> &B) const
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

double MySolver::dotprod(const Tensor<1,3> &A, const Vector<double> &B) const
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void MySolver::assemble_0 ()
{
  const QGauss<3>  quadrature_formula(fe.degree+2);
  FEValues<3> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  FEValuesExtractors::Vector A(0);
  
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  RightHandSide right_hand_side;
  vector<Vector<double> > rhs_values (n_q_points,Vector<double>(3));
  vector<vector<Tensor<1,3>>> rhs_values_grads(n_q_points, vector<Tensor<1,3>>(3));  
  Tensor<1,3> value_i, value_j, rot_rhs;
  
  system_rhs=0;
  system_matrix=0;

  typename DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    if( cell->is_locally_owned() )
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);
      right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          value_i[0] = fe_values.shape_value_component(i,q_point,0);
          value_i[1] = fe_values.shape_value_component(i,q_point,1);
          value_i[2] = fe_values.shape_value_component(i,q_point,2);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
	    value_j[0] = fe_values.shape_value_component(j,q_point,0);
	    value_j[1] = fe_values.shape_value_component(j,q_point,1);
	    value_j[2] = fe_values.shape_value_component(j,q_point,2);
	    cell_matrix(i,j) +=  dotprod(value_i,value_j)*fe_values.JxW(q_point);
          }
	  cell_rhs(i) += dotprod(value_i,rhs_values[q_point])*fe_values.JxW(q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void MySolver::assemble_B ()
{
  const QGauss<3>  quadrature_formula(fe.degree+2);
  FEValues<3> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  FEValuesExtractors::Vector A(0);
  
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  RightHandSide right_hand_side;
  vector<Vector<double> > rhs_values (n_q_points,Vector<double>(3));
  vector<vector<Tensor<1,3>>> rhs_values_grads(n_q_points, vector<Tensor<1,3>>(3));  
  Tensor<1,3> value_i, rot_rhs;

  system_rhs=0;
  system_matrix=0;  
  
  typename DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    if( cell->is_locally_owned() )
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);
      right_hand_side.vector_gradient_list (fe_values.get_quadrature_points(), rhs_values_grads);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          value_i[0] = fe_values.shape_value_component(i,q_point,0);
          value_i[1] = fe_values.shape_value_component(i,q_point,1);
          value_i[2] = fe_values.shape_value_component(i,q_point,2);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
	    cell_matrix(i,j) += fe_values.JxW(q_point)*(fe_values[A].curl(i,q_point)*fe_values[A].curl(j,q_point)+fe_values[A].value(i,q_point)*fe_values[A].value(j,q_point));
          }
          rot_rhs[0] = rhs_values_grads[q_point][2][1]-rhs_values_grads[q_point][1][2];
	  rot_rhs[1] = rhs_values_grads[q_point][0][2]-rhs_values_grads[q_point][2][0];
	  rot_rhs[2] = rhs_values_grads[q_point][1][0]-rhs_values_grads[q_point][0][1];
	  
          cell_rhs(i) += dotprod(value_i,rot_rhs)*fe_values.JxW(q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void MySolver::assemble_E ()
{
  const QGauss<3>  quadrature_formula(fe.degree+2);
  FEValues<3> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  FEValuesExtractors::Vector A(0);
  
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  RightHandSide right_hand_side;
  vector<Vector<double> > rhs_values (n_q_points,Vector<double>(3));
  Tensor<1,3> value_i;

  system_rhs=0;
  system_matrix=0;  
  
  typename DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    if( cell->is_locally_owned() )
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);
      right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          value_i[0] = fe_values.shape_value_component(i,q_point,0);
          value_i[1] = fe_values.shape_value_component(i,q_point,1);
          value_i[2] = fe_values.shape_value_component(i,q_point,2);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
	    cell_matrix(i,j) += fe_values.JxW(q_point)*(fe_values[A].curl(i,q_point)*fe_values[A].curl(j,q_point)+fe_values[A].value(i,q_point)*fe_values[A].value(j,q_point));
          }
          cell_rhs(i) -= dotprod(value_i,rhs_values[q_point])*fe_values.JxW(q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void MySolver::solve ()
{
  sol=0;
//   SolverControl solver_control;
//   PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
//   solver.set_symmetric_mode(false);
//   solver.solve(system_matrix, sol, system_rhs);

  PETScWrappers::PreconditionSOR preconditioner;
  PETScWrappers::PreconditionSOR::AdditionalData data;
  preconditioner.initialize(system_matrix, data);
  SolverControl solver_control (dof_handler.n_dofs(), 1e-15);
  PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
  solver.solve(system_matrix, sol, system_rhs, preconditioner);
  
  
  constraints.distribute (sol);
}

void MySolver::run()
{
  pcout << "1" << endl;
  make_grid();
  setup_system();
  assemble_0();
  solve();
  m_J=sol;

  pcout << "2" << endl;
  assemble_E();
  solve();
  m_E=sol;
  
  pcout << "3" << endl;
  assemble_B();
  solve();
  m_B=sol;

  pcout << "4" << endl;
  output_results("res.vtu");
}

void MySolver::output_results ( string filename )
{
  vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(3, DataComponentInterpretation::component_is_part_of_vector);

  DataOut<3> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (m_B, "B", DataOut<3>::type_dof_data, data_component_interpretation);
  data_out.add_data_vector (m_E, "E", DataOut<3>::type_dof_data, data_component_interpretation);
  data_out.add_data_vector (m_J, "J", DataOut<3>::type_dof_data, data_component_interpretation);
  data_out.build_patches ();
  data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);
}
} // end of namespace

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  using namespace realtime_propagation;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MySolver solver;
    solver.run();
  }

  return EXIT_SUCCESS;
}
