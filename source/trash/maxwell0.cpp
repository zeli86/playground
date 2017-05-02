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
#include <deal.II/grid/manifold_lib.h>
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
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "mpi.h"

namespace realtime_propagation
{
enum Status { SUCCESS, FAILED };

using namespace std;
using namespace dealii;

class Crho : public Function<3>
{
public:
  Crho () : Function<3>(){}
  virtual double value ( const Point<3> &p, const unsigned component = 0) const;
};

double Crho::value( const Point<3> &p, const unsigned component ) const
{
  if( sqrt(p(0)*p(0)+p(1)*p(1)) < 0.5 ) return 1.0;
  return 0.0;
}

class ComputeCurl : public DataPostprocessorVector<3>
{
public:
  ComputeCurl ();
  virtual void compute_derived_quantities_vector (const vector<Vector<double>> &uh,
      const vector<vector<Tensor<1,3>>> &duh,
      const vector<vector<Tensor<2,3>>> &dduh,
      const vector<Point<3>>            &normals,
      const vector<Point<3>>            &evaluation_points,
      vector<Vector<double>>            &computed_quantities) const;
};

ComputeCurl::ComputeCurl () : DataPostprocessorVector<3> ( "B", update_gradients )
{}

void  ComputeCurl::compute_derived_quantities_vector (
  const vector<Vector<double>>      & /*uh*/,
  const vector<vector<Tensor<1,3>>> & duh,
  const vector<vector<Tensor<2,3>>> & /*dduh*/,
  const vector<Point<3>>            & /*normals*/,
  const vector<Point<3>>            & /*evaluation_points*/,
  vector<Vector<double>>            &computed_quantities
) const
{
  for (unsigned int q=0; q<computed_quantities.size(); q++)
  {
    computed_quantities[q](0) = (duh[q][2][1]-duh[q][1][2]);
    computed_quantities[q](1) = (duh[q][0][2]-duh[q][2][0]);
    computed_quantities[q](2) = (duh[q][1][0]-duh[q][0][1]);
  }
}
// RIGHT HAND SIDE CLASS
template <int dim>
class RightHandSide :  public Function<dim>
{
public:
  RightHandSide ();
  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &values) const;
  virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                  std::vector<Vector<double> >   &value_list) const;
private:
  const double PI = dealii::numbers::PI;
  const double bc_constant = 0.1;
};
template <int dim>
RightHandSide<dim>::RightHandSide () :
  Function<dim> (dim)
{}
template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
  Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());

  values(0) = 0;
  values(1) = 0;
  values(2) = 0; //1/(p(0)*p(0)+p(1)*p(1));
}
template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >   &value_list) const
{
  Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));
  const unsigned int n_points = points.size();
  for (unsigned int p=0; p<n_points; ++p)
  {
    RightHandSide<dim>::vector_value (points[p], value_list[p]);
  }
}
// END RIGHT HAND SIDE MEMBERS


class MySolver
{
public:
  MySolver();
  ~MySolver();

  void run ();

protected:
  void make_grid();
  void setup_system();
  void assemble();
  void assemble_2();
  void assemble_3();
  void assemble_4();

  double dotprod(const Tensor<1,3> &A, const Tensor<1,3> &B) const;
  double dotprod(const Tensor<1,3> &A, const Vector<double> &B) const;

  void solve();
  void solve_2();
  void output_results ( string );
  void output_results_2 ( string );

  MPI_Comm mpi_communicator;
  parallel::distributed::Triangulation<3> triangulation;
  //FESystem<3> fe;
  FE_Q<3> fe;
  FE_Nedelec<3> fe_2;
  DoFHandler<3> dof_handler;
  DoFHandler<3> dof_handler_2;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  IndexSet locally_owned_dofs_2;
  IndexSet locally_relevant_dofs_2;
  ConstraintMatrix constraints;
  ConstraintMatrix constraints_2;

  LA::MPI::SparseMatrix system_matrix;
  LA::MPI::SparseMatrix system_matrix_2;
  LA::MPI::Vector system_rhs;
  LA::MPI::Vector system_rhs_2;
  LA::MPI::Vector sol;
  LA::MPI::Vector sol_2;
  LA::MPI::Vector m_A; // Vektorpotential
  LA::MPI::Vector m_E; // Vektorpotential
  LA::MPI::Vector m_phi; // Vektorpotential
  LA::MPI::Vector m_rho; // Vektorpotential
  Vector<double> m_error_per_cell;

  ConditionalOStream pcout;

  bool m_root;
};

MySolver::MySolver ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator, typename Triangulation<3>::MeshSmoothing(Triangulation<3>::smoothing_on_refinement|Triangulation<3>::smoothing_on_coarsening)),
  fe (1),
  dof_handler (triangulation),
  fe_2(0),
  dof_handler_2 (triangulation),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
}

MySolver::~MySolver ()
{
  dof_handler.clear ();
  dof_handler_2.clear ();
}

void MySolver::make_grid ()
{
  Point<3,double> pt1( -1, -1, -1 );
  Point<3,double> pt2( 1, 1, 1 );

  GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
  triangulation.refine_global(5);

//     unsigned int tmp1[2], tmp2[2];
//     tmp1[0] = triangulation.n_cells();
//     tmp1[1] = triangulation.n_active_cells();
//     MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);
}

void MySolver::setup_system()
{
  dof_handler.distribute_dofs (fe);
  dof_handler_2.distribute_dofs (fe_2);

  m_error_per_cell.reinit(triangulation.n_active_cells());
  
  /// 2
  locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);
  sol_2.reinit(locally_owned_dofs_2, mpi_communicator);
  system_rhs_2.reinit(locally_owned_dofs_2, mpi_communicator);
  m_A.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);
  m_E.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

  constraints_2.clear ();
  constraints_2.reinit (locally_relevant_dofs_2);
  DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
  VectorTools::project_boundary_values_curl_conforming(dof_handler_2, 0, ZeroFunction<3>(fe_2.n_components()), 0, constraints_2);
  constraints_2.close ();

  CompressedSimpleSparsityPattern csp_2 (locally_relevant_dofs_2);
  DoFTools::make_sparsity_pattern (dof_handler_2, csp_2, constraints_2, false);
  SparsityTools::distribute_sparsity_pattern (csp_2, dof_handler_2.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs_2);
  system_matrix_2.reinit (locally_owned_dofs_2, locally_owned_dofs_2, csp_2, mpi_communicator);

  ///
  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

  sol.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  m_phi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  m_rho.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<3>(), constraints);  
  //VectorTools::interpolate_boundary_values (dof_handler, 2, ConstantFunction<3>(1.0), constraints);  
  //VectorTools::interpolate_boundary_values (dof_handler, 3, ConstantFunction<3>(1.0), constraints);  
  constraints.close ();

  CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
  system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);
  
  printf( "dof_handler.n_dofs() = %d\n", dof_handler.n_dofs());
  printf( "dof_handler_2.n_dofs() = %d\n", dof_handler_2.n_dofs());
  
  printf( "system_matrix: (%d,%d)\n", system_matrix.m(), system_matrix.n() );
  printf( "system_matrix_2: (%d,%d)\n", system_matrix_2.m(), system_matrix_2.n() );
}

double MySolver::dotprod(const Tensor<1,3> &A, const Tensor<1,3> &B) const
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

double MySolver::dotprod(const Tensor<1,3> &A, const Vector<double> &B) const
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void MySolver::assemble ()
{
  const QGauss<3>  quadrature_formula(fe.degree+2);
  FEValues<3> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs (dofs_per_cell);
  vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  Crho right_hand_side;
  double rhs;
  
  double JxW;
  
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
      
      int id = cell->material_id();
      if( id == 1 ) rhs=1.0;
      else rhs=0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
	JxW = fe_values.JxW(q_point);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
	    cell_matrix(i,j) += JxW*fe_values.shape_grad(i,q_point)*fe_values.shape_grad(j,q_point);
          }
          //cell_rhs(i) += JxW*right_hand_side.value(fe_values.quadrature_point(q_point))*fe_values.shape_value(i,q_point);
          cell_rhs(i) += JxW*rhs*fe_values.shape_value(i,q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void MySolver::assemble_2 ()
{
  const QGauss<3>  quadrature_formula(fe_2.degree+2);
  FEValues<3> fe_values (fe_2, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  Functions::FEFieldFunction<3,DoFHandler<3>,LA::MPI::Vector> funcE(dof_handler,m_phi);
  
  FEValuesExtractors::Vector A(0);
  
  const unsigned int   dofs_per_cell = fe_2.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  vector<Tensor<1,3>> E(n_q_points);
  Tensor<1,3> value_i;
  
  system_rhs_2=0;
  system_matrix_2=0;

  typename DoFHandler<3>::active_cell_iterator cell = dof_handler_2.begin_active(), endc = dof_handler_2.end();
  typename DoFHandler<3>::active_cell_iterator cellE = dof_handler.begin_active();
  for (; cell!=endc; ++cell, ++cellE)
  {
    if( cell->is_locally_owned() )
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);
      
      funcE.set_active_cell(cellE);
      funcE.gradient_list(fe_values.get_quadrature_points(),E);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          value_i[0] = fe_values.shape_value_component(i,q_point,0);
          value_i[1] = fe_values.shape_value_component(i,q_point,1);
          value_i[2] = fe_values.shape_value_component(i,q_point,2);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            cell_matrix(i,j) += fe_values.JxW(q_point)*fe_values[A].curl(i,q_point)*fe_values[A].curl(j,q_point);
          }
          cell_rhs(i) += dotprod(value_i,E[q_point])*fe_values.JxW(q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints_2.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix_2, system_rhs_2);
    }
  }
  system_rhs_2.compress(VectorOperation::add);
  system_matrix_2.compress(VectorOperation::add);
}

void MySolver::assemble_3 ()
{
  const QGauss<3>  quadrature_formula(fe_2.degree+2);
  FEValues<3> fe_values (fe_2, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  Functions::FEFieldFunction<3,DoFHandler<3>,LA::MPI::Vector> funcE(dof_handler,m_phi);
  
  FEValuesExtractors::Vector A(0);
  
  const unsigned int   dofs_per_cell = fe_2.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  vector<Tensor<1,3>> E(n_q_points);
  Tensor<1,3> value_i, value_j;

  system_rhs_2=0;
  system_matrix_2=0;
  
  typename DoFHandler<3>::active_cell_iterator cell = dof_handler_2.begin_active(), endc = dof_handler_2.end();
  typename DoFHandler<3>::active_cell_iterator cellE = dof_handler.begin_active();
  for (; cell!=endc; ++cell, ++cellE)
  {
    if( cell->is_locally_owned() )
    {
      cell_matrix = 0;
      cell_rhs = 0;
      fe_values.reinit (cell);
      
      funcE.set_active_cell(cellE);
      funcE.gradient_list(fe_values.get_quadrature_points(),E);

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
            cell_matrix(i,j) += fe_values.JxW(q_point)*dotprod(value_i,value_j);	    
          }
          cell_rhs(i) += dotprod(value_i,E[q_point])*fe_values.JxW(q_point);
        }
      }
      cell->get_dof_indices (local_dof_indices);
      constraints_2.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix_2, system_rhs_2);
    }
  }
  system_rhs_2.compress(VectorOperation::add);
  system_matrix_2.compress(VectorOperation::add);
}

void MySolver::assemble_4 ()
{
  const QGauss<3>  quadrature_formula(fe.degree+2);
  FEValues<3> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs (dofs_per_cell);
  vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  Crho right_hand_side;
  double rhs;
  
  double JxW;
  
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
      
      int id = cell->material_id();
      if( id == 1 ) rhs=1.0; else rhs=0; 

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
	JxW = fe_values.JxW(q_point);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
	    cell_matrix(i,j) -= JxW*fe_values.shape_value(i,q_point)*fe_values.shape_value(j,q_point);
          }
          cell_rhs(i) += JxW*rhs*fe_values.shape_value(i,q_point);
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
  SolverControl solver_control;
  PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
  solver.set_symmetric_mode(true);
  solver.solve(system_matrix, sol, system_rhs);

  constraints.distribute (sol);
  m_phi = sol;
}

void MySolver::solve_2 ()
{
  SolverControl solver_control;
  PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
  solver.set_symmetric_mode(false);
  solver.solve(system_matrix_2, sol_2, system_rhs_2);

  constraints_2.distribute (sol_2);
  m_A = sol_2;
}

void MySolver::run()
{
//   GridIn<3> gridin;
//   gridin.attach_triangulation(triangulation);
//   std::ifstream f("cube.msh");
//   gridin.read_msh(f);

//   if( m_root )
//   {
//     GridOut gridout;
//     std::ofstream fout("grid.msh");
//     gridout.write_msh(triangulation, fout);
//   }
  
  make_grid();

  setup_system();
  assemble_4();
  solve();
  m_rho=sol;
  assemble();
  solve();
  output_results("phi.vtu");
  
  assemble_3();
  solve_2();
  m_E=sol_2;

  assemble_2();
  solve_2();
  output_results_2("B.vtu");
}

void MySolver::output_results ( string filename )
{
  Crho rho;
  VectorTools::interpolate(dof_handler, rho, m_rho );
  
  DataOut<3> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (m_phi, "phi");
  data_out.add_data_vector (m_rho, "rho");
  data_out.build_patches ();
  data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);
}

void MySolver::output_results_2 ( string filename )
{
  ComputeCurl cc;

  vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(3, DataComponentInterpretation::component_is_part_of_vector);

  DataOut<3> data_out;
  data_out.attach_dof_handler (dof_handler_2);
  data_out.add_data_vector (m_A, "A", DataOut<3>::type_dof_data, data_component_interpretation);
  data_out.add_data_vector (m_E, "E", DataOut<3>::type_dof_data, data_component_interpretation);
  data_out.add_data_vector (m_A, cc);
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
